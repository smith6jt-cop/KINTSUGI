"""
KRONOS Embedding Extraction from KINTSUGI Projects

Loads registered, marker-named 2D TIFFs from a KINTSUGI project,
tiles them into 256x256 patches, normalizes per-marker, and runs
KRONOS inference to produce embeddings.

The pipeline:
1. Load registered EDF images (marker-named TIFFs per cycle)
2. Stack into (markers, H, W) tensor
3. Tile into overlapping 256x256 patches
4. Z-score normalize per marker
5. Run KRONOS inference → patch/marker/token embeddings
6. Save to HDF5 for downstream analysis
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

import numpy as np

from .markers import MarkerMapper, MarkerMapping
from .model import EmbeddingResult, KronosConfig, KronosModel

logger = logging.getLogger(__name__)

# Default tiling parameters
DEFAULT_PATCH_SIZE = 256
DEFAULT_OVERLAP = 64  # pixels of overlap between adjacent patches


class KronosEmbedder:
    """Extract KRONOS embeddings from a KINTSUGI project.

    Parameters
    ----------
    config : KronosConfig, optional
        Model configuration.
    model : KronosModel, optional
        Pre-loaded model instance. If None, creates and loads one.

    Example
    -------
    >>> embedder = KronosEmbedder()
    >>> mapper = MarkerMapper.from_project("./my_project")
    >>> result = embedder.embed_project("./my_project", mapper)
    >>> result.save("./my_project/data/processed/embeddings/kronos.h5")
    """

    def __init__(
        self,
        config: KronosConfig | None = None,
        model: KronosModel | None = None,
    ):
        self.config = config or KronosConfig()
        self.model = model
        self._ensure_model()

    def _ensure_model(self) -> None:
        """Load model if not already loaded."""
        if self.model is None:
            self.model = KronosModel(self.config)
        if not self.model.is_loaded:
            self.model.load()

    def embed_project(
        self,
        project_dir: str | Path,
        mapper: MarkerMapper,
        patch_size: int = DEFAULT_PATCH_SIZE,
        overlap: int = DEFAULT_OVERLAP,
        output_path: str | Path | None = None,
    ) -> EmbeddingResult:
        """Extract embeddings from an entire KINTSUGI project.

        Parameters
        ----------
        project_dir : Path
            Path to KINTSUGI project root.
        mapper : MarkerMapper
            Configured marker mapper with channel-to-KRONOS mappings.
        patch_size : int
            Patch size for tiling (default 224 for KRONOS ViT-S/16).
        overlap : int
            Overlap between adjacent patches in pixels.
        output_path : Path, optional
            If provided, save embeddings to this HDF5 path.

        Returns
        -------
        EmbeddingResult
        """
        t0 = time.perf_counter()

        # Load and stack marker images
        logger.info("Loading registered images from %s", project_dir)
        image_stack, used_mappings = self._load_marker_images(project_dir, mapper)
        logger.info(
            "Loaded image stack: shape=%s, markers=%d",
            image_stack.shape,
            len(used_mappings),
        )

        # Tile into patches
        patches, coords = self._tile_image(image_stack, patch_size, overlap)
        logger.info("Created %d patches of size %dx%d", len(patches), patch_size, patch_size)

        # Compute dataset-level stats for any unmatched markers
        unmatched = [m for m in used_mappings if not m.matched]
        if unmatched:
            logger.info("Computing dataset-level normalization for %d unmatched markers", len(unmatched))
            self.compute_dataset_normalization(patches, used_mappings)

        # Get normalization parameters and marker IDs
        mean, std = mapper.get_normalization_arrays(used_mappings)
        marker_ids_np = mapper.get_marker_ids(used_mappings)

        # Run inference
        import torch

        patches_tensor = torch.from_numpy(patches).float()
        mean_tensor = torch.from_numpy(mean).float()
        std_tensor = torch.from_numpy(std).float()
        marker_ids_tensor = (
            torch.from_numpy(marker_ids_np).long() if marker_ids_np is not None else None
        )

        logger.info("Running KRONOS inference on %d patches...", len(patches))
        patch_emb, marker_emb, token_emb = self.model.embed_batched(
            patches_tensor, mean_tensor, std_tensor,
            marker_ids=marker_ids_tensor,
            is_raw_16bit=True,  # KINTSUGI images are 16-bit
        )

        elapsed = time.perf_counter() - t0
        marker_names = [m.kintsugi_name for m in used_mappings]

        result = EmbeddingResult(
            patch_embeddings=patch_emb,
            marker_embeddings=marker_emb,
            token_embeddings=token_emb,
            patch_coords=coords,
            marker_names=marker_names,
            metadata={
                "project_dir": str(project_dir),
                "n_markers": len(used_mappings),
                "patch_size": patch_size,
                "overlap": overlap,
                "image_shape": list(image_stack.shape),
                "model_type": self.config.model_type,
                "elapsed_seconds": elapsed,
            },
        )

        if output_path is not None:
            result.save(output_path)

        logger.info("Embedding extraction complete in %.1fs", elapsed)
        return result

    def embed_region(
        self,
        images: dict[str, np.ndarray],
        mapper: MarkerMapper,
        patch_size: int = DEFAULT_PATCH_SIZE,
        overlap: int = DEFAULT_OVERLAP,
    ) -> EmbeddingResult:
        """Extract embeddings from a set of marker images (arbitrary region).

        Parameters
        ----------
        images : dict
            {marker_name: 2D np.ndarray} mapping.
        mapper : MarkerMapper
            Configured marker mapper.
        patch_size : int
            Patch size for tiling.
        overlap : int
            Overlap between patches.

        Returns
        -------
        EmbeddingResult
        """
        # Build stack from provided images
        used_mappings = []
        arrays = []
        for m in mapper.mappings:
            if m.kintsugi_name in images:
                arrays.append(images[m.kintsugi_name].astype(np.float32))
                used_mappings.append(m)

        if not arrays:
            raise ValueError("No matching marker images found.")

        image_stack = np.stack(arrays, axis=0)  # (M, H, W)
        patches, coords = self._tile_image(image_stack, patch_size, overlap)

        mean, std = mapper.get_normalization_arrays(used_mappings)
        marker_ids_np = mapper.get_marker_ids(used_mappings)

        import torch

        patches_tensor = torch.from_numpy(patches).float()
        mean_tensor = torch.from_numpy(mean).float()
        std_tensor = torch.from_numpy(std).float()
        marker_ids_tensor = (
            torch.from_numpy(marker_ids_np).long() if marker_ids_np is not None else None
        )

        patch_emb, marker_emb, token_emb = self.model.embed_batched(
            patches_tensor, mean_tensor, std_tensor,
            marker_ids=marker_ids_tensor,
            is_raw_16bit=True,
        )

        return EmbeddingResult(
            patch_embeddings=patch_emb,
            marker_embeddings=marker_emb,
            token_embeddings=token_emb,
            patch_coords=coords,
            marker_names=[m.kintsugi_name for m in used_mappings],
        )

    def _load_marker_images(
        self,
        project_dir: str | Path,
        mapper: MarkerMapper,
    ) -> tuple[np.ndarray, list[MarkerMapping]]:
        """Load registered images and stack into (M, H, W) array.

        Searches for marker-named TIFFs in the registered/ and edf/ directories.

        Returns
        -------
        tuple
            (image_stack, used_mappings) — stack shape is (M, H, W).
        """
        from kintsugi.project import KintsugiProject

        project = KintsugiProject.load(project_dir)

        # Look for registered images first, then EDF
        search_dirs = []
        if project.registered_dir.exists():
            search_dirs.append(project.registered_dir)
        if project.edf_dir.exists():
            search_dirs.append(project.edf_dir)

        if not search_dirs:
            raise FileNotFoundError(
                f"No registered or EDF images found in {project_dir}. "
                "Run the processing pipeline first."
            )

        # Build index of available image files
        available_files: dict[str, Path] = {}
        for search_dir in search_dirs:
            for tif in sorted(search_dir.rglob("*.tif")):
                stem = tif.stem
                if stem not in available_files:
                    available_files[stem] = tif
            for tif in sorted(search_dir.rglob("*.tiff")):
                stem = tif.stem
                if stem not in available_files:
                    available_files[stem] = tif

        # Match mappings to available files
        used_mappings = []
        image_paths = []

        for m in mapper.mappings:
            # Try exact name match
            if m.kintsugi_name in available_files:
                image_paths.append(available_files[m.kintsugi_name])
                used_mappings.append(m)
                continue

            # Try case-insensitive match
            lower_name = m.kintsugi_name.lower()
            for stem, path in available_files.items():
                if stem.lower() == lower_name:
                    image_paths.append(path)
                    used_mappings.append(m)
                    break

        if not used_mappings:
            available = list(available_files.keys())[:20]
            raise FileNotFoundError(
                f"No matching marker images found. "
                f"Available files: {available}...\n"
                f"Expected markers: {[m.kintsugi_name for m in mapper.mappings[:10]]}..."
            )

        logger.info(
            "Found %d / %d marker images",
            len(used_mappings),
            len(mapper.mappings),
        )

        # Load images
        try:
            import tifffile
        except ImportError:
            raise ImportError("tifffile required for loading images: pip install tifffile")

        images = []
        ref_shape = None
        for path in image_paths:
            img = tifffile.imread(str(path)).astype(np.float32)
            # Handle multi-page TIFFs (take first page / max projection)
            if img.ndim == 3:
                img = img.max(axis=0)
            if ref_shape is None:
                ref_shape = img.shape
            elif img.shape != ref_shape:
                logger.warning(
                    "Image %s shape %s != reference %s, resizing",
                    path.name,
                    img.shape,
                    ref_shape,
                )
                from skimage.transform import resize

                img = resize(img, ref_shape, preserve_range=True).astype(np.float32)
            images.append(img)

        image_stack = np.stack(images, axis=0)  # (M, H, W)
        return image_stack, used_mappings

    def _tile_image(
        self,
        image_stack: np.ndarray,
        patch_size: int,
        overlap: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Tile a (M, H, W) image stack into (N, M, patch_size, patch_size) patches.

        Parameters
        ----------
        image_stack : np.ndarray
            Shape (M, H, W).
        patch_size : int
            Size of each square patch.
        overlap : int
            Overlap between adjacent patches.

        Returns
        -------
        tuple
            (patches, coords) where patches has shape (N, M, H, W) and
            coords has shape (N, 4) with (row, col, height, width).
        """
        _, h, w = image_stack.shape
        stride = patch_size - overlap

        # Calculate grid positions
        row_positions = list(range(0, max(h - patch_size + 1, 1), stride))
        col_positions = list(range(0, max(w - patch_size + 1, 1), stride))

        # Ensure we cover the edges
        if row_positions[-1] + patch_size < h:
            row_positions.append(h - patch_size)
        if col_positions[-1] + patch_size < w:
            col_positions.append(w - patch_size)

        patches = []
        coords = []

        for r in row_positions:
            for c in col_positions:
                patch = image_stack[:, r : r + patch_size, c : c + patch_size]
                # Pad if needed (edge patches)
                if patch.shape[1] != patch_size or patch.shape[2] != patch_size:
                    padded = np.zeros(
                        (image_stack.shape[0], patch_size, patch_size),
                        dtype=image_stack.dtype,
                    )
                    padded[:, : patch.shape[1], : patch.shape[2]] = patch
                    patch = padded
                patches.append(patch)
                coords.append([r, c, patch_size, patch_size])

        return np.stack(patches, axis=0), np.array(coords)

    @staticmethod
    def compute_dataset_normalization(
        patches: np.ndarray,
        mappings: list[MarkerMapping],
        max_intensity: float = 65535.0,
    ) -> None:
        """Compute dataset-level normalization stats for unmatched markers.

        Updates the MarkerMapping objects in-place with mean/std computed
        on [0, 1]-normalized data (after dividing by max_intensity).
        Only updates unmatched markers; matched markers keep their
        KRONOS metadata values.

        Parameters
        ----------
        patches : np.ndarray
            Shape (N, M, H, W). Raw 16-bit values.
        mappings : list of MarkerMapping
            Mappings to update (modified in-place).
        max_intensity : float
            Maximum intensity for [0, 1] normalization.
        """
        for i, m in enumerate(mappings):
            if not m.matched:
                # Compute stats on [0, 1]-normalized data to match KRONOS convention
                marker_data = patches[:, i, :, :].astype(np.float32) / max_intensity
                m.mean = float(np.mean(marker_data))
                m.std = float(np.std(marker_data))
                if m.std < 1e-8:
                    m.std = 1.0
                logger.debug(
                    "Dataset-level stats for %s: mean=%.4f, std=%.4f",
                    m.kintsugi_name,
                    m.mean,
                    m.std,
                )
