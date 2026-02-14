"""
KRONOS Model Loading and Inference

Wraps the KRONOS Vision Transformer for use within KINTSUGI.
Handles model loading from HuggingFace Hub or local cache,
device management, and batch inference.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import torch

logger = logging.getLogger(__name__)

# Default model configuration
DEFAULT_MODEL_TYPE = "vits16"
DEFAULT_PATCH_SIZE = 224  # Model's configured img_size (ViT-S/16 = 224)
DEFAULT_HF_REPO = "MahmoodLab/kronos"
MAX_INTENSITY_16BIT = 65535.0  # For [0,1] normalization of 16-bit images


@dataclass
class KronosConfig:
    """Configuration for KRONOS model loading and inference."""

    model_type: str = DEFAULT_MODEL_TYPE
    patch_size: int = DEFAULT_PATCH_SIZE
    token_overlap: bool = False
    cache_dir: str | Path = "./model_assets"
    hf_auth_token: str | None = None
    device: str = "auto"
    batch_size: int = 32

    @property
    def device_resolved(self) -> str:
        """Resolve 'auto' to actual device."""
        if self.device != "auto":
            return self.device
        try:
            import torch

            return "cuda" if torch.cuda.is_available() else "cpu"
        except ImportError:
            return "cpu"


@dataclass
class EmbeddingResult:
    """Container for KRONOS embedding outputs.

    Attributes
    ----------
    patch_embeddings : np.ndarray
        Shape (N, embedding_dim). CLS token — aggregated representation per patch.
        Best for clustering, search, and region classification.
    marker_embeddings : np.ndarray
        Shape (N, n_markers, embedding_dim). Per-marker feature vectors,
        spatially averaged. Best for marker-specific analyses.
    token_embeddings : np.ndarray
        Shape (N, n_markers, tokens_row, tokens_col, embedding_dim).
        Full spatial-molecular resolution. Best for fine-grained spatial analysis.
        With patch_size=224 and stride=16: tokens_row=tokens_col=14.
    patch_coords : np.ndarray
        Shape (N, 4). Patch coordinates as (row, col, height, width) in pixels.
    marker_names : list[str]
        Names of markers used for embedding extraction.
    metadata : dict
        Additional metadata (model type, device, timing, etc.).
    """

    patch_embeddings: np.ndarray
    marker_embeddings: np.ndarray
    token_embeddings: np.ndarray
    patch_coords: np.ndarray
    marker_names: list[str] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)

    @property
    def n_patches(self) -> int:
        return self.patch_embeddings.shape[0]

    @property
    def embedding_dim(self) -> int:
        return self.patch_embeddings.shape[1]

    def save(self, path: str | Path) -> None:
        """Save embeddings to HDF5 file."""
        import h5py

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(path, "w") as f:
            f.create_dataset("patch_embeddings", data=self.patch_embeddings)
            f.create_dataset("marker_embeddings", data=self.marker_embeddings)
            f.create_dataset("token_embeddings", data=self.token_embeddings)
            f.create_dataset("patch_coords", data=self.patch_coords)
            f.attrs["marker_names"] = self.marker_names
            for k, v in self.metadata.items():
                try:
                    f.attrs[k] = v
                except TypeError:
                    f.attrs[k] = str(v)
        logger.info("Saved embeddings to %s (%d patches)", path, self.n_patches)

    @classmethod
    def load(cls, path: str | Path) -> EmbeddingResult:
        """Load embeddings from HDF5 file."""
        import h5py

        with h5py.File(path, "r") as f:
            result = cls(
                patch_embeddings=f["patch_embeddings"][:],
                marker_embeddings=f["marker_embeddings"][:],
                token_embeddings=f["token_embeddings"][:],
                patch_coords=f["patch_coords"][:],
                marker_names=list(f.attrs.get("marker_names", [])),
                metadata={k: v for k, v in f.attrs.items() if k != "marker_names"},
            )
        logger.info("Loaded embeddings from %s (%d patches)", path, result.n_patches)
        return result


class KronosModel:
    """Wrapper for the KRONOS Vision Transformer model.

    Handles model loading, device management, and inference.

    Parameters
    ----------
    config : KronosConfig, optional
        Model configuration. Uses defaults if not provided.

    Example
    -------
    >>> model = KronosModel()
    >>> model.load()
    >>> embeddings = model.forward(batch_tensor, mean, std)
    """

    def __init__(self, config: KronosConfig | None = None):
        self.config = config or KronosConfig()
        self._model = None
        self._precision = None
        self._embedding_dim = None
        self._device = None

    @property
    def is_loaded(self) -> bool:
        return self._model is not None

    @property
    def precision(self) -> "torch.dtype | None":
        return self._precision

    @property
    def embedding_dim(self) -> int | None:
        return self._embedding_dim

    @property
    def device(self) -> "torch.device":
        if self._device is None:
            import torch

            self._device = torch.device(self.config.device_resolved)
        return self._device

    def load(self) -> None:
        """Load the KRONOS model from HuggingFace Hub or local cache.

        Raises
        ------
        ImportError
            If the kronos package is not installed.
        RuntimeError
            If model loading fails.
        """
        try:
            from kronos import create_model_from_pretrained
        except ImportError:
            raise ImportError(
                "KRONOS package not installed. Install with:\n"
                "  git clone https://github.com/mahmoodlab/KRONOS.git\n"
                "  cd KRONOS && pip install -e .\n"
                "Model weights require HuggingFace access: MahmoodLab/KRONOS"
            )

        import torch

        cache_dir = Path(self.config.cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)

        logger.info(
            "Loading KRONOS model (type=%s, device=%s)",
            self.config.model_type,
            self.config.device_resolved,
        )

        model, precision, embedding_dim = create_model_from_pretrained(
            checkpoint_path=f"hf_hub:{DEFAULT_HF_REPO}",
            cfg_path=None,
            hf_auth_token=self.config.hf_auth_token,
            cache_dir=str(cache_dir),
            cfg={
                "model_type": self.config.model_type,
                "token_overlap": self.config.token_overlap,
            },
        )

        self._device = torch.device(self.config.device_resolved)
        model.to(self._device)
        model.eval()

        self._model = model
        self._precision = precision
        self._embedding_dim = embedding_dim

        logger.info(
            "KRONOS model loaded: precision=%s, embedding_dim=%d, device=%s",
            precision,
            embedding_dim,
            self._device,
        )

    def forward(
        self,
        batch: "torch.Tensor",
        mean: "torch.Tensor",
        std: "torch.Tensor",
        marker_ids: "torch.Tensor | None" = None,
        is_raw_16bit: bool = True,
    ) -> tuple["torch.Tensor", "torch.Tensor", "torch.Tensor"]:
        """Run inference on a batch.

        Applies two-step normalization:
        1. Intensity normalization: raw 16-bit values → [0, 1] (divide by 65535)
        2. Z-score normalization: (x - mean) / std using marker_metadata.csv stats

        Parameters
        ----------
        batch : torch.Tensor
            Shape (B, M, H, W) where M is marker count.
        mean : torch.Tensor
            Per-marker mean for z-score normalization. Shape (M,).
        std : torch.Tensor
            Per-marker std for z-score normalization. Shape (M,).
        marker_ids : torch.Tensor, optional
            Per-marker IDs for sinusoidal positional encoding. Shape (M,).
            If None, KRONOS uses sequential IDs with +4 offset.
            For best results, pass actual IDs from marker_metadata.csv.
        is_raw_16bit : bool
            If True, divides by 65535 before z-score normalization.
            Set False if data is already in [0, 1] range.

        Returns
        -------
        tuple of torch.Tensor
            (patch_embeddings, marker_embeddings, token_embeddings)
        """
        import torch

        if not self.is_loaded:
            raise RuntimeError("Model not loaded. Call model.load() first.")

        batch = batch.to(device=self.device, dtype=self._precision)
        mean = mean.to(device=self.device, dtype=self._precision)
        std = std.to(device=self.device, dtype=self._precision)

        # Step 1: Intensity normalization to [0, 1]
        if is_raw_16bit:
            batch = batch / MAX_INTENSITY_16BIT

        # Step 2: Z-score normalization using KRONOS marker statistics
        batch = (batch - mean[None, :, None, None]) / (std[None, :, None, None] + 1e-8)

        # Build marker_ids list (one tensor per batch item) if provided
        batch_marker_ids = None
        if marker_ids is not None:
            marker_ids = marker_ids.to(device=self.device)
            batch_marker_ids = [marker_ids for _ in range(batch.shape[0])]

        with torch.no_grad():
            if batch_marker_ids is not None:
                patch_emb, marker_emb, token_emb = self._model(
                    batch, marker_ids=batch_marker_ids
                )
            else:
                patch_emb, marker_emb, token_emb = self._model(batch)

        return patch_emb, marker_emb, token_emb

    def embed_batched(
        self,
        patches: "torch.Tensor",
        mean: "torch.Tensor",
        std: "torch.Tensor",
        marker_ids: "torch.Tensor | None" = None,
        batch_size: int | None = None,
        is_raw_16bit: bool = True,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Run batched inference, returning numpy arrays.

        Parameters
        ----------
        patches : torch.Tensor
            Shape (N, M, H, W). Raw 16-bit values if is_raw_16bit=True.
        mean : torch.Tensor
            Per-marker mean. Shape (M,).
        std : torch.Tensor
            Per-marker std. Shape (M,).
        marker_ids : torch.Tensor, optional
            Per-marker IDs from marker_metadata.csv. Shape (M,).
        batch_size : int, optional
            Override config batch_size.
        is_raw_16bit : bool
            If True, divides by 65535 before z-score normalization.

        Returns
        -------
        tuple of np.ndarray
            (patch_embeddings, marker_embeddings, token_embeddings)
        """
        import torch

        bs = batch_size or self.config.batch_size
        n_patches = patches.shape[0]

        all_patch = []
        all_marker = []
        all_token = []

        for i in range(0, n_patches, bs):
            batch = patches[i : i + bs]
            p_emb, m_emb, t_emb = self.forward(
                batch, mean, std,
                marker_ids=marker_ids,
                is_raw_16bit=is_raw_16bit,
            )
            all_patch.append(p_emb.cpu().numpy())
            all_marker.append(m_emb.cpu().numpy())
            # Flatten token embeddings to 3D for concatenation:
            # (B, M, R, C, D) -> (B, M*R*C, D)
            t_flat = t_emb.reshape(t_emb.shape[0], -1, t_emb.shape[-1])
            all_token.append(t_flat.cpu().numpy())

            if (i + bs) % (bs * 10) == 0:
                logger.info("  Processed %d / %d patches", min(i + bs, n_patches), n_patches)

        return (
            np.concatenate(all_patch, axis=0),
            np.concatenate(all_marker, axis=0),
            np.concatenate(all_token, axis=0),
        )
