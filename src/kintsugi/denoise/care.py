"""
CARE-style Denoising

Pure Python implementation of Content-Aware Image Restoration (CARE)
for fluorescence microscopy denoising.

CARE uses supervised learning with paired noisy-clean training data,
or pre-trained models for common fluorescence microscopy scenarios.

Based on: Weigert et al., "Content-aware image restoration: pushing the limits
of fluorescence microscopy" Nature Methods, 2018.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Union
from urllib.request import urlretrieve

import numpy as np

if TYPE_CHECKING:
    import dask.array

logger = logging.getLogger("kintsugi.denoise.care")

# Type alias
ArrayLike = Union[np.ndarray, "dask.array.Array"]

# Pre-trained model registry
PRETRAINED_MODELS = {
    "fluorescence_2d": {
        "description": "General 2D fluorescence microscopy denoising",
        "url": None,  # Placeholder - would host models
        "input_shape": (None, None, 1),
        "trained_on": "Diverse fluorescence images",
    },
    "confocal_nuclei": {
        "description": "Confocal microscopy nuclear staining (DAPI/Hoechst)",
        "url": None,
        "input_shape": (None, None, 1),
        "trained_on": "Nuclear stains from confocal microscopy",
    },
    "widefield_membrane": {
        "description": "Widefield microscopy membrane markers",
        "url": None,
        "input_shape": (None, None, 1),
        "trained_on": "Membrane markers from widefield microscopy",
    },
}


@dataclass
class CAREConfig:
    """Configuration for CARE training and inference."""

    # Network architecture
    n_filters_base: int = 32
    n_depth: int = 4
    kernel_size: int = 3
    batch_norm: bool = True
    residual: bool = True  # Learn residual (noise) instead of clean image

    # Training
    patch_size: int = 64
    batch_size: int = 32
    n_epochs: int = 100
    learning_rate: float = 1e-4
    augment: bool = True

    # Loss
    loss_type: str = "mae"  # "mae", "mse", "ssim", "mixed"

    # Device
    device: str = "auto"


class CAREDenoiser:
    """
    CARE denoiser with training and inference.

    CARE uses a U-Net architecture trained on paired noisy-clean images.
    Can also use pre-trained models or be trained on synthetic noise.

    Example
    -------
    >>> denoiser = CAREDenoiser()
    >>> denoiser.train(noisy_images, clean_images, epochs=100)
    >>> denoised = denoiser.predict(noisy_image)

    >>> # Or use pre-trained model
    >>> denoiser = CAREDenoiser.from_pretrained("fluorescence_2d")
    >>> denoised = denoiser.predict(noisy_image)
    """

    def __init__(self, config: CAREConfig | None = None):
        """
        Initialize CARE denoiser.

        Parameters
        ----------
        config : CAREConfig, optional
            Configuration for the denoiser
        """
        self.config = config or CAREConfig()
        self.model = None
        self.device = None
        self.trained = False
        self._mean = 0.0
        self._std = 1.0

    @classmethod
    def from_pretrained(
        cls,
        model_name: str,
        cache_dir: str | Path | None = None,
    ) -> CAREDenoiser:
        """
        Load a pre-trained CARE model.

        Parameters
        ----------
        model_name : str
            Name of the pre-trained model (see list_pretrained_models())
        cache_dir : str or Path, optional
            Directory to cache downloaded models

        Returns
        -------
        CAREDenoiser
            Loaded denoiser ready for inference
        """
        if model_name not in PRETRAINED_MODELS:
            raise ValueError(
                f"Unknown model: {model_name}. " f"Available: {list(PRETRAINED_MODELS.keys())}"
            )

        model_info = PRETRAINED_MODELS[model_name]

        if model_info["url"] is None:
            raise NotImplementedError(
                f"Pre-trained model '{model_name}' is not yet available. "
                "Train a custom model or use another method."
            )

        # Download model if needed
        if cache_dir is None:
            cache_dir = Path.home() / ".kintsugi" / "models" / "care"
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)

        model_path = cache_dir / f"{model_name}.pth"

        if not model_path.exists():
            logger.info(f"Downloading {model_name}...")
            urlretrieve(model_info["url"], model_path)

        denoiser = cls()
        denoiser.load(model_path)

        return denoiser

    def _setup_device(self):
        """Set up PyTorch device with multi-GPU support."""
        try:
            import torch

            from kintsugi.gpu import get_gpu_manager

            self._gpu_manager = get_gpu_manager()

            if self.config.device == "auto":
                self.device = self._gpu_manager.get_torch_device()
            else:
                self.device = torch.device(self.config.device)

            if self._gpu_manager.device_count > 1:
                logger.info(
                    f"CARE using {self._gpu_manager.device_count} GPUs: "
                    f"{self._gpu_manager.device_ids}"
                )
            else:
                logger.info(f"CARE using device: {self.device}")

        except ImportError:
            raise ImportError(
                "PyTorch is required for CARE denoising. " "Install with: pip install torch"
            )

    def _build_model(self, n_channels: int = 1):
        """Build the CARE U-Net model."""
        import torch
        import torch.nn as nn

        class ConvBlock(nn.Module):
            def __init__(self, in_ch, out_ch, kernel_size=3, batch_norm=True):
                super().__init__()
                padding = kernel_size // 2
                layers = [
                    nn.Conv2d(in_ch, out_ch, kernel_size, padding=padding),
                    nn.LeakyReLU(0.1, inplace=True),
                ]
                if batch_norm:
                    layers.insert(1, nn.BatchNorm2d(out_ch))
                layers.extend(
                    [
                        nn.Conv2d(out_ch, out_ch, kernel_size, padding=padding),
                        nn.LeakyReLU(0.1, inplace=True),
                    ]
                )
                if batch_norm:
                    layers.insert(-1, nn.BatchNorm2d(out_ch))
                self.block = nn.Sequential(*layers)

            def forward(self, x):
                return self.block(x)

        class CAREUNet(nn.Module):
            def __init__(
                self, n_channels, n_filters_base, n_depth, kernel_size, batch_norm, residual
            ):
                super().__init__()
                self.residual = residual
                self.n_depth = n_depth

                # Encoder
                self.encoders = nn.ModuleList()
                self.pools = nn.ModuleList()
                in_ch = n_channels
                out_ch = n_filters_base
                for _i in range(n_depth):
                    self.encoders.append(ConvBlock(in_ch, out_ch, kernel_size, batch_norm))
                    self.pools.append(nn.MaxPool2d(2))
                    in_ch = out_ch
                    out_ch = min(out_ch * 2, 512)

                # Bottleneck
                self.bottleneck = ConvBlock(in_ch, out_ch, kernel_size, batch_norm)

                # Decoder
                self.upconvs = nn.ModuleList()
                self.decoders = nn.ModuleList()
                for i in range(n_depth):
                    self.upconvs.append(nn.ConvTranspose2d(out_ch, out_ch // 2, 2, stride=2))
                    # Concatenate with skip connection
                    decoder_in_ch = (
                        out_ch // 2 + self.encoders[n_depth - 1 - i].block[0].out_channels
                    )
                    self.decoders.append(
                        ConvBlock(decoder_in_ch, out_ch // 2, kernel_size, batch_norm)
                    )
                    out_ch //= 2

                # Output - predict either clean image or residual
                self.output = nn.Conv2d(out_ch, n_channels, 1)

            def forward(self, x):
                input_x = x

                # Encoder path
                enc_features = []
                for encoder, pool in zip(self.encoders, self.pools, strict=False):
                    x = encoder(x)
                    enc_features.append(x)
                    x = pool(x)

                # Bottleneck
                x = self.bottleneck(x)

                # Decoder path
                for i, (upconv, decoder) in enumerate(
                    zip(self.upconvs, self.decoders, strict=False)
                ):
                    x = upconv(x)
                    enc_feat = enc_features[self.n_depth - 1 - i]
                    # Handle size mismatch
                    if x.shape[2:] != enc_feat.shape[2:]:
                        x = nn.functional.interpolate(
                            x, size=enc_feat.shape[2:], mode="bilinear", align_corners=False
                        )
                    x = torch.cat([x, enc_feat], dim=1)
                    x = decoder(x)

                output = self.output(x)

                if self.residual:
                    return input_x - output  # Learn noise, output clean
                return output

        model = CAREUNet(
            n_channels=n_channels,
            n_filters_base=self.config.n_filters_base,
            n_depth=self.config.n_depth,
            kernel_size=self.config.kernel_size,
            batch_norm=self.config.batch_norm,
            residual=self.config.residual,
        )

        # Use GPUManager for multi-GPU wrapping
        if hasattr(self, "_gpu_manager"):
            return self._gpu_manager.wrap_model(model)
        return model.to(self.device)

    def _augment_batch(
        self,
        noisy: np.ndarray,
        clean: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Apply random augmentations to training batch."""
        # Random flips
        if np.random.rand() > 0.5:
            noisy = noisy[:, :, ::-1, :]
            clean = clean[:, :, ::-1, :]
        if np.random.rand() > 0.5:
            noisy = noisy[:, :, :, ::-1]
            clean = clean[:, :, :, ::-1]

        # Random 90-degree rotations
        k = np.random.randint(4)
        if k > 0:
            noisy = np.rot90(noisy, k, axes=(2, 3))
            clean = np.rot90(clean, k, axes=(2, 3))

        return noisy.copy(), clean.copy()

    def _extract_paired_patches(
        self,
        noisy_images: list[np.ndarray],
        clean_images: list[np.ndarray],
        n_patches_per_image: int = 128,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Extract paired patches from noisy and clean images."""
        noisy_patches = []
        clean_patches = []
        patch_size = self.config.patch_size

        for noisy, clean in zip(noisy_images, clean_images, strict=False):
            if noisy.ndim == 2:
                noisy = noisy[np.newaxis, :, :]
                clean = clean[np.newaxis, :, :]

            c, h, w = noisy.shape

            for _ in range(n_patches_per_image):
                if h > patch_size and w > patch_size:
                    y = np.random.randint(0, h - patch_size)
                    x = np.random.randint(0, w - patch_size)
                    noisy_patch = noisy[:, y : y + patch_size, x : x + patch_size]
                    clean_patch = clean[:, y : y + patch_size, x : x + patch_size]
                else:
                    # Pad if needed
                    pad_h = max(0, patch_size - h)
                    pad_w = max(0, patch_size - w)
                    noisy_padded = np.pad(noisy, ((0, 0), (0, pad_h), (0, pad_w)), mode="reflect")
                    clean_padded = np.pad(clean, ((0, 0), (0, pad_h), (0, pad_w)), mode="reflect")
                    noisy_patch = noisy_padded[:, :patch_size, :patch_size]
                    clean_patch = clean_padded[:, :patch_size, :patch_size]

                noisy_patches.append(noisy_patch)
                clean_patches.append(clean_patch)

        return np.stack(noisy_patches), np.stack(clean_patches)

    def train(
        self,
        noisy_images: list[np.ndarray],
        clean_images: list[np.ndarray],
        epochs: int | None = None,
        verbose: bool = True,
    ) -> dict[str, Any]:
        """
        Train the CARE model on paired noisy-clean images.

        Parameters
        ----------
        noisy_images : list of np.ndarray
            List of noisy training images
        clean_images : list of np.ndarray
            List of corresponding clean training images
        epochs : int, optional
            Number of training epochs
        verbose : bool
            Print training progress

        Returns
        -------
        dict
            Training history
        """
        import torch

        if len(noisy_images) != len(clean_images):
            raise ValueError("Number of noisy and clean images must match")

        self._setup_device()

        # Normalize
        all_pixels = np.concatenate([img.ravel() for img in noisy_images + clean_images])
        self._mean = float(np.mean(all_pixels))
        self._std = float(np.std(all_pixels))
        if self._std < 1e-8:
            self._std = 1.0

        noisy_norm = [(img - self._mean) / self._std for img in noisy_images]
        clean_norm = [(img - self._mean) / self._std for img in clean_images]

        # Determine channels
        n_channels = 1 if noisy_images[0].ndim == 2 else noisy_images[0].shape[0]

        # Build model
        self.model = self._build_model(n_channels)

        # Extract patches
        if verbose:
            logger.info("Extracting training patches...")
        noisy_patches, clean_patches = self._extract_paired_patches(noisy_norm, clean_norm)

        if verbose:
            logger.info(f"Extracted {len(noisy_patches)} patch pairs")

        # Training setup
        n_epochs = epochs or self.config.n_epochs
        optimizer = torch.optim.Adam(self.model.parameters(), lr=self.config.learning_rate)
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode="min", factor=0.5, patience=10
        )

        # Loss function
        if self.config.loss_type == "mae":
            criterion = torch.nn.L1Loss()
        elif self.config.loss_type == "mse":
            criterion = torch.nn.MSELoss()
        else:
            criterion = torch.nn.L1Loss()

        history = {"loss": [], "val_loss": []}

        # Split for validation
        n_val = max(1, len(noisy_patches) // 10)
        indices = np.random.permutation(len(noisy_patches))
        train_idx = indices[n_val:]
        val_idx = indices[:n_val]

        for epoch in range(n_epochs):
            # Training
            self.model.train()
            np.random.shuffle(train_idx)
            epoch_loss = 0.0
            n_batches = 0

            for i in range(0, len(train_idx), self.config.batch_size):
                batch_idx = train_idx[i : i + self.config.batch_size]
                noisy_batch = noisy_patches[batch_idx]
                clean_batch = clean_patches[batch_idx]

                if self.config.augment:
                    noisy_batch, clean_batch = self._augment_batch(noisy_batch, clean_batch)

                noisy_tensor = torch.from_numpy(noisy_batch).float().to(self.device)
                clean_tensor = torch.from_numpy(clean_batch).float().to(self.device)

                optimizer.zero_grad()
                output = self.model(noisy_tensor)
                loss = criterion(output, clean_tensor)
                loss.backward()
                optimizer.step()

                epoch_loss += loss.item()
                n_batches += 1

            avg_train_loss = epoch_loss / n_batches
            history["loss"].append(avg_train_loss)

            # Validation
            self.model.eval()
            with torch.no_grad():
                val_noisy = torch.from_numpy(noisy_patches[val_idx]).float().to(self.device)
                val_clean = torch.from_numpy(clean_patches[val_idx]).float().to(self.device)
                val_output = self.model(val_noisy)
                val_loss = criterion(val_output, val_clean).item()

            history["val_loss"].append(val_loss)
            scheduler.step(val_loss)

            if verbose and (epoch + 1) % 10 == 0:
                lr = optimizer.param_groups[0]["lr"]
                logger.info(
                    f"Epoch {epoch + 1}/{n_epochs}, "
                    f"Loss: {avg_train_loss:.6f}, Val: {val_loss:.6f}, LR: {lr:.2e}"
                )

        self.trained = True

        if verbose:
            logger.info("Training complete!")

        return history

    def train_synthetic(
        self,
        clean_images: list[np.ndarray],
        noise_model: str = "gaussian",
        noise_level: float = 0.1,
        epochs: int | None = None,
        verbose: bool = True,
    ) -> dict[str, Any]:
        """
        Train using synthetic noise added to clean images.

        Useful when clean images are available but not paired noisy images.

        Parameters
        ----------
        clean_images : list of np.ndarray
            List of clean training images
        noise_model : str
            Type of synthetic noise: "gaussian", "poisson", "mixed"
        noise_level : float
            Noise level (interpretation depends on noise_model)
        epochs : int, optional
            Number of training epochs
        verbose : bool
            Print training progress

        Returns
        -------
        dict
            Training history
        """
        noisy_images = []

        for clean in clean_images:
            clean_float = clean.astype(np.float64)

            if noise_model == "gaussian":
                noise = np.random.normal(0, noise_level * clean_float.std(), clean.shape)
                noisy = clean_float + noise

            elif noise_model == "poisson":
                # Scale to have reasonable count rate
                scale = 1.0 / (noise_level + 1e-8)
                noisy = np.random.poisson(clean_float * scale) / scale

            elif noise_model == "mixed":
                # Poisson + Gaussian (common in fluorescence)
                scale = 1.0 / (noise_level + 1e-8)
                poisson = np.random.poisson(clean_float * scale) / scale
                gaussian = np.random.normal(0, noise_level * clean_float.std() * 0.5, clean.shape)
                noisy = poisson + gaussian

            else:
                raise ValueError(f"Unknown noise model: {noise_model}")

            noisy_images.append(np.clip(noisy, clean.min(), clean.max()).astype(clean.dtype))

        return self.train(noisy_images, clean_images, epochs=epochs, verbose=verbose)

    def predict(
        self,
        image: ArrayLike,
        tile_size: int = 256,
        tile_overlap: int = 32,
    ) -> np.ndarray:
        """
        Denoise an image using the trained model.

        Parameters
        ----------
        image : array-like
            Noisy input image
        tile_size : int
            Size of tiles for processing large images
        tile_overlap : int
            Overlap between tiles

        Returns
        -------
        np.ndarray
            Denoised image
        """
        import torch

        if not self.trained:
            raise RuntimeError("Model not trained. Call train() first.")

        # Convert to numpy
        if hasattr(image, "compute"):
            img = image.compute()
        else:
            img = np.asarray(image)

        original_dtype = img.dtype
        original_shape = img.shape

        # Handle 2D images
        if img.ndim == 2:
            img = img[np.newaxis, :, :]

        # Normalize
        img_norm = (img.astype(np.float32) - self._mean) / self._std

        # Process with tiling
        c, h, w = img_norm.shape
        output = np.zeros_like(img_norm)
        weights = np.zeros((h, w), dtype=np.float32)

        self.model.eval()

        with torch.no_grad():
            for y in range(0, h, tile_size - tile_overlap):
                for x in range(0, w, tile_size - tile_overlap):
                    y_end = min(y + tile_size, h)
                    x_end = min(x + tile_size, w)
                    tile = img_norm[:, y:y_end, x:x_end]

                    # Pad if necessary
                    pad_h = tile_size - tile.shape[1]
                    pad_w = tile_size - tile.shape[2]
                    if pad_h > 0 or pad_w > 0:
                        tile = np.pad(tile, ((0, 0), (0, pad_h), (0, pad_w)), mode="reflect")

                    # Predict
                    tile_tensor = torch.from_numpy(tile[np.newaxis]).float().to(self.device)
                    pred = self.model(tile_tensor)[0].cpu().numpy()

                    # Remove padding
                    pred = pred[:, : y_end - y, : x_end - x]

                    # Blend weights
                    tile_h, tile_w = pred.shape[1], pred.shape[2]
                    blend_y = np.linspace(0, 1, min(tile_overlap, tile_h))
                    blend_x = np.linspace(0, 1, min(tile_overlap, tile_w))

                    weight = np.ones((tile_h, tile_w))
                    if y > 0:
                        weight[: len(blend_y), :] *= blend_y[:, np.newaxis]
                    if x > 0:
                        weight[:, : len(blend_x)] *= blend_x[np.newaxis, :]

                    output[:, y:y_end, x:x_end] += pred * weight
                    weights[y:y_end, x:x_end] += weight

        # Normalize
        weights = np.maximum(weights, 1e-8)
        output = output / weights[np.newaxis, :, :]

        # Denormalize
        output = output * self._std + self._mean

        # Restore shape
        if len(original_shape) == 2:
            output = output[0]

        # Clip and restore dtype
        if np.issubdtype(original_dtype, np.integer):
            info = np.iinfo(original_dtype)
            output = np.clip(output, info.min, info.max)

        return output.astype(original_dtype)

    def save(self, path: str | Path):
        """Save trained model."""
        import torch

        if not self.trained:
            raise RuntimeError("Model not trained.")

        torch.save(
            {
                "model_state": self.model.state_dict(),
                "config": self.config,
                "mean": self._mean,
                "std": self._std,
            },
            path,
        )

        logger.info(f"Model saved to {path}")

    def load(self, path: str | Path):
        """Load trained model."""
        import torch

        checkpoint = torch.load(path, map_location="cpu")
        self.config = checkpoint["config"]
        self._mean = checkpoint["mean"]
        self._std = checkpoint["std"]

        self._setup_device()
        self.model = self._build_model(n_channels=1)
        self.model.load_state_dict(checkpoint["model_state"])
        self.trained = True

        logger.info(f"Model loaded from {path}")


def denoise_care(
    image: ArrayLike,
    model: str | CAREDenoiser | None = None,
    **kwargs,
) -> np.ndarray:
    """
    Denoise an image using CARE.

    Parameters
    ----------
    image : array-like
        Noisy input image
    model : str or CAREDenoiser, optional
        Either a model name for pre-trained models, or a CAREDenoiser instance.
        If None, uses default pre-trained model.
    **kwargs
        Additional arguments passed to predict()

    Returns
    -------
    np.ndarray
        Denoised image
    """
    if model is None:
        model = "fluorescence_2d"

    if isinstance(model, str):
        denoiser = CAREDenoiser.from_pretrained(model)
    elif isinstance(model, CAREDenoiser):
        denoiser = model
    else:
        raise TypeError(f"model must be str or CAREDenoiser, got {type(model)}")

    return denoiser.predict(image, **kwargs)


def list_pretrained_models() -> dict[str, dict[str, Any]]:
    """
    List available pre-trained CARE models.

    Returns
    -------
    dict
        Dictionary mapping model names to their descriptions
    """
    return {
        name: {
            "description": info["description"],
            "trained_on": info["trained_on"],
            "available": info["url"] is not None,
        }
        for name, info in PRETRAINED_MODELS.items()
    }
