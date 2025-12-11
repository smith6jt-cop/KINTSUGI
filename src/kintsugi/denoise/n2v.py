"""
Noise2Void (N2V) Denoising

Pure Python implementation of self-supervised denoising using
blind-spot networks. Trains on noisy images without clean references.

Based on: Krull et al., "Noise2Void - Learning Denoising from Single Noisy Images"
IEEE/CVF Conference on Computer Vision and Pattern Recognition, 2019.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Union

import numpy as np

if TYPE_CHECKING:
    import dask.array

logger = logging.getLogger("kintsugi.denoise.n2v")

# Type alias
ArrayLike = Union[np.ndarray, "dask.array.Array"]


@dataclass
class N2VConfig:
    """Configuration for N2V training and inference."""

    # Network architecture
    n_filters: int = 32
    n_depth: int = 3
    kernel_size: int = 3
    batch_norm: bool = True

    # Training
    patch_size: int = 64
    batch_size: int = 32
    n_epochs: int = 100
    learning_rate: float = 1e-3
    n_patches_per_image: int = 128

    # N2V specific
    blind_spot_radius: int = 1
    masking_method: str = "uniform"  # "uniform", "median", "gaussian"
    mask_percentage: float = 0.2  # Percentage of pixels to mask per patch

    # Device
    device: str = "auto"  # "auto", "cuda", "cpu"


class N2VDenoiser:
    """
    Noise2Void denoiser with training and inference.

    N2V is a self-supervised denoising method that trains on noisy
    images alone, without requiring clean ground truth data. It uses
    a blind-spot network that predicts each pixel from its surroundings,
    preventing the network from learning the identity.

    Example
    -------
    >>> denoiser = N2VDenoiser()
    >>> denoiser.train([noisy_image1, noisy_image2], epochs=50)
    >>> denoised = denoiser.predict(noisy_image)
    """

    def __init__(self, config: N2VConfig | None = None):
        """
        Initialize N2V denoiser.

        Parameters
        ----------
        config : N2VConfig, optional
            Configuration for the denoiser
        """
        self.config = config or N2VConfig()
        self.model = None
        self.device = None
        self.trained = False
        self._mean = 0.0
        self._std = 1.0

    def _setup_device(self):
        """Set up PyTorch device."""
        try:
            import torch

            if self.config.device == "auto":
                self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
            else:
                self.device = torch.device(self.config.device)

            logger.info(f"N2V using device: {self.device}")

        except ImportError:
            raise ImportError(
                "PyTorch is required for N2V denoising. " "Install with: pip install torch"
            )

    def _build_model(self, n_channels: int = 1):
        """Build the U-Net model for N2V."""
        import torch
        import torch.nn as nn

        class ConvBlock(nn.Module):
            def __init__(self, in_ch, out_ch, kernel_size=3, batch_norm=True):
                super().__init__()
                padding = kernel_size // 2
                layers = [
                    nn.Conv2d(in_ch, out_ch, kernel_size, padding=padding),
                    nn.ReLU(inplace=True),
                ]
                if batch_norm:
                    layers.insert(1, nn.BatchNorm2d(out_ch))
                layers.extend(
                    [
                        nn.Conv2d(out_ch, out_ch, kernel_size, padding=padding),
                        nn.ReLU(inplace=True),
                    ]
                )
                if batch_norm:
                    layers.insert(-1, nn.BatchNorm2d(out_ch))
                self.block = nn.Sequential(*layers)

            def forward(self, x):
                return self.block(x)

        class N2VUNet(nn.Module):
            def __init__(self, n_channels, n_filters, n_depth, kernel_size, batch_norm):
                super().__init__()
                self.n_depth = n_depth

                # Encoder
                self.encoders = nn.ModuleList()
                self.pools = nn.ModuleList()
                in_ch = n_channels
                out_ch = n_filters
                for _ in range(n_depth):
                    self.encoders.append(ConvBlock(in_ch, out_ch, kernel_size, batch_norm))
                    self.pools.append(nn.MaxPool2d(2))
                    in_ch = out_ch
                    out_ch *= 2

                # Bottleneck
                self.bottleneck = ConvBlock(in_ch, out_ch, kernel_size, batch_norm)

                # Decoder
                self.upconvs = nn.ModuleList()
                self.decoders = nn.ModuleList()
                for _ in range(n_depth):
                    self.upconvs.append(nn.ConvTranspose2d(out_ch, out_ch // 2, 2, stride=2))
                    self.decoders.append(ConvBlock(out_ch, out_ch // 2, kernel_size, batch_norm))
                    out_ch //= 2

                # Output
                self.output = nn.Conv2d(out_ch, n_channels, 1)

            def forward(self, x):
                # Encoder path
                enc_features = []
                for encoder, pool in zip(self.encoders, self.pools, strict=False):
                    x = encoder(x)
                    enc_features.append(x)
                    x = pool(x)

                # Bottleneck
                x = self.bottleneck(x)

                # Decoder path
                for upconv, decoder, enc_feat in zip(
                    self.upconvs, self.decoders, reversed(enc_features), strict=False
                ):
                    x = upconv(x)
                    # Handle size mismatch
                    if x.shape != enc_feat.shape:
                        x = nn.functional.interpolate(x, size=enc_feat.shape[2:])
                    x = torch.cat([x, enc_feat], dim=1)
                    x = decoder(x)

                return self.output(x)

        model = N2VUNet(
            n_channels=n_channels,
            n_filters=self.config.n_filters,
            n_depth=self.config.n_depth,
            kernel_size=self.config.kernel_size,
            batch_norm=self.config.batch_norm,
        )

        return model.to(self.device)

    def _extract_patches(
        self,
        images: list[np.ndarray],
    ) -> np.ndarray:
        """Extract random patches from images for training."""
        patches = []
        patch_size = self.config.patch_size
        n_patches = self.config.n_patches_per_image

        for img in images:
            if img.ndim == 2:
                img = img[np.newaxis, :, :]  # Add channel dim

            c, h, w = img.shape

            for _ in range(n_patches):
                if h > patch_size and w > patch_size:
                    y = np.random.randint(0, h - patch_size)
                    x = np.random.randint(0, w - patch_size)
                    patch = img[:, y : y + patch_size, x : x + patch_size]
                else:
                    # Pad if image smaller than patch
                    pad_h = max(0, patch_size - h)
                    pad_w = max(0, patch_size - w)
                    padded = np.pad(img, ((0, 0), (0, pad_h), (0, pad_w)), mode="reflect")
                    patch = padded[:, :patch_size, :patch_size]

                patches.append(patch)

        return np.stack(patches, axis=0)

    def _create_n2v_mask(
        self,
        patches: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Create N2V blind-spot masks.

        Returns input patches, target values, and mask locations.
        """
        batch_size = patches.shape[0]
        n_channels = patches.shape[1]
        h, w = patches.shape[2], patches.shape[3]
        n_masked = int(h * w * self.config.mask_percentage)

        # Output arrays
        masked_patches = patches.copy()
        target_values = np.zeros((batch_size, n_channels, n_masked), dtype=patches.dtype)
        mask_coords = np.zeros((batch_size, 2, n_masked), dtype=np.int64)

        radius = self.config.blind_spot_radius

        for b in range(batch_size):
            # Random positions to mask
            y_coords = np.random.randint(radius, h - radius, size=n_masked)
            x_coords = np.random.randint(radius, w - radius, size=n_masked)

            for i, (y, x) in enumerate(zip(y_coords, x_coords, strict=False)):
                # Store original value as target
                target_values[b, :, i] = patches[b, :, y, x]
                mask_coords[b, 0, i] = y
                mask_coords[b, 1, i] = x

                # Replace with value from blind spot neighborhood
                if self.config.masking_method == "uniform":
                    # Random neighbor
                    dy, dx = np.random.randint(-radius, radius + 1, size=2)
                    if dy == 0 and dx == 0:
                        dx = 1
                    replacement = patches[b, :, y + dy, x + dx]

                elif self.config.masking_method == "median":
                    # Median of neighborhood
                    neighborhood = patches[
                        b, :, y - radius : y + radius + 1, x - radius : x + radius + 1
                    ]
                    replacement = np.median(neighborhood, axis=(1, 2))

                else:  # gaussian
                    # Gaussian-weighted average of neighborhood
                    neighborhood = patches[
                        b, :, y - radius : y + radius + 1, x - radius : x + radius + 1
                    ]
                    size = 2 * radius + 1
                    y_grid, x_grid = np.mgrid[:size, :size]
                    center = radius
                    weights = np.exp(
                        -((y_grid - center) ** 2 + (x_grid - center) ** 2) / (2 * radius**2)
                    )
                    weights[center, center] = 0
                    weights /= weights.sum()
                    replacement = np.sum(neighborhood * weights, axis=(1, 2))

                masked_patches[b, :, y, x] = replacement

        return masked_patches, target_values, mask_coords

    def train(
        self,
        images: list[np.ndarray],
        epochs: int | None = None,
        verbose: bool = True,
    ) -> dict[str, Any]:
        """
        Train the N2V model on noisy images.

        Parameters
        ----------
        images : list of np.ndarray
            List of noisy training images (2D or 3D with channel first)
        epochs : int, optional
            Number of training epochs (overrides config)
        verbose : bool
            Print training progress

        Returns
        -------
        dict
            Training history with loss values
        """
        import torch

        self._setup_device()

        # Normalize images
        all_pixels = np.concatenate([img.ravel() for img in images])
        self._mean = float(np.mean(all_pixels))
        self._std = float(np.std(all_pixels))
        if self._std < 1e-8:
            self._std = 1.0

        normalized_images = [(img - self._mean) / self._std for img in images]

        # Determine number of channels
        if images[0].ndim == 2:
            n_channels = 1
        else:
            n_channels = images[0].shape[0]

        # Build model
        self.model = self._build_model(n_channels)

        # Extract patches
        if verbose:
            logger.info("Extracting training patches...")
        patches = self._extract_patches(normalized_images)

        if verbose:
            logger.info(f"Extracted {len(patches)} patches of size {self.config.patch_size}")

        # Training setup
        n_epochs = epochs or self.config.n_epochs
        optimizer = torch.optim.Adam(self.model.parameters(), lr=self.config.learning_rate)
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode="min", factor=0.5, patience=10
        )

        history = {"loss": []}

        for epoch in range(n_epochs):
            # Shuffle patches
            np.random.shuffle(patches)

            epoch_loss = 0.0
            n_batches = 0

            for i in range(0, len(patches), self.config.batch_size):
                batch = patches[i : i + self.config.batch_size]

                # Create N2V masks
                masked, targets, coords = self._create_n2v_mask(batch)

                # Convert to tensors
                masked_tensor = torch.from_numpy(masked).float().to(self.device)
                targets_tensor = torch.from_numpy(targets).float().to(self.device)
                coords_tensor = torch.from_numpy(coords).long().to(self.device)

                # Forward pass
                optimizer.zero_grad()
                output = self.model(masked_tensor)

                # Compute loss only at masked positions
                loss = 0.0
                for b in range(len(batch)):
                    for c in range(n_channels):
                        for j in range(coords_tensor.shape[2]):
                            y, x = coords_tensor[b, 0, j], coords_tensor[b, 1, j]
                            pred = output[b, c, y, x]
                            target = targets_tensor[b, c, j]
                            loss += (pred - target) ** 2

                loss = loss / (len(batch) * n_channels * coords_tensor.shape[2])

                # Backward pass
                loss.backward()
                optimizer.step()

                epoch_loss += loss.item()
                n_batches += 1

            avg_loss = epoch_loss / n_batches
            history["loss"].append(avg_loss)
            scheduler.step(avg_loss)

            if verbose and (epoch + 1) % 10 == 0:
                lr = optimizer.param_groups[0]["lr"]
                logger.info(f"Epoch {epoch + 1}/{n_epochs}, Loss: {avg_loss:.6f}, LR: {lr:.2e}")

        self.trained = True

        if verbose:
            logger.info("Training complete!")

        return history

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

        # Process with tiling for large images
        c, h, w = img_norm.shape
        output = np.zeros_like(img_norm)
        weights = np.zeros((h, w), dtype=np.float32)

        self.model.eval()

        with torch.no_grad():
            for y in range(0, h, tile_size - tile_overlap):
                for x in range(0, w, tile_size - tile_overlap):
                    # Extract tile
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

                    # Blend with existing output using linear ramp
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

        # Normalize by weights
        weights = np.maximum(weights, 1e-8)
        output = output / weights[np.newaxis, :, :]

        # Denormalize
        output = output * self._std + self._mean

        # Restore original shape
        if len(original_shape) == 2:
            output = output[0]

        # Clip to valid range and restore dtype
        if np.issubdtype(original_dtype, np.integer):
            info = np.iinfo(original_dtype)
            output = np.clip(output, info.min, info.max)

        return output.astype(original_dtype)

    def save(self, path: str | Path):
        """Save trained model to file."""
        import torch

        if not self.trained:
            raise RuntimeError("Model not trained.")

        path = Path(path)
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
        """Load trained model from file."""
        import torch

        path = Path(path)
        checkpoint = torch.load(path, map_location="cpu")

        self.config = checkpoint["config"]
        self._mean = checkpoint["mean"]
        self._std = checkpoint["std"]

        self._setup_device()
        self.model = self._build_model(n_channels=1)  # Assume single channel
        self.model.load_state_dict(checkpoint["model_state"])
        self.trained = True

        logger.info(f"Model loaded from {path}")


def denoise_n2v(
    image: ArrayLike,
    patch_size: int = 64,
    n_epochs: int = 50,
    device: str = "auto",
    model_path: str | Path | None = None,
    save_model: bool = False,
    **kwargs,
) -> np.ndarray:
    """
    Denoise an image using Noise2Void.

    This is a convenience function that handles training and inference.
    For repeated use on similar images, consider using N2VDenoiser directly.

    Parameters
    ----------
    image : array-like
        Noisy input image
    patch_size : int
        Patch size for training
    n_epochs : int
        Number of training epochs
    device : str
        Device to use ("auto", "cuda", "cpu")
    model_path : str or Path, optional
        Path to pre-trained model (skip training if provided)
    save_model : bool
        Save trained model to model_path after training
    **kwargs
        Additional arguments passed to N2VConfig

    Returns
    -------
    np.ndarray
        Denoised image
    """
    config = N2VConfig(
        patch_size=patch_size,
        n_epochs=n_epochs,
        device=device,
        **kwargs,
    )

    denoiser = N2VDenoiser(config)

    if model_path and Path(model_path).exists():
        denoiser.load(model_path)
    else:
        # Convert to numpy if needed
        if hasattr(image, "compute"):
            img = image.compute()
        else:
            img = np.asarray(image)

        denoiser.train([img], verbose=True)

        if save_model and model_path:
            denoiser.save(model_path)

    return denoiser.predict(image)


def train_n2v(
    images: list[np.ndarray],
    output_path: str | Path,
    patch_size: int = 64,
    n_epochs: int = 100,
    device: str = "auto",
    **kwargs,
) -> N2VDenoiser:
    """
    Train an N2V model on multiple images and save it.

    Parameters
    ----------
    images : list of np.ndarray
        List of noisy training images
    output_path : str or Path
        Path to save the trained model
    patch_size : int
        Patch size for training
    n_epochs : int
        Number of training epochs
    device : str
        Device to use
    **kwargs
        Additional arguments for N2VConfig

    Returns
    -------
    N2VDenoiser
        Trained denoiser
    """
    config = N2VConfig(
        patch_size=patch_size,
        n_epochs=n_epochs,
        device=device,
        **kwargs,
    )

    denoiser = N2VDenoiser(config)
    denoiser.train(images, verbose=True)
    denoiser.save(output_path)

    return denoiser
