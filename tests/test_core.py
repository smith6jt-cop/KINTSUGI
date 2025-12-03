"""
Tests for core functionality without heavy dependencies.
"""

import numpy as np


class TestImageProcessing:
    """Test basic image processing functionality."""

    def test_numpy_array_creation(self, sample_image):
        """Test sample image fixture."""
        assert sample_image.shape == (256, 256)
        assert sample_image.dtype == np.uint8

    def test_multichannel_image(self, sample_multichannel_image):
        """Test multichannel image fixture."""
        assert sample_multichannel_image.shape == (256, 256, 3)
        assert sample_multichannel_image.dtype == np.uint8

    def test_image_stack(self, sample_stack):
        """Test image stack fixture."""
        assert sample_stack.shape == (10, 256, 256)
        assert sample_stack.dtype == np.uint8


class TestTiffIO:
    """Test TIFF I/O functionality."""

    def test_tiff_write_read(self, sample_tiff, sample_image):
        """Test writing and reading TIFF files."""
        import tifffile

        # Read the written file
        loaded = tifffile.imread(sample_tiff)

        assert loaded.shape == sample_image.shape
        assert loaded.dtype == sample_image.dtype
        np.testing.assert_array_equal(loaded, sample_image)

    def test_tiff_metadata(self, temp_dir, sample_image):
        """Test TIFF with metadata."""
        import tifffile

        tiff_path = temp_dir / "with_metadata.tif"

        # Write with metadata
        tifffile.imwrite(
            tiff_path,
            sample_image,
            metadata={"test_key": "test_value"},
        )

        assert tiff_path.exists()


class TestConfiguration:
    """Test configuration handling."""

    def test_config_structure(self, sample_config):
        """Test sample config fixture structure."""
        assert "src_dir" in sample_config
        assert "dst_dir" in sample_config
        assert "image_type" in sample_config
        assert "max_image_dim_px" in sample_config

    def test_config_defaults(self, sample_config):
        """Test config default values."""
        assert sample_config["image_type"] == "tif"
        assert sample_config["series"] == 0
        assert sample_config["align_to_reference"] is True

    def test_config_json_serialization(self, sample_config):
        """Test config can be serialized to JSON."""
        import json

        json_str = json.dumps(sample_config)
        loaded = json.loads(json_str)

        assert loaded == sample_config


class TestPathHandling:
    """Test path handling utilities."""

    def test_temp_dir_fixture(self, temp_dir):
        """Test temp_dir fixture."""
        assert temp_dir.exists()
        assert temp_dir.is_dir()

    def test_create_subdirectories(self, temp_dir):
        """Test creating subdirectories."""
        subdir = temp_dir / "level1" / "level2"
        subdir.mkdir(parents=True)

        assert subdir.exists()
        assert subdir.is_dir()


class TestNumpyOperations:
    """Test NumPy operations used in image processing."""

    def test_image_normalization(self, sample_image):
        """Test image normalization."""
        normalized = sample_image.astype(np.float32) / 255.0

        assert normalized.min() >= 0.0
        assert normalized.max() <= 1.0
        assert normalized.dtype == np.float32

    def test_image_statistics(self, sample_image):
        """Test computing image statistics."""
        mean = np.mean(sample_image)
        std = np.std(sample_image)
        min_val = np.min(sample_image)
        max_val = np.max(sample_image)

        assert 0 <= mean <= 255
        assert std >= 0
        assert 0 <= min_val <= max_val <= 255

    def test_image_resize_basic(self, sample_image):
        """Test basic image resizing with numpy."""
        # Simple downsampling by factor of 2
        downsampled = sample_image[::2, ::2]

        assert downsampled.shape == (128, 128)

    def test_image_channel_operations(self, sample_multichannel_image):
        """Test channel operations."""
        # Split channels
        r = sample_multichannel_image[:, :, 0]
        g = sample_multichannel_image[:, :, 1]
        b = sample_multichannel_image[:, :, 2]

        assert r.shape == (256, 256)
        assert g.shape == (256, 256)
        assert b.shape == (256, 256)

        # Stack channels
        stacked = np.stack([r, g, b], axis=-1)
        np.testing.assert_array_equal(stacked, sample_multichannel_image)


class TestSciPyOperations:
    """Test SciPy operations used in image processing."""

    def test_gaussian_filter(self, sample_image):
        """Test Gaussian filtering."""
        from scipy.ndimage import gaussian_filter

        smoothed = gaussian_filter(sample_image.astype(np.float32), sigma=2)

        assert smoothed.shape == sample_image.shape
        # Smoothed image should have less variance
        assert np.std(smoothed) < np.std(sample_image.astype(np.float32))

    def test_median_filter(self, sample_image):
        """Test median filtering."""
        from scipy.ndimage import median_filter

        filtered = median_filter(sample_image, size=3)

        assert filtered.shape == sample_image.shape
        assert filtered.dtype == sample_image.dtype

    def test_affine_transform(self, sample_image):
        """Test affine transformation."""
        from scipy.ndimage import affine_transform

        # Identity transform
        identity_matrix = np.eye(2)
        transformed = affine_transform(
            sample_image.astype(np.float32),
            identity_matrix,
        )

        assert transformed.shape == sample_image.shape
