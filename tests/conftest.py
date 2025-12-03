"""
Pytest configuration and fixtures for KINTSUGI tests.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

# Add project paths
PROJECT_ROOT = Path(__file__).parent.parent
NOTEBOOKS_PATH = PROJECT_ROOT / "notebooks"
SRC_PATH = PROJECT_ROOT / "src"

if str(NOTEBOOKS_PATH) not in sys.path:
    sys.path.insert(0, str(NOTEBOOKS_PATH))
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))


@pytest.fixture
def sample_image():
    """Generate a sample test image."""
    np.random.seed(42)
    return np.random.randint(0, 255, (256, 256), dtype=np.uint8)


@pytest.fixture
def sample_multichannel_image():
    """Generate a sample multi-channel test image."""
    np.random.seed(42)
    return np.random.randint(0, 255, (256, 256, 3), dtype=np.uint8)


@pytest.fixture
def sample_stack():
    """Generate a sample image stack."""
    np.random.seed(42)
    return np.random.randint(0, 255, (10, 256, 256), dtype=np.uint8)


@pytest.fixture
def temp_dir(tmp_path):
    """Create a temporary directory for test outputs."""
    return tmp_path


@pytest.fixture
def sample_tiff(temp_dir, sample_image):
    """Create a sample TIFF file for testing."""
    import tifffile

    tiff_path = temp_dir / "sample.tif"
    tifffile.imwrite(tiff_path, sample_image)
    return tiff_path


@pytest.fixture
def sample_config():
    """Generate a sample configuration dictionary."""
    return {
        "src_dir": "/path/to/source",
        "dst_dir": "/path/to/output",
        "reference_image": "cycle1.tif",
        "image_type": "tif",
        "series": 0,
        "max_image_dim_px": 2048,
        "max_processed_image_dim_px": 2048,
        "align_to_reference": True,
        "create_masks": True,
        "resolution_xyu": [0.325, 0.325, "um"],
        "channel_names": ["DAPI", "CD31", "CD3e"],
    }
