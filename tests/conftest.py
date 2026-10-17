"""
Pytest configuration and fixtures for KINTSUGI tests.
"""

import os
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


@pytest.fixture(autouse=True)
def _no_configured_mcp_project(monkeypatch):
    """Start every test without a configured MCP project.

    The MCP server's default project comes from ``kintsugi mcp start
    --project-dir`` (module state) or ``KINTSUGI_PROJECT_DIR``; neither may leak
    from the developer's shell or from one test into the next.
    """
    monkeypatch.delenv("KINTSUGI_PROJECT_DIR", raising=False)
    try:
        from kintsugi.mcp.session import set_project_root
    except ImportError:  # pragma: no cover - package not importable
        yield
        return
    set_project_root(None)
    yield
    set_project_root(None)


@pytest.fixture(scope="session", autouse=True)
def _require_mcp_when_demanded():
    """Turn MCP skips into a hard failure when CI demands the SDK.

    CI sets KINTSUGI_REQUIRE_MCP=1 after installing the ``mcp`` extra, so a
    broken SDK import can never degrade into silently skipped protocol tests.
    """
    if os.environ.get("KINTSUGI_REQUIRE_MCP"):
        from kintsugi.mcp._sdk import MCP_AVAILABLE, unavailable_message

        if not MCP_AVAILABLE:
            pytest.exit(f"KINTSUGI_REQUIRE_MCP is set but {unavailable_message()}", returncode=1)


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


@pytest.fixture
def sample_zstack():
    """Generate a sample 3D z-stack for testing."""
    np.random.seed(42)
    return np.random.randint(0, 65535, (10, 256, 256), dtype=np.uint16)


@pytest.fixture
def sample_tile_stack():
    """Generate a stack of tile images for illumination correction tests."""
    np.random.seed(42)
    return np.random.randint(1000, 10000, (25, 128, 128), dtype=np.float64)


@pytest.fixture
def sample_signal_blank_pair():
    """Generate a signal/blank image pair for autofluorescence tests."""
    np.random.seed(42)
    # Signal with true signal + autofluorescence
    signal = np.zeros((256, 256), dtype=np.uint16)
    signal[100:150, 100:150] = 25000  # True signal
    af_component = np.random.randint(1000, 5000, (256, 256), dtype=np.uint16)
    signal = (signal + af_component).astype(np.uint16)

    # Blank captures autofluorescence
    blank = np.random.randint(1000, 5000, (256, 256), dtype=np.uint16)

    return signal, blank
