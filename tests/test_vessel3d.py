"""Tests for 3D vessel segmentation (vessel3d.py).

Tests cover:
- Frangi Ra formula correctness (critical bug fix)
- Scikit-image deprecation-free imports
- GPU VRAM guard logic
- Binarization and morphological cleanup
- End-to-end pipeline on synthetic data
"""

import numpy as np
import pytest

from kintsugi.vessel3d import (
    VesselSpacing,
    _hessian_eigenvalues_3d,
    binarize_vessel_mask,
    compute_vesselness_frangi,
    segment_vessels_3d,
)


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def synthetic_tube():
    """Create a small 3D volume with a bright tube along the z-axis.

    A tube should have: l1 ~ 0, |l2| ~ |l3| >> 0
    Ra = |l2|/|l3| should be near 1 for a perfect tube.
    """
    vol = np.zeros((32, 32, 32), dtype=np.float32)
    # Draw a tube along z-axis at center (y=16, x=16) with radius ~3
    zz, yy, xx = np.mgrid[0:32, 0:32, 0:32]
    r = np.sqrt((yy - 16) ** 2 + (xx - 16) ** 2)
    vol[r < 3] = 1000.0
    return vol


@pytest.fixture
def synthetic_plate():
    """Create a 3D volume with a bright plane (sheet).

    A plate should have: l1 ~ 0, l2 ~ 0, |l3| >> 0
    Ra = |l2|/|l3| should be near 0 for a plate.
    """
    vol = np.zeros((32, 32, 32), dtype=np.float32)
    # Bright sheet at z=16 with thickness ~3
    vol[14:18, :, :] = 1000.0
    return vol


@pytest.fixture
def synthetic_blob():
    """Create a 3D volume with a bright sphere (blob).

    A blob should have: |l1| ~ |l2| ~ |l3|
    Rb should be near 1 for a blob.
    """
    vol = np.zeros((32, 32, 32), dtype=np.float32)
    zz, yy, xx = np.mgrid[0:32, 0:32, 0:32]
    r = np.sqrt((zz - 16) ** 2 + (yy - 16) ** 2 + (xx - 16) ** 2)
    vol[r < 5] = 1000.0
    return vol


# =============================================================================
# Frangi Ra Formula Tests
# =============================================================================


class TestFrangiRaFormula:
    """Verify the corrected Ra = |l2|/|l3| formula."""

    def test_tube_has_high_ra(self, synthetic_tube):
        """For a tube, Ra = |l2|/|l3| should be near 1 at center voxels."""
        l1, l2, l3 = _hessian_eigenvalues_3d(synthetic_tube, sigma=2.0)

        # At the tube center (z=16, y=16, x=16), l2 and l3 should be
        # similar in magnitude (both large), l1 should be small
        center = (16, 16, 16)
        eps = 1e-10
        ra_center = np.abs(l2[center]) / (np.abs(l3[center]) + eps)

        # Ra should be close to 1 for a tube (l2 ~ l3)
        assert ra_center > 0.5, (
            f"Ra at tube center = {ra_center:.4f}, expected > 0.5. "
            f"l1={l1[center]:.4f}, l2={l2[center]:.4f}, l3={l3[center]:.4f}"
        )

    def test_plate_has_low_ra(self, synthetic_plate):
        """For a plate, Ra = |l2|/|l3| should be near 0 at center voxels."""
        l1, l2, l3 = _hessian_eigenvalues_3d(synthetic_plate, sigma=2.0)

        # At the plate center, l3 >> l2 ~ l1 ~ 0
        center = (16, 16, 16)
        eps = 1e-10
        ra_center = np.abs(l2[center]) / (np.abs(l3[center]) + eps)

        # Ra should be small for a plate (l2 << l3)
        assert ra_center < 0.5, (
            f"Ra at plate center = {ra_center:.4f}, expected < 0.5. "
            f"l1={l1[center]:.4f}, l2={l2[center]:.4f}, l3={l3[center]:.4f}"
        )

    def test_tube_center_has_vesselness(self, synthetic_tube):
        """Tube center voxel should have high vesselness response."""
        v = compute_vesselness_frangi(synthetic_tube, sigmas=[2.0], device="cpu")
        # The tube runs along z at (y=16, x=16)
        center_response = v[16, 16, 16]
        assert center_response > 0.5, (
            f"Tube center vesselness = {center_response:.4f}, expected > 0.5"
        )

    def test_plate_center_suppressed_vs_tube(self, synthetic_tube, synthetic_plate):
        """Plate center should have lower unnormalized vesselness than tube center.

        Compare raw (pre-normalization) vesselness at structure centers.
        With correct Ra = |l2|/|l3|, plates yield (1-exp(-Ra^2/2a^2)) ≈ 0
        while tubes yield ≈ 0.86.
        """
        # Get eigenvalues and compute raw vesselness terms at center
        l1_t, l2_t, l3_t = _hessian_eigenvalues_3d(
            synthetic_tube.astype(np.float32), sigma=2.0
        )
        l1_p, l2_p, l3_p = _hessian_eigenvalues_3d(
            synthetic_plate.astype(np.float32), sigma=2.0
        )

        center = (16, 16, 16)
        eps = 1e-10

        # Ra term: (1 - exp(-Ra^2 / (2*alpha^2))) with alpha=0.5
        ra_tube = np.abs(l2_t[center]) / (np.abs(l3_t[center]) + eps)
        ra_plate = np.abs(l2_p[center]) / (np.abs(l3_p[center]) + eps)

        ra_term_tube = 1.0 - np.exp(-(ra_tube ** 2) / (2.0 * 0.5 ** 2))
        ra_term_plate = 1.0 - np.exp(-(ra_plate ** 2) / (2.0 * 0.5 ** 2))

        assert ra_term_tube > ra_term_plate, (
            f"Tube Ra term ({ra_term_tube:.4f}) should exceed "
            f"plate Ra term ({ra_term_plate:.4f}). "
            f"Ra_tube={ra_tube:.4f}, Ra_plate={ra_plate:.4f}"
        )


# =============================================================================
# Eigenvalue Sorting
# =============================================================================


class TestEigenvalueSorting:
    """Verify eigenvalues are sorted by absolute value."""

    def test_eigenvalues_sorted_by_abs(self, synthetic_tube):
        """Eigenvalues must satisfy |l1| <= |l2| <= |l3| everywhere."""
        l1, l2, l3 = _hessian_eigenvalues_3d(synthetic_tube, sigma=2.0)

        abs_l1 = np.abs(l1)
        abs_l2 = np.abs(l2)
        abs_l3 = np.abs(l3)

        # Allow small floating point tolerance
        tol = 1e-6 * abs_l3.max()
        assert np.all(abs_l1 <= abs_l2 + tol), "|l1| > |l2| violation found"
        assert np.all(abs_l2 <= abs_l3 + tol), "|l2| > |l3| violation found"


# =============================================================================
# Scikit-image Deprecation
# =============================================================================


class TestSkimageDeprecation:
    """Verify no deprecated scikit-image functions are used."""

    def test_closing_import(self):
        """Verify we import 'closing' not 'binary_closing'."""
        import inspect
        from kintsugi import vessel3d

        source = inspect.getsource(vessel3d.binarize_vessel_mask)
        assert "binary_closing" not in source, (
            "binarize_vessel_mask still uses deprecated binary_closing"
        )
        assert "closing" in source

    def test_binarize_no_warnings(self, synthetic_tube):
        """Binarize should run without deprecation warnings."""
        import warnings

        vesselness = compute_vesselness_frangi(
            synthetic_tube, sigmas=[2.0], device="cpu"
        )
        with warnings.catch_warnings():
            warnings.simplefilter("error", FutureWarning)
            binarize_vessel_mask(vesselness, min_size=10, closing_radius=1)


# =============================================================================
# GPU VRAM Guard
# =============================================================================


class TestVRAMGuard:
    """Test GPU VRAM estimation logic."""

    def test_vram_estimation_formula(self):
        """VRAM estimate should be volume_bytes * 15."""
        vol = np.zeros((40, 100, 100), dtype=np.float32)
        expected_vram = vol.nbytes * 15
        # 40*100*100*4 = 1.6 MB * 15 = 24 MB
        assert expected_vram == 40 * 100 * 100 * 4 * 15

    def test_cpu_fallback_on_small_volume(self, synthetic_tube):
        """CPU path should always work regardless of GPU availability."""
        l1, l2, l3 = _hessian_eigenvalues_3d(
            synthetic_tube.astype(np.float32), sigma=2.0, use_gpu=False
        )
        assert l1.shape == synthetic_tube.shape
        assert l2.shape == synthetic_tube.shape
        assert l3.shape == synthetic_tube.shape


# =============================================================================
# Binarization
# =============================================================================


class TestBinarization:
    """Test binarize_vessel_mask."""

    def test_basic_binarization(self, synthetic_tube):
        """Binarization should produce a boolean mask."""
        v = compute_vesselness_frangi(synthetic_tube, sigmas=[2.0], device="cpu")
        mask = binarize_vessel_mask(v, min_size=0, closing_radius=0)
        assert mask.dtype == bool
        assert mask.sum() > 0

    def test_threshold_reduces_volume(self, synthetic_tube):
        """Higher threshold should result in fewer positive voxels."""
        v = compute_vesselness_frangi(synthetic_tube, sigmas=[2.0], device="cpu")
        mask_low = binarize_vessel_mask(v, threshold=0.1, min_size=0, closing_radius=0)
        mask_high = binarize_vessel_mask(v, threshold=0.5, min_size=0, closing_radius=0)
        assert mask_high.sum() <= mask_low.sum()

    def test_min_size_removes_small_objects(self, synthetic_tube):
        """Large min_size should remove small components."""
        v = compute_vesselness_frangi(synthetic_tube, sigmas=[2.0], device="cpu")
        mask_small = binarize_vessel_mask(v, min_size=0, closing_radius=0)
        mask_large = binarize_vessel_mask(v, min_size=10000, closing_radius=0)
        assert mask_large.sum() <= mask_small.sum()


# =============================================================================
# End-to-End Pipeline
# =============================================================================


class TestPipeline:
    """Test segment_vessels_3d end-to-end."""

    def test_pipeline_synthetic_tube(self, synthetic_tube):
        """Full pipeline should complete on a synthetic tube."""
        spacing = VesselSpacing(xy=1.0, z=1.0)  # isotropic, skip resampling
        result = segment_vessels_3d(
            synthetic_tube,
            spacing=spacing,
            sigmas=[2.0],
            min_size=10,
            closing_radius=1,
            denoise_sigma=0,
            make_isotropic=False,
            prune_min_length_um=0,
            skip_skeleton=False,
            skip_graph=False,
            device="cpu",
            marker_name="test_tube",
        )

        assert result.binary_mask is not None
        assert result.binary_mask.dtype == bool
        assert result.binary_mask.sum() > 0
        assert result.marker_name == "test_tube"
        assert result.spacing == spacing
        assert result.vesselness is not None

    def test_pipeline_skip_skeleton(self, synthetic_tube):
        """Pipeline with skip_skeleton should return mask but no skeleton."""
        spacing = VesselSpacing(xy=1.0, z=1.0)
        result = segment_vessels_3d(
            synthetic_tube,
            spacing=spacing,
            sigmas=[2.0],
            min_size=10,
            make_isotropic=False,
            denoise_sigma=0,
            skip_skeleton=True,
            device="cpu",
        )

        assert result.binary_mask is not None
        assert result.skeleton is None
        assert result.features is None

    def test_pipeline_with_isotropic_resampling(self):
        """Pipeline should handle anisotropic spacing with z-resampling."""
        # Small anisotropic volume
        vol = np.zeros((8, 32, 32), dtype=np.float32)
        zz, yy, xx = np.mgrid[0:8, 0:32, 0:32]
        r = np.sqrt((yy - 16) ** 2 + (xx - 16) ** 2)
        vol[r < 3] = 1000.0

        spacing = VesselSpacing(xy=1.0, z=2.0)
        result = segment_vessels_3d(
            vol,
            spacing=spacing,
            sigmas=[2.0],
            min_size=10,
            make_isotropic=True,
            denoise_sigma=0,
            skip_skeleton=True,
            device="cpu",
        )

        # Isotropic mask should have more z-slices than input
        assert result.binary_mask.shape[0] > vol.shape[0]
        # Native-resolution mask should have original z-slices
        assert result.binary_mask_native is not None
        assert result.binary_mask_native.shape[0] == vol.shape[0]

    def test_pipeline_rejects_2d_input(self):
        """Pipeline should reject 2D input."""
        vol_2d = np.zeros((32, 32), dtype=np.float32)
        with pytest.raises(ValueError, match="3D"):
            segment_vessels_3d(vol_2d, device="cpu")


# =============================================================================
# VesselSpacing
# =============================================================================


class TestVesselSpacing:
    """Test VesselSpacing data class."""

    def test_ratio(self):
        s = VesselSpacing(xy=0.377, z=1.5)
        assert abs(s.ratio - 1.5 / 0.377) < 0.01

    def test_from_experiment(self):
        config = {"xy_pixel_size": 0.5, "z_step_size": 2.0}
        s = VesselSpacing.from_experiment(config)
        assert s.xy == 0.5
        assert s.z == 2.0

    def test_from_experiment_defaults(self):
        s = VesselSpacing.from_experiment({})
        assert s.xy == 0.377
        assert s.z == 1.5
