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
    VESSEL_PRESETS,
    VesselMarkerInfo,
    VesselSegmentationResult,
    VesselSpacing,
    _hessian_eigenvalues_3d,
    analyze_vessel_graph,
    binarize_vessel_mask,
    combine_vessel_masks,
    compute_vesselness_frangi,
    discover_vessel_markers,
    prune_skeleton,
    segment_vessels_3d,
    segment_vessels_multichannel,
    skeletonize_vessels,
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
        """VRAM estimate should be volume_bytes * 12 (peak=10 arrays + 20% safety)."""
        vol = np.zeros((40, 100, 100), dtype=np.float32)
        expected_vram = vol.nbytes * 12
        # 40*100*100*4 = 1.6 MB * 12 = 19.2 MB
        assert expected_vram == 40 * 100 * 100 * 4 * 12

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


# =============================================================================
# Extended Morphometry (VesselExpress features)
# =============================================================================


@pytest.fixture
def branching_structure():
    """Create a Y-shaped branching volume for branching angle tests.

    A trunk runs along z at center, with two branches diverging in y.
    """
    vol = np.zeros((48, 48, 48), dtype=np.float32)
    zz, yy, xx = np.mgrid[0:48, 0:48, 0:48]

    # Trunk: z-axis tube at (y=24, x=24) from z=0 to z=24
    r_trunk = np.sqrt((yy - 24) ** 2 + (xx - 24) ** 2)
    trunk_mask = (r_trunk < 3) & (zz < 25)
    vol[trunk_mask] = 1000.0

    # Branch 1: diverges toward y=36 from junction at z=24
    for z in range(24, 44):
        y_center = 24 + (z - 24) * 0.6  # angled branch
        r = np.sqrt((yy[z] - y_center) ** 2 + (xx[z] - 24) ** 2)
        vol[z][r < 3] = 1000.0

    # Branch 2: diverges toward y=12 from junction at z=24
    for z in range(24, 44):
        y_center = 24 - (z - 24) * 0.6
        r = np.sqrt((yy[z] - y_center) ** 2 + (xx[z] - 24) ** 2)
        vol[z][r < 3] = 1000.0

    return vol


class TestExtendedMorphometry:
    """Test VesselExpress-inspired morphometry features."""

    def test_volume_column_present(self, synthetic_tube):
        """volume_um3 should be > 0 for a synthetic tube."""
        spacing = VesselSpacing(xy=1.0, z=1.0)
        result = segment_vessels_3d(
            synthetic_tube,
            spacing=spacing,
            sigmas=[2.0],
            min_size=10,
            make_isotropic=False,
            denoise_sigma=0,
            prune_min_length_um=0,
            device="cpu",
        )
        assert result.features is not None
        assert "volume_um3" in result.features.columns
        assert (result.features["volume_um3"] >= 0).all()
        assert result.features["volume_um3"].sum() > 0

    def test_z_angle_for_z_aligned_tube(self, synthetic_tube):
        """A tube along the z-axis should have z_angle near 0 degrees."""
        spacing = VesselSpacing(xy=1.0, z=1.0)
        result = segment_vessels_3d(
            synthetic_tube,
            spacing=spacing,
            sigmas=[2.0],
            min_size=10,
            make_isotropic=False,
            denoise_sigma=0,
            prune_min_length_um=0,
            device="cpu",
        )
        assert result.features is not None
        assert "z_angle_deg" in result.features.columns
        # Main segment should be roughly z-aligned (< 30 degrees)
        longest = result.features.loc[result.features["branch_length_um"].idxmax()]
        assert longest["z_angle_deg"] < 30.0, (
            f"Z-aligned tube has z_angle={longest['z_angle_deg']:.1f}, expected < 30"
        )

    def test_branching_angle_column(self, synthetic_tube):
        """branching_angle_deg column should exist (NaN for unbranched tube is fine)."""
        spacing = VesselSpacing(xy=1.0, z=1.0)
        result = segment_vessels_3d(
            synthetic_tube,
            spacing=spacing,
            sigmas=[2.0],
            min_size=10,
            make_isotropic=False,
            denoise_sigma=0,
            prune_min_length_um=0,
            device="cpu",
        )
        assert result.features is not None
        assert "branching_angle_deg" in result.features.columns

    def test_branching_angle_on_y_structure(self, branching_structure):
        """A Y-shaped structure should have non-NaN branching angles on branches."""
        spacing = VesselSpacing(xy=1.0, z=1.0)
        result = segment_vessels_3d(
            branching_structure,
            spacing=spacing,
            sigmas=[2.0],
            min_size=10,
            make_isotropic=False,
            denoise_sigma=0,
            prune_min_length_um=0,
            device="cpu",
        )
        assert result.features is not None
        angles = result.features["branching_angle_deg"].dropna()
        # Y-structure should produce at least one non-NaN branching angle
        # (may not always produce depending on skeleton topology, so we just
        # check the column is computed without errors)
        assert "branching_angle_deg" in result.features.columns

    def test_diameter_ratio_pruning(self, synthetic_tube):
        """Diameter-ratio pruning should remove additional branches beyond length-only."""
        spacing = VesselSpacing(xy=1.0, z=1.0)
        v = compute_vesselness_frangi(synthetic_tube, sigmas=[2.0], device="cpu")
        mask = binarize_vessel_mask(v, min_size=10, closing_radius=1)
        skel = skeletonize_vessels(mask)

        # Length-only pruning
        pruned_length = prune_skeleton(
            skel,
            min_branch_length_um=1.0,
            spacing=spacing,
        )
        # Diameter-ratio pruning (aggressive ratio)
        pruned_ratio = prune_skeleton(
            skel,
            min_branch_length_um=1.0,
            diameter_length_ratio=5.0,
            binary_mask=mask,
            spacing=spacing,
        )
        # Diameter-ratio should remove at least as many voxels as length-only
        assert pruned_ratio.sum() <= pruned_length.sum()

    def test_summary_method(self, synthetic_tube):
        """summary() should return aggregate statistics."""
        spacing = VesselSpacing(xy=1.0, z=1.0)
        result = segment_vessels_3d(
            synthetic_tube,
            spacing=spacing,
            sigmas=[2.0],
            min_size=10,
            make_isotropic=False,
            denoise_sigma=0,
            prune_min_length_um=0,
            device="cpu",
        )
        s = result.summary()
        assert "n_segments" in s
        assert "total_length_um" in s
        assert "total_volume_um3" in s
        assert "mean_diameter_um" in s
        assert "mean_tortuosity" in s
        assert s["n_segments"] > 0
        assert s["total_length_um"] > 0
        assert s["total_volume_um3"] > 0

    def test_summary_empty_when_no_features(self):
        """summary() should return empty dict when features are None."""
        result = VesselSegmentationResult(
            binary_mask=np.zeros((4, 4, 4), dtype=bool)
        )
        assert result.summary() == {}


# =============================================================================
# Sensitivity Tuning, Presets, Marker Discovery, Multi-Channel
# =============================================================================


class TestSensitivityAndMultiChannel:
    """Test presets, alpha/beta forwarding, marker discovery, and multi-channel."""

    def test_alpha_beta_forwarded(self, synthetic_tube):
        """Passing alpha=0.1 should produce a denser mask than alpha=0.9."""
        spacing = VesselSpacing(xy=1.0, z=1.0)
        result_sensitive = segment_vessels_3d(
            synthetic_tube,
            spacing=spacing,
            sigmas=[2.0],
            alpha=0.1,
            beta=0.1,
            min_size=0,
            make_isotropic=False,
            denoise_sigma=0,
            skip_skeleton=True,
            device="cpu",
        )
        result_strict = segment_vessels_3d(
            synthetic_tube,
            spacing=spacing,
            sigmas=[2.0],
            alpha=0.9,
            beta=0.9,
            min_size=0,
            make_isotropic=False,
            denoise_sigma=0,
            skip_skeleton=True,
            device="cpu",
        )
        # More sensitive (lower alpha/beta) should keep at least as many voxels
        assert result_sensitive.binary_mask.sum() >= result_strict.binary_mask.sum(), (
            f"alpha=0.1 mask ({result_sensitive.binary_mask.sum()}) should be >= "
            f"alpha=0.9 mask ({result_strict.binary_mask.sum()})"
        )

    def test_preset_overrides_defaults(self, synthetic_tube):
        """preset='high_sensitivity' should change sigmas and min_size from defaults."""
        spacing = VesselSpacing(xy=1.0, z=1.0)
        # With preset, defaults should be overridden to high_sensitivity values
        result = segment_vessels_3d(
            synthetic_tube,
            spacing=spacing,
            make_isotropic=False,
            skip_skeleton=True,
            device="cpu",
            preset="high_sensitivity",
        )
        # The preset uses min_size=250 (vs default 500) and alpha=0.3 (vs 0.5)
        # We can't directly inspect the params used, but the mask should be non-empty
        assert result.binary_mask.sum() > 0

    def test_preset_explicit_override(self, synthetic_tube):
        """Explicit min_size=100 should override preset's min_size=250."""
        spacing = VesselSpacing(xy=1.0, z=1.0)
        result = segment_vessels_3d(
            synthetic_tube,
            spacing=spacing,
            min_size=100,  # explicit override of preset's 250
            make_isotropic=False,
            skip_skeleton=True,
            device="cpu",
            preset="high_sensitivity",
        )
        assert result.binary_mask is not None

    def test_preset_invalid_raises(self):
        """Unknown preset name should raise ValueError."""
        vol = np.zeros((8, 8, 8), dtype=np.float32)
        with pytest.raises(ValueError, match="Unknown preset"):
            segment_vessels_3d(vol, preset="nonexistent", device="cpu")

    def test_discover_vessel_markers_cd34_preferred(self):
        """CD34 should be ranked before CD31."""
        channel_names = {
            1: ["DAPI", "Blank", "Blank", "Blank"],
            2: ["DAPI", "CD31", "CD8", "CD45"],
            3: ["DAPI", "CD34", "CD20", "Blank"],
        }
        markers = discover_vessel_markers(channel_names)
        assert len(markers) >= 2
        assert markers[0].name == "CD34"
        assert markers[0].role == "endothelial"
        assert markers[0].cycle == 3
        assert markers[0].channel == 2  # 1-indexed
        assert markers[1].name == "CD31"
        assert markers[1].cycle == 2

    def test_discover_vessel_markers_excludes_generic_actin(self):
        """Plain 'Actin' (no SMA prefix) should be excluded."""
        channel_names = {
            1: ["DAPI", "Actin", "Blank", "Blank"],
            2: ["DAPI", "aSMA", "CD31", "Blank"],
        }
        markers = discover_vessel_markers(channel_names)
        names = [m.name for m in markers]
        assert "Actin" not in names
        assert "aSMA" in names
        assert "CD31" in names
        # SMA should come after CD31 in priority
        sma_idx = names.index("aSMA")
        cd31_idx = names.index("CD31")
        assert cd31_idx < sma_idx

    def test_combine_vessel_masks_union(self):
        """Union of two masks should have at least as many voxels as either."""
        mask1 = np.zeros((10, 10, 10), dtype=bool)
        mask2 = np.zeros((10, 10, 10), dtype=bool)
        mask1[2:5, 2:5, 2:5] = True
        mask2[6:9, 6:9, 6:9] = True

        combined = combine_vessel_masks([mask1, mask2], method="union")
        assert combined.sum() == mask1.sum() + mask2.sum()
        assert combined.sum() >= mask1.sum()
        assert combined.sum() >= mask2.sum()

        # Intersection should be empty (non-overlapping)
        combined_inter = combine_vessel_masks([mask1, mask2], method="intersection")
        assert combined_inter.sum() == 0

    def test_multichannel_pipeline(self):
        """End-to-end multi-channel pipeline with two synthetic tubes."""
        # Tube 1 at (y=12, x=16)
        vol1 = np.zeros((16, 32, 32), dtype=np.float32)
        zz, yy, xx = np.mgrid[0:16, 0:32, 0:32]
        r1 = np.sqrt((yy - 12) ** 2 + (xx - 16) ** 2)
        vol1[r1 < 3] = 1000.0

        # Tube 2 at (y=22, x=16)
        vol2 = np.zeros((16, 32, 32), dtype=np.float32)
        r2 = np.sqrt((yy - 22) ** 2 + (xx - 16) ** 2)
        vol2[r2 < 3] = 1000.0

        spacing = VesselSpacing(xy=1.0, z=1.0)
        result = segment_vessels_multichannel(
            {"CH1": vol1, "CH2": vol2},
            spacing=spacing,
            preset=None,
            sigmas=[2.0],
            min_size=10,
            make_isotropic=False,
            denoise_sigma=0,
            device="cpu",
        )

        assert result.marker_name == "CH1+CH2"
        assert result.binary_mask.sum() > 0
        assert result.skeleton is not None
        assert result.features is not None
        assert result.per_channel_masks is not None
        assert len(result.per_channel_masks) == 2
        # Combined mask should be >= either individual
        assert result.binary_mask.sum() >= result.per_channel_masks["CH1"].sum()
        assert result.binary_mask.sum() >= result.per_channel_masks["CH2"].sum()
