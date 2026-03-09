"""Tests for channel clustering and parameter propagation."""

import numpy as np
import pytest

from kintsugi.signal.clustering import (
    cluster_channels,
    generate_cluster_qc,
    get_cluster_representatives,
    plot_cluster_summary,
    propagate_cluster_parameters,
)
from kintsugi.signal.features import extract_channel_features


def _make_channel_image(intensity, noise_level=100, size=64, seed=None):
    """Create a synthetic channel image with controllable intensity."""
    rng = np.random.RandomState(seed)
    img = np.full((size, size), intensity, dtype=np.float64)
    # Add some structure
    img[size // 4 : 3 * size // 4, size // 4 : 3 * size // 4] += intensity * 0.5
    img += rng.normal(0, noise_level, img.shape)
    return np.clip(img, 0, 65535).astype(np.uint16)


def _make_test_features(n_channels=10):
    """Create test channel features for clustering."""
    channels = {}
    features = {}
    for i in range(n_channels):
        # Vary intensity and noise to create natural clusters
        if i < 4:  # bright group
            intensity = 5000 + i * 200
            name = f"bright_ch{i}"
        elif i < 7:  # dim group
            intensity = 500 + i * 50
            name = f"dim_ch{i}"
        else:  # medium group
            intensity = 2000 + i * 100
            name = f"med_ch{i}"

        img = _make_channel_image(intensity, seed=i)
        channels[name] = img
        features[name] = extract_channel_features(img, channel_name=name)

    return channels, features


class TestClusterChannels:
    def test_auto_n_clusters(self):
        _, features = _make_test_features(10)
        assignments = cluster_channels(features, wavelength_aware=False)
        assert len(assignments) == 10
        assert all(isinstance(v, int) for v in assignments.values())
        # Should find at least 2 clusters with varied features
        assert len(set(assignments.values())) >= 2

    def test_explicit_n_clusters(self):
        _, features = _make_test_features(10)
        assignments = cluster_channels(features, n_clusters=3, wavelength_aware=False)
        assert len(set(assignments.values())) == 3

    def test_fewer_than_3_channels(self):
        """Fewer than 3 channels should return single cluster."""
        channels = {}
        features = {}
        for i in range(2):
            img = _make_channel_image(3000, seed=i)
            name = f"ch{i}"
            channels[name] = img
            features[name] = extract_channel_features(img)
        assignments = cluster_channels(features)
        assert all(v == 0 for v in assignments.values())

    def test_kmeans_method(self):
        _, features = _make_test_features(10)
        assignments = cluster_channels(
            features, n_clusters=3, method="kmeans", wavelength_aware=False
        )
        assert len(assignments) == 10
        assert len(set(assignments.values())) == 3

    def test_wavelength_aware(self):
        """Channels with different wavelength groups should not cluster together."""
        features = {}
        # Create channels with specific wavelength groups
        for i, name in enumerate(["DAPI-01", "DAPI-02", "DAPI-03", "AF488-CD3", "AF488-CD20", "AF488-CD45"]):
            img = _make_channel_image(3000, seed=i)
            features[name] = extract_channel_features(img, channel_name=name)

        assignments = cluster_channels(features, wavelength_aware=True)

        # DAPI channels (405 group) should share cluster(s)
        dapi_clusters = {assignments[n] for n in ["DAPI-01", "DAPI-02", "DAPI-03"]}
        af488_clusters = {assignments[n] for n in ["AF488-CD3", "AF488-CD20", "AF488-CD45"]}
        # No overlap between wavelength groups
        assert dapi_clusters.isdisjoint(af488_clusters)

    def test_all_channels_assigned(self):
        _, features = _make_test_features(8)
        assignments = cluster_channels(features, wavelength_aware=False)
        assert set(assignments.keys()) == set(features.keys())


class TestGetClusterRepresentatives:
    def test_one_representative_per_cluster(self):
        _, features = _make_test_features(10)
        assignments = cluster_channels(features, n_clusters=3, wavelength_aware=False)
        reps = get_cluster_representatives(features, assignments)

        assert len(reps) == 3
        for cid, (name, count) in reps.items():
            assert name in features
            assert count > 0

    def test_representative_is_member(self):
        _, features = _make_test_features(10)
        assignments = cluster_channels(features, n_clusters=2, wavelength_aware=False)
        reps = get_cluster_representatives(features, assignments)

        for cid, (rep_name, _) in reps.items():
            assert assignments[rep_name] == cid

    def test_member_counts_sum_to_total(self):
        _, features = _make_test_features(10)
        assignments = cluster_channels(features, n_clusters=3, wavelength_aware=False)
        reps = get_cluster_representatives(features, assignments)

        total = sum(count for _, count in reps.values())
        assert total == 10


class TestPropagateClusterParameters:
    def test_output_files_created(self, tmp_path):
        rng = np.random.RandomState(42)
        marker_dict = {
            "CD3": rng.randint(100, 5000, (64, 64), dtype=np.uint16),
            "CD20": rng.randint(100, 5000, (64, 64), dtype=np.uint16),
            "Blank1": rng.randint(100, 800, (64, 64), dtype=np.uint16),
        }
        params = {
            "blank_params": {
                "blank_clip_factor": 200,
                "blank_scale_factor": 1.0,
            }
        }
        blank_map = {"CD3": "Blank1", "CD20": "Blank1"}

        results = propagate_cluster_parameters(
            marker_dict, ["CD3", "CD20"], params, blank_map, tmp_path
        )

        assert len(results) == 2
        for ch, res in results.items():
            assert res["status"] == "success"
            assert res["output_path"].exists()
            assert 0.0 <= res["quality_score"] <= 1.0

    def test_quality_scores_returned(self, tmp_path):
        rng = np.random.RandomState(42)
        marker_dict = {
            "CD3": rng.randint(100, 5000, (64, 64), dtype=np.uint16),
            "Blank1": rng.randint(100, 800, (64, 64), dtype=np.uint16),
        }
        params = {"blank_params": {"blank_clip_factor": 0, "blank_scale_factor": 1.0}}
        blank_map = {"CD3": "Blank1"}

        results = propagate_cluster_parameters(
            marker_dict, ["CD3"], params, blank_map, tmp_path
        )
        assert "quality_score" in results["CD3"]
        assert "quality_details" in results["CD3"]

    def test_error_handling(self, tmp_path):
        """Missing channel should not crash, should report error."""
        marker_dict = {}
        params = {"blank_params": {"blank_clip_factor": 0}}
        results = propagate_cluster_parameters(
            marker_dict, ["missing_ch"], params, {}, tmp_path
        )
        assert results["missing_ch"]["status"] == "error"

    def test_custom_process_fn(self, tmp_path):
        rng = np.random.RandomState(42)
        marker_dict = {
            "CD3": rng.randint(100, 5000, (64, 64), dtype=np.uint16),
        }

        def custom_fn(sig, blank, params):
            return sig  # identity

        results = propagate_cluster_parameters(
            marker_dict, ["CD3"], {}, {}, tmp_path, process_fn=custom_fn
        )
        assert results["CD3"]["status"] == "success"

    def test_no_blank_available(self, tmp_path):
        """Channel without blank should still process (copy signal)."""
        rng = np.random.RandomState(42)
        marker_dict = {
            "CD3": rng.randint(100, 5000, (64, 64), dtype=np.uint16),
        }
        params = {"blank_params": {"blank_clip_factor": 0, "blank_scale_factor": 1.0}}
        results = propagate_cluster_parameters(
            marker_dict, ["CD3"], params, {}, tmp_path
        )
        assert results["CD3"]["status"] == "success"

    def test_with_clean_params(self, tmp_path):
        rng = np.random.RandomState(42)
        marker_dict = {
            "CD3": rng.randint(100, 5000, (64, 64), dtype=np.uint16),
            "Blank1": rng.randint(100, 800, (64, 64), dtype=np.uint16),
        }
        params = {
            "blank_params": {"blank_clip_factor": 200, "blank_scale_factor": 1.0},
            "clean_params": {"background_threshold": 100, "smooth": False, "remove_small": True, "small_size": 30},
        }
        blank_map = {"CD3": "Blank1"}
        results = propagate_cluster_parameters(
            marker_dict, ["CD3"], params, blank_map, tmp_path
        )
        assert results["CD3"]["status"] == "success"


class TestPlotClusterSummary:
    def test_returns_figure(self):
        _, features = _make_test_features(6)
        assignments = cluster_channels(features, n_clusters=2, wavelength_aware=False)
        fig = plot_cluster_summary(features, assignments)

        import matplotlib.figure

        assert isinstance(fig, matplotlib.figure.Figure)
        import matplotlib.pyplot as plt

        plt.close(fig)

    def test_save_path(self, tmp_path):
        _, features = _make_test_features(6)
        assignments = cluster_channels(features, n_clusters=2, wavelength_aware=False)
        save_path = tmp_path / "clusters.png"
        fig = plot_cluster_summary(features, assignments, save_path=save_path)
        assert save_path.exists()
        import matplotlib.pyplot as plt

        plt.close(fig)


class TestGenerateClusterQC:
    def test_returns_figure(self, tmp_path):
        rng = np.random.RandomState(42)
        marker_dict = {
            "CD3": rng.randint(100, 5000, (64, 64), dtype=np.uint16),
            "CD20": rng.randint(100, 5000, (64, 64), dtype=np.uint16),
        }
        # Create fake processing results (no actual files)
        processing_results = {
            "CD3": {"output_path": None, "quality_score": 0.0, "status": "error"},
            "CD20": {"output_path": None, "quality_score": 0.0, "status": "error"},
        }
        fig = generate_cluster_qc(
            marker_dict, processing_results, ["CD3", "CD20"], cluster_id=0
        )

        import matplotlib.figure

        assert isinstance(fig, matplotlib.figure.Figure)
        import matplotlib.pyplot as plt

        plt.close(fig)

    def test_with_saved_files(self, tmp_path):
        import tifffile

        rng = np.random.RandomState(42)
        img = rng.randint(100, 5000, (64, 64), dtype=np.uint16)
        marker_dict = {"CD3": img}
        out_path = tmp_path / "CD3.tif"
        tifffile.imwrite(str(out_path), img)

        processing_results = {
            "CD3": {
                "output_path": out_path,
                "quality_score": 0.7,
                "status": "success",
            }
        }
        fig = generate_cluster_qc(
            marker_dict, processing_results, ["CD3"], cluster_id=1
        )

        import matplotlib.figure

        assert isinstance(fig, matplotlib.figure.Figure)
        import matplotlib.pyplot as plt

        plt.close(fig)
