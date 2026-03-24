"""Tests for _load_channel_names_lightweight and importlib-based tissue type parsing."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from kintsugi.cli import _load_channel_names_lightweight


class TestLoadChannelNamesLightweight:
    """Tests for the stdlib-only channel name parser."""

    def test_simple_list_format(self, tmp_path):
        """Parse simple DAPI-XX delimited format."""
        meta = tmp_path / "meta"
        meta.mkdir()
        (meta / "CHANNELNAMES.txt").write_text(
            "DAPI-01\nBlank\nBlank\nBlank\nDAPI-02\nCD31\nCD8\nCD45\n"
        )
        result = _load_channel_names_lightweight(meta, channels_per_cycle=4)
        assert sorted(result.keys()) == [1, 2]
        assert result[1][0] == "DAPI-01"
        assert result[2][0] == "DAPI_cyc02"
        assert result[2][1] == "CD31"

    def test_cycle_prefixed_format(self, tmp_path):
        """Parse colon-separated cycle-prefixed format."""
        meta = tmp_path / "meta"
        meta.mkdir()
        (meta / "CHANNELNAMES.txt").write_text(
            "1: DAPI, Blank, Blank, Blank\n2: DAPI, CD31, CD8, CD45\n"
        )
        result = _load_channel_names_lightweight(meta, channels_per_cycle=4)
        assert sorted(result.keys()) == [1, 2]
        assert result[1][0] == "DAPI"
        assert result[2][0] == "DAPI_cyc02"

    def test_tab_separated_format(self, tmp_path):
        """Parse tab-separated format."""
        meta = tmp_path / "meta"
        meta.mkdir()
        (meta / "CHANNELNAMES.txt").write_text(
            "1\tDAPI\tBlank\tBlank\tBlank\n2\tDAPI\tCD31\tCD8\tCD45\n"
        )
        result = _load_channel_names_lightweight(meta, channels_per_cycle=5)
        assert sorted(result.keys()) == [1, 2]

    def test_uniquification_blanks(self, tmp_path):
        """Duplicate names get cycle+channel suffixes."""
        meta = tmp_path / "meta"
        meta.mkdir()
        (meta / "CHANNELNAMES.txt").write_text(
            "DAPI-01\nBlank\nBlank\nBlank\nDAPI-02\nCD31\nCD8\nBlank\n"
        )
        result = _load_channel_names_lightweight(meta, channels_per_cycle=4)
        # Blanks should all be unique
        all_names = [n for names in result.values() for n in names]
        non_dapi = [n for n in all_names if not n.upper().startswith("DAPI")]
        assert len(non_dapi) == len(set(non_dapi)), f"Duplicate names found: {non_dapi}"

    def test_uniquification_dapi_cycles(self, tmp_path):
        """DAPI after cycle 1 gets _cycNN suffix."""
        meta = tmp_path / "meta"
        meta.mkdir()
        (meta / "CHANNELNAMES.txt").write_text(
            "1: DAPI, CD3, CD4, CD8\n2: DAPI, CD20, Ki67, FOXP3\n3: DAPI, CD31, SMA, Blank\n"
        )
        result = _load_channel_names_lightweight(meta, channels_per_cycle=4)
        assert result[1][0] == "DAPI"
        assert result[2][0] == "DAPI_cyc02"
        assert result[3][0] == "DAPI_cyc03"

    def test_file_not_found(self, tmp_path):
        """Raises FileNotFoundError when no channel file exists."""
        meta = tmp_path / "meta"
        meta.mkdir()
        with pytest.raises(FileNotFoundError):
            _load_channel_names_lightweight(meta)

    def test_empty_file(self, tmp_path):
        """Raises ValueError when file is empty."""
        meta = tmp_path / "meta"
        meta.mkdir()
        (meta / "CHANNELNAMES.txt").write_text("# just a comment\n\n")
        with pytest.raises(ValueError, match="empty"):
            _load_channel_names_lightweight(meta)

    def test_comments_skipped(self, tmp_path):
        """Lines starting with # are ignored."""
        meta = tmp_path / "meta"
        meta.mkdir()
        (meta / "CHANNELNAMES.txt").write_text(
            "# Header comment\n1: DAPI, CD3, CD4, CD8\n# Another comment\n"
        )
        result = _load_channel_names_lightweight(meta, channels_per_cycle=4)
        assert 1 in result

    def test_alternate_filenames(self, tmp_path):
        """Finds channel names under alternate filenames."""
        meta = tmp_path / "meta"
        meta.mkdir()
        (meta / "markers.txt").write_text("1: DAPI, CD3, CD4, CD8\n")
        result = _load_channel_names_lightweight(meta, channels_per_cycle=4)
        assert 1 in result


class TestImportlibTissueType:
    """Test that parse_tissue_type works via importlib without importing kintsugi.signal."""

    def test_tissue_type_detection(self):
        """parse_tissue_type detects tissue from project name via importlib."""
        import importlib.util

        _mod_name = "_test_batch_multi"
        _path = Path(__file__).parent.parent / "src" / "kintsugi" / "signal" / "batch_multi.py"
        _spec = importlib.util.spec_from_file_location(_mod_name, _path)
        assert _spec is not None and _spec.loader is not None
        _mod = importlib.util.module_from_spec(_spec)
        sys.modules[_mod_name] = _mod
        try:
            _spec.loader.exec_module(_mod)
            assert _mod.parse_tissue_type("CX_19-002_spleen_CC2-A", None) == "spleen"
            assert _mod.parse_tissue_type("CX_19-003_LN_CC1-B", None) == "lymph_node"
            assert _mod.parse_tissue_type("some_unknown_project", None) == "unknown"
        finally:
            sys.modules.pop(_mod_name, None)

    def test_no_signal_package_imported(self):
        """Importing batch_multi via importlib should NOT load kintsugi.signal."""
        import importlib.util

        # Remove any cached signal import
        signal_keys = [k for k in sys.modules if k.startswith("kintsugi.signal")]
        saved = {k: sys.modules.pop(k) for k in signal_keys}

        _mod_name = "_test_batch_multi_isolated"
        _path = Path(__file__).parent.parent / "src" / "kintsugi" / "signal" / "batch_multi.py"
        _spec = importlib.util.spec_from_file_location(_mod_name, _path)
        _mod = importlib.util.module_from_spec(_spec)
        sys.modules[_mod_name] = _mod
        try:
            _spec.loader.exec_module(_mod)
            # kintsugi.signal should NOT have been imported
            assert "kintsugi.signal" not in sys.modules
        finally:
            sys.modules.pop(_mod_name, None)
            # Restore saved modules
            sys.modules.update(saved)
