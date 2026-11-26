"""Staging public HuBMAP CODEX downloads into KINTSUGI projects (no network)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from kintsugi.hubmap import (
    CYCLE_DIR_RE,
    HUBMAP_ACKNOWLEDGEMENT,
    DatasetInfo,
    HubmapError,
    fetch_dataset,
    find_source_dir,
    manifest_text,
    stage_dataset,
    validate_hubmap_id,
)

HBM = "HBM283.GNGR.785"
UUID = "0a22e0b91026c2164b238ff0e735077b"

ENTITY = {
    "hubmap_id": HBM,
    "uuid": UUID,
    "title": "CODEX data from the lymph node of a 17-year-old unknown male",
    "dataset_type": "CODEX",
    "status": "Published",
    "data_access_level": "public",
    "doi_url": "https://doi.org/10.35079/HBM283.GNGR.785",
    "group_name": "University of Florida TMC",
    "contains_human_genetic_sequences": False,
    "metadata": {
        "acquisition_instrument_vendor": "Keyence",
        "number_of_biomarker_imaging_rounds": 13,
        "resolution_x_value": 377.4038462,
        "operator_email": "someone@example.org",
    },
}

CODEX_EXPERIMENT = {
    "numCycles": 3,
    "numZPlanes": 2,
    "regionHeight": 2,
    "regionWidth": 2,
    "numChannels": 4,
    "xyResolution": 377.4038462,
    "zPitch": 1500.0,
    "aperture": 0.75,
    "tileOverlapX": 0.3,
    "tileOverlapY": 0.3,
    "wavelengths": [358, 488, 550, 650],
}


def _download(root: Path, cycles=(1, 2, 3), region: int = 1, nest: str = "src_LN") -> Path:
    """A minimal raw CODEX download: Globus nesting, cycle dirs, tiles, metadata."""
    source = root / HBM / nest
    for c in cycles:
        cyc = source / f"Cyc{c}_reg{region}"
        cyc.mkdir(parents=True)
        for tile in (1, 2, 3, 4):
            for z in (1, 2):
                for ch in (1, 2, 3, 4):
                    (cyc / f"{region}_{tile:05d}_Z{z:03d}_CH{ch}.tif").write_bytes(b"")
    (source / "experiment.json").write_text(json.dumps(CODEX_EXPERIMENT))
    (source / "channelnames.txt").write_text(
        "\n".join(
            ["DAPI-01", "Blank", "Blank", "Blank", "DAPI-02", "CD31", "CD8", "CD45"]
            + ["DAPI-03", "CD3", "CD20", "Ki67"]
        )
    )
    return root / HBM


@pytest.fixture
def project(tmp_path):
    from kintsugi.project import KintsugiProject

    root = tmp_path / "project"
    KintsugiProject.create(root, name="HBM283_LN")
    return root


def _fake_fetch(responses: dict):
    def fetch(url: str, as_json: bool = True):
        for key, value in responses.items():
            if key in url:
                if isinstance(value, Exception):
                    raise value
                return value if as_json else str(value)
        raise OSError(f"unexpected url {url}")

    return fetch


def test_cycle_directory_names():
    for name, cycle, region in [
        ("Cyc1_reg1", 1, "1"),
        ("cyc001_reg001_211118", 1, "1"),
        ("Cyc13_reg2", 13, "2"),
        ("cyc04", 4, None),
        ("Cycle3", 3, None),
    ]:
        m = CYCLE_DIR_RE.match(name)
        assert m and int(m["cycle"]) == cycle and m["region"] == region
    for name in ("processed", "cycles", "src_LN", "Cyc_reg1"):
        assert not CYCLE_DIR_RE.match(name)


def test_validate_hubmap_id():
    assert validate_hubmap_id(" hbm283.gngr.785 ") == HBM
    with pytest.raises(HubmapError):
        validate_hubmap_id("HBM283")


def test_fetch_dataset_keeps_public_metadata_only():
    info = fetch_dataset(HBM, fetch=_fake_fetch({"entities": ENTITY}))
    assert info.uuid == UUID
    assert info.doi_url.endswith(HBM)
    assert info.portal_url.endswith(UUID)
    assert info.acquisition["acquisition_instrument_vendor"] == "Keyence"
    assert "operator_email" not in info.acquisition


def test_fetch_dataset_refuses_protected_data():
    protected = {**ENTITY, "data_access_level": "protected"}
    with pytest.raises(HubmapError, match="not public"):
        fetch_dataset(HBM, fetch=_fake_fetch({"entities": protected}))


def test_manifest_uses_the_search_api_and_falls_back():
    assert manifest_text(HBM, fetch=_fake_fetch({"param-search": f"{HBM} /"})) == f"{HBM} /\n"
    offline = _fake_fetch({"param-search": OSError("offline")})
    assert manifest_text(HBM, fetch=offline) == f"{HBM} /\n"


def test_find_source_dir_descends_globus_nesting(tmp_path):
    download = _download(tmp_path / "dl")
    assert find_source_dir(download).name == "src_LN"
    assert find_source_dir(tmp_path).name == "src_LN"  # three levels down
    with pytest.raises(HubmapError, match="No CODEX cycle directories"):
        find_source_dir(tmp_path, max_depth=1)


def test_find_source_dir_without_cycles(tmp_path):
    (tmp_path / "empty" / "nested").mkdir(parents=True)
    with pytest.raises(HubmapError, match="No CODEX cycle directories"):
        find_source_dir(tmp_path / "empty")


def test_stage_links_cycles_and_writes_metadata(tmp_path, project):
    download = _download(tmp_path / "dl")
    dataset = DatasetInfo.from_entity(ENTITY)
    report = stage_dataset(download, project, hubmap_id=HBM, dataset=dataset)

    raw = project / "data" / "raw"
    for c in (1, 2, 3):
        link = raw / f"cyc{c:03d}"
        assert link.is_symlink() and link.resolve().name == f"Cyc{c}_reg1"
        assert len(list(link.glob("*_Z*_CH1.tif"))) == 8
    assert report.tiles_per_cycle == {1: 32, 2: 32, 3: 32}
    assert (raw / ".staged").read_text().startswith(f"staged from HuBMAP {HBM}")

    exp = json.loads((project / "meta" / "experiment.json").read_text())
    assert exp["n_cycles"] == 3 and exp["n_zplanes"] == 2
    assert exp["tile_rows"] == 2 and exp["tile_cols"] == 2
    assert exp["xy_pixel_size"] == pytest.approx(377.4038462)
    assert json.loads((project / "meta" / "experiment.source.json").read_text()) == (
        CODEX_EXPERIMENT
    )
    assert (project / "meta" / "CHANNELNAMES.txt").read_text().startswith("DAPI-01\nBlank")

    prov = json.loads((project / "meta" / "hubmap_provenance.json").read_text())
    assert prov["dataset"]["hubmap_id"] == HBM
    assert prov["citation"]["acknowledgement"] == HUBMAP_ACKNOWLEDGEMENT
    assert prov["citation"]["dataset_doi"] == ENTITY["doi_url"]
    assert prov["portal_url"].endswith(UUID)
    assert prov["mode"] == "symlink"
    assert report.warnings == []


def test_stage_refuses_to_overwrite_without_force(tmp_path, project):
    download = _download(tmp_path / "dl")
    stage_dataset(download, project)
    with pytest.raises(HubmapError, match="Already staged"):
        stage_dataset(download, project)
    report = stage_dataset(download, project, force=True)
    assert sorted(report.cycles) == [1, 2, 3]


def test_stage_copy_mode_and_region_selection(tmp_path, project):
    download = _download(tmp_path / "dl", cycles=(1, 2), region=2)
    with pytest.raises(HubmapError, match="region 1"):
        stage_dataset(download, project)
    report = stage_dataset(download, project, region=2, copy=True)
    link = project / "data" / "raw" / "cyc001"
    assert link.is_dir() and not link.is_symlink()
    assert report.region == 2


def test_stage_warns_about_gaps_and_missing_metadata(tmp_path, project):
    download = _download(tmp_path / "dl", cycles=(1, 3))
    source = download / "src_LN"
    (source / "experiment.json").unlink()
    (source / "channelnames.txt").unlink()
    report = stage_dataset(download, project)
    text = " ".join(report.warnings)
    assert "not contiguous" in text
    assert "No experiment.json" in text
    assert "No channelnames.txt" in text


def test_stage_requires_an_initialized_project(tmp_path):
    download = _download(tmp_path / "dl")
    (tmp_path / "plain").mkdir()
    with pytest.raises(HubmapError, match="kintsugi init"):
        stage_dataset(download, tmp_path / "plain")


def test_cli_stage_offline(tmp_path, project):
    from click.testing import CliRunner

    from kintsugi.cli import main

    download = _download(tmp_path / "dl")
    result = CliRunner().invoke(
        main,
        ["hubmap", "stage", str(download), str(project), "--hubmap-id", HBM, "--offline"],
    )
    assert result.exit_code == 0, result.output
    assert "Staged 3 cycles" in result.output
    prov = json.loads((project / "meta" / "hubmap_provenance.json").read_text())
    assert prov["dataset"]["hubmap_id"] == HBM


def test_stage_refuses_the_projects_own_raw_directory(tmp_path, project):
    download = _download(tmp_path / "dl")
    stage_dataset(download, project, copy=True)
    raw = project / "data" / "raw"
    before = sorted(p.name for p in (raw / "cyc001").iterdir())
    with pytest.raises(HubmapError, match="inside the project's raw data"):
        stage_dataset(raw, project, force=True, copy=True)
    assert sorted(p.name for p in (raw / "cyc001").iterdir()) == before


def test_stage_refuses_a_source_inside_a_replaced_target(tmp_path, project):
    from kintsugi.hubmap import _cycle_dirs, _preflight

    raw = project / "data" / "raw"
    inner = raw / "cyc001" / "nested"
    _download(inner, cycles=(1,))
    source = inner / HBM / "src_LN"
    cycles = _cycle_dirs(source)
    with pytest.raises(HubmapError, match="overlaps the staging target"):
        _preflight(source, tmp_path / "elsewhere", cycles, {1: raw / "cyc001"}, force=True)


def test_stage_conflict_changes_nothing(tmp_path, project):
    download = _download(tmp_path / "dl")
    raw = project / "data" / "raw"
    (raw / "cyc002").mkdir(parents=True)
    with pytest.raises(HubmapError, match="Already staged"):
        stage_dataset(download, project)
    assert not (raw / "cyc001").exists() and not (raw / "cyc001").is_symlink()
    assert not (raw / ".staged").exists()
    assert not (project / "meta" / "hubmap_provenance.json").exists()
