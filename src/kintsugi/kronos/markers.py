"""
Marker Name Mapping Between KINTSUGI and KRONOS

KRONOS was trained on 175 protein markers. This module handles mapping
KINTSUGI channel names (from CHANNELNAMES.txt) to KRONOS marker IDs
with normalization statistics (mean/std) from marker_metadata.csv.

Markers not in KRONOS's vocabulary are handled gracefully — they can
still be included with dataset-level normalization, or excluded.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)

# Common aliases: KINTSUGI name -> canonical KRONOS name
# Extend as needed when processing new panels
_MARKER_ALIASES = {
    # DAPI variants
    "DAPI": "DAPI",
    "DAPI-01": "DAPI",
    "DAPI-02": "DAPI",
    "DAPI-03": "DAPI",
    "Hoechst": "DAPI",
    "Hoechst 33342": "DAPI",
    # Common CD markers
    "CD3e": "CD3e",
    "CD3": "CD3e",
    "CD4": "CD4",
    "CD8": "CD8",
    "CD8a": "CD8",
    "CD11b": "CD11b",
    "CD11c": "CD11c",
    "CD14": "CD14",
    "CD15": "CD15",
    "CD16": "CD16",
    "CD19": "CD19",
    "CD20": "CD20",
    "CD21": "CD21",
    "CD25": "CD25",
    "CD27": "CD27",
    "CD31": "CD31",
    "CD34": "CD34",
    "CD38": "CD38",
    "CD44": "CD44",
    "CD45": "CD45",
    "CD45RA": "CD45RA",
    "CD45RO": "CD45RO",
    "CD56": "CD56",
    "CD57": "CD57",
    "CD68": "CD68",
    "CD69": "CD69",
    "CD103": "CD103",
    "CD107a": "CD107a",
    "CD127": "CD127",
    "CD138": "CD138",
    "CD163": "CD163",
    "CD206": "CD206",
    "CD209": "CD209",
    # Functional markers
    "FOXP3": "FOXP3",
    "FoxP3": "FOXP3",
    "Ki67": "Ki67",
    "Ki-67": "Ki67",
    "MKI67": "Ki67",
    "Granzyme B": "Granzyme B",
    "GranzymeB": "Granzyme B",
    "GZMB": "Granzyme B",
    "PD-1": "PD-1",
    "PDCD1": "PD-1",
    "PD-L1": "PD-L1",
    "CD274": "PD-L1",
    "CTLA-4": "CTLA-4",
    "CTLA4": "CTLA-4",
    "LAG-3": "LAG-3",
    "LAG3": "LAG-3",
    "TIM-3": "TIM-3",
    "HAVCR2": "TIM-3",
    # Structural markers
    "Vimentin": "Vimentin",
    "VIM": "Vimentin",
    "E-Cadherin": "E-cadherin",
    "E-cadherin": "E-cadherin",
    "CDH1": "E-cadherin",
    "Pan-Cytokeratin": "Pan-CK",
    "PanCK": "Pan-CK",
    "Pan-CK": "Pan-CK",
    "CK": "Pan-CK",
    "Cytokeratin": "Pan-CK",
    "SMA": "aSMA",
    "aSMA": "aSMA",
    "alpha-SMA": "aSMA",
    "Collagen I": "Collagen I",
    "Collagen IV": "Collagen IV",
    # Other common markers
    "HLA-DR": "HLA-DR",
    "HLADR": "HLA-DR",
    "Podoplanin": "Podoplanin",
    "PDPN": "Podoplanin",
    "ICOS": "ICOS",
    "BCL6": "BCL6",
    "BCL-6": "BCL6",
    "IRF4": "IRF4",
    "Tbet": "T-bet",
    "T-bet": "T-bet",
    "TBX21": "T-bet",
}


@dataclass
class MarkerMapping:
    """Mapping result for a single KINTSUGI channel.

    Attributes
    ----------
    kintsugi_name : str
        Original channel name from CHANNELNAMES.txt.
    kronos_name : str or None
        Matched KRONOS marker name, or None if not found.
    kronos_id : int or None
        KRONOS marker index (row in marker_metadata.csv).
    mean : float
        Normalization mean (from KRONOS metadata or computed from data).
    std : float
        Normalization std (from KRONOS metadata or computed from data).
    matched : bool
        True if marker was found in KRONOS vocabulary.
    cycle : int
        Cycle number (1-indexed).
    channel_idx : int
        Channel index within cycle (0-indexed).
    """

    kintsugi_name: str
    kronos_name: str | None = None
    kronos_id: int | None = None
    mean: float = 0.0
    std: float = 1.0
    matched: bool = False
    cycle: int = 0
    channel_idx: int = 0


@dataclass
class MarkerMapper:
    """Maps KINTSUGI channel names to KRONOS marker IDs.

    Loads marker_metadata.csv from the KRONOS model assets and matches
    against KINTSUGI's CHANNELNAMES.txt.

    Parameters
    ----------
    kronos_metadata_path : Path, optional
        Path to marker_metadata.csv. If None, loads from HF cache.
    """

    kronos_markers: dict[str, dict[str, Any]] = field(default_factory=dict)
    mappings: list[MarkerMapping] = field(default_factory=list)
    _loaded: bool = False

    def load_kronos_metadata(
        self,
        metadata_path: str | Path | None = None,
        cache_dir: str | Path = "./model_assets",
    ) -> None:
        """Load KRONOS marker metadata (names, mean, std).

        Parameters
        ----------
        metadata_path : Path, optional
            Direct path to marker_metadata.csv.
        cache_dir : Path
            Directory containing cached model assets.
        """
        import csv

        if metadata_path is not None:
            csv_path = Path(metadata_path)
        else:
            # Try common locations
            candidates = [
                Path(cache_dir) / "marker_metadata.csv",
                Path(cache_dir) / "models--MahmoodLab--KRONOS" / "marker_metadata.csv",
            ]
            # Also check HF cache structure
            hf_cache = Path(cache_dir) / "models--MahmoodLab--KRONOS"
            if hf_cache.exists():
                # HF cache uses snapshot directories
                for snapshot_dir in sorted(hf_cache.glob("snapshots/*/marker_metadata.csv")):
                    candidates.insert(0, snapshot_dir)

            csv_path = None
            for c in candidates:
                if c.exists():
                    csv_path = c
                    break

            if csv_path is None:
                raise FileNotFoundError(
                    f"marker_metadata.csv not found. Searched: {candidates}\n"
                    "Download the KRONOS model first:\n"
                    "  from kronos import create_model_from_pretrained\n"
                    "  create_model_from_pretrained("
                    'checkpoint_path="hf_hub:MahmoodLab/kronos", '
                    f'cache_dir="{cache_dir}")'
                )

        logger.info("Loading KRONOS marker metadata from %s", csv_path)
        self.kronos_markers.clear()

        with open(csv_path, newline="") as f:
            reader = csv.DictReader(f)
            for i, row in enumerate(reader):
                name = row.get("marker_name", row.get("name", row.get("marker", "")))
                name = name.strip()
                if not name:
                    continue
                # marker_id is the sinusoidal positional encoding ID from KRONOS
                # Falls back to row index if not present in CSV
                marker_id = int(row.get("marker_id", row.get("id", i)))
                self.kronos_markers[name.lower()] = {
                    "name": name,
                    "id": marker_id,
                    "mean": float(row.get("marker_mean", row.get("mean", 0.0))),
                    "std": float(row.get("marker_std", row.get("std", 1.0))),
                }

        self._loaded = True
        logger.info("Loaded %d KRONOS markers", len(self.kronos_markers))

    def map_channels(
        self,
        channel_dict: dict[int, list[str]],
        skip_blank: bool = True,
        skip_dapi: bool = False,
    ) -> list[MarkerMapping]:
        """Map KINTSUGI channel names to KRONOS markers.

        Parameters
        ----------
        channel_dict : dict
            {cycle_number: [channel_names]} from load_channel_names().
        skip_blank : bool
            Skip channels named 'Blank' or 'Empty'.
        skip_dapi : bool
            Skip DAPI channels (nuclear stain, not protein marker).

        Returns
        -------
        list of MarkerMapping
        """
        if not self._loaded:
            raise RuntimeError("Call load_kronos_metadata() first.")

        self.mappings = []
        skip_names = {"blank", "empty", "none", "na", "n/a"}
        if skip_dapi:
            skip_names.update({"dapi", "hoechst", "hoechst 33342"})

        for cycle, channels in sorted(channel_dict.items()):
            for ch_idx, ch_name in enumerate(channels):
                clean = ch_name.strip()
                if skip_blank and clean.lower() in skip_names:
                    continue
                if skip_dapi and clean.lower().startswith("dapi"):
                    continue

                mapping = self._match_single(clean, cycle, ch_idx)
                self.mappings.append(mapping)

        matched = sum(1 for m in self.mappings if m.matched)
        logger.info(
            "Mapped %d / %d channels to KRONOS markers",
            matched,
            len(self.mappings),
        )
        return self.mappings

    def _match_single(self, name: str, cycle: int, ch_idx: int) -> MarkerMapping:
        """Match a single channel name to KRONOS vocabulary."""
        # Try direct match (case-insensitive)
        lower = name.lower()
        if lower in self.kronos_markers:
            info = self.kronos_markers[lower]
            return MarkerMapping(
                kintsugi_name=name,
                kronos_name=info["name"],
                kronos_id=info["id"],
                mean=info["mean"],
                std=info["std"],
                matched=True,
                cycle=cycle,
                channel_idx=ch_idx,
            )

        # Try alias lookup
        canonical = _MARKER_ALIASES.get(name) or _MARKER_ALIASES.get(name.upper())
        if canonical and canonical.lower() in self.kronos_markers:
            info = self.kronos_markers[canonical.lower()]
            return MarkerMapping(
                kintsugi_name=name,
                kronos_name=info["name"],
                kronos_id=info["id"],
                mean=info["mean"],
                std=info["std"],
                matched=True,
                cycle=cycle,
                channel_idx=ch_idx,
            )

        # Try fuzzy: strip suffixes like _b, _c (unique-ified names)
        stripped = re.sub(r"_[a-z]$", "", name)
        if stripped != name:
            return self._match_single(stripped, cycle, ch_idx)

        # Try fuzzy: strip cycle prefixes like "DAPI-02"
        dash_stripped = re.sub(r"-\d+$", "", name)
        if dash_stripped != name and dash_stripped.lower() in self.kronos_markers:
            info = self.kronos_markers[dash_stripped.lower()]
            return MarkerMapping(
                kintsugi_name=name,
                kronos_name=info["name"],
                kronos_id=info["id"],
                mean=info["mean"],
                std=info["std"],
                matched=True,
                cycle=cycle,
                channel_idx=ch_idx,
            )

        # No match — will use dataset-level normalization
        return MarkerMapping(
            kintsugi_name=name,
            kronos_name=None,
            kronos_id=None,
            mean=0.0,
            std=1.0,
            matched=False,
            cycle=cycle,
            channel_idx=ch_idx,
        )

    def get_normalization_arrays(
        self,
        selected: list[MarkerMapping] | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Get mean/std arrays for the mapped markers.

        Parameters
        ----------
        selected : list, optional
            Subset of mappings. Uses all mappings if None.

        Returns
        -------
        tuple of np.ndarray
            (mean_array, std_array), each shape (n_markers,).
        """
        mappings = selected or self.mappings
        means = np.array([m.mean for m in mappings], dtype=np.float32)
        stds = np.array([m.std for m in mappings], dtype=np.float32)
        return means, stds

    def get_marker_ids(
        self,
        selected: list[MarkerMapping] | None = None,
    ) -> np.ndarray | None:
        """Get KRONOS marker IDs for sinusoidal positional encoding.

        Parameters
        ----------
        selected : list, optional
            Subset of mappings. Uses all mappings if None.

        Returns
        -------
        np.ndarray or None
            Integer marker IDs shape (n_markers,), or None if no markers
            have KRONOS IDs (all unmatched).
        """
        mappings = selected or self.mappings
        ids = []
        has_any = False
        for m in mappings:
            if m.kronos_id is not None:
                ids.append(m.kronos_id)
                has_any = True
            else:
                # Default: sequential with +4 offset (KRONOS convention)
                ids.append(len(ids) + 4)
        return np.array(ids, dtype=np.int64) if has_any else None

    def summary(self) -> str:
        """Human-readable mapping summary."""
        lines = ["KRONOS Marker Mapping Summary", "=" * 40]
        matched = [m for m in self.mappings if m.matched]
        unmatched = [m for m in self.mappings if not m.matched]

        lines.append(f"Total channels: {len(self.mappings)}")
        lines.append(f"Matched to KRONOS: {len(matched)}")
        lines.append(f"Unmatched: {len(unmatched)}")

        if matched:
            lines.append("\nMatched markers:")
            for m in matched:
                lines.append(f"  {m.kintsugi_name} -> {m.kronos_name} (id={m.kronos_id})")

        if unmatched:
            lines.append("\nUnmatched markers (will use dataset-level normalization):")
            for m in unmatched:
                lines.append(f"  {m.kintsugi_name}")

        return "\n".join(lines)

    @classmethod
    def from_project(
        cls,
        project_dir: str | Path,
        cache_dir: str | Path = "./model_assets",
        metadata_path: str | Path | None = None,
        skip_blank: bool = True,
        skip_dapi: bool = False,
    ) -> MarkerMapper:
        """Create a mapper from a KINTSUGI project directory.

        Loads channel names from meta/CHANNELNAMES.txt and matches
        against KRONOS marker metadata.

        Parameters
        ----------
        project_dir : Path
            Path to KINTSUGI project root.
        cache_dir : Path
            Directory for KRONOS model cache.
        metadata_path : Path, optional
            Direct path to marker_metadata.csv.
        skip_blank : bool
            Skip Blank channels.
        skip_dapi : bool
            Skip DAPI channels.

        Returns
        -------
        MarkerMapper
        """
        from kintsugi.project import KintsugiProject

        project = KintsugiProject.load(project_dir)
        channel_dict = project.load_channel_names()
        if channel_dict is None:
            raise FileNotFoundError(
                f"No CHANNELNAMES.txt found in {project.paths.meta}. "
                "Create one with marker names for each channel."
            )

        mapper = cls()
        mapper.load_kronos_metadata(metadata_path=metadata_path, cache_dir=cache_dir)
        mapper.map_channels(channel_dict, skip_blank=skip_blank, skip_dapi=skip_dapi)
        return mapper
