"""
Parameter Learning System for KINTSUGI

Learns optimal processing parameters based on tissue type, marker name,
and image characteristics. Provides recommendations for new images based
on historical successes.

Architecture:
    - SQLite database stores all parameter records
    - Fuzzy matching for tissue/marker names
    - Weighted scoring based on quality outcomes and recency
    - Confidence intervals for recommendations
"""

from __future__ import annotations

import json
import re
import sqlite3
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

# Common tissue type aliases for normalization
TISSUE_ALIASES = {
    # Lymphoid
    "tonsil": ["tonsil", "palatine tonsil", "tonsillar"],
    "lymph_node": ["lymph node", "ln", "lymphnode", "lymph"],
    "spleen": ["spleen", "splenic"],
    "thymus": ["thymus", "thymic"],
    # GI tract
    "colon": ["colon", "colonic", "large intestine", "large bowel"],
    "small_intestine": ["small intestine", "small bowel", "ileum", "jejunum", "duodenum"],
    "stomach": ["stomach", "gastric"],
    "esophagus": ["esophagus", "esophageal"],
    # Solid organs
    "liver": ["liver", "hepatic"],
    "kidney": ["kidney", "renal"],
    "pancreas": ["pancreas", "pancreatic"],
    "lung": ["lung", "pulmonary"],
    "heart": ["heart", "cardiac", "myocardium"],
    # Skin
    "skin": ["skin", "dermis", "epidermis", "cutaneous"],
    # Brain
    "brain": ["brain", "cerebral", "cerebrum", "cortex", "hippocampus"],
    # Reproductive
    "breast": ["breast", "mammary"],
    "prostate": ["prostate", "prostatic"],
    "ovary": ["ovary", "ovarian"],
    "uterus": ["uterus", "uterine", "endometrium"],
    # Tumors
    "tumor": ["tumor", "tumour", "cancer", "carcinoma", "adenocarcinoma", "malignant"],
    "metastasis": ["metastasis", "metastatic", "met", "mets"],
}

# Common marker aliases for normalization
MARKER_ALIASES = {
    # Nuclear
    "dapi": ["dapi", "hoechst", "nuclear", "nuclei"],
    # T cells
    "cd3": ["cd3", "cd3e", "t cell", "t-cell", "pan-t"],
    "cd4": ["cd4", "helper t", "th"],
    "cd8": ["cd8", "cd8a", "cytotoxic t", "ctl"],
    "foxp3": ["foxp3", "treg", "regulatory t"],
    # B cells
    "cd20": ["cd20", "b cell", "b-cell", "pan-b"],
    "cd19": ["cd19", "b cell", "b-cell"],
    # Myeloid
    "cd68": ["cd68", "macrophage", "mac"],
    "cd163": ["cd163", "m2 macrophage"],
    "cd11c": ["cd11c", "dendritic", "dc"],
    "cd14": ["cd14", "monocyte"],
    # Other immune
    "cd45": ["cd45", "lca", "leukocyte common", "pan-leukocyte"],
    "pd1": ["pd1", "pd-1", "pdcd1"],
    "pdl1": ["pdl1", "pd-l1", "cd274"],
    # Structural
    "panck": ["panck", "pan-ck", "pan-cytokeratin", "cytokeratin", "ck", "epithelial"],
    "ecad": ["ecad", "e-cadherin", "cdh1"],
    "vimentin": ["vimentin", "vim", "mesenchymal"],
    "sma": ["sma", "α-sma", "alpha-sma", "smooth muscle actin", "acta2"],
    "cd31": ["cd31", "pecam", "endothelial"],
    # Proliferation
    "ki67": ["ki67", "ki-67", "mki67", "proliferation"],
    # Functional
    "granzyme_b": ["granzyme b", "gzmb", "grb"],
    "perforin": ["perforin", "prf1"],
    # Autofluorescence/blank
    "blank": ["blank", "af", "autofluorescence", "auto", "empty"],
}


def normalize_tissue_type(tissue: str) -> str:
    """Normalize tissue type to canonical form."""
    tissue_lower = tissue.lower().strip()

    for canonical, aliases in TISSUE_ALIASES.items():
        for alias in aliases:
            if alias in tissue_lower or tissue_lower in alias:
                return canonical

    # Return cleaned version if no match
    return re.sub(r"[^a-z0-9]", "_", tissue_lower)


def normalize_marker_name(marker: str) -> str:
    """Normalize marker name to canonical form."""
    marker_lower = marker.lower().strip()

    for canonical, aliases in MARKER_ALIASES.items():
        for alias in aliases:
            if alias in marker_lower or marker_lower in alias:
                return canonical

    # Return cleaned version if no match
    return re.sub(r"[^a-z0-9]", "_", marker_lower)


@dataclass
class ImageCharacteristics:
    """Measured characteristics of an image for matching."""

    snr: float = 0.0
    dynamic_range: float = 0.0
    mean_intensity: float = 0.0
    background_level: float = 0.0
    saturation_ratio: float = 0.0
    sparsity: float = 0.0  # Ratio of non-zero pixels

    def to_dict(self) -> dict[str, float]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ImageCharacteristics:
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})

    @classmethod
    def from_image(cls, image: np.ndarray) -> ImageCharacteristics:
        """Compute characteristics from an image array."""
        # Sample for large images
        if image.size > 4_000_000:
            step = int(np.sqrt(image.size / 1_000_000))
            sample = image[::step, ::step]
        else:
            sample = image

        sample = sample.astype(np.float64)

        # Basic stats
        mean_intensity = float(np.mean(sample))
        p10 = float(np.percentile(sample, 10))
        p90 = float(np.percentile(sample, 90))

        # SNR
        noise_region = sample[sample < p10]
        noise_std = float(np.std(noise_region)) if len(noise_region) > 0 else 1.0
        snr = p90 / noise_std if noise_std > 0 else 0.0

        # Dynamic range
        dynamic_range = float(np.max(sample) - np.min(sample))

        # Saturation
        max_val = 65535 if sample.max() > 255 else 255
        saturation_ratio = float(np.sum(sample >= max_val * 0.99) / sample.size)

        # Sparsity
        sparsity = float(np.sum(sample > p10) / sample.size)

        return cls(
            snr=snr,
            dynamic_range=dynamic_range,
            mean_intensity=mean_intensity,
            background_level=p10,
            saturation_ratio=saturation_ratio,
            sparsity=sparsity,
        )

    def similarity_score(self, other: ImageCharacteristics) -> float:
        """Compute similarity score between two image characteristics (0-1)."""
        scores = []

        # SNR similarity (log scale)
        if self.snr > 0 and other.snr > 0:
            snr_diff = abs(np.log10(self.snr + 1) - np.log10(other.snr + 1))
            scores.append(max(0, 1 - snr_diff / 2))

        # Dynamic range similarity (normalized)
        if self.dynamic_range > 0 and other.dynamic_range > 0:
            dr_ratio = min(self.dynamic_range, other.dynamic_range) / max(
                self.dynamic_range, other.dynamic_range
            )
            scores.append(dr_ratio)

        # Sparsity similarity
        sparsity_diff = abs(self.sparsity - other.sparsity)
        scores.append(max(0, 1 - sparsity_diff * 2))

        return np.mean(scores) if scores else 0.5


@dataclass
class ParameterRecord:
    """A record of parameters used for processing."""

    id: int | None = None
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())

    # Context
    tissue_type: str = ""
    tissue_type_normalized: str = ""
    marker_name: str = ""
    marker_name_normalized: str = ""

    # Image characteristics
    characteristics: ImageCharacteristics = field(default_factory=ImageCharacteristics)

    # Processing parameters
    operation: str = ""  # e.g., "blank_subtraction", "denoise", "full_pipeline"
    parameters: dict[str, Any] = field(default_factory=dict)

    # Outcome
    quality_before: float = 0.0
    quality_after: float = 0.0
    quality_improvement: float = 0.0
    user_approved: bool = False
    user_notes: str = ""

    # Metadata
    project_path: str = ""
    channel_name: str = ""

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["characteristics"] = self.characteristics.to_dict()
        return d


class ParameterLearningEngine:
    """
    Learns and recommends processing parameters based on tissue type and marker.

    Usage:
        engine = ParameterLearningEngine("/path/to/project")

        # Record a successful parameter set
        engine.record_success(
            tissue_type="tonsil",
            marker_name="CD3",
            operation="blank_subtraction",
            parameters={"blank_scale_factor": 0.8, "erosion": 1},
            quality_before=0.6,
            quality_after=0.92,
        )

        # Get recommendations for a new image
        recommendations = engine.recommend_parameters(
            tissue_type="tonsil",
            marker_name="CD4",  # Similar T-cell marker
            operation="blank_subtraction",
        )
    """

    def __init__(self, project_path: str | None = None, global_db: bool = True):
        """
        Initialize the learning engine.

        Args:
            project_path: Path to project (for project-specific learning)
            global_db: Whether to also use global learning database
        """
        self.project_path = Path(project_path) if project_path else None
        self.global_db = global_db

        # Project-specific database
        if self.project_path:
            self.project_db_path = self.project_path / ".claude_workflow" / "parameter_learning.db"
            self.project_db_path.parent.mkdir(parents=True, exist_ok=True)
            self._init_db(self.project_db_path)
        else:
            self.project_db_path = None

        # Global database (user-wide learning)
        if global_db:
            self.global_db_path = Path.home() / ".kintsugi" / "parameter_learning.db"
            self.global_db_path.parent.mkdir(parents=True, exist_ok=True)
            self._init_db(self.global_db_path)
        else:
            self.global_db_path = None

    def _init_db(self, db_path: Path):
        """Initialize SQLite database schema."""
        conn = sqlite3.connect(str(db_path))
        cursor = conn.cursor()

        cursor.execute("""
            CREATE TABLE IF NOT EXISTS parameter_records (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                timestamp TEXT NOT NULL,

                tissue_type TEXT,
                tissue_type_normalized TEXT,
                marker_name TEXT,
                marker_name_normalized TEXT,

                characteristics_json TEXT,

                operation TEXT NOT NULL,
                parameters_json TEXT NOT NULL,

                quality_before REAL,
                quality_after REAL,
                quality_improvement REAL,
                user_approved INTEGER DEFAULT 0,
                user_notes TEXT,

                project_path TEXT,
                channel_name TEXT
            )
        """)

        # Index for fast lookups
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_tissue_marker
            ON parameter_records(tissue_type_normalized, marker_name_normalized, operation)
        """)

        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_approved
            ON parameter_records(user_approved, quality_improvement)
        """)

        conn.commit()
        conn.close()

    def record_parameters(
        self,
        tissue_type: str,
        marker_name: str,
        operation: str,
        parameters: dict[str, Any],
        quality_before: float = 0.0,
        quality_after: float = 0.0,
        user_approved: bool = False,
        user_notes: str = "",
        characteristics: ImageCharacteristics | None = None,
        channel_name: str = "",
        to_global: bool = True,
    ) -> int:
        """
        Record a parameter set with its outcome.

        Args:
            tissue_type: Tissue type (will be normalized)
            marker_name: Marker name (will be normalized)
            operation: Operation name (e.g., "blank_subtraction")
            parameters: Parameter dictionary
            quality_before: Quality score before processing
            quality_after: Quality score after processing
            user_approved: Whether user approved the result
            user_notes: Optional user notes
            characteristics: Image characteristics (optional)
            channel_name: Channel identifier
            to_global: Also save to global database

        Returns:
            Record ID
        """
        record = ParameterRecord(
            tissue_type=tissue_type,
            tissue_type_normalized=normalize_tissue_type(tissue_type),
            marker_name=marker_name,
            marker_name_normalized=normalize_marker_name(marker_name),
            operation=operation,
            parameters=parameters,
            quality_before=quality_before,
            quality_after=quality_after,
            quality_improvement=quality_after - quality_before,
            user_approved=user_approved,
            user_notes=user_notes,
            characteristics=characteristics or ImageCharacteristics(),
            project_path=str(self.project_path) if self.project_path else "",
            channel_name=channel_name,
        )

        record_id = None

        # Save to project database
        if self.project_db_path:
            record_id = self._insert_record(self.project_db_path, record)

        # Save to global database
        if to_global and self.global_db_path:
            self._insert_record(self.global_db_path, record)

        return record_id or 0

    def _insert_record(self, db_path: Path, record: ParameterRecord) -> int:
        """Insert a record into the database."""
        conn = sqlite3.connect(str(db_path))
        cursor = conn.cursor()

        cursor.execute(
            """
            INSERT INTO parameter_records (
                timestamp, tissue_type, tissue_type_normalized,
                marker_name, marker_name_normalized,
                characteristics_json, operation, parameters_json,
                quality_before, quality_after, quality_improvement,
                user_approved, user_notes, project_path, channel_name
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
            (
                record.timestamp,
                record.tissue_type,
                record.tissue_type_normalized,
                record.marker_name,
                record.marker_name_normalized,
                json.dumps(record.characteristics.to_dict()),
                record.operation,
                json.dumps(record.parameters),
                record.quality_before,
                record.quality_after,
                record.quality_improvement,
                1 if record.user_approved else 0,
                record.user_notes,
                record.project_path,
                record.channel_name,
            ),
        )

        record_id = cursor.lastrowid
        conn.commit()
        conn.close()

        return record_id

    def recommend_parameters(
        self,
        tissue_type: str,
        marker_name: str,
        operation: str,
        characteristics: ImageCharacteristics | None = None,
        min_quality_improvement: float = 0.0,
        approved_only: bool = True,
        max_results: int = 5,
    ) -> dict[str, Any]:
        """
        Recommend parameters based on similar tissue/marker combinations.

        Args:
            tissue_type: Target tissue type
            marker_name: Target marker name
            operation: Operation to get parameters for
            characteristics: Current image characteristics for similarity matching
            min_quality_improvement: Minimum quality improvement to consider
            approved_only: Only use user-approved records
            max_results: Maximum number of similar records to consider

        Returns:
            Dictionary with recommended parameters and confidence
        """
        tissue_norm = normalize_tissue_type(tissue_type)
        marker_norm = normalize_marker_name(marker_name)

        # Query both databases
        all_records = []

        if self.project_db_path and self.project_db_path.exists():
            records = self._query_similar_records(
                self.project_db_path,
                tissue_norm,
                marker_norm,
                operation,
                min_quality_improvement,
                approved_only,
            )
            # Project records get higher weight
            for r in records:
                r["source"] = "project"
                r["weight"] = 1.5
            all_records.extend(records)

        if self.global_db_path and self.global_db_path.exists():
            records = self._query_similar_records(
                self.global_db_path,
                tissue_norm,
                marker_norm,
                operation,
                min_quality_improvement,
                approved_only,
            )
            for r in records:
                r["source"] = "global"
                r["weight"] = 1.0
            all_records.extend(records)

        if not all_records:
            return {
                "status": "no_history",
                "message": f"No parameter history for {tissue_type}/{marker_name}/{operation}",
                "recommended_parameters": None,
                "confidence": 0.0,
                "similar_records": [],
            }

        # Score and rank records
        scored_records = self._score_records(
            all_records,
            tissue_norm,
            marker_norm,
            characteristics,
        )

        # Take top records
        top_records = sorted(scored_records, key=lambda x: x["score"], reverse=True)[:max_results]

        # Aggregate parameters with weighted voting
        recommended_params = self._aggregate_parameters(top_records)

        # Calculate confidence
        confidence = self._calculate_confidence(top_records, tissue_norm, marker_norm)

        return {
            "status": "success",
            "tissue_type": tissue_type,
            "tissue_type_normalized": tissue_norm,
            "marker_name": marker_name,
            "marker_name_normalized": marker_norm,
            "operation": operation,
            "recommended_parameters": recommended_params,
            "confidence": round(confidence, 3),
            "confidence_level": self._confidence_level(confidence),
            "based_on_records": len(top_records),
            "similar_records": [
                {
                    "tissue": r["tissue_type"],
                    "marker": r["marker_name"],
                    "quality_improvement": r["quality_improvement"],
                    "score": round(r["score"], 3),
                    "parameters": r["parameters"],
                }
                for r in top_records[:3]
            ],
        }

    def _query_similar_records(
        self,
        db_path: Path,
        tissue_norm: str,
        marker_norm: str,
        operation: str,
        min_improvement: float,
        approved_only: bool,
    ) -> list[dict[str, Any]]:
        """Query records similar to the target tissue/marker."""
        conn = sqlite3.connect(str(db_path))
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        # Build query with fuzzy matching
        query = """
            SELECT * FROM parameter_records
            WHERE operation = ?
            AND quality_improvement >= ?
        """
        params = [operation, min_improvement]

        if approved_only:
            query += " AND user_approved = 1"

        # Order by relevance
        query += """
            ORDER BY
                CASE
                    WHEN tissue_type_normalized = ? AND marker_name_normalized = ? THEN 4
                    WHEN tissue_type_normalized = ? THEN 3
                    WHEN marker_name_normalized = ? THEN 2
                    ELSE 1
                END DESC,
                quality_improvement DESC,
                timestamp DESC
            LIMIT 50
        """
        params.extend([tissue_norm, marker_norm, tissue_norm, marker_norm])

        cursor.execute(query, params)
        rows = cursor.fetchall()
        conn.close()

        records = []
        for row in rows:
            records.append(
                {
                    "id": row["id"],
                    "timestamp": row["timestamp"],
                    "tissue_type": row["tissue_type"],
                    "tissue_type_normalized": row["tissue_type_normalized"],
                    "marker_name": row["marker_name"],
                    "marker_name_normalized": row["marker_name_normalized"],
                    "characteristics": (
                        json.loads(row["characteristics_json"])
                        if row["characteristics_json"]
                        else {}
                    ),
                    "parameters": json.loads(row["parameters_json"]),
                    "quality_before": row["quality_before"],
                    "quality_after": row["quality_after"],
                    "quality_improvement": row["quality_improvement"],
                    "user_approved": bool(row["user_approved"]),
                }
            )

        return records

    def _score_records(
        self,
        records: list[dict[str, Any]],
        tissue_norm: str,
        marker_norm: str,
        characteristics: ImageCharacteristics | None,
    ) -> list[dict[str, Any]]:
        """Score records based on similarity to target."""
        for record in records:
            score = 0.0

            # Exact tissue match: +3
            # Same tissue family: +1.5
            if record["tissue_type_normalized"] == tissue_norm:
                score += 3.0
            elif self._same_tissue_family(record["tissue_type_normalized"], tissue_norm):
                score += 1.5

            # Exact marker match: +3
            # Same marker family: +2
            if record["marker_name_normalized"] == marker_norm:
                score += 3.0
            elif self._same_marker_family(record["marker_name_normalized"], marker_norm):
                score += 2.0

            # Quality improvement bonus
            score += min(2.0, record["quality_improvement"] * 2)

            # User approved bonus
            if record.get("user_approved"):
                score += 1.0

            # Image characteristics similarity
            if characteristics and record.get("characteristics"):
                char_record = ImageCharacteristics.from_dict(record["characteristics"])
                char_similarity = characteristics.similarity_score(char_record)
                score += char_similarity * 1.5

            # Recency bonus (decay over 30 days)
            try:
                record_time = datetime.fromisoformat(record["timestamp"])
                days_old = (datetime.now() - record_time).days
                recency_factor = max(0, 1 - days_old / 30)
                score += recency_factor * 0.5
            except Exception:
                pass

            # Apply source weight
            score *= record.get("weight", 1.0)

            record["score"] = score

        return records

    def _same_tissue_family(self, tissue1: str, tissue2: str) -> bool:
        """Check if two tissues are in the same family."""
        families = {
            "lymphoid": ["tonsil", "lymph_node", "spleen", "thymus"],
            "gi": ["colon", "small_intestine", "stomach", "esophagus"],
            "solid_organ": ["liver", "kidney", "pancreas"],
            "neural": ["brain"],
        }

        for _family, members in families.items():
            if tissue1 in members and tissue2 in members:
                return True
        return False

    def _same_marker_family(self, marker1: str, marker2: str) -> bool:
        """Check if two markers are in the same family."""
        families = {
            "t_cell": ["cd3", "cd4", "cd8", "foxp3"],
            "b_cell": ["cd20", "cd19"],
            "myeloid": ["cd68", "cd163", "cd11c", "cd14"],
            "structural": ["panck", "vimentin", "sma", "ecad"],
            "checkpoint": ["pd1", "pdl1"],
        }

        for _family, members in families.items():
            if marker1 in members and marker2 in members:
                return True
        return False

    def _aggregate_parameters(self, records: list[dict[str, Any]]) -> dict[str, Any]:
        """Aggregate parameters from multiple records with weighted voting."""
        if not records:
            return {}

        # Collect all parameter values with weights
        param_values: dict[str, list[tuple[Any, float]]] = {}

        for record in records:
            weight = record.get("score", 1.0)
            for key, value in record["parameters"].items():
                if key not in param_values:
                    param_values[key] = []
                param_values[key].append((value, weight))

        # Aggregate each parameter
        aggregated = {}
        for key, values_weights in param_values.items():
            values = [v for v, w in values_weights]
            weights = [w for v, w in values_weights]

            # Check type of first value
            first_value = values[0]

            if isinstance(first_value, bool):
                # Boolean: weighted majority vote
                true_weight = sum(w for v, w in values_weights if v)
                false_weight = sum(w for v, w in values_weights if not v)
                aggregated[key] = true_weight > false_weight

            elif isinstance(first_value, (int, float)):
                # Numeric: weighted average
                total_weight = sum(weights)
                if total_weight > 0:
                    weighted_avg = sum(v * w for v, w in values_weights) / total_weight
                    # Round to appropriate precision
                    if isinstance(first_value, int):
                        aggregated[key] = int(round(weighted_avg))
                    else:
                        aggregated[key] = round(weighted_avg, 3)
                else:
                    aggregated[key] = first_value

            else:
                # Other types: take most common weighted
                value_weights: dict[Any, float] = {}
                for v, w in values_weights:
                    v_key = str(v)  # Convert to string for dict key
                    value_weights[v_key] = value_weights.get(v_key, 0) + w
                best_value_str = max(value_weights, key=value_weights.get)
                # Find original value
                for v, _w in values_weights:
                    if str(v) == best_value_str:
                        aggregated[key] = v
                        break

        return aggregated

    def _calculate_confidence(
        self,
        records: list[dict[str, Any]],
        tissue_norm: str,
        marker_norm: str,
    ) -> float:
        """Calculate confidence score for the recommendation."""
        if not records:
            return 0.0

        factors = []

        # Number of supporting records
        n_records = len(records)
        factors.append(min(1.0, n_records / 5))  # Max confidence at 5+ records

        # Exact matches
        exact_matches = sum(
            1
            for r in records
            if r["tissue_type_normalized"] == tissue_norm
            and r["marker_name_normalized"] == marker_norm
        )
        factors.append(min(1.0, exact_matches / 3))

        # Average quality improvement
        avg_improvement = np.mean([r["quality_improvement"] for r in records])
        factors.append(min(1.0, avg_improvement * 2))

        # Consistency of parameters (low variance = high confidence)
        if len(records) > 1:
            param_keys = set()
            for r in records:
                param_keys.update(r["parameters"].keys())

            consistency_scores = []
            for key in param_keys:
                values = [r["parameters"].get(key) for r in records if key in r["parameters"]]
                if len(values) > 1 and all(isinstance(v, (int, float)) for v in values):
                    cv = np.std(values) / (np.mean(values) + 1e-6)  # Coefficient of variation
                    consistency_scores.append(max(0, 1 - cv))

            if consistency_scores:
                factors.append(np.mean(consistency_scores))

        return np.mean(factors) if factors else 0.0

    def _confidence_level(self, confidence: float) -> str:
        """Convert confidence score to human-readable level."""
        if confidence >= 0.8:
            return "high"
        elif confidence >= 0.6:
            return "medium"
        elif confidence >= 0.4:
            return "low"
        else:
            return "very_low"

    def get_statistics(self) -> dict[str, Any]:
        """Get statistics about the learning database."""
        stats = {
            "project_records": 0,
            "global_records": 0,
            "unique_tissues": set(),
            "unique_markers": set(),
            "operations": {},
        }

        for db_path, key in [(self.project_db_path, "project"), (self.global_db_path, "global")]:
            if db_path and db_path.exists():
                conn = sqlite3.connect(str(db_path))
                cursor = conn.cursor()

                cursor.execute("SELECT COUNT(*) FROM parameter_records")
                stats[f"{key}_records"] = cursor.fetchone()[0]

                cursor.execute("SELECT DISTINCT tissue_type_normalized FROM parameter_records")
                stats["unique_tissues"].update(r[0] for r in cursor.fetchall() if r[0])

                cursor.execute("SELECT DISTINCT marker_name_normalized FROM parameter_records")
                stats["unique_markers"].update(r[0] for r in cursor.fetchall() if r[0])

                cursor.execute("""
                    SELECT operation, COUNT(*), AVG(quality_improvement)
                    FROM parameter_records
                    GROUP BY operation
                """)
                for op, count, avg_improvement in cursor.fetchall():
                    if op not in stats["operations"]:
                        stats["operations"][op] = {"count": 0, "avg_improvement": 0}
                    stats["operations"][op]["count"] += count
                    stats["operations"][op]["avg_improvement"] = avg_improvement or 0

                conn.close()

        stats["unique_tissues"] = list(stats["unique_tissues"])
        stats["unique_markers"] = list(stats["unique_markers"])

        return stats

    def export_learned_parameters(self, output_path: Path | None = None) -> Path:
        """Export all learned parameters to a JSON file."""
        if output_path is None:
            output_path = (self.project_path or Path.cwd()) / "learned_parameters.json"

        all_records = []

        for db_path, source in [(self.project_db_path, "project"), (self.global_db_path, "global")]:
            if db_path and db_path.exists():
                conn = sqlite3.connect(str(db_path))
                conn.row_factory = sqlite3.Row
                cursor = conn.cursor()

                cursor.execute("SELECT * FROM parameter_records WHERE user_approved = 1")
                for row in cursor.fetchall():
                    all_records.append(
                        {
                            "source": source,
                            "tissue_type": row["tissue_type"],
                            "marker_name": row["marker_name"],
                            "operation": row["operation"],
                            "parameters": json.loads(row["parameters_json"]),
                            "quality_improvement": row["quality_improvement"],
                            "timestamp": row["timestamp"],
                        }
                    )

                conn.close()

        export_data = {
            "exported": datetime.now().isoformat(),
            "total_records": len(all_records),
            "records": all_records,
        }

        with open(output_path, "w") as f:
            json.dump(export_data, f, indent=2)

        return output_path
