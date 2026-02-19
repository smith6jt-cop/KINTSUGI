# Signal Isolation

Signal isolation removes autofluorescence (AF) from multiplex immunofluorescence images, separating true marker signal from tissue-intrinsic background. KINTSUGI provides two subtraction methods — **global** and **weighted multi-range** — plus integrated denoising, contrast enhancement, and quality assessment.

## What is Autofluorescence?

Biological tissues emit autofluorescence from endogenous fluorophores such as collagen, elastin, lipofuscin, and NADH. In multiplex IF experiments (e.g., CODEX), a **blank channel** — exposed to the same excitation light but with no antibody conjugate — captures this autofluorescence pattern. Subtracting the blank from each signal channel removes the AF contribution.

The challenge: a single global subtraction factor works well for bright markers (CD3, CD20) but can destroy dim markers (FOXP3, CD163, CD25) where AF intensity approaches or exceeds the true signal.

## Methods

### Global Subtraction

The original method applies a single scale factor to the entire blank channel:

```
result = signal - min(signal, blank * scale_factor)
```

The `min()` operation prevents over-subtraction — no pixel can go below zero.

**When to use:** Bright markers with good signal-to-AF ratio, quick exploratory analysis, or when weighted parameters haven't been established.

**Key parameters:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `blank_clip_factor` | Zero out blank pixels below this value (removes noise) | 0 |
| `blank_scale_factor` | Multiply blank before subtraction (>1 = more aggressive) | 1.0 |
| `smooth_low` | Apply uniform filter to dim regions | False |
| `smooth_high` | Apply uniform filter to bright regions | False |
| `erosion` | Erode signal mask to clean edges (disk radius) | 0 |

### Weighted Multi-Range Subtraction

The weighted method segments the signal histogram into intensity ranges (default: 5) and computes a per-range subtraction weight based on the signal-to-AF ratio:

```
weight_map = f(signal_intensity)     # per-pixel weight from range membership
result = signal - min(signal, blank * base_scale * weight_map)
```

Each intensity range gets a different weight:

| Range | AF vs Signal | Weight | Behavior |
|-------|-------------|--------|----------|
| **Background** | Both near zero | 0.0 | No subtraction |
| **Very dim** | AF dominates (ratio > 1.5) | 0.3-0.5 | Gentle — protect dim signal |
| **Dim** | Mixed (ratio 0.8-1.5) | 0.5-0.8 | Moderate subtraction |
| **Medium** | Signal moderate (ratio 0.3-0.8) | 0.8-1.0 | Near-full subtraction |
| **Bright** | Signal dominant (ratio < 0.3) | 1.0-1.15 | Aggressive AF removal |

Transitions between ranges use cosine blending to prevent discontinuities.

**When to use:**
- Dim markers where global subtraction destroys signal (FOXP3, CD163, CD25, CD11c)
- Tissues with heterogeneous AF (collagen-rich stroma adjacent to low-AF lymphoid regions)
- When different intensity regions need different treatment

**Key parameters:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `base_scale_factor` | Base scale applied before per-range weighting | 1.0 |
| `blank_clip_factor` | Zero out blank pixels below this value | 0 |
| `n_ranges` | Number of intensity ranges | 5 |
| `range_method` | Boundary computation: `"percentile"` or `"otsu"` | `"percentile"` |
| `transition_width` | Fraction of range width for cosine blending (0-0.5) | 0.1 |
| `ranges` | Pre-computed ranges (overrides `n_ranges`/`range_method`) | None |

All global parameters (`smooth_low`, `smooth_high`, `erosion`, etc.) are also available as post-processing options.

## Usage

### Claude-Guided (MCP) — Recommended

The Claude Code MCP integration provides an interactive, AI-assisted workflow. Claude analyzes the image, suggests parameters, applies subtraction, and records successful parameters for future use.

**Setup:**
```bash
pip install kintsugi[claude]
kintsugi mcp config /path/to/project  # If not using kintsugi init
```

**Workflow:**

1. **Preview range analysis** — Claude calls `analyze_weighted_subtraction()` to show the intensity ranges, per-range weights, and expected behavior before applying anything:

   ```
   User: "Load CD163 from cycle 3 and analyze for weighted blank subtraction"
   Claude: [Uses load_channel, then analyze_weighted_subtraction]
   Claude: "Analysis shows 5 ranges. Very_dim range (weight 0.35) protects
            the dim CD163+ cells. Bright range (weight 1.12) aggressively
            removes collagen AF."
   ```

2. **Apply weighted subtraction** — Once the user approves the preview, Claude applies the subtraction:

   ```
   User: "Apply those parameters"
   Claude: [Uses subtract_blank with method="weighted"]
   Claude: "Subtraction complete. Quality score: 0.82. SNR improved 1.4x."
   ```

3. **Approve and learn** — Recording the successful parameters improves future recommendations:

   ```
   User: "That looks good, approve it"
   Claude: [Uses approve_and_learn to record parameters]
   Claude: "Parameters recorded for CD163/tonsil. Future recommendations
            will use this as a baseline."
   ```

**Available MCP tools:**

| Tool | Purpose |
|------|---------|
| `load_channel` | Load a channel image from the project |
| `analyze_weighted_subtraction` | Preview per-range weights without applying |
| `subtract_blank` | Apply subtraction (`method="global"` or `"weighted"`) |
| `suggest_with_learning` | Get recommendations combining analysis + learned history |
| `approve_and_learn` | Record approved parameters for future learning |
| `assess_quality` | Run quality assessment on the result |
| `get_learned_parameters` | Query the learning database directly |

### Python API

#### AutofluorescenceSubtractor class

The high-level class manages parameter suggestion, subtraction, quality assessment, and learning in one interface:

```python
from kintsugi.signal.subtractor import AutofluorescenceSubtractor

# Global subtraction (default)
subtractor = AutofluorescenceSubtractor(
    project_dir="./my_project",
    tissue_type="tonsil",
)
result = subtractor.process(signal, blank, marker="CD3")
print(f"Quality: {result.quality_metrics['quality_score']:.3f}")

# Weighted subtraction
subtractor = AutofluorescenceSubtractor(
    project_dir="./my_project",
    tissue_type="tonsil",
    method="weighted",
)
result = subtractor.process(signal, blank, marker="CD163")
print(f"Quality: {result.quality_metrics['quality_score']:.3f}")
print(f"Per-range metrics: {result.range_metrics}")
```

The `method` parameter can also be overridden per-call:

```python
subtractor = AutofluorescenceSubtractor(project_dir="./my_project")

# Global for bright markers
result_cd3 = subtractor.process(signal_cd3, blank, marker="CD3", method="global")

# Weighted for dim markers
result_foxp3 = subtractor.process(signal_foxp3, blank, marker="FOXP3", method="weighted")
```

#### Direct functions

For lower-level control, use the functions directly:

```python
from kintsugi.signal.autofluorescence import (
    subtract_autofluorescence,           # Global method
    subtract_autofluorescence_weighted,   # Weighted method
    analyze_for_weighted_subtraction,     # Preview analysis
    compute_intensity_ranges,             # Compute ranges
    build_weight_map,                     # Build per-pixel weights
)

# Analyze first
analysis = analyze_for_weighted_subtraction(
    signal, blank, tissue_type="tonsil", marker_name="CD163"
)
print(f"Ranges: {len(analysis['ranges'])}")
for r in analysis["ranges"]:
    print(f"  {r['label']}: weight={r['weight']:.2f}, pixels={r['pixel_fraction']:.1%}")

# Apply
result = subtract_autofluorescence_weighted(
    signal, blank,
    blank_clip_factor=analysis["blank_clip_factor"],
    base_scale_factor=analysis["base_scale_factor"],
    ranges=analysis["ranges"],
    transition_width=0.1,
)
```

#### Batch processing

Process multiple channels with the same subtractor instance:

```python
subtractor = AutofluorescenceSubtractor(
    project_dir="./my_project",
    tissue_type="tonsil",
    method="weighted",
)

results = subtractor.process_batch(
    channels={"CD3": cd3_img, "CD163": cd163_img, "FOXP3": foxp3_img},
    blank_channels={"Blank1b": blank_img},
)

for marker, result in results.items():
    print(f"{marker}: quality={result.quality_metrics['quality_score']:.3f}")
```

### Batch Bootstrap (CLI)

Pre-populate the learning database from existing processed data without interactive sessions:

```bash
# Scan EDF outputs for signal/blank pairs and record weighted parameters
kintsugi mcp pretrain /path/to/project --tissue-type tonsil

# Preview what would be processed
kintsugi mcp pretrain /path/to/project --tissue-type tonsil --dry-run
```

The bootstrap scanner:
1. Finds EDF outputs in `data/processed/edf/cyc##/`
2. Identifies signal/blank pairs per cycle (skips DAPI)
3. Computes weighted subtraction parameters for each pair
4. Records to the learning database with `algorithm_version="weighted_v1"`
5. Marks as `user_approved=False` (auto-computed, not human-verified)

This is useful when setting up a new project — run bootstrap first, then use Claude-guided or Python API for fine-tuning. The learning system weights human-approved parameters higher than bootstrapped ones.

## Parameter Reference

### Global Subtraction Parameters

| Parameter | Type | Range | Default | Description |
|-----------|------|-------|---------|-------------|
| `blank_clip_factor` | int | 0-65535 | 0 | Zero out blank pixels below this threshold |
| `blank_scale_factor` | float | 0.1-3.0 | 1.0 | Scale blank intensity before subtraction |
| `smooth_low` | bool | — | False | Smooth dim regions with uniform filter |
| `low_size` | int | 1-10 | 2 | Uniform filter kernel size for dim regions |
| `low_percentile` | int | 0-100 | 60 | Percentile threshold defining "dim" |
| `smooth_high` | bool | — | False | Smooth bright regions with uniform filter |
| `high_size` | int | 1-10 | 2 | Uniform filter kernel size for bright regions |
| `high_percentile` | int | 0-100 | 90 | Percentile threshold defining "bright" |
| `erosion` | int | 0-10 | 0 | Morphological erosion disk radius |

### Weighted Subtraction Parameters

All global parameters above, plus:

| Parameter | Type | Range | Default | Description |
|-----------|------|-------|---------|-------------|
| `base_scale_factor` | float | 0.1-3.0 | 1.0 | Base scale applied before per-range weighting |
| `n_ranges` | int | 3-10 | 5 | Number of intensity ranges |
| `range_method` | str | — | `"percentile"` | Boundary method: `"percentile"` or `"otsu"` |
| `transition_width` | float | 0.0-0.5 | 0.1 | Fraction of range width for cosine blending |
| `ranges` | list | — | None | Pre-computed ranges (overrides auto-computation) |

### Quality Metrics

Both methods compute quality metrics after subtraction:

| Metric | Description | Good Range |
|--------|-------------|------------|
| `quality_score` | Overall quality (0-1) | > 0.5 |
| `snr_improvement` | SNR change relative to original | > 0 |
| `af_removal` | Fraction of AF correlation removed | > 0.5 |
| `signal_preservation` | Fraction of non-zero pixels retained | > 0.3 |
| `residual_correlation` | Remaining correlation with blank | < 0.3 |

The weighted method additionally reports per-range metrics (SNR, AF removal, signal preservation per intensity range).

## Parameter Learning

KINTSUGI's parameter learning system builds a database of successful parameters indexed by tissue type, marker name, and algorithm version. Over time, this enables increasingly accurate recommendations.

### How it works

1. **Record**: After the user approves a subtraction result (via MCP `approve_and_learn` or the Python API with `auto_learn=True`), the parameters and quality score are stored in a SQLite database.

2. **Recommend**: When processing a new channel, the system queries the database for similar tissue/marker combinations. If matches exist with sufficient confidence, it recommends those parameters as a starting point.

3. **Merge**: When both image analysis and learned parameters are available, the system merges them — weighting learned parameters more heavily when confidence is high.

### Database location

The learning database is stored per-project at:
```
<project_dir>/data/.kintsugi_learning.db
```

### Checking learning statistics

```python
from kintsugi.mcp.tools.learning import ParameterLearningEngine

engine = ParameterLearningEngine(str(project_dir))
stats = engine.get_statistics()
print(f"Total records: {stats['total_records']}")
print(f"Tissue types: {stats['tissue_types']}")
print(f"Markers: {stats['markers']}")
```

Or via the MCP tool:
```
User: "Show me the learning database statistics"
Claude: [Uses get_learning_statistics]
```

### Backwards compatibility

The weighted method is fully backwards compatible:
- `subtract_blank()` defaults to `method="global"` (original behavior)
- Uniform weights of 1.0 produce identical output to the global method
- Existing parameter databases continue to work; weighted parameters are stored with `algorithm_version="weighted_v1"`
