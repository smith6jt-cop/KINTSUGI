# Signal Processing — CLAUDE.md

See also: Root `CLAUDE.md` for project overview, `../CLAUDE.md` for KRONOS/vessel3d/cleanup, `../../workflow/CLAUDE.md` for Snakemake/registration.

## Weighted Autofluorescence Subtraction (Feb 2026)

Replaces the single global `blank_scale_factor` with per-intensity-range weights. Protects dim signal (FOXP3, CD163) while aggressively removing bright autofluorescence (collagen, lipofuscin).

**Algorithm**: Segment signal histogram into N ranges (default 5: background, very_dim, dim, medium, bright). Per-range weight is computed from signal-to-AF ratio — dim regions where AF > signal get weight 0.3-0.5 (gentle subtraction), bright signal regions get weight 1.0+ (aggressive removal). Cosine transitions between ranges prevent discontinuities.

**Key files:**
- `autofluorescence.py` — `subtract_autofluorescence_weighted()`, `compute_intensity_ranges()`, `build_weight_map()`, `analyze_for_weighted_subtraction()`, `compute_weighted_subtraction_quality()`, Dask variant
- `utils.py` — `adaptive_range_boundaries()`, `smooth_membership()`
- `subtractor.py` — `IntensityRange`, `WeightedSubtractionParameters` dataclasses, `method="weighted"` in `AutofluorescenceSubtractor`
- `bootstrap.py` — `bootstrap_from_project()` for batch pretraining
- `../claude/parameter_learning.py` — `algorithm_version` column, range aggregation
- `../mcp/tools/signal_isolation.py` — `analyze_weighted_subtraction()` tool, `method` param on `subtract_blank()`

**CLI**: `kintsugi mcp pretrain <project_dir> [--tissue-type TYPE] [--dry-run]` — scans EDF outputs for signal/blank pairs, computes weighted parameters, records to learning database with `algorithm_version="weighted_v1"`.

**Backwards compatible**: `subtract_blank()` defaults to `method="global"` (original behavior). Uniform weights=1.0 produces identical output to global method.

## Batch Signal Isolation (Feb 2026)

Recipe-driven and auto-analyzed batch signal isolation with multi-step subtraction, background cleaning, parameter learning, and self-normalized QC.

**Two processing paths:**
1. **Recipe path** (`--recipe-dir`): Loads legacy Notebook 3 `*_param.txt` files → primary subtraction → optional second subtraction → optional background cleaning. ~30 sec/channel.
2. **Auto-analysis path** (default): `analyze_for_subtraction()` suggests parameters, `select_method()` chooses global vs weighted. ~2 min/channel with sigma=0.

**CLI usage:**
```bash
kintsugi workflow isolate plan .                                          # Auto-analysis dry-run
kintsugi workflow isolate plan . --recipe-dir .../Processing_parameters   # Recipe preview
kintsugi workflow isolate run . --recipe-dir .../Processing_parameters --tissue-type spleen  # Recipe processing
kintsugi workflow isolate run . --method weighted                         # Force weighted for all
kintsugi workflow isolate run . --channels CD45,CD1c --force              # Specific markers
kintsugi workflow isolate run . --no-learn                                # Skip learning DB recording
kintsugi workflow isolate qc .                                            # Generate QC pages
kintsugi workflow isolate status .                                        # Summary table from manifest
```

CLI options: `--method {auto,global,weighted}`, `--tissue-type`, `--output-dir`, `--force`, `--channels`, `--tile-smooth-sigma` (default 0), `--recipe-dir`, `--learn/--no-learn`.

**Key files:**
- `batch.py` — `SubtractionParams`, `CleanParams`, `MarkerRecipe` dataclasses; `load_legacy_recipes()`, `clean_background()`, blank resolution helpers; `discover_channels()`, `select_method()`, `process_channel()`, `process_batch()`
- `isolation_qc.py` — `generate_qc_pages()` (3-column self-normalized layout), `generate_summary_table()`
- `../../cli.py` — `@workflow.group("isolate")` with `plan`, `run`, `qc`, `status` subcommands
- `../claude/parameter_learning.py` — `ParameterLearningEngine` for recording outcomes
- `../../../tests/test_batch_signal_isolation.py` — 66 tests

### Recipe System

**`load_legacy_recipes(recipe_dir)`** parses Notebook 3 `*_param.txt` files into `MarkerRecipe` objects:
- `SubtractionParams`: blank_name, blank_clip_factor, blank_scale_factor, smooth_low/high, erosion
- `CleanParams`: background_threshold, smooth, smooth_threshold, footprint, remove_small, small_size
- `MarkerRecipe`: marker_name, primary (required), second (optional), clean (optional)
- Skips `dask.array<...>` and `datetime.datetime(...)` artifact lines via `_parse_param_value()`

**Blank name resolution** (`_resolve_blank_path()`):
1. Normalize: `Blank1b` → `Blank_1b`, `Blank13c` → `Blank_13c` (regex)
2. Exact match in location map (all registered channel names → paths)
3. Fuzzy match: strip hyphens/underscores, case-insensitive (`HLADR` → `HLA-DR`)

**`clean_background(image, params)`** — pure numpy reimplementation of `Kutils.clean()`:
1. Zero pixels below `background_threshold`
2. Median filter in transition zone (between 0 and `smooth_threshold`)
3. `remove_small_objects` + morphological closing

### Method Selection (auto-analysis path)

`select_method()` decision logic:
- `blank_p99 / signal_p99 > 1.2` → weighted (blank dominates, protects dim signal)
- `af_contribution > 0.3 AND correlation > 0.4` → weighted (structured AF)
- `correlation > 0.5 AND dynamic_range > 5000` → weighted (correlated AF)
- Otherwise → global

### Parameter Learning Integration

When `--learn` (default) and `--tissue-type` are provided:
- Records recipe outcomes to `ParameterLearningEngine` after successful processing
- Two operations per channel: `blank_subtraction` and `clean_background`
- `algorithm_version="recipe_v1"` distinguishes from auto-analyzed params
- Cross-dataset transfer: tuned params from CX_19-001 spleen propagate to other spleen datasets

### Tile Smoothing

`smooth_blank_for_subtraction(blank, sigma)` — large-sigma Gaussian blur on blank (not signal) before subtraction. Removes tile-to-tile intensity variation from BaSiC correction while preserving gradual AF spatial gradient. Default sigma changed from 500 to **0** (recipes don't need smoothing). Signal is NOT smoothed.

### QC

Three-column layout per channel — Before (self-normalized p1-p99, gray), After (self-normalized p1-p99, gray), Difference (inferno). Self-normalization ensures dim cellular signal is visible. Pages of 6 channels.

**Manifest**: `signal_isolation_manifest.json` records per-channel method, parameters, quality metrics, analysis, and warnings. Includes `recipe_dir` when recipes are used.

**Channel discovery**: `discover_channels()` parses `workflow/config.yaml` channel_names. With recipes, uses recipe blank names via resolution chain; without recipes, uses positional mapping (CH2→Blank_1a, CH3→Blank_1b, CH4→Blank_1c).

**Key gotcha**: `compute_weighted_subtraction_quality()` returns `{"global": {...}, "per_range": [...]}` — must extract `["global"]` for flat quality_score access.

### Validated Results

CX_19-001_SP_CC2-A28: 28 channels (24 recipe + 4 auto-weighted), mean quality 0.803, 0 errors, ~30 sec/channel with recipes.
