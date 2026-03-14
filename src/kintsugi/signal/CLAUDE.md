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
- `../../../tests/test_batch_signal_isolation.py` — 94 tests

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
4. Positional fallback: extract suffix (a/b/c) from blank name, try `Blank_1{suffix}` — handles recipe blanks from higher cycles (e.g., `Blank_13b` → `Blank_1b`) when only cycle 1 blanks exist

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

Three-column layout per channel — Before (normalized p1-p99, gray), After (normalized p1-p99, gray), Difference (inferno). Two normalization modes:
- `independent` (default): each image uses its own p1-p99 range — ensures dim signal is visible
- `matched`: both Before and After use Before's p1-p99 range — prevents tile-boundary artifacts (1-5% BaSiC stitching residuals) from being amplified to full contrast. CLI: `--normalize-mode matched`

**Manifest**: `signal_isolation_manifest.json` records per-channel method, parameters, quality metrics, analysis, and warnings. Includes `recipe_dir` when recipes are used.

**Channel discovery**: `discover_channels()` parses `workflow/config.yaml` channel_names. With recipes, uses recipe blank names via resolution chain; without recipes, uses positional mapping (CH2→Blank_1a, CH3→Blank_1b, CH4→Blank_1c).

**Key gotcha**: `compute_weighted_subtraction_quality()` returns `{"global": {...}, "per_range": [...]}` — must extract `["global"]` for flat quality_score access.

### Safety Mechanisms (Mar 2026)

**Self-referential secondary blank validation** (`_validate_secondary_blank()`): Prevents a marker from using itself as its own secondary blank (e.g., CD20 recipe lists CD20 as `blankID2`). Case-insensitive comparison. Validation runs before AND after name resolution (resolved name may differ from recipe name).

**Duplicate blank validation** (`_validate_duplicate_blank()`): Prevents using the same blank for both primary and secondary subtraction, which doubles AF removal and destroys signal. Checks both raw recipe names (e.g., `Blank_1b` == `Blank_1b`) and resolved file names (e.g., `Blank_13b` resolves to `Blank_1b` via positional fallback). When duplicate detected, second subtraction is skipped with warning `"secondary_blank_duplicate_blank"` or `"secondary_blank_duplicate_blank_resolved"`.

**Post-subtraction quality gate** (`_check_over_subtraction()`): Detects over-subtracted images:
- `zero_percent >= 95%` → always fails (image essentially destroyed)
- `quality_score < quality_gate AND zero_percent > 70%` → fails (poor quality + significant loss)
- `quality_gate <= 0` disables the check entirely

When a recipe produces over-subtracted results, `_auto_fallback()` attempts auto-analysis without the recipe. If the fallback improves both zero% AND quality score, it replaces the recipe result (`recipe_source="auto_fallback"`). **The fallback result is re-validated against the quality gate** — if still over-subtracted, `status="failed"` and warning includes `"auto_fallback_also_failed"`. Otherwise the result is saved with `status="failed"`.

Config: `signal_isolation.quality_gate` (default 0.6) in `workflow/config.yaml`.

**Outlier clipping** (`clip_outliers()`): Hot pixels survive subtraction unchanged (signal=55000, blank=5000 → result=50000). `clip_outliers(image, percentile=99.5)` clips values above the given percentile of the **non-zero** pixel distribution. Applied after all subtractions + background cleaning, before quality metrics. Config: `signal_isolation.clip_percentile` (default 99.5), set to 0 to disable.

### Validated Results

CX_19-001_SP_CC2-A28: 28 channels (24 recipe + 4 auto-weighted), mean quality 0.803, 0 errors, ~30 sec/channel with recipes.

## Multi-Project Batch Signal Isolation (Feb 2026)

Orchestrates autofluorescence subtraction across multiple KINTSUGI projects with cascading recipe resolution.

**CLI:**
```bash
kintsugi workflow isolate batch /path/to/KINTSUGI_Projects --dry-run          # Preview all projects
kintsugi workflow isolate batch . -d CX_19-003_spleen_CC1-A                   # Single project
kintsugi workflow isolate batch . --learn --method auto                        # Full batch
kintsugi workflow isolate batch . -f --template-recipe-dir /path/to/recipes   # Force + custom template
```

**Cascading recipe resolution** (per marker, in priority order):

| Tier | Source | How |
|------|--------|-----|
| 1 | Own recipes | `{project}/data/processed/Processing_parameters/` etc. |
| 2 | Learned DB | `ParameterLearningEngine.recommend_parameters()` with confidence >= 0.6 |
| 3 | Template | Default: CX_19-001's 26 validated recipes (match by marker name) |
| 4 | Auto | `select_method()` picks global vs weighted per marker |

**Key files:**
- `batch_multi.py` — `discover_projects()`, `parse_tissue_type()`, `resolve_recipes_for_project()`, `learned_params_to_recipe()`, `process_all_projects()`
- `../../cli.py` — `@isolate_group.command("batch")` with `--dry-run`, `--dataset`, `--force`, `--template-recipe-dir`, `--min-confidence`
- `../../../tests/test_batch_multi_project.py` — 34 tests

**Tissue type parsing** from project name (regex, case-insensitive):
- `_SP_` or `spleen` → `"spleen"`, `_LN_` or `lymph-node` → `"lymph_node"`, `_TH_` or `thymus` → `"thymus"`
- Falls back to `experiment.json` name field, then `"unknown"`

**Output:** Per-project `signal_isolation_manifest.json` + cross-project `batch_isolation_report_{timestamp}.json` in projects dir. Sequential processing with error isolation per project.
