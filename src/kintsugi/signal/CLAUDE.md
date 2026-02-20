# Signal Processing — CLAUDE.md

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

Automated per-marker parameter optimization and batch autofluorescence subtraction. Replaces the naive `blank_scale_factor=1.0` script with intelligent method selection, blank smoothing, and QC reporting.

**CLI usage:**
```bash
kintsugi workflow isolate plan .                    # Dry-run: per-channel parameter preview
kintsugi workflow isolate run .                     # Process all channels
kintsugi workflow isolate run . --method weighted   # Force weighted method for all
kintsugi workflow isolate run . --channels CD45,CD1c  # Process specific markers
kintsugi workflow isolate run . --force             # Overwrite existing outputs
kintsugi workflow isolate qc .                      # Generate QC pages
kintsugi workflow isolate status .                  # Summary table from manifest
```

CLI options: `--method {auto,global,weighted}`, `--tissue-type`, `--output-dir`, `--force`, `--channels`, `--tile-smooth-sigma` (default 500, 0 to disable).

**Key files:**
- `batch.py` — `ChannelSpec`, `ChannelResult`, `BatchResult` dataclasses; `discover_channels()`, `select_method()`, `process_channel()`, `process_batch()`, `smooth_blank_for_subtraction()`
- `isolation_qc.py` — `generate_qc_pages()` (3-column self-normalized layout), `generate_summary_table()` (Rich-formatted)
- `../../cli.py` — `@workflow.group("isolate")` with `plan`, `run`, `qc`, `status` subcommands
- `../../../tests/test_batch_signal_isolation.py` — 29 tests

**Method selection algorithm** (`select_method()`):
- `blank_p99 / signal_p99 > 1.2` → weighted (blank dominates, protects dim signal)
- `af_contribution > 0.3 AND correlation > 0.4` → weighted (structured AF)
- `correlation > 0.5 AND dynamic_range > 5000` → weighted (correlated AF)
- Otherwise → global

**Tile compensation — blank smoothing**: Stitched images are mosaics of 63 tiles (7x9 grid, 30% overlap). Blank (cycle 1) and signal (cycle N) have different tile intensity patterns from BaSiC correction. Pixel-wise subtraction creates a visible grid artifact. `smooth_blank_for_subtraction(blank, sigma=500)` applies large-sigma Gaussian blur to the blank before subtraction, removing tile-to-tile variation while preserving the gradual AF spatial gradient. Signal is NOT smoothed.

**QC pages**: Three-column layout per channel — Before (self-normalized p1-p99, gray), After (self-normalized p1-p99, gray), Difference (inferno). Self-normalization ensures dim cellular signal is visible instead of appearing all-black. Pages of 6 channels, output to `qc_plots/signal_isolation_qc_p{N}.png`.

**Manifest**: `data/processed/signal_isolated/signal_isolation_manifest.json` records per-channel method, parameters, quality metrics, blank/signal ratio, correlation, and warnings. Used by QC and status commands.

**Channel discovery**: Parses `workflow/config.yaml` channel_names dict. Maps channel position to blank: CH2→Blank_1a, CH3→Blank_1b, CH4→Blank_1c. Skips DAPI, Blank, Empty channels.

**Key gotcha**: `compute_weighted_subtraction_quality()` returns `{"global": {...}, "per_range": [...]}` — must extract `["global"]` for flat quality_score access.
