# Core Package — CLAUDE.md

Feature-specific documentation for `src/kintsugi/` modules. See also:
- Root `CLAUDE.md` — project overview, architecture, build commands, testing
- `signal/CLAUDE.md` — weighted AF subtraction, batch signal isolation
- `../../workflow/CLAUDE.md` — Snakemake workflow, registration, batch processing

## KRONOS Foundation Model Integration (Feb 2026)

[KRONOS](https://github.com/mahmoodlab/KRONOS) is a panel-agnostic foundation model for spatial proteomics (Harvard/BWH, CC-BY-NC-ND-4.0). Trained via self-supervised learning on 47M single-marker patches across 175 protein markers.

**Pipeline position**: Post-registration analysis — runs after EDF/registration produces aligned, marker-named 2D TIFFs.

**Key files:**
- `kronos/__init__.py` — Lazy loading, dependency checks
- `kronos/model.py` — `KronosModel`, `KronosConfig`, `EmbeddingResult` (HDF5 save/load)
- `kronos/markers.py` — `MarkerMapper` maps CHANNELNAMES.txt to KRONOS's 175 markers
- `kronos/embeddings.py` — `KronosEmbedder` loads registered images, tiles, normalizes, runs inference
- `kronos/analysis.py` — Clustering (Leiden/KMeans), UMAP/PCA, spatial search, cross-dataset comparison
- `../../notebooks/5_KRONOS_Analysis.ipynb` — Standalone analysis notebook

**Usage:**
```python
from kintsugi.kronos import KronosEmbedder, MarkerMapper
mapper = MarkerMapper.from_project("./my_project", cache_dir="./model_assets")
embedder = KronosEmbedder()
result = embedder.embed_project("./my_project", mapper)
result.save("./embeddings.h5")
```

**Prerequisites**: `pip install -e KRONOS` (from GitHub clone), HuggingFace access for model weights.

**Model specs**: ViT-S/16, 175 MB weights, input `(batch, markers, 256, 256)`, output 3 embedding types (patch/marker/token). Per-marker z-score normalization via `marker_metadata.csv`.

**KRONOS does NOT replace**: stitching, deconvolution, EDF, registration, segmentation masks, or signal isolation — it consumes their outputs.

## 3D Vessel Segmentation (Feb 2026)

Segments blood and lymphatic vessels from deconvolved 3D z-stacks using Frangi vesselness filtering, morphological cleanup, skeletonization, and graph-based morphometry.

**Pipeline position**: Step 3.5 — after deconvolution, before EDF. Uses the full 3D z-stack (EDF collapses z). Runs in parallel with EDF; does not block downstream steps.

**Key files:**
- `vessel3d.py` — Core processing: `compute_vesselness_frangi()`, `binarize_vessel_mask()`, `skeletonize_vessels()`, `prune_skeleton()`, `analyze_vessel_graph()`, `export_vessel_results()`, `segment_vessels_3d()` (high-level pipeline)
- `vessel3d_viz.py` — Visualization: `ortho_view()`, `overlay_mask_on_raw()`, `render_skeleton_mip()`, `plot_vessel_features()`
- `../../notebooks/2.5_Vessel_3D_Segmentation.ipynb` — Interactive notebook (7 sections with verification checkpoints)
- `../../slurm/jobs/03b_vessel3d.sh` — HPC batch job script (Frangi + skeleton, automated)

**Data classes:**
- `VesselSpacing(xy, z)` — Physical voxel spacing; `from_experiment()` loads from experiment.json
- `VesselSegmentationResult` — Container for mask, skeleton, features, vesselness map

**Usage:**
```python
from kintsugi.vessel3d import segment_vessels_3d, VesselSpacing, export_vessel_results

spacing = VesselSpacing(xy=0.377, z=1.5)  # or VesselSpacing.from_experiment(config)
result = segment_vessels_3d(
    volume, spacing=spacing, sigmas=[1, 2, 4, 8],
    device='auto', marker_name='CD31',
)
export_vessel_results(result, "data/processed/vessel_3d/")
```

**Pipeline steps:**
1. Preprocess: isotropic z-resampling (scipy/CuPy zoom) + light Gaussian denoising
2. Frangi vesselness: multi-scale Hessian eigenvalue analysis (GPU-accelerated, Frangi 1998)
3. Binarize: Otsu threshold + `remove_small_objects` + `closing` + hole filling
4. Skeletonize: Lee (1994) 3D thinning via `skimage.morphology.skeletonize`
5. Graph analysis: `skan.Skeleton` + `skan.summarize` + distance transform for radius

**Output files** (in `data/processed/vessel_3d/`):
- `binary_mask_{marker}.tif` — 3D vessel mask (isotropic)
- `binary_mask_{marker}_native.tif` — 3D mask at native resolution
- `skeleton_{marker}.tif` — 3D skeleton
- `vessel_features_{marker}.csv` — Per-segment morphometry (length, diameter, tortuosity)
- `vessel_graph_{marker}.graphml` — NetworkX graph
- `vessel_2d_projection_{marker}.tif` — MIP for Notebook 4 spatial analysis

**Dependencies**: `skan>=0.11.0` and `networkx>=3.0` in the `analysis` optional group. GPU acceleration via CuPy (existing `gpu` group). Skeletonization is CPU-only.

**Bug fixes (Feb 2026):**

| Bug | Impact | Fix |
|-----|--------|-----|
| Ra = \|λ₁\|/\|λ₂\| (wrong) | Tube response suppressed — small/faint vessels missed | Ra = \|λ₂\|/\|λ₃\| per Frangi 1998 eq. 11 |
| `binary_closing` deprecated | FutureWarning in skimage 0.26, removed in 0.28 | Use `closing` from `skimage.morphology` |
| `remove_small_objects(min_size=)` deprecated | FutureWarning in skimage 0.26 | Use `max_size=min_size-1` parameter |
| L4 GPU OOM on isotropic volumes | `gaussian_filter` crashes on 23 GB VRAM GPUs | VRAM guard: estimate 15x volume size, CPU fallback |

**GPU memory**: Hessian eigenvalues use Cardano's analytical formula for symmetric 3x3 matrices — operates element-wise on 6 Hessian arrays. The naive `(N, 3, 3)` + `eigvalsh` approach OOMs on stitched volumes (108 GB allocation for 3 billion voxels). Sorting uses a 3-element sorting network (compare-swap) to avoid index array allocation. VRAM guard queries `cp.cuda.Device().mem_info` before GPU path and falls back to CPU if insufficient.

**SLURM**: `03b_vessel3d.sh` runs between `03_deconvolution.sh` and `04_edf.sh`. Environment variables: `VESSEL_CYCLE`, `VESSEL_CHANNEL`, `VESSEL_MARKER`, `VESSEL_SIGMAS`, `VESSEL_MIN_SIZE`. Recommended: 256 GB RAM, 4 h wall time. Auto-detects GPU VRAM < 40 GB and forces CPU mode (prevents L4 OOM).

**Tests**: `tests/test_vessel3d.py` — 19 tests covering Ra formula correctness (tube vs plate eigenvalues), eigenvalue sorting invariants, scikit-image deprecation-free imports, GPU VRAM guard, binarization, end-to-end pipeline, and VesselSpacing.

## Pipeline-Aware Cleanup (Feb 2026)

Safe deletion of intermediate data (`stitched/`, `deconvolved/`, `raw/`) with dependency graph validation, QC gating, and staged deletion with recovery.

**Root cause**: 14 datasets lost deconvolved z-stacks when the old `cleanup_datasets.sh` ran before vessel3d completed — it only checked EDF, missing that `deconvolved/` has two consumers (edf + vessel3d).

**Key files:**
- `cleanup.py` — Core module: dependency graph, assessment, trash staging, recovery
- `../../tests/test_cleanup.py` — 42 tests covering all safety scenarios
- `cli.py` — CLI commands under `@workflow.group("cleanup")`
- `../../workflow/Snakefile` — `cleanup_safe` rule aggregates all consumer sentinels

**CLI usage:**
```bash
kintsugi workflow cleanup status .                    # Show what's safe/blocked and why
kintsugi workflow cleanup plan .                      # Dry-run: what would be deleted
kintsugi workflow cleanup execute .                   # Delete safe intermediates (trash mode)
kintsugi workflow cleanup execute . --no-trash        # Permanent deletion
kintsugi workflow cleanup execute . --skip-vessel3d   # Override vessel3d blocking (emergency)
kintsugi workflow cleanup execute . --force           # Skip interactive confirmation
kintsugi workflow cleanup recover .                   # List trash entries
kintsugi workflow cleanup recover . --entry NAME      # Restore from trash
kintsugi workflow cleanup purge .                     # Delete trash entries >7 days old
kintsugi workflow cleanup purge . --all               # Delete all trash entries
```

**Safety mechanisms:**
1. **Dependency graph** — `deconvolved/` requires BOTH edf AND vessel3d (if enabled) to complete
2. **Conservative default** — if `optional_stages` absent from config, blocks deletion with guidance
3. **QC gate** — all cleanup blocked until stitch/decon/edf QC sentinels present
4. **Staged deletion** — `data/.trash/{dir}_{timestamp}/` with JSON receipt for recovery
5. **Per-cycle verification** — vessel3d can require only specific cycles (e.g., `cycles: [2, 3]`)

**Config integration** (`workflow/config.yaml`):
```yaml
optional_stages:
  vessel3d:
    enabled: false    # Set true if vessel3d is planned for this dataset
    cycles: []        # Empty = all cycles; or specify [2, 3] for subset
  spillover:
    enabled: false    # Placeholder for future spillover correction stage
```

**Decision logic:**
- If `optional_stages` section missing → BLOCK deconvolved deletion (conservative)
- If `vessel3d.enabled: false` → deconvolved SAFE (only edf needs it)
- If `vessel3d.enabled: true` + all sentinels present → deconvolved SAFE
- If `vessel3d.enabled: true` + sentinels missing → deconvolved BLOCKED
- If any QC sentinel missing → EVERYTHING BLOCKED

**Snakemake integration:** `snakemake cleanup_safe --profile profiles/slurm` verifies safety via DAG (preferred method is the Python CLI).

**Batch processing:** `cleanup_datasets.sh` is now a thin wrapper calling `kintsugi workflow cleanup` for each staged dataset.
