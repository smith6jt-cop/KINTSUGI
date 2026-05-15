# Notebooks CLAUDE.md

Context for Claude Code when working in the `notebooks/` directory. See root `CLAUDE.md` for project overview, CLI usage, and architecture.

## Jupyter Notebook Workflow

**CRITICAL: Autoreload is enabled in all KINTSUGI notebooks:**

```python
%load_ext autoreload
%autoreload 2
```

This means:
- **NEVER tell the user to restart the kernel** after editing Python modules
- **NEVER tell the user to reload/reopen the notebook**
- After editing `.py` files (Kio.py, kintsugi/*.py, etc.), changes take effect automatically on next cell execution
- Just re-run the relevant cell - no restart needed

Only suggest kernel restart for: new package installation, environment variable changes, or C extension recompilation.

### Troubleshooting "Function Not Found" Errors

When users report `NameError: name 'function_name' is not defined`:

1. **Check if function is defined in a notebook cell** (not a module):
   - Some functions like `run_deconvolution` are defined inside the notebook, not imported from modules
   - These require running the definition cell before the cell that uses them

2. **Verify cell execution order**:
   ```python
   # Use this to analyze notebook cell order
   import json
   with open('notebook.ipynb') as f:
       nb = json.load(f)
   for i, cell in enumerate(nb['cells']):
       source = ''.join(cell.get('source', []))
       if 'function_name' in source:
           is_def = 'def function_name' in source
           print(f"Cell {i}: {'DEFINES' if is_def else 'CALLS'} function_name")
   ```

3. **Common patterns in 2_Cycle_Processing.ipynb**:
   | Function | Defined In | Called In |
   |----------|------------|-----------|
   | `run_deconvolution` | Cell 24 | Cell 25 |
   | `process_edf_tiff` | Cell 31 | Cell 32 |
   | `visualize_deconvolution` | Cell 27 | (manual) |

4. **Solution**: Run cells sequentially from the top, or at minimum run the definition cell before the calling cell.

## Shared Notebook Utilities (`Kio.py`)

Canonical home for helpers used by more than one notebook. Add functions here instead of duplicating in notebook cells.

| Helper | Purpose | Used by |
|--------|---------|---------|
| `load_channel_names()` | Parse `CHANNELNAMES.txt` (multiple formats) | Notebooks 1, 2, 3 |
| `make_channel_names_unique()` | Disambiguate duplicate Blank/Empty names | Notebooks 1, 2 |
| `save_channel_names_template()` | Write a `CHANNELNAMES.txt` starter | Notebook 1 |
| `extract_channels_from_ome()` | OME-TIFF → per-channel TIFFs | Notebook 4 |
| `normalize_cycle_dirs()` | Rename long-form CODEX cycle folders (e.g. `cyc004_reg001_211206_201615` → `cyc004`) | Notebooks 1 cell 5, 2 cell 6 |
| `load_tiles_parallel()` / `load_stack_parallel()` | Threaded tile/z-stack loading | Notebook 2 |
| `ProgressCounter`, `ProcessingConfig`, `zarr_write_lock` | Shared dataclasses | Notebook 2 |

### Cycle Directory Normalization

Raw CODEX folders can arrive as `cyc004_reg001_211206_201615`, but Notebook 1 cell 14 and Notebook 2 cells 13/16/42 all hard-code short-form paths like `f"cyc{cycle:03d}"`. The `normalize_cycle_dirs(raw_dir, width=3, dry_run=False)` helper renames in place, is idempotent, and raises `RuntimeError` on collisions (never overwrites).

This was previously a cell in Notebook 2 that got dropped. Both notebooks now call it automatically before `find_raw_cycles()`. The SLURM/Snakemake path uses a different mechanism (`workflow/scripts/stitch.py:find_cycle_dir()` resolves long-form at glob time) and is unaffected.

## Notebook 4: Segmentation & Spatial Analysis

Notebook 4 provides a comprehensive workflow for cell segmentation, feature extraction, and spatial analysis:

### Current Capabilities

**Segmentation (InstanSeg):**
- Nuclear segmentation from DAPI
- Combined cell + nucleus segmentation
- ECM (extracellular matrix) region segmentation via watershed from cell boundaries
- Matched compartment masks (cytoplasm, nuclear, ECM) with shared labels

**Feature Extraction:**
- Uses `napari-simpleitk-image-processing.label_statistics` for:
  - Morphological features (size, perimeter, shape descriptors)
  - Position features (centroids)
  - Intensity statistics (mean, median, min, max, sigma, variance, sum)
  - Moments for texture analysis
- Separate feature DataFrames for cell, nuclear, and ECM compartments

**SOM Clustering (pyFlowSOM):**
- Self-organizing map clustering with configurable grid size
- Learning rate scheduling (start -> end)
- Iterative refinement with visualization
- Separate clustering for: cells, nuclei, ECM, combined features

**Spatial Analysis (scanpy + scimap):**
- AnnData integration for single-cell analysis framework
- PCA, UMAP dimensionality reduction
- PAGA graph-based analysis
- Leiden community detection clustering
- Combined SOM + Leiden consensus phenotyping
- Spatial scatter plots with cluster overlays
- Hierarchical clustering for merged phenotypes
- Differential marker analysis (Wilcoxon rank-sum)

### Key Dependencies
- `instanseg` - Deep learning instance segmentation
- `pyFlowSOM` - Self-organizing maps for clustering
- `napari-simpleitk-image-processing` - Feature extraction
- `scanpy` - Single-cell analysis framework
- `scimap` - Spatial analysis toolkit
- `napari` - Interactive visualization

## Migration Notes

The old Notebooks 3 (Signal Isolation) and 5 (DL Channel Refinement) have been replaced with a unified `3_Signal_Isolation_QC.ipynb` that supports both Claude-guided and interactive workflows. See `notebooks/MIGRATION_GUIDE.md` for detailed transition guidance.
