<!-- Written 2026-09-02 from a read-only survey of the repository (notebook inventories with adversarial re-reads, infra touchpoint survey, source-verified marimo 0.24.0 research, and a three-way design judge panel). Cell indices refer to nb["cells"] positions in the .ipynb files at commit 716d1eb. -->

# Plan: Migrate KINTSUGI workflow notebooks from Jupyter to marimo

## Context

KINTSUGI's interactive workflow lives in seven Jupyter notebooks under `notebooks/`. Four of them (1, 2, 3_QC, 4) are copied into every project by `KintsugiProject.setup_notebooks()` and re-synced by `scripts/sync_to_projects.py`; 2.5 (vessel 3D) and 5 (KRONOS) are repo-only; `3_Signal_Isolation.ipynb` is a deprecated legacy notebook. The user wants all workflow steps moved to marimo (plain `.py` reactive notebooks, git-diffable, runnable as scripts, no ipykernel/ipywidgets/JupyterLab).

Exploration (four read-only survey workflows, ~25 agents including adversarial re-reads, plus source-verified marimo research and a three-way design judge panel) established that this is not a file-format conversion. The blockers are:

- **The notebooks' interactive layer is ipywidgets.** `notebooks/Kview/` is a vendored stackview fork built on ipycanvas + ipywidgets + ipyevents; the notebooks also call the pip `stackview` package directly (a core dependency). marimo does not render classic ipywidgets (verified in marimo source `repr_formatters.py` and open issue marimo-team/marimo#5099). 18 call sites also reference a `Kview2` module that no longer exists, so those cells are already broken.
- **Notebook 2 holds ~1,500 lines of processing logic inside cells** (z-plane stitch loop, deconvolution wrapper, EDF wrapper, registration driver, five `visualize_*` helpers). Several duplicate `workflow/scripts/*.py` and have drifted (registration `rigid_only=True` in the notebook is documented as wrong in `workflow/CLAUDE.md`).
- **Every notebook redefines globals across cells** (nb1 33 names, nb2 30, nb3_QC 44, nb4 51), which marimo forbids. `marimo convert` fixes neither this nor the magics.
- **~78 places outside the notebooks assume Jupyter**: `project.py` (`setup_notebooks`, VS Code settings, IPython probe in `prompt_for_existing_directory`), `sync_to_projects.py`, `deps.py` `NOTEBOOK_REQUIREMENTS`, four conda env files, CI `test-notebooks` job, README/docs/CLAUDE.md files, `hpc.py` README template, MCP tool `generate_jupyter_cell`, `param_tuner.py`/`artifact_tuner.py`.
- **`init_project()` blocks on `input()`** when a project directory is non-empty; a marimo cell cannot answer it.

Intended outcome: six marimo notebooks (`1_Single_Channel_Eval.py`, `2_Cycle_Processing.py`, `3_Signal_Isolation_QC.py`, `4_Segmentation_Analysis.py`, `2.5_Vessel_3D_Segmentation.py`, `5_KRONOS_Analysis.py`) that are thin reactive DAGs over module code shared with the Snakemake scripts, a marimo-native viewer package replacing Kview/stackview, distribution/sync/docs/CI updated, Jupyter packages removed from the environments, and a documented HPC launch recipe.

## Decisions taken (defaults; flag if any should change)

| # | Decision | Choice | Reason |
|---|---|---|---|
| 1 | Directory name | Keep `notebooks/` | ~20 `sys.path`/PYTHONPATH sites hard-code it (conftest, cli.py:179, bridges, all `slurm/jobs/*.sh`, all `workflow/scripts/*.py`, hpc.py:876) |
| 2 | File names | Keep numeric prefixes (`1_Single_Channel_Eval.py`) | Not valid identifiers, so they can never shadow a module while `notebooks/` is first on `sys.path` |
| 3 | Scope | Port 1, 2, 3_QC, 4, 2.5, 5. Retire legacy `3_Signal_Isolation.ipynb` after folding its six uncovered capabilities into 3_QC | MIGRATION_GUIDE already deprecates it; `setup_notebooks()` never distributed it |
| 4 | Widget replacement | New `src/kintsugi/marimo/` package (`mo.ui.*` + `mo.image`/`mo.mpl.interactive`); no anywidget ports; `Kview.interact`'s `globals()` scraping not ported | ipywidgets unsupported; anywidget path has open reliability issue #10494 |
| 5 | Dependencies | `marimo>=0.24,<0.25` and `plotly` become core deps; jupyterlab/ipywidgets/ipykernel leave `envs/*.yml`; `stackview` leaves core in the last PR | End state has no Jupyter dependency; marimo's wheel is light (starlette/uvicorn, no numpy) |
| 6 | Coexistence | `.ipynb` stays beside its `.py` until that notebook passes HPC acceptance, then one PR flips the manifest and `git rm`s it | Sandbox cannot run GPU notebooks; per-notebook rollback is `git revert` |
| 7 | Runtime config | Per-notebook PEP 723 `[tool.marimo.runtime]` block (travels with project copies) + repo `pyproject.toml` `[tool.marimo]` | Script metadata has highest precedence and is copied with the file; project dirs have no pyproject |
| 8 | Module homes | Shared drivers go into existing modules (`kintsugi.stitch_blend`, `kintsugi.edf`, `kintsugi.kreg`, `kintsugi.signal.batch`), never a parallel package | Avoids a third copy of each driver |
| 9 | Interactive vs batch pixels | Interactive signal isolation previews and applies through `kintsugi.signal.batch` (same code path as Snakemake) | Today 3_QC's apply loop uses `Kutils.ini_params`/`clean` while Snakemake uses `process_batch`, so users tune one algorithm and ship another |
| 10 | HPC access | `marimo edit --headless --port` inside an sbatch GPU job + SSH jump-host tunnel, token auth kept, never `--session-ttl`/`--timeout`/`--sandbox`/`--host 0.0.0.0` | marimo's documented remote recipe; `edit` sessions are not reaped on websocket disconnect by default |
| 11 | Zarr output branch in nb2 | Dropped (TIFF only) | `KintsugiZarr(mode='w')` truncates the store (`zarr_io.py:130-133`); Snakemake is TIFF-only |
| 12 | Legacy blend in nb1/nb2 | Dropped in favour of `stitch_blend.stitch_with_blending` | Module docstring already declares it the replacement |
| 13 | Registration defaults | `rigid_only=False` with tuned non-rigid params | `workflow/CLAUDE.md` documents the notebook's `rigid_only=True` as wrong |
| 14 | 2.5 output layout | `data/processed/vessel_3d/cyc{NN}/` (workflow layout) | Notebook's flat layout collides across cycles |
| 15 | napari in nb4 | Out-of-process launcher behind a run button + in-cell `klabels_overlay`; behaviour recorded by an HPC spike | napari/Qt undocumented by marimo (zero issues) |
| 16 | Notebook sync | Opt-in `--notebooks` flag with `.bak` on MD5 mismatch; module sync unchanged | User-edited `.py` notebooks must never be clobbered silently |
| 17 | `generate_jupyter_cell` MCP tool | Renamed `generate_marimo_cell`; deprecated alias kept one minor release | Public tool name |
| 18 | Post-commit hook claim | Docs corrected to "manual sync" | No hook exists in the repo despite README/CLAUDE.md claims |
| 19 | `kintsugi workspace` + out-of-repo `.code-workspace` | Kept as a VS Code source-editing convenience; CLAUDE.md reworded | Marimo notebooks do not depend on it |
| 20 | Wheel shipping `notebooks/` | Unchanged (not shipped); recorded as debt | Out of scope |

## Verified marimo facts the plan relies on (marimo 0.24.0, 2026-08-17)

- Python >= 3.10, conda-forge `marimo` on linux-64 (not ppc64le/aarch64), light dependencies.
- Multiple-definition rule; underscore-prefixed names are cell-local and freed after the cell runs; imports count as definitions.
- `marimo convert nb.ipynb -o nb.py` then `marimo check --fix/--strict` (breaking rules are the worklist). No magic translation.
- Runtime keys (`marimo/_config/config.py`): `on_cell_change = "lazy"`, `auto_instantiate = false`, `auto_reload = "lazy"|"autorun"` (the `%autoreload` replacement), `watcher_on_save`, `output_max_bytes` (default 8 MB), `std_stream_max_bytes` (default 1 MB), `pythonpath`, `reactive_tests`. Settings pinned in pyproject/script metadata cannot be changed from the editor UI.
- Gating: `mo.ui.run_button()` + `mo.stop()`; `mo.ui.dictionary(...).form()` (value is `None` until submitted; a form-wrapped slider cannot drive a live preview); `mo.cache` (memory); `mo.persistent_cache(save_path=...)` pickles and hashes array arguments by value (avoid for multi-GB arrays).
- Display: `mo.image(arr, vmin, vmax)` renders arrays; `mo.mpl.interactive(fig)` pan/zoom (SubFigure falls back to PNG); `mo.ui.matplotlib(ax)` box/lasso masks; return `fig` as the last expression, no `plt.show()`.
- `mo.status.progress_bar`/`spinner` replace tqdm; `mo.output.replace/append` replace `display`/`clear_output`.
- `with app.setup:` + `@app.function` make importable top-level functions (may only reference setup-cell names).
- `mo.notebook_dir()`, `mo.cli_args()`, `mo.app_meta().mode` (`edit|run|script|test`); `python nb.py -- --project X` runs the DAG as a script, but a `run_button` is `False` and a form is `None` without a frontend.
- Remote: `marimo edit --headless --port N` on the node; `--token` on by default; `--session-ttl` default `None` for `edit` (never reaped), `120` for `run`.
- pytest collects cells containing only `test_*` functions; `marimo export ipynb --include-outputs` executes a notebook.
- napari/Qt/CuPy undocumented by marimo; the only CUDA issue (#806, memory not freed on rerun) is closed but warrants an empirical check.

## Module readiness (drives port order)

| Tier | Stages | Meaning for the port |
|---|---|---|
| 1. Thin already | vessel3d (nb 2.5), `Kprocess.run_*_qc`, `kintsugi.signal` features/clustering/optimizer/predictor, spillover, kronos (nb 5), `Kdecon.decon`, `EDFProcessor`, `KCorrectGPU.fit/transform` | Marimo cell = one module call |
| 2. Duplicated driver glue | per-z-plane stitch loop (`workflow/scripts/stitch.py:292-397` vs nb2 cell 16), EDF file wrapper + quality gate (`workflow/scripts/edf.py:160-300` vs nb2 cell 31), registration driver (`workflow/scripts/registration.py:224-430` vs nb2 cell 38, drifted) | Extract one shared function into the existing module, re-point both callers, prove equivalence |
| 3. Notebook-only | nb3_QC apply loop (bypasses `signal.batch.process_batch`), CLAHE/gaussian chain (only in `Kutils`), nb4 quantification / ECM watershed / label matching / SOM+Leiden merge / ~55 scanpy cells, `Kanalysis.visualize_som` (fit fused with frame rendering, imports `tqdm.notebook`) | Extract to modules where load-bearing; translate the rest cell by cell |

## Per-notebook findings that shape the port

Cell indices are `nb['cells']` indices (markdown included). Each notebook was inventoried by one agent and adversarially re-read by another.

| Notebook | Code cells | Globals defined in >1 cell | Notebook-defined logic | Widgets | Already broken today |
|---|---|---|---|---|---|
| 1_Single_Channel_Eval | 34 | 33 | `smooth_tile_borders_2d` (dup of nb2), `optimize_for_display`, `load_channel_names`/`save_channel_names_template` (dups of `Kio:39/:203`), `_resize_images_list` (dup of `Kio.resize_images_list:548`) | `stackview` + `tqdm.notebook` imported, unused; `Kview2.imshow` cell 76 | Never executed (only cell 9 has an execution count); cell 76 NameError; cells 30/70 parameter blocks are dead (nothing reads them); cells 60/62 read `stitch_dir` that nothing in this notebook writes |
| 2_Cycle_Processing | 29 | 34 | ~1,500 lines: `process_zplane`/`process_channel_zplanes_parallel` (cell 16), stitch driver + multi-GPU queue (19), `run_deconvolution` (24), `process_edf_tiff` (31), registration driver (38), 5 `visualize_*` (13/21/27/34/40), `smooth_tile_borders_2d` (18) | none (tqdm.notebook, `clear_output`) | `KintsugiZarr(mode='w')` truncates the store; `Kio.cleanup_gpu_memory` shadows `kintsugi.gpu.cleanup_gpu_memory` with a different positional signature; cell 42 deletes TIFFs; `stage_times` mutated by 4 cells; `process_zplane` returns `None` on two paths that callers unpack |
| 3_Signal_Isolation_QC | 33 | 44 | `scan_registered_channels` (4), `extract_channels_from_registered` (7), `generate_reference_histogram` (43) | `Kview2.interact` ×4, `Kview2.switch`, `stackview.switch/crop/curtain/insight/imshow` | All 5 tuner cells NameError (`Kview2`); `clahe_params={'clip':…}` and `clean_params={'radius','amount',…}` do not match `Kutils.CLAHE`/`Kutils.clean` signatures; `from kintsugi.mcp.tools.learning import ParameterLearningEngine` fails (real home `kintsugi.claude.parameter_learning`); 29 of 33 cells use `'x' in dir()` guards; outputs written into the directory the loader globs (feedback loop); `result` bound to three types |
| 4_Segmentation_Analysis | 86 | 46 | none as functions, but ~55 inline scanpy cells, 3 copy-pasted `label_statistics` loops, 4 copy-pasted SOM blocks, 3 copy-pasted AnnData builders | `Kview2.curtain` ×5, `stackview.slice` ×4, napari ×2 | `proc_dir` never defined (cells 7/15/18/26/28/45); cell 26 reads `cells_noResolve`, a path Kseg never writes (it writes `<stem>_instanseg_prediction.tiff`); cells 109/110 use `leiden_Res1` but cell 106 writes `leiden_Res0.6`; cell 87 reads `_ecm.h5ad` after cell 86 writes `_ECM.h5ad`; cells 82/84/86 read `*_clusters` columns `visualize_som` never emits; `np.isnan` on uint16 `X`; cell 53 overwrites `cell_stats.csv` in place |
| 2.5_Vessel_3D_Segmentation | 18 | 2 | none (thin over `kintsugi.vessel3d`) | none | cell 6 writes into a directory created only at cell 24; flat output layout differs from `workflow/scripts/vessel3d.py:70` |
| 5_KRONOS_Analysis | 17 | 4 | none (thin over `kintsugi.kronos`) | none | never executed; `total_mem` should be `total_memory` (cell 3); `plt.cm.get_cmap` (removed in matplotlib 3.9) at `kronos/analysis.py:221` and `:387`; UMAP fit three times |
| 3_Signal_Isolation (legacy) | 37 | 23 | byte-identical copies of 6 `Kutils` functions in cell 4 | `Kview2.interact` ×8, stackview ×13 | `import Kview2` fails; does not run top-to-bottom. Capabilities absent from 3_QC: gaussian chain (`find_gauss`/`clipped_gauss`/`gauss_factor`), second blank subtraction (`SubtractionParams.second` slot exists), OME-TIFF channel splitting (`Kio.extract_channels_from_ome:442` exists), animated-curtain GIF, shape-mismatch audit, legacy `*_param.txt` recipe writer (the only format `signal.batch.load_legacy_recipes:184` reads; `batch_multi._write_recipes_as_param_files:760` already produces it) |

Cross-cutting facts:

- **Kview usage is 9 of 24 exports**: `interact`, `switch`, `curtain`, `slice`, `crop`, `imshow`, `insight`, `animate_curtain`, `create_colormap`. The rest (`annotate`, `picker`, `scatterplot`, …) are unused and will not be ported.
- **`Kutils.py` is already the pure numeric layer** (`ini_params:206`, `denoise:238`, `find_gauss:262`, `clipped_gauss:267`, `gauss_factor:278`, `CLAHE:292`, `clean:48`); `Kview.interact` only auto-generates sliders from their signatures via `Context(globals())`.
- **Pure, portable Kview parts**: `_colormaps.py` (`create_colormap:38`, `_labels_lut:108`), `_image_widget._img_to_rgb:69-103`, `_utilities.py` (`merge_rgb`, `numpy_to_gif_bytestream:66`), `_imshow.py` (matplotlib), the curtain compositing closure `_curtain.py:99-104`, `_switch._image_stack_to_rgb:162`.
- **`Kview_qc.py` returns Plotly `go.Figure` from 6 of 9 viewers** and is imported by no notebook; `plotly` is undeclared in `pyproject.toml`. Only `view_all_edf_interactive` and `QCPanel.show` are ipywidgets.
- **`param_tuner.py` / `artifact_tuner.py`** are used by no notebook; only `generate_jupyter_cell` emits code calling them; its `full_pipeline` branch imports a nonexistent `kintsugi.claude.pipeline`. `kintsugi.claude.__init__` imports them eagerly.
- **Zero tests** cover Kview, Kutils, Kview_qc, param_tuner, artifact_tuner.
- **Two path conventions**: notebooks use `KintsugiProject.paths`; `workflow/scripts/*.py` rebuild paths by hand from `snakemake.params.project_dir`; `workflow/scripts/stitch.py:34-37` dereferences `snakemake` at import, so those scripts cannot be imported in tests (textual checks only).
- **Console cap**: nb2's drivers print per-tile progress; marimo truncates console output at 1 MB by default.

## Approach

Risk-first, thin-notebook migration in 15 PRs on `claude/kintsugi-jupyter-marimo-migration-a9gff1`:

1. A safety net first (synthetic project + hash-frozen copies of the legacy cells) so CPU equivalence of extracted drivers can be proven without a GPU.
2. Shared enablers (dependencies, runtime config, a notebook manifest as the single source of truth, the `kintsugi.marimo` helper package, safe init/sync, agent-guidance rewrite).
3. Tier-1 notebooks (2.5, 5) ported first as low-risk proofs, plus the HPC launch story.
4. Driver extraction into existing modules (signal isolation, stitch, EDF, registration, nb4 quantification/analysis), each with equivalence tests, with `workflow/scripts/*.py` re-pointed at the same functions.
5. Ports of 3_QC, 1, 2, 4 in that order; `.ipynb` retired per notebook after HPC acceptance; Jupyter stack removed last.

Sandbox note: this environment has no numpy/marimo/kintsugi installed, no GPU, and no `test_data/mini_project`. Every executing session starts with `python -m venv .venv && source .venv/bin/activate && pip install -e '.[dev]'` (and `pip install pyvips-binary` if libvips is missing), all via the configured proxy.

## Phase A — Safety net (PR1)

`test(notebooks): synthetic project fixture + frozen legacy cells`

- `tests/fixtures/synthetic_project.py::make_synthetic_project(root, n_rows=2, n_cols=2, n_z=3, n_ch=2, tile=(96,96), overlap=0.3, seed=0)` writes overlapping crops of one smooth random field as `data/raw/cyc001/1_{tile:05d}_Z{z:03d}_CH{c}.tif`, `meta/CHANNELNAMES.txt`, `meta/experiment.json`, `kintsugi_project.json`.
- `tests/legacy_reference/` holds verbatim copies between `# BEGIN FROZEN`/`# END FROZEN` markers of: nb2 cell 16 (`process_zplane`), nb2 cell 31 (`process_edf_tiff`), nb3_QC cell 55 (apply loop), nb4 cells 30/32/47/49/51 (compartments + quantification), nb1 cell 41 (`smooth_tile_borders_2d`). `tests/test_legacy_reference_frozen.py` asserts `sha256(cell source) == sha256(frozen text)` and skips once the `.ipynb` is gone.
- `tests/fixtures/legacy_loader.py::load_frozen(name, **namespace)` execs the frozen text into a namespace whose `KCorrectGPU` subclass forces `use_gpu=False`, whose `stitch_images` wrapper forces `use_gpu=False`, and whose `cleanup_gpu_memory` is a no-op (nb2 cell 16 hard-codes `use_gpu=True` at two sites; text is never edited, only names are shimmed).

## Phase B — Enablers (PR2, PR3)

### PR2 `build(marimo): deps, runtime config, notebook manifest, spike scripts`

Dependency/config table:

| Site | Add now | Remove in PR14 |
|---|---|---|
| `pyproject.toml` core deps (lines ~59-60) | `marimo>=0.24,<0.25`, `plotly>=5.18` | `stackview>=0.18.0` |
| `requirements.txt` | marimo, plotly | `stackview` (:24), Snyk `jupyter-core` (:71), `ipython` (:75) |
| `envs/env-hpc.yml` (90-93, 143), `env-linux.yml`/`env-macos.yml` (50-53, 83), `env-windows.yml` (51-54, 84) | `- marimo>=0.24`, `- plotly` | jupyterlab, ipywidgets, ipykernel, stackview |
| `src/kintsugi/deps.py:678` core check | — | `("stackview","0.18.0")` → `("marimo","0.24.0")`; docstring :264 |
| `scripts/verify_hpc_env.sh:83`, `envs/activate.d/env_vars.sh:7` | — | stackview references |
| `.gitignore:91-96` | `__marimo__/`, `.marimo.toml` | `.ipynb_checkpoints`, IPython lines |
| `constraints.txt`, `tests/test_deps.py` | unchanged (marimo is core, no new group) | — |

`pyproject.toml` additions:

```toml
[tool.marimo.runtime]
on_cell_change = "lazy"
auto_instantiate = false
auto_reload = "lazy"          # %autoreload 2 replacement (config.py key, not the UI label)
watcher_on_save = "lazy"
output_max_bytes = 33554432
std_stream_max_bytes = 4000000
reactive_tests = false
pythonpath = ["src", "notebooks"]

[tool.marimo.display]
default_width = "full"

[tool.marimo.formatting]
line_length = 100

[tool.ruff.lint.per-file-ignores]
"notebooks/[0-9]*.py" = ["B018", "E402", "E731", "F811", "F401"]
"scripts/marimo_spikes/*.py" = ["B018", "E402", "F811"]
```

Every notebook starts with the same PEP 723 block (project copies have no pyproject); `bootstrap.SCRIPT_METADATA_TEMPLATE` is the single copy and a test asserts each notebook begins with it. No `dependencies` key, never `--sandbox` on HPC (uv would isolate from CuPy).

```python
# /// script
# requires-python = ">=3.10"
# [tool.marimo.runtime]
# on_cell_change = "lazy"
# auto_instantiate = false
# auto_reload = "lazy"
# output_max_bytes = 33554432
# std_stream_max_bytes = 4000000
# [tool.marimo.display]
# default_width = "full"
# ///
```

Notebook manifest `src/kintsugi/notebook_manifest.py` replaces the four divergent lists (`project.py:1783-1841`, `scripts/sync_to_projects.py:57-79` incl. stale `Kview2`/`Kpipeline.py`/`Kvis.py`, `deps.py:157-168`, README/docs tables):

```python
@dataclass(frozen=True)
class NotebookSpec:
    stem: str; title: str; format: Literal["ipynb", "py"]; distributed: bool
    requires: tuple[str, ...]; steps: tuple[str, ...] = (); legacy: bool = False
    @property
    def filename(self) -> str: return f"{self.stem}.{self.format}"

NOTEBOOKS = (...)  # 7 entries; format flips "ipynb" -> "py" per notebook as it is retired
SUPPORT_DIRS = ("Kreg", "Kstitch", "Kview", "Kdecon", "Kseg")          # Kview removed in PR14
SUPPORT_FILES = ("Kio.py", "Kprocess.py", "Kutils.py", "Kanalysis.py", "Kview_qc.py",
                 "config_example.json", "MIGRATION_GUIDE.md")
def distributed_notebooks() -> list[NotebookSpec]
def notebook_requirements() -> dict[str, list[str]]    # deps.NOTEBOOK_REQUIREMENTS = this
def notebook_table_markdown() -> str                   # README:170-176, docs/DEPENDENCY_GUIDE.md:121-133
```

Spike scripts under `scripts/marimo_spikes/` with results recorded in `docs/marimo/SPIKES.md` (soft gates: a PR ships the fallback behind a flag and records `UNKNOWN` if the user has not run the HPC spike):

| ID | Question | Where | PASS | Fallback |
|---|---|---|---|---|
| S1 | CuPy memory across gated re-runs; threaded `cp.cuda.Device(i).use()` affinity | HPC GPU node | 5 re-runs of a 4 GB alloc, `mem_info` within 5% of baseline | nb1/nb2 GPU stages become submit-and-poll launchers over `kintsugi workflow run` |
| S2 | In-flight cell survives a tunnel drop under headless `marimo edit` in sbatch | HPC | 5 s heartbeat log has no gap >10 s across a 60 s outage | Long stages run as `python nb.py -- ...` in the same allocation |
| S3 | 10k×10k array display cost | sandbox | `mo.image` at `max_size=1024` <300 ms, <2 MB | Fixed 1024 previews + crop controls for full-res regions |
| S4 | `init_project` without blocking | sandbox pytest | `open_project()` on a non-empty dir with `builtins.input` patched to raise returns a project | code fix |
| S5 | sync/init never clobber user notebooks or leave stale `.ipynb` | sandbox pytest | `test_sync_to_projects.py`, `test_project_init.py` green | code fix |
| S6 | napari from a marimo cell | HPC/workstation | kernel responsive while `Popen` runs | nb4 writes a napari session script for Open OnDemand desktop |
| S7 | CI validates notebooks without a GPU | GitHub Actions | `marimo check --strict`, compileall, ruff, script-mode smokes pass | optional-dep skips |
| S8 | `auto_reload="lazy"` watcher cost on `/blue` with cupy/torch loaded | HPC | dependent cell marked stale within 5 s of editing `Kio.py`; idle CPU <5% | header sets `auto_reload="off"` for project copies; docs say restart `marimo edit` |

CI in PR2 also adds `pip show ipycanvas ipywidgets` to the `test-notebooks` job so the current `from Kview import ...` check's dependence on stackview's transitive installs is explicit.

### PR3 `feat(marimo): helper package, safe init/sync, CLAUDE.md rewrite`

New package `src/kintsugi/marimo/` (marimo imported inside functions; `import kintsugi` never requires marimo):

- `_render.py` (pure numpy, no marimo import): verbatim `create_colormap`, `labels_lut`, `to_rgb` (= `_img_to_rgb`), `merge_rgb`, `numpy_to_gif_bytes`; new `composite_curtain` (from `_curtain.py:99-104`), `blend_stack` (from `_switch.py:92-110`), `downsample_for_display(arr, max_size=1024)` (stride-slices dask arrays lazily, then computes the small view; test asserts no full-resolution compute), `encode_png`, `image_stats`.
- `viewers.py` (control/render pairs, because reactivity only crosses cell boundaries):

| Old call (sites) | New cells |
|---|---|
| `stackview.slice` (nb4 c60/63/66/70) | `z = kv.zslider(stack)` ; `kv.kslice(stack, z.value)` |
| `stackview.curtain`/`Kview2.curtain` (nb3 c39; nb4 c27/28/35-37) | `cc = kv.curtain_controls(w)` ; `kv.kcurtain(a, b, cc.value, 'gray', 'turbo')` |
| `stackview.switch`/`Kview2.switch` (nb3 c11/40) | `sel = kv.switch_controls(names, toggleable)` ; `kv.kswitch(images, sel.value, colormaps)` |
| `stackview.crop` (nb3 c12-13) | `crop = kv.crop_controls(shape)` ; `sl = kv.crop_slices(crop.value)` ; `kv.kcrop_preview(img, sl)`; `sl` feeds downstream cells |
| `Kview2.interact(fn, context=globals())` (nb3 c18/28/30/33) | `ctrl = tuners.*_controls().form()` ; preview cell calls the module function on `downsample_for_display(...)` |
| `stackview.imshow`/`Kview2.imshow` (nb1 c76) | `kv.kimage(arr, colormap='turbo')` |
| `stackview.insight` (nb3 c43) | `kv.kinsight(img)` |
| `stackview.animate_curtain` (legacy c77) | `kv.kgif(frames, path)` |
| `Kseg.utils.show_images`, napari (nb4 c21/55/126) | `kv.klabels_overlay(image, labels)`; napari via out-of-process launcher |
| `%matplotlib widget`, `plt.show()`, `display()` | cell returns `fig`; `kv.kfigure(fig, interactive=True)` = `mo.mpl.interactive`; helpers call `plt.close` |
| `tqdm`, `tqdm.notebook`, `ProgressCounter`, `log()` | drivers take `on_progress`/`on_log`; notebooks pass `bootstrap.progress`; verbose driver logs go to `<project>/slurm/logs/marimo/<stage>.log` |

- `bootstrap.py` (the piece every notebook and every CI smoke depends on):

```python
SCRIPT_METADATA_TEMPLATE: str
def add_module_paths(notebook_dir: Path | None) -> list[Path]
    # <project>/notebooks, <KINTSUGI>/notebooks, <KINTSUGI>/src via KintsugiProject._detect_kintsugi_path (project.py:1185)
def resolve_project_dir(default: Path | None = None) -> Path | None
    # mo.cli_args()['project'] > $KINTSUGI_PROJECT_DIR > mo.notebook_dir().parent if it holds kintsugi_project.json or meta/experiment.json > default
def project_input(default: Path | None) -> mo.ui.text
def open_project(project_dir, name=None) -> KintsugiProject
    # init_project(project_dir, name=name, mode="auto", interactive=False, on_existing="keep", check_processed=False)  (project.py:2314-2323); raises RuntimeError on None
def mode() -> str                      # mo.app_meta().mode or "script"
def steps_requested() -> set[str]      # --steps a,b | all | none (default none)
def run_button(step, label=None) -> mo.ui.run_button
def gate(button, step, output=None) -> None
    # mo.stop(not (button.value or (mode() == "script" and step in steps_requested())), output)
def resolve_params(form, defaults: dict, *, cli_prefix=None) -> dict
    # edit/run: mo.stop(form.value is None, form); return form.value
    # script: defaults | {k: parsed(v) for --<prefix>.<k>=v in mo.cli_args()} | json from --params <file>
def progress(iterable, title="")       # mo.status.progress_bar in edit/run, tqdm in script
@dataclass class GPUContext(use_gpu, multi_gpu, device_ids)
def gpu_context(force_cpu=False) -> GPUContext   # GPU required, fail loud (workflow/CLAUDE.md policy) unless force_cpu
```

- `tuners.py`: `subtraction_controls`, `clean_controls`, `denoise_controls`, `clahe_controls`, `gauss_controls`, `artifact_controls` return `mo.ui.dictionary` mirroring the `Kutils` / `param_tuner.py` knobs; `to_subtraction_params`/`to_clean_params`/`to_recipe` map onto `SubtractionParams` (`batch.py:126`) / `CleanParams` (`:142`) with `backgrnd_thresh → background_threshold`, `smooth_thresh → smooth_threshold` (read at `:277/:279`).

Standard cell layout every port follows:

```python
with app.setup:
    import marimo as mo
    import numpy as np
    from pathlib import Path
    from kintsugi.marimo import bootstrap as kb, viewers as kv
    kb.add_module_paths(mo.notebook_dir())

@app.cell
def _():
    project_dir_ui = kb.project_input(kb.resolve_project_dir())
    project_dir_ui
    return (project_dir_ui,)

@app.cell
def _(project_dir_ui):
    mo.stop(not project_dir_ui.value, mo.md("Enter a project directory"))
    project = kb.open_project(project_dir_ui.value, name=Path(project_dir_ui.value).name)
    return (project,)

@app.cell
def _(project):
    stitch_form = mo.ui.dictionary({...defaults from project.load_experiment_config()...}).form()
    stitch_btn = kb.run_button("stitch")
    mo.vstack([stitch_form, stitch_btn])
    return stitch_form, stitch_btn

@app.cell
def _(project, stitch_form, stitch_btn):
    _p = kb.resolve_params(stitch_form, STITCH_DEFAULTS, cli_prefix="stitch")
    kb.gate(stitch_btn, "stitch")
    stitch_result = stitch_cycle(..., StitchConfig(**_p), on_progress=kb.progress)
    return (stitch_result,)
```

Rules enforced by `tests/test_marimo_notebooks.py`: one setup cell for imports; unique output names; temporaries and loop targets underscore-prefixed; big arrays `del`-ed after saving; stage results are small dataclasses pointing at files; no `def` longer than 25 lines inside a cell; no magics; no `/blue/` literals; `warnings.filterwarnings` narrowed.

Also in PR3:

- `project.setup_notebooks()` (`project.py:1761-1843`) iterates `distributed_notebooks()`; when copying `<stem>.py` it moves an existing `<stem>.ipynb` to `notebooks/legacy_ipynb/`, backs up a differing destination `.py` to `notebooks/_backup/<stem>.<ts>.py`, and (once `Kview` leaves `SUPPORT_DIRS`) archives a stale `notebooks/Kview/`; returns the `copied` list. `kintsugi init --force` (`cli.py:4574`) inherits this. Printed next steps (`project.py:1370-1406`) say `kintsugi notebook launch .`.
- `scripts/sync_to_projects.py`: lists come from the manifest; notebook files are opt-in via `--notebooks` with `.bak` on MD5 mismatch (the clobber path is the explicit `SYNC_FILES` list fed to `sync_file:192-226`, not the `**/*.py` glob, which only walks `SYNC_DIRECTORIES`); module sync unchanged; `tests/test_sync_to_projects.py` added.
- `_create_vscode_config` (`project.py:1960-1981`): drop `jupyter.notebookFileRoot` and `.ipynb_checkpoints`; add `**/__marimo__` exclude; write `.vscode/extensions.json` recommending `marimo-team.vscode-marimo`.
- `prompt_for_existing_directory` (`project.py:2515+`): delete the IPython probe; `if not sys.stdin.isatty(): return "keep"`.
- `notebooks/CLAUDE.md` rewritten (reactive DAG semantics, one-definition rule, gating, `auto_reload` replaces `%autoreload`, header escape hatch, port status table) and root `CLAUDE.md` "Jupyter Notebook Workflow" section (lines ~346-350) replaced.

## Phase C — Tier-1 ports and HPC launch (PR4, PR5)

### PR4 `feat(notebooks): port 2.5 + 5; CI notebooks job`

Common porting recipe (`docs/marimo/PORTING_RECIPE.md`):

1. `python scripts/convert_notebook_to_marimo.py notebooks/N.ipynb` (runs `marimo convert -o`, strips `%`/`!` lines, injects the PEP 723 header, sets `app = marimo.App(width="full", app_title=...)`, prints `marimo check --strict`); commit as the first commit of the PR so hand edits are the reviewable diff.
2. Fix breaking rules with the cheapest correct edit, in order: underscore locals > merge cells > rename genuinely different values > `def _run(): ...; x = _run()` > move to module (only when the code already exists in `src/` or `workflow/scripts/`).
3. Setup/project cells per the standard layout; delete tkinter, OpenCL probes, `'X' in dir()` guards.
4. Replace widget sites per the viewers table; gate expensive cells with `kb.run_button` + `kb.gate`; parameter blocks become `mo.ui.dictionary(...).form()` read through `kb.resolve_params`; live previews use plain sliders.
5. Verify: `marimo check --strict`, `ruff check`, `python -m compileall`, `pytest tests/test_marimo_notebooks.py tests/test_notebook_script_mode.py`, `python N.py -- --project <synthetic> --steps none` plus CPU steps where they exist.

| Notebook | DAG (form → run_button → module call → viewer) | Deleted / fixed |
|---|---|---|
| 2.5_Vessel_3D | target form (cycle, channel, marker dropdown from `discover_vessel_markers`, device) → gated load (`Kio.load_stack_parallel(max_workers=4)`) / `preprocess_volume` / `compute_vesselness_frangi` / otsu hist + threshold slider / `binarize_vessel_mask` / `skeletonize_vessels`+`prune_skeleton` / `analyze_vessel_graph` / `export_vessel_results`; radio for the one-call `segment_vessels_3d` path | mkdir before first savefig; `plt.close` after each viz; per-cycle output layout; `fig`/`axes` underscored |
| 5_KRONOS | config form → dependency status via `kintsugi.gpu` → `MarkerMapper` → gated `KronosModel.load` → embed source radio (compute/load) → Leiden/KMeans sliders → single `umap_coords` cell → plots → query slider → spatial search → cross-dataset form (no `/blue` literals) → `marker_variance` → gated export | `total_mem` bug; mkdir OUTPUT_DIR; `get_cmap` ×2 → `matplotlib.colormaps[...].resampled(n)`; new `kronos.analysis.marker_variance`/`plot_marker_variance` |

Both `.ipynb` files are removed and the manifest flipped in this PR (repo-only, CPU-verifiable). CI `test-notebooks` job (`ci.yml:81-131`) becomes `notebooks`: `pip install -e '.[dev]'`; `marimo check --strict notebooks/[0-9]*.py`; `python -m compileall -q notebooks/[0-9]*.py`; `pytest tests/test_marimo_notebooks.py tests/test_notebook_script_mode.py`; Kreg file-existence and `Kstitch` import checks kept. `lint` job extends `ruff check` to `notebooks/[0-9]*.py`.

### PR5 `feat(hpc,cli): kintsugi notebook launch + SLURM session job`

`slurm/jobs/00_marimo_session.sh` (lives in `slurm/jobs/` so `project.py:2061-2073` symlinks it into projects):

```bash
#SBATCH --gres=gpu:1 --time=08:00:00 --mem=64G --cpus-per-task=8
module load conda && conda activate KINTSUGI
export KINTSUGI_PROJECT_DIR="${PROJECT_DIR}"
export MARIMO_OUTPUT_MAX_BYTES=33554432 MARIMO_STD_STREAM_MAX_BYTES=4000000
PORT=$(( 20000 + RANDOM % 20000 ))
echo "ssh -J ${USER}@<login-node> -N -L 2718:127.0.0.1:${PORT} ${USER}@$(hostname)" > "${PROJECT_DIR}/slurm/marimo_session.txt"
exec marimo edit "${NOTEBOOK:-${PROJECT_DIR}/notebooks}" --headless --host 127.0.0.1 --port "${PORT}"
# token auth ON by default; never --no-token, --host 0.0.0.0, --session-ttl, --timeout, --sandbox
```

marimo binds loopback on the compute node, so the tunnel must jump through the login node (`-J`) to the node's `127.0.0.1`; fallback if compute nodes refuse `ssh -J` is `--host 0.0.0.0` behind token auth with a single-hop `-L 2718:<node>:<port>`. CLI: `kintsugi notebook launch [--project DIR] [--notebook N] [--dry-run]`, `kintsugi notebook status` (prints the tunnel line and token URL), `kintsugi notebook list`. `tests/test_dashboard.py:219` gains a `marimo-session` case. Docs: `docs/marimo/HPC.md`. Multi-hour runs are steered to `kintsugi workflow run` or `python N.py -- --project P --steps ...` for logs and restartability.

## Phase D — Driver extraction (PR6, PR8, PR11)

### PR6 `refactor(signal): preview/writer/enhance; MCP generate_marimo_cell`

- `src/kintsugi/signal/batch.py`: factor `run_recipe_pipeline(signal, blank, recipe, *, second_blank=None, clip_percentile)` out of `_process_channel_recipe` (:854); `process_channel` (:1123) calls it; add `preview_recipe(signal, blank, primary, clean=None, second=None, second_blank=None)`; promote `batch_multi._write_recipes_as_param_files` (:760) to public `write_legacy_recipes(recipes, recipe_dir)` with a round-trip test through `load_legacy_recipes` (:184) asserting `background_threshold`/`smooth_threshold` survive.
- New `src/kintsugi/signal/enhance.py`: `apply_blank_subtraction`, `apply_denoise`, `clahe` (from `Kutils.CLAHE:292`), `apply_clean` (delegates to `clean_background:398`), `apply_gaussian_subtract` (from `Kutils.find_gauss/clipped_gauss/gauss_factor:262-290`); `qc.stripe_artifact.get_fft_visualization` from `artifact_tuner.py:85-96`.
- MCP: rename `generate_jupyter_cell` → `generate_marimo_cell` in `tool_specs.py:483`, `tools/workflow.py:390`, `tests/test_mcp_tools.py:55,853-869`, `docs/cli.md:123`, `CLAUDE.md`; keep a deprecated alias spec; result key `marimo_cell`; snippets become two-cell `tuners.*_controls().form()` + `enhance.*`/`preview_recipe`; `full_pipeline` emits `process_channel`; test `compile()`s the snippet.
- `kintsugi.claude.__init__:36-42` eager tuner imports → PEP 562 `__getattr__` shim with `DeprecationWarning`.
- `tests/test_signal_parity.py`: `preview_recipe` vs `Kutils.ini_params + Kutils.clean` on synthetic arrays; prints max abs diff (decision 9 evidence).

### PR8 `refactor(stitch,edf,kreg): shared drivers; workflow scripts re-pointed`

- `src/kintsugi/stitch_blend.py`: `StitchConfig` dataclass (grid, overlap, `pou`, `initial_ncc_threshold`, BaSiC knobs, `tile_order` replacing the hard-coded `corrected[[0,1,3,2]]` 2×2 reorder, blend params), `correct_tiles` (`KCorrectGPU.fit:432` + `.transform:592` with the `flatfield_min=0.1` clamp), `compute_stitch_model`, `stitch_zplane`, `stitch_cycle(config, on_progress, on_log)` with the multi-GPU queue from nb2 cell 19. Model application extends `stitch_corrected_tiles` (:327) rather than adding a third path. `Kio.ProcessingConfig` (`Kio.py:849`) becomes a deprecated alias. `workflow/scripts/stitch.py` deletes its private `alphanumeric_key/load_tiles_parallel/_load_and_correct_tiles/_stitch_with_model/_compute_stitch_model/process_zplane` and calls the module; boundary-QC fallback (:399-552), config hash and skip logic stay. Equivalence: frozen nb2 `process_zplane` (CPU-forced) vs new on the synthetic project → `np.array_equal` per z-plane TIFF, `assert_frame_equal(result_df)`.
- `src/kintsugi/edf.py`: move `quality_gate_zstack` from `workflow/scripts/edf.py:160` (byte-identical to nb2 cell 31's gate) into the module; `EDFParams(z_margin=1, ...)` where `z_margin=1` reproduces the workflow bounds (`workflow/scripts/edf.py` z_start=2/z_end=n-1) and `z_margin=0` the notebook's `1..n`; `process_edf_channel(decon_channel_dir, output_path, params, *, device_id=None, processor=None, on_log=None)`; `Kio.channel_output_name` replaces `workflow/scripts/edf.py:107-115`. Equivalence vs frozen `process_edf_tiff` with `z_margin=0`, `backend='numpy'`.
- `src/kintsugi/kreg.py`: `RegistrationParams` (+ `from_default_parameters(project)` reading `configs/default_parameters.json`, written by `project.py:1845`) and `run_registration(...)` whose body is `workflow/scripts/registration.py:224-430` (Valis with `non_rigid_registrar_cls=None`, `register()`, `register_micro(...)`, warp, save, staging cleanup). The notebook's `rigid_only=True` and its copy-on-exception fallback are not ported. Kreg imported lazily inside the function; VALIS-free helpers unit-tested; VALIS itself is HPC acceptance A8.
- Deconvolution: `Kdecon.decon` (`notebooks/Kdecon/main.py:378`) is already shared; notebooks pass kwargs derived from `project.load_experiment_config()` (`project.py:1497`).
- `tests/test_workflow.py` gains a textual check that `workflow/scripts/stitch.py`/`edf.py`/`registration.py` import the shared drivers (the scripts dereference `snakemake` at import and cannot be imported).

### PR11 `refactor(segment,analysis): quantify, compartments, SOM split, phenotype`

- `src/kintsugi/segment/quantify.py::quantify_compartments(image, labels_by_compartment, marker_names, ...)` — the three `label_statistics` loops (nb4 cells 47/49/51) with the load-bearing `cell_/nuclear_/ecm_` + `{marker}_{stat}` + `centroid_0/1` column contract asserted against the frozen reference.
- `src/kintsugi/segment/compartments.py::derive_compartment_labels(cell_labels, nuclear_labels)` — ECM watershed + matched cytoplasm/nuclear labels (cells 30/32).
- `src/kintsugi/analysis/som.py::fit_som` + `plot_som_pass`; `Kanalysis.visualize_som` (:341) becomes a wrapper and loses `from tqdm.notebook import tqdm` (`Kanalysis.py:396`).
- `src/kintsugi/analysis/phenotype.py`: `merge_compartment_tables`, `build_anndata`, `run_standard_pipeline` (filter/pca/normalize/log1p/scale/neighbors/paga/umap/leiden), `merge_clusters_hierarchical(n_final)` (cell 113).

## Phase E — Ports of the distributed notebooks (PR7, PR9, PR10, PR12)

Each PR: convert-first commit, then hand edits; `.py` lands beside the `.ipynb` (manifest still `ipynb`) except where noted.

| Notebook | PR | DAG | Deleted / fixed |
|---|---|---|---|
| 3_Signal_Isolation_QC (+ legacy fold-in) | PR7 (`feat(notebooks)!`) | project → channels (`discover_channels:511` + `Kio` copy button) → load → `kswitch` browse → crop controls → pair form (signal, blank, tissue promoted to config) → `subtraction_controls().form()` → `preview_recipe` + `kcurtain` → **second blank** → weighted (`AutofluorescenceSubtractor`) → denoise (`kintsugi.denoise`) → **gaussian chain** (`enhance.apply_gaussian_subtract`) → CLAHE → clean → `ImageQC` → compare (`kcurtain`/`kswitch`+DAPI) → recipe (`write_legacy_recipes` + `mo.download` JSON) → apply one (`process_channel`) → **GIF export** → learn (`kintsugi.claude.parameter_learning`) → 7B HITL (`mo.ui.table` selection) → batch (`process_batch:1405`, `recipe_dir=configs/recipes`) → QC pages (`generate_qc_pages` + `mo.pdf`) → **merged OME export** with correct `Channel` metadata; **OME-TIFF splitting** via `Kio.extract_channels_from_ome:442`; **shape audit** table | `clip`→`clip_limit`; clean kwargs → real signature; `result` split three ways; broken learning import; outputs to `signal_isolated/` not the input dir; hand-rolled cell 55 loop deleted; `3_Signal_Isolation.ipynb` removed |
| 1_Single_Channel_Eval | PR9 | project → hw → detect (`normalize_cycle_dirs` dry-run + confirm) → select form → tiles → `zslider`/`kslice` → BaSiC form (cell 30 literals become live) → `correct_tiles` → flatfield/darkfield views → `kcurtain` before/after → stitch form (`tile_order`) → `compute_stitch_model` → `meta/result.pkl` → stitched views → save → decon form/run → compare → EDF form → `process_edf_channel` → `kimage` | `_resize_images_list`, `smooth_tile_borders_2d`, 144-line channel-name paste, OpenCL cell, tkinter, `optimize_for_display`, `Kview2.imshow`; `SELECTED_GPU` threaded through all four GPU stages |
| 2_Cycle_Processing | PR10 | project form → hw (`gpu_context`, force-CPU switch) → experiment form (defaults from `experiment.json`) → `StitchConfig` → names → raw preview (`Kview_qc.view_raw_tiles`) → `run_raw_qc` → `stitch_cycle` per cycle → view → `run_stitched_qc` → decon form → `decon` per cycle × channel → `kcurtain` → `run_decon_qc` → EDF form → `process_edf_channel` → view → `run_edf_qc` → registration form (`RegistrationParams.from_default_parameters`) → `run_registration` → `run_registration_qc` → overlay → summary | zarr branch, legacy blend + cell 18, destructive cell 42, `log()`/`ProgressCounter`/tqdm/IPython, dead shims; script mode `--steps stitch,edf --force-cpu` |
| 4_Segmentation_Analysis | PR12 | project → markers (`CHANNELNAMES.txt` canonical) → inputs (dropdowns over files) → InstanSeg form/run (`Kseg`) → `klabels_overlay` → save → `derive_compartment_labels` → spillover form/run → `quantify_compartments` → `filter_cells_pipeline` (→ `cell_stats_qc.csv`) → SOM form/`fit_som` → pass slider/`plot_som_pass` → `merge_compartment_tables` → `build_anndata` → `run_standard_pipeline` → umap/spatial → `merge_clusters_hierarchical` (n_final slider) → dotplot/rank_genes → h5ad; a "load existing" switch replaces ~10 read-back cells | `proc_dir`; `cells_noResolve` path; `ECM_`/`ecm_` prefix; `leiden_Res1` key; h5ad case; dotplot on the right AnnData; cell 53 self-overwrite; napari → `launch_napari` run button |

## Phase F — Retirement and cleanup (PR13, PR14, PR15)

- PR13 `feat(notebooks)!: retire accepted .ipynb` (may be split per notebook): manifest flip + `git rm` after each notebook's HPC acceptance item; `git grep '\.ipynb'` sweep; project migration verified by `kintsugi init --force`.
- PR14 `build(deps)!: remove Jupyter stack`: delete `notebooks/Kview/`, `src/kintsugi/kview.py`, `param_tuner.py`, `artifact_tuner.py`, `notebooks/check_setup.py`, `test_multi-root.ini`; drop stackview/jupyterlab/ipywidgets/ipykernel per the table; rewrite `Kview_qc.view_all_edf_interactive`/`QCPanel.show` with `mo.ui`; remove `from Kview import ...` from CI; fresh-venv full pytest.
- PR15 `docs: coherence sweep, debt register, CHANGELOG`: README (`Notebooks` heading kept for the CI substring check, `.ipynb` links → `.py`, "Running Notebooks" → `marimo edit`/`kintsugi notebook launch`, autoreload note removed, `Kview2`→`Kview`, hook claim corrected), `docs/quickstart.md`, `docs/workflows.md:136-138`, `docs/installation.md:304-305`, `docs/TROUBLESHOOTING.md:176`, `docs/DEPENDENCY_GUIDE.md`, `notebooks/MIGRATION_GUIDE.md`, `workflow/CLAUDE.md:5-16`, `src/kintsugi/CLAUDE.md:20,46`, `hpc.py:1205-1207`, `deprecation.py:108-110`, `dl_refinement.py:37-38`, `hitl.py:10,36`, `artifact_scanner.py:621`; new `docs/marimo/MIGRATION.md` (Jupyter → marimo for users), `docs/marimo/DEBT.md`; CHANGELOG via commitizen (breaking: MCP tool rename, dropped ipywidgets).

PR dependency graph: 1 → 2 → 3 → {4, 5, 6, 8, 11 in parallel}; 7 needs 6; 9 and 10 need 8; 12 needs 11; 13 needs acceptance; 14 needs 13 and A11; 15 last.

## Critical files

- `src/kintsugi/project.py` (`setup_notebooks` 1761-1843, `_create_vscode_config` 1960-1981, `init_project` 2314-2323, `prompt_for_existing_directory` 2515+, `_detect_kintsugi_path` 1185, `load_experiment_config` 1497, `_create_default_params_file` 1845)
- `src/kintsugi/marimo/bootstrap.py` (new; `gate`/`resolve_params`/`open_project` — every notebook and CI smoke depends on it)
- `src/kintsugi/stitch_blend.py` with `workflow/scripts/stitch.py:292-397` (largest driver unification)
- `src/kintsugi/edf.py` with `workflow/scripts/edf.py:107-300`; `src/kintsugi/kreg.py` with `workflow/scripts/registration.py:224-430`
- `src/kintsugi/signal/batch.py` (`process_channel` 1123, `_process_channel_recipe` 854, `load_legacy_recipes` 184, `process_batch` 1405) and `batch_multi.py:760`
- `scripts/sync_to_projects.py` (57-79, 142-226), `.github/workflows/ci.yml` (30-37, 81-131), `pyproject.toml` (59-60, 313-342), `envs/*.yml`
- `notebooks/Kview/{_colormaps,_image_widget,_utilities,_curtain,_switch}.py` (sources of the pure render helpers), `notebooks/Kutils.py`, `notebooks/Kanalysis.py:341-420`
- `notebooks/CLAUDE.md`, root `CLAUDE.md`

## Verification

Static and unit (every PR, sandbox and CI):

```bash
python -m venv .venv && source .venv/bin/activate && pip install -e '.[dev]'
ruff check src/ tests/ notebooks/[0-9]*.py && ruff format --check src/ tests/
marimo check --strict notebooks/[0-9]*.py && python -m compileall -q notebooks/[0-9]*.py
pytest tests/ -v
```

Key new tests: `test_marimo_notebooks.py` (PEP 723 header == template; banned tokens `%`, `get_ipython`, `ipywidgets`, `stackview`, `Kview`, `Kview2`, `tqdm.notebook`, `IPython`, `input(`; no `/blue/`; thinness rule; parity guards such as `process_batch(` present and `Kutils.ini_params(` absent in nb3; importing the module defines `app` without importing cupy/torch/scanpy), `test_notebook_script_mode.py` (2.5 full CPU run; nb1 `correct,stitch --force-cpu`; nb2 `stitch,edf --force-cpu` and `none`; nb3 `extract,load,apply`; nb4 `quantify`; nb5 `none`), `test_stitch_driver.py`/`test_edf_driver.py`/`test_registration_driver.py` (frozen-cell equivalence), `test_recipe_preview.py`, `test_signal_parity.py`, `test_quantify.py`, `test_compartments.py`, `test_som_split.py`, `test_marimo_render.py`/`_bootstrap.py`/`_viewers.py`/`_tuners.py`, `test_sync_to_projects.py`, `test_notebook_manifest.py`, `test_notebook_launch.py`, `test_legacy_reference_frozen.py`; updates to `test_project_init.py` (fixture filename from the manifest; assert `copied`), `test_mcp_tools.py`, `test_dashboard.py`, `test_deps.py`. Heavy deps behind `skipif`; GPU behind `check_gpu()`.

HPC acceptance checklist (user, real data, GPU node; each item gates the matching retirement flip):

- A1 `kintsugi notebook launch <repo> --notebook 2.5`; run the printed `ssh -J` line; browser opens with the token URL; no cell runs on load; settings menu shows runtime options as overridden. Record whether `-J` to the compute node works.
- A2 S2: start the session-survival spike, kill the tunnel for 60 s, reconnect; log has no gap >10 s.
- A3 S1 + S8: CuPy memory returns to baseline after 5 gated re-runs; editing `notebooks/Kio.py` on `/blue` marks the dependent cell stale within 5 s with idle CPU <5%.
- A4 2.5 end to end on a real deconvolved cycle (1904CC1-1L CD31); output file set matches the previous run; features CSV rows within 1%.
- A5 5_KRONOS: model load + `embed_project` on one small project.
- A6 nb1 on `test_data/mini_project`: BaSiC/stitch/decon/EDF via buttons; `meta/result.pkl` positions equal pre-migration.
- A7 nb3_QC: tune CD3e/Blank (slider re-render <1 s at 1024), apply one, then Snakemake signal isolation with `recipe_dir=configs/recipes`: md5 of `CD3e.tif` identical; `_param.txt` loads via `load_legacy_recipes`; `kintsugi mcp tools` lists `generate_marimo_cell`.
- A8 After PR8: `kintsugi workflow run test_data/mini_project` — stitched/EDF TIFF md5s identical to a pre-migration run; `boundary_qc.json` still written; re-register one dataset: sentinel `method=rigid+nonrigid`. Redeploy scripts to all projects and confirm 0 stale.
- A9 nb2 interactively on mini_project through all gates; then `kintsugi workflow run` on the same project reports every channel SKIPPED; repeat one stitch with a 2-minute tunnel drop.
- A10 nb4 on 1904CC1-1L: segmentation → quantify → SOM → phenotype; `cell_stats.csv` header identical to the old file; h5ad written; napari launcher behaviour recorded.
- A11 `kintsugi init --force <project>`: `.py` notebooks present, `.ipynb` and stale `Kview/` under `notebooks/legacy_ipynb/`, `.vscode/settings.json` has no `jupyter.*` key; `sync_to_projects.py --dry-run` lists no notebook files without `--notebooks`; with `--notebooks` a `.bak` is written for an edited notebook.
- A12 After PR14: `conda env update -f envs/env-hpc.yml --prune`; `kintsugi check --strict` passes; `import stackview` fails (expected); all six notebooks open with `marimo edit`; `marimo check --strict notebooks/[0-9]*.py` clean.

## Risks

| Risk | Mitigation |
|---|---|
| PR8 rewrites the production Snakemake drivers used on 25+ projects; sandbox proves only CPU equivalence | Hash-frozen legacy cells + `np.array_equal` tests before merge; A8 md5 parity is a hard gate; boundary QC and skip logic untouched |
| CuPy memory or threaded device affinity misbehaves under reactive re-runs (historical #806) | S1 spike; fallback: GPU stages become submit-and-poll launchers over `kintsugi workflow run` |
| In-flight execution across a websocket drop unverified | S2 spike; never pass `--session-ttl`/`--timeout`; docs steer multi-hour runs to script mode / Snakemake |
| Script-mode smokes stop at the first form/button | `bootstrap.gate()`/`resolve_params()` honour `mode()=='script'` + `--steps`/`--params`; unit-tested before any port |
| `preview_recipe` differs numerically from the legacy `Kutils` path users tuned with | `test_signal_parity.py` reports max abs diff before PR7; recipes stay `_param.txt` so old outputs can be re-run either way |
| 10k×10k images make slider drags slow or exceed `output_max_bytes` | S3 sets `max_size=1024` defaults; dask-safe `downsample_for_display`; header raises the caps |
| Existing project copies keep stale `Kview/` and `.ipynb`; `import Kview` fails once ipywidgets is gone | `setup_notebooks` archives on `--force`; A11 verified before PR14 |
| `auto_reload="lazy"` watcher slow on `/blue` | S8; fallback `auto_reload="off"` in project copies, restart `marimo edit` |
| Pinned header settings surprise users | Documented escape hatch (edit the header) in `notebooks/CLAUDE.md` and `docs/marimo/MIGRATION.md` |
| Unattended sessions stall on HPC-only gates | Spike gates are soft (ship fallback behind a flag, record `UNKNOWN`); only A8 and per-notebook retirement are hard user gates |

## Carried debt (recorded in `docs/marimo/DEBT.md`)

Wheel does not ship `notebooks/`; `Kutils.py`/`Kio.py`/`Kanalysis.py` stay under `notebooks/` on `sys.path` via `add_module_paths`; `kintsugi.segment` classical module still unused by nb4 (Kseg/InstanSeg is the path); boundary-QC fallback and quick-metrics z-plane repair remain one-sided (workflow-only vs notebook-only) and are surfaced as explicit switches.
