# Changelog

All notable changes to KINTSUGI will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **evals**: the agentic tool-selection eval is provider-pluggable.
  `evals/mcp/providers.py` (standard library only) speaks the Anthropic Messages API
  (`--provider anthropic`) and the OpenAI Chat Completions API (`--provider openai`),
  which also reaches other platforms and local servers such as Ollama and vLLM through
  `OPENAI_BASE_URL` or `--base-url`. The report records the provider, model and
  endpoint host. There is no default provider or model: the harness takes
  `--provider` and `--model`, asks for them in a terminal, and otherwise exits with a
  usage error. The `mcp-eval` workflow runs the pass only when a provider and model are
  chosen (dispatch inputs `provider`, `model`, `base_url`; for the schedule, variables
  `EVAL_PROVIDER`, `EVAL_MODEL`, `EVAL_BASE_URL`) and writes
  `mcp_eval_report_<provider>.json`.
- **mcp**: persistent per-project workflow state. `load_channel`, the processing tools,
  `save_processed`, `propagate_parameters` and `approve_and_learn` record channel status
  and steps in `<project>/.kintsugi/workflow_state.json` (with an `audit_log` in
  `workflow_audit.db`), which survives server restarts. New tools `get_workflow_status`
  (typed output) and `set_review_decision`, resource `kintsugi://project/workflow`, CLI
  `kintsugi mcp status`; `get_session_state` gains a `workflow` summary and
  `get_processing_history` a `source` key, falling back to persisted steps for a channel
  that is no longer loaded. Derived images inherit `project`, `cycle` and `root` from
  the image they were computed from. 34 tools.
- **mcp**: client registry `kintsugi.mcp.clients` (11 clients), `kintsugi mcp config
  --client NAME` (repeatable, or `all`), `kintsugi mcp clients`, `kintsugi init
  --mcp-client NAME` (repeatable, `all` or `none`). No client is configured unless
  chosen: without `--client`, `mcp config` asks in a terminal and otherwise exits 2;
  `init` asks (default: none) or creates the project without MCP files, and `init
  --force` refreshes only the clients already configured. KINTSUGI works without the
  `mcp` package and without any assistant. Existing client files are merged and only
  the `kintsugi` entry is replaced, with `--force`. Generated projects get an `AGENTS.md` (an existing
  one only has its tool list refreshed); the repository gains a root `AGENTS.md`, and
  `kintsugi mcp tools --check` covers it in CI.
- **mcp**: every tool is advertised with a title and MCP behaviour annotations
  (readOnly / destructive / idempotent / openWorld); the server publishes
  instructions, a version and a website URL.
- **mcp**: resources `kintsugi://session/images`,
  `kintsugi://session/history/{channel}`, `kintsugi://learning/stats` and
  `kintsugi://project/channels`, and the `signal_isolation_workflow` prompt.
- **mcp**: progress notifications on long-running tools (weighted
  `subtract_blank`, `denoise_advanced`, `optimize_parameters`, `cluster_channels`,
  `propagate_parameters`); `get_thumbnail` returns the preview as an image content
  block.
- **mcp**: typed output schemas for `get_image_stats`, `assess_quality`,
  `suggest_with_learning` and `get_session_state`; structured output on all tools.
- **mcp**: new tools `get_session_state` and `unload_channel`; 25 previously hidden
  parameters are now advertised (including `subtract_blank` `method="weighted"`).
- **cli**: `kintsugi mcp tools --markdown`, `--write FILE` and `--check FILE`
  (backed by the new `kintsugi.mcp.docs_sync`) generate and verify the tool
  reference embedded in `docs/mcp.md`, `docs/cli.md`, `README.md` and `CLAUDE.md`;
  `kintsugi mcp config --force` overwrites an existing `kintsugi` entry.
- **learning**: `LEARNING_OPERATIONS` and `HISTORY_TO_LEARNING_OP` constants used by
  the MCP learning tools to validate operation names; recommendations report
  `approved_records` and `bootstrap_records`.
- **tests**: `test_mcp_protocol.py` (registry parity, metadata on the wire,
  schema/signature parity, validation errors as `isError`, structured round trips),
  `test_mcp_config.py`, `test_mcp_stdio.py` (raw JSON-RPC over a real subprocess)
  and `test_mcp_docs_sync.py`; `KINTSUGI_REQUIRE_MCP=1` turns MCP skips into
  failures in CI.
- **evals**: the structural pass reads tools through the client and, under
  `--strict`, also fails on missing titles, missing annotations and untargeted
  tools; 12 new tasks (29 in total) cover every tool.
- **ci**: the test job installs the `mcp` extra; the release workflow
  smoke-tests the built wheel over stdio; the MCP eval workflow also triggers on
  pin and CLI changes.
- **docs**: `docs/mcp.md` is the canonical MCP page; the other pages link to it.
- **mcp**: `compare_channels`, `get_intensity_profile` and `batch_assess_quality`
  (implemented and tested but never registered) are now on the tool surface with
  titles, annotations and eval tasks; 32 tools in total.
- **project**: `kintsugi init` and `kintsugi mcp config` add
  `.claude/settings.local.json` (a per-user file) to the project's `.gitignore`.
- **docs**: `CITATION.cff` with the STAR Protocols article as the preferred
  citation; the README BibTeX entry now carries the full author list and DOI.

### Changed
- **docs**: the documentation, notebooks, docstrings and package metadata describe the
  MCP server as optional and usable with any MCP client; client names appear only
  where a sentence is about that client. `docs/mcp.md` is titled "MCP server" and has
  a "Connect your client" section whose table is generated from the client registry
  (new `mcp-tools:clients` block, checked like the tool tables) plus setup for any
  other client (stdio entry, MCP Inspector, Python SDK). `pyproject.toml` and
  `CITATION.cff` no longer describe KINTSUGI as "for Claude Code" (keywords
  `mcp-server`, `ai-assistant`).
- **deprecation**: the Notebook 3/5 notices list the Python API first and the MCP
  route as optional for any client; `DEPRECATION_MAP` points `ini_params*` to
  `kintsugi.signal.subtract_autofluorescence` / `subtract_autofluorescence_dask`
  instead of an MCP tool function.
- **docs**: guidance for AI coding assistants lives in `AGENTS.md` files (repository
  root, `src/kintsugi/`, `src/kintsugi/signal/`, `workflow/`, `notebooks/`), read by
  any assistant; each `CLAUDE.md` is now the one line `@AGENTS.md`. The generated tool
  list is checked in `AGENTS.md` only.
- **BREAKING (evals)**: `evals/mcp/harness.py --agentic` requires `--provider` and
  `--model` (no default model, no provider inferred from API keys) and no longer uses
  the `anthropic` SDK; the `mcp-eval` workflow has no default model and runs the
  agentic pass only when a provider and model are configured. `ruff` in CI now also
  checks `evals/`.
- **BREAKING**: `kintsugi.project.write_claude_config` is replaced by
  `write_mcp_config(root, clients, *, force=False, require_project=True)`, where
  `clients` is required (no default client). `kintsugi init` no longer writes Claude
  Code files by default: it and `kintsugi mcp config` write the configuration of the
  client(s) the user chooses among those that keep it in the project (Claude Code `.mcp.json` +
  `.claude/settings.local.json`, VS Code `.vscode/mcp.json`, Cursor `.cursor/mcp.json`,
  Gemini CLI `.gemini/settings.json`, Antigravity `.agents/mcp_config.json`, Codex
  `.codex/config.toml`, Zed `.zed/settings.json`), an `AGENTS.md`, and `.gitignore`
  lines for the generated files, and print snippets for Claude Desktop, Windsurf, Cline
  and JetBrains. `kintsugi mcp config` refuses a directory without
  `kintsugi_project.json` unless `--force` is given.
- **BREAKING**: project state moved from `<project>/.claude_workflow/` to
  `<project>/.kintsugi/` (`kintsugi.assist.state_dir.project_state_dir`). An existing
  `.claude_workflow/` is read as is and renamed to `.kintsugi/` on the first write, so
  learned parameters carry over; nothing writes to `.claude_workflow/` any more.
  `kintsugi init` adds `.kintsugi/` to the project `.gitignore`.
- **BREAKING**: the MCP server no longer reads `CLAUDE_PROJECT_DIR`. Its default
  project comes from `kintsugi mcp start --project-dir PATH` (the directory must contain
  `kintsugi_project.json`), else `KINTSUGI_PROJECT_DIR`, else the working
  directory only if it is a KINTSUGI project. Without one, `load_channel`,
  `list_channels`, `save_processed` and the `kintsugi://project/channels` resource
  return an error asking for `project_path`, and `save_processed` never creates
  directories under an arbitrary working directory. Learning tools default to the
  configured project's database plus the user-wide one (previously user-wide
  only), and `get_session_state` reports `project_root: null` when no project is
  configured. `kintsugi mcp config` / `kintsugi init` now write
  `"args": ["mcp", "start", "--project-dir", "<project>"]`; regenerate existing
  entries with `kintsugi mcp config PROJECT --force`.
- **BREAKING**: the package `kintsugi.claude` is now `kintsugi.assist` (parameter
  learning, Jupyter tuners, workflow state), and the optional dependency group /
  pip extra `claude` is now `mcp` (`pip install -e ".[mcp]"` or `kintsugi install
  mcp`). The MCP server works with any MCP client, so neither name should point at
  one vendor. There are no aliases: `import kintsugi.claude` fails and `kintsugi
  install claude` exits with a pointer to the new name. The package's exports now
  resolve lazily, so importing the learning engine no longer imports ipywidgets or
  stackview; the unused `httpx` dependency left the extra.
- **notebooks**: the parameter-learning cells of `3_Signal_Isolation_QC.ipynb`
  imported `ParameterLearningEngine` from a module that never defined it, read
  result keys the engine does not return (`found`, `sample_count`) and recorded
  cleaning parameters under an unknown operation; they now import from
  `kintsugi.assist`, check `status`, read `based_on_records` and record
  `clean_background`.
- **mcp**: the server now requires the 2.x Python SDK (`mcp>=2.0,<3`,
  `mcp.server.mcpserver.MCPServer`) through the import gate `kintsugi.mcp._sdk`,
  which distinguishes "mcp is not installed" from "an unsupported mcp 1.x is
  installed". The 1.x `mcp.server.fastmcp` API is no longer supported. Install
  with `pip install -e ".[mcp]"` or `kintsugi install mcp`.
- **mcp**: generated Claude Code configuration changed shape. `.mcp.json` holds
  `{"mcpServers": {"kintsugi": {"type": "stdio", "command": "<absolute path to
  kintsugi>", "args": ["mcp", "start"]}}}` with no `cwd` (Claude Code starts
  project servers in the project directory and sets `CLAUDE_PROJECT_DIR`), and
  `.claude/settings.local.json` approves only that server
  (`{"enabledMcpjsonServers": ["kintsugi"]}`, no `enableAllProjectMcpServers`).
  `kintsugi mcp config` delegates to `project.write_claude_config()`.
- **mcp**: input schemas are derived from the handler signatures
  (`Annotated[..., Field(...)]`, `Literal[...]`) and enforced — a bad enum or
  out-of-range value is rejected with `isError` before the handler runs — and
  in-band errors are reported as `isError`. The hand-written schema copy and the
  private `server._tool_manager` override are gone; `tool_specs.py` is a thin
  registry of name, handler, title, description and annotations.
- **mcp**: `save_processed` and `propagate_parameters` are marked destructive
  (they overwrite files); `get_thumbnail.max_size` is bounded at 4096; two
  never-read parameters were removed (`subtract_blank.record_success` /
  `quality_threshold`, `assess_quality.model_path` / `tile_size`).
- **mcp**: `approve_and_learn` maps session history names (`subtract_blank`) onto
  learning operations (`blank_subtraction`) and rejects unknown keys in-band;
  learning tools take a `Literal` operation.
- **mcp**: `project_path` arguments default to `CLAUDE_PROJECT_DIR` or the
  server's working directory.
- **cli**: `kintsugi mcp start` writes its banner and errors to stderr; stdout
  belongs to the JSON-RPC stream. `KintsugiProject.load()` logs instead of
  printing, since it runs inside the server.
- **learning**: unapproved (bootstrap) records are included in recommendations by
  default at `BOOTSTRAP_WEIGHT` (0.6) of an approved record, and confidence is
  scaled by `0.5 + 0.5 * approved_fraction` so bootstrap-only evidence stays below
  the "learned" tier.
- **deps**: `pydantic` is a core dependency; the `mcp` install group matches
  the pyproject extra (adds `optuna`, `smo`); `mcp>=2.0,<3` in `pyproject.toml`,
  `deps.py` and `envs/env-hpc.yml`.
- **docs**: install text corrected everywhere — `pip install kintsugi[mcp]`
  installed an unrelated PyPI package; it now reads `pip install -e ".[mcp]"`.
  Tool tables are generated from the registry; stale file references
  (`kintsugi-dev.code-workspace`, `test_data/mini_project`, `scripts/release.py`),
  install-group counts, test counts and the denoise method lists were corrected.
- **mcp**: the in-memory session (`_loaded_images`, `_processing_history` and
  their accessors) lives in `kintsugi.mcp.session`; the tool modules import it from
  there instead of from each other, and the three private copies of `_get_image`
  are gone. `signal_isolation` re-exports the names for existing imports.
- **mcp**: `suggest_parameters` and `suggest_with_learning` share one rule set
  (`kintsugi.mcp.heuristics`) instead of two diverging copies; the blank channel is
  sampled at the same stride as the signal, and a blank channel that is named but
  not loaded is an error rather than silently ignored.
- **learning**: `ParameterLearningEngine.recommend_parameters_many()` answers all
  operations with one query per database (a window function ranks the top 50 per
  operation) and `record_parameters_many()` writes every operation of an approval
  in one transaction per database. `suggest_with_learning` and `approve_and_learn`
  use them; `approve_and_learn` now reports each `record_id`.

### Removed
- The repository no longer ships MCP client configuration: the root `.mcp.json` and
  `.claude/commands/` are untracked and gitignored, together with the other clients'
  directories, so each contributor keeps their own. A clone that pulls this change
  loses its tracked copies; restore them as local, ignored files with
  `git restore --source=a0645c3 --worktree -- .mcp.json .claude/commands`.
- **mcp**: the unregistered, untested `set_project_context` handler and the
  `_basic_parameter_suggestion` fallback (its only trigger, `kintsugi.signal`
  failing to import, cannot happen with the core dependencies).
- **docs**: four zero-byte files (`docs/-r`, `docs/2.0.0`, `docs/2024.0.0`,
  `docs/7.0.0`) and the unreferenced `docs/CD31_curtain.gif`; the stale Java/JVM
  instructions in `TROUBLESHOOTING.md`.

### Fixed
- **release**: `cz bump` failed (`TypeError: unhashable type: 'Table'`) because
  `[tool.commitizen.changelog_format]` was a table, and `major_version_zero = true`
  would have kept a breaking change below 2.0.0; `changelog_file` is now a plain key
  and `major_version_zero` is removed, so `cz bump --dry-run` reports 1.2.3 → 2.0.0.
- **mcp**: numpy values in a tool result are converted before serialisation; a
  default (global) `subtract_blank` call returned a `numpy.bool_` and failed with
  "Error executing tool" over the protocol.
- **mcp**: `assess_quality` and `predict_parameters` (read-only tools) no longer append to
  the processing history; `assess_quality` stores its score in the image metadata.
  `estimate_background` records its history on the corrected output, so saving it
  writes a non-empty history file. `save_processed` serialises numpy values in the
  history file.
- **assist**: `WorkflowState` creates nothing until the first change, merges concurrent
  writers under a lock with atomic writes, no longer replaces an existing channel in
  `add_channel`, and reports unknown channels instead of ignoring them.
- **cli**: `kintsugi mcp config --print-only` printed `args` without `--project-dir`.
- **learning**: constructing `ParameterLearningEngine` no longer creates directories or
  databases; they are created on the first recorded parameters, so read-only MCP tools
  (`get_learned_parameters`, `suggest_with_learning`, `get_learning_statistics`, the
  `kintsugi://learning/stats` resource) leave no files behind. A database written before
  the `algorithm_version` column existed is upgraded when it is first read.
- **learning**: a record mirrored into both the project and the user-wide
  database is counted once in recommendations (it was counted twice, inflating
  confidence).
- **mcp**: the server imported `mcp.server.fastmcp.FastMCP`, which does not exist
  in the 2.x SDK that `pip install mcp` resolves, while the project pinned
  `mcp<2`; it now targets `mcp>=2.0,<3` and starts.
- **learning**: marker and tissue aliases matched by substring, so CD31, CD34 and
  CD38 collapsed onto `cd3`, CD45 (and "smooth muscle actin") onto `cd4`, and the
  tissue "T" became "tonsil", silently merging unrelated records. Matching is now
  exact on the cleaned name or on a complete token sequence, longest alias first.
- **learning**: records written by `kintsugi mcp pretrain` (`user_approved=False`)
  never reached `get_learned_parameters` or `suggest_with_learning` because
  `recommend_parameters()` defaulted to `approved_only=True`.
- **mcp**: a single cached learning engine was reused whenever `project_path` was
  omitted, so a call without a project after a call with one kept reading and
  writing the previous project's database; engines are now keyed by the resolved
  project path.
- **mcp**: `generate_jupyter_cell` interpolated values with `str()` (a quote in a
  channel name broke the code), offered a `full_pipeline` template that imported
  a module that does not exist, and told the tuner to look up MCP session channels
  that live in another process. Cells now interpolate with `repr()`, the
  `full_pipeline` branch is gone, and the cell loads the channel from its file on
  disk (`param_tuner._load_image()` accepts a path).
- **mcp**: `list_channels` joined the `cycle` argument onto the search directory
  unchecked, so `cycle="../.."` listed files outside the project; path components
  are now rejected.
- **mcp**: `load_channel` took the first unsorted glob hit for `*CD3*`, so whether
  `CD3` loaded `CD3.tif` or `CD31.tif` depended on filesystem order; it now ranks
  matches (exact stem, whole token, substring) and reports a tie as an in-band
  `candidates` error. `list_channels` output is sorted for the same reason.
- **mcp**: `denoise_advanced` passed any existing `model_path` to `torch.load` and
  `predict_parameters` passed any path to `joblib.load`; checkpoints are now
  confined and loaded safely (see Security).
- **mcp**: the learned-parameter merge in `subtract_blank` never applied — the
  auto-suggest path imported `ParameterLearningEngine` from
  `kintsugi.mcp.tools.learning` (which only imports it lazily) and tested
  `recommendation["found"]` while the engine returns `status == "success"`. Both
  the MCP tool and `signal.subtractor` now import from
  `kintsugi.assist.parameter_learning` and test the status key.
- **mcp**: `cluster_channels` never stored the loaded images in the session or set
  `metadata.cluster_id`, which is exactly what `propagate_parameters` filters on,
  so the pair always failed with "No channels found". Clustered channels are now
  kept in the session with `cluster_id`, `source_path` and `wavelength_group`, and
  the result lists `loaded_channels`.
- **mcp**: `get_intensity_profile` with `position=1.0` indexed one row/column past
  the image; the last row/column is now returned.

### Security
- **mcp**: model checkpoint paths given to `denoise_advanced` (`model_path`) and
  `predict_parameters` must lie under the project, the Claude Code project
  (`CLAUDE_PROJECT_DIR`) or `~/.kintsugi/models`
  (`path_safety.allowed_model_roots` / `ensure_within_any`) and carry the expected
  extension; N2V and CARE checkpoints are loaded with `torch.load(...,
  weights_only=True)` (their config is saved as a plain dict so that works), and
  the Random Forest predictor validates the loaded object's shape.
- **mcp**: `list_channels` rejects path components in `cycle`, and
  `cluster_channels` confines the files it loads to the project (`rglob` follows
  symlinks, so the resolved path is checked).

## [1.2.3] - 2025-12-12

### Changed
- Merge branch 'main' of https://github.com/smith6jt-cop/KINTSUGI
- Refactor MCP configuration process and enhance project setup

## [1.2.2] - 2025-12-12

### Changed
- Fix badge formatting in README.md

## [1.2.1] - 2025-12-12

### Changed
- Merge branch 'main' of https://github.com/smith6jt-cop/KINTSUGI
- Update project setup and configuration files

## [1.2.0] - 2025-12-12

### Added
- **release**: add automatic release versioning system
- integrate Claude Code with KINTSUGI for AI-assisted image processing; add MCP server commands and update documentation
- add KINTSUGI project management framework

### Changed
- Merge pull request #18 from smith6jt-cop/claude/auto-release-versioning-012gaU2Hg5kvPfBGLnxffbev
- Add configuration setup for Claude Code and VS Code
- Clean up etc
- Restore install command for optional dependencies
- Remove settings.local.json from git tracking
- Update Skills Registry and fix slash commands
- Add Autofluorescence Subtractor and Utility Functions
- Add new notebook for automated channel quality assessment and update migration guide
- Fix linting errors: ruff and black formatting compliance
- Fix PyPI collision: use direct package names in install commands
- Add segmentation module with SAM, classical methods, and post-processing utilities
- Streamline installation: base env + optional feature groups
- Deprecate Java, Maven, and PyImageJ/CLIJ2 dependencies; transition to pure Python implementations for all processing
- Add GPU acceleration for phase correlation matrix computation and optimize image processing
- Optimize SVD initialization in KCorrectGPU for performance
- Add OME-Zarr I/O utilities for KINTSUGI pipeline
- Merge pull request #17 from smith6jt-cop/claude/fix-cupy-fallback-01WEBXNQsFDLDDE2wg64Mfpz
- Add cuda-nvrtc package for CuPy JIT compilation support
- Merge pull request #16 from smith6jt-cop/claude/fix-cupy-fallback-01WEBXNQsFDLDDE2wg64Mfpz
- Auto-add conda CUDA bin to PATH on Windows for stitching
- Pin numpy<2 after cupy-cuda12x install to prevent version conflicts
- Merge pull request #15 from smith6jt-cop/claude/fix-cupy-fallback-01WEBXNQsFDLDDE2wg64Mfpz
- Fix CuPy on Windows: use pip cupy-cuda12x instead of conda cupy
- Add CuPy and CUDA runtime to base environments
- Remove CPU fallback when use_gpu=True for stitching
- Merge pull request #14 from smith6jt-cop/claude/resolve-branch-conflicts-01EeMNK6c4qPysVzFXvSiNcV
- Use correct libvips Windows builds repository
- Use GitHub API to find correct libvips Windows asset
- Fix Windows libvips installation in CI
- Add ome-types and install package in docs build
- Add beautifulsoup4 to docs dependencies
- Add matplotlib and cycler to docs dependencies
- Add colour-science to docs requirements
- Add base environment files and Read the Docs configuration
- Merge remote-tracking branch 'origin/claude/evaluate-kstitch-alternatives-01HEYPnuanYGqDdvw7LF8t6Q' into claude/resolve-branch-conflicts-01EeMNK6c4qPysVzFXvSiNcV
- Add weightedstats to autodoc_mock_imports
- Add pqdm to autodoc_mock_imports
- Fix ruff linting issues in edf.py
- Add matplotlib and visualization libs to autodoc_mock_imports
- Add shapely and geospatial deps to autodoc_mock_imports
- Add SimpleITK and other heavy deps to autodoc_mock_imports
- Add autodoc_mock_imports to prevent doc build import errors
- Add scikit-image and other deps to docs requirements
- Fix documentation build CI failures
- Fix GPU compatibility issues in Kstitch
- Add Kstitch alternatives evaluation document
- Add PyImageJ/CLIJ2 vs Python EDF evaluation and pure Python implementation
- Update documentation for improved clarity and consistency in command syntax
- Merge branch 'claude/fix-docs-build-01NbA31uDHtnALSW5tj9HcWU' of https://github.com/smith6jt-cop/KINTSUGI
- Fix documentation build warnings
- Merge branch 'main' of https://github.com/smith6jt-cop/KINTSUGI
- Remove outdated environment configuration files for KINTSUGI_new and KINTSUGI_streamlined, streamlining dependencies and channels for improved package management and performance.
- Merge pull request #10 from smith6jt-cop/claude/fix-docs-build-01NbA31uDHtnALSW5tj9HcWU
- Add complete Sphinx documentation with Windows support and GitHub Pages
- Merge pull request #9 from smith6jt-cop/claude/remove-matlab-deconvolution-01S3U9sduDVJtHUcVBiC1zMh
- Replace MATLAB deconvolution with pure Python implementation
- Merge pull request #8 from smith6jt-cop/claude/analyze-project-structure-01CUT7TkESRjzbTJSR8En179
- Make cupy import optional in Kstitch stitching module
- Fix Windows CI exit code from previous choco failure
- Add GitHub releases fallback for Windows libvips installation
- Make Windows libvips CI installation more robust
- Add Windows libvips installation to CI workflow
- Fix recursion error in lazy module loading
- Move packages with native dependencies to optional checks
- Fix linting errors and formatting issues
- Implement comprehensive package management improvements
- Add comprehensive dependency analysis for KINTSUGI project
- Merge pull request #7 from smith6jt-cop/claude/document-kintsugi-folder-012vCftHqEtCxWEUT6JiMafn
- Rename kintsugi module to dl_refinement and move to notebooks folder
- Merge pull request #6 from smith6jt-cop/claude/merge-branches-to-main-01CA88zSLyXHMqRGHhHdJYyF
- Merge Development branch: Vessel Analysis, enhanced Cropper layout, and code refactoring
- Merge pull request #4 from smith6jt-cop/copilot/fix-7f321af5-455e-431c-a392-bed5b8680d99
- Merge pull request #5 from smith6jt-cop/claude/repo-familiarization-013erQkhT9C5EKnA5javu1ga
- Add PR description documentation
- Add deep learning channel refinement module
- Revise README for citation and installation clarity
- Refactor code structure for improved readability and maintainability
- Testing improvements and adding Vessel Analysis
- Add comprehensive implementation summary documentation
- Complete registration workflow implementation for PANCREAS and THYMUS images
- Initial plan
- Update README.md
- Add files via upload
- Enhance _Cropper layout by adding custom width and description width to range sliders
- Merge branch 'main' of https://github.com/smith6jt-cop/KINTSUGI
- Update instructions for signal subtraction parameters in Signal Isolation notebook
- Add KDecon v4 log output and refactor SimpleElastix groupwise registration
- Update env.yml
- Update environment configuration and notebooks
- Merge branch 'main' of https://github.com/smith6jt-cop/KINTSUGI
- For revisions
- Update README.md
- Update README.md
- Merge clustering
- Refactor code structure for improved readability and maintainability
- Refactor Kutils and Kview2: enhance array handling and improve layout for better performance
- Refactor Kutils and Kview2: streamline array handling and optimize layout for enhanced performance
- Refactor Kutils and Kview2: update array copy method and simplify crop view layout for improved clarity
- ``` Refactor code structure for improved readability and maintainability ```
- Refactor clean function parameters in Kutils: rename smooth_threshold to smooth_thresh for consistency
- Refactor clean function parameters in Kutils: rename background_threshold to backgrnd_thresh for consistency
- Refactor layout dimensions in Kview2: enhance crop view and grid sizes for better UI alignment
- Refactor CPU detection and processing in MicroRigidRegistrar: improve batch processing and error handling
- Update README.md
- Update README.md

## [1.1.0] - 2024-12-01

### Added
- Initial public release
- Multi-cycle immunofluorescence image registration
- 2D/3D GPU/CPU illumination correction
- Image stitching with GPU acceleration
- Deconvolution support
- Extended depth of focus processing
- Autofluorescence removal
- Segmentation workflows
- Clustering and spatial analysis
- MCP server for Claude Code integration
- Comprehensive documentation

### Changed
- Updated to numpy<2.0 for compatibility
- Improved dependency management with optional groups

[Unreleased]: https://github.com/smith6jt-cop/KINTSUGI/compare/v1.2.3...HEAD
[1.1.0]: https://github.com/smith6jt-cop/KINTSUGI/releases/tag/v1.1.0
[1.2.0]: https://github.com/smith6jt-cop/KINTSUGI/compare/v1.1.0...v1.2.0
[1.2.1]: https://github.com/smith6jt-cop/KINTSUGI/compare/v1.1.0...v1.2.1
[1.2.2]: https://github.com/smith6jt-cop/KINTSUGI/compare/v1.1.0...v1.2.2
[1.2.3]: https://github.com/smith6jt-cop/KINTSUGI/compare/v1.1.0...v1.2.3
