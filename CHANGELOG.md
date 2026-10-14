# Changelog

All notable changes to KINTSUGI will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
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
  tools; 9 new tasks (26 in total) cover every tool.
- **ci**: the test job installs the `claude` extra; the release workflow
  smoke-tests the built wheel over stdio; the MCP eval workflow also triggers on
  pin and CLI changes.
- **docs**: `docs/mcp.md` is the canonical MCP page; the other pages link to it.

### Changed
- **mcp**: the server now requires the 2.x Python SDK (`mcp>=2.0,<3`,
  `mcp.server.mcpserver.MCPServer`) through the import gate `kintsugi.mcp._sdk`,
  which distinguishes "mcp is not installed" from "an unsupported mcp 1.x is
  installed". The 1.x `mcp.server.fastmcp` API is no longer supported. Install
  with `pip install -e ".[claude]"` or `kintsugi install claude`.
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
- **deps**: `pydantic` is a core dependency; the `claude` install group matches
  the pyproject extra (adds `optuna`, `smo`); `mcp>=2.0,<3` in `pyproject.toml`,
  `deps.py` and `envs/env-hpc.yml`.
- **docs**: install text corrected everywhere — `pip install kintsugi[claude]`
  installed an unrelated PyPI package; it now reads `pip install -e ".[claude]"`.
  Tool tables are generated from the registry; stale file references
  (`kintsugi-dev.code-workspace`, `test_data/mini_project`, `scripts/release.py`),
  install-group counts, test counts and the denoise method lists were corrected.

### Fixed
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
- **mcp**: `denoise_advanced` passed any existing `model_path` to `torch.load` and
  `predict_parameters` passed any path to `joblib.load`; checkpoints are now
  confined and loaded safely (see Security).
- **mcp**: the learned-parameter merge in `subtract_blank` never applied — the
  auto-suggest path imported `ParameterLearningEngine` from
  `kintsugi.mcp.tools.learning` (which only imports it lazily) and tested
  `recommendation["found"]` while the engine returns `status == "success"`. Both
  the MCP tool and `signal.subtractor` now import from
  `kintsugi.claude.parameter_learning` and test the status key.
- **mcp**: `cluster_channels` never stored the loaded images in the session or set
  `metadata.cluster_id`, which is exactly what `propagate_parameters` filters on,
  so the pair always failed with "No channels found". Clustered channels are now
  kept in the session with `cluster_id`, `source_path` and `wavelength_group`, and
  the result lists `loaded_channels`.

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
