# Changelog

All notable changes to KINTSUGI will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/smith6jt-cop/KINTSUGI/compare/v1.2.0...HEAD
[1.1.0]: https://github.com/smith6jt-cop/KINTSUGI/releases/tag/v1.1.0
[1.2.0]: https://github.com/smith6jt-cop/KINTSUGI/compare/v1.1.0...v1.2.0
