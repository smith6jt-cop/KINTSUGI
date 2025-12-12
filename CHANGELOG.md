# Changelog

All notable changes to KINTSUGI will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Claude Code configuration setup for VS Code integration
- Skills Registry for AI-assisted workflows
- Autofluorescence Subtractor and utility functions
- Automated channel quality assessment notebook
- Claude Code MCP integration for AI-assisted image processing
- Segmentation module with SAM, classical methods, and post-processing utilities
- KINTSUGI project management framework
- GPU acceleration for phase correlation matrix computation
- OME-Zarr I/O utilities

### Changed
- Streamlined installation with base env + optional feature groups
- Optimized SVD initialization in KCorrectGPU for performance
- Updated migration guide for notebook transitions

### Deprecated
- Java, Maven, and PyImageJ/CLIJ2 dependencies (transitioned to pure Python)

### Fixed
- Linting errors: ruff and black formatting compliance
- PyPI collision with direct package names in install commands
- CuPy fallback issues on Windows
- libvips installation in CI for Windows

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

[Unreleased]: https://github.com/smith6jt-cop/KINTSUGI/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/smith6jt-cop/KINTSUGI/releases/tag/v1.1.0
