# KINTSUGI Dependency Analysis Report

> **⚠️ DEPRECATION NOTICE (2025):** This document is historical. Java, Maven, FIJI, and CLIJ2
> dependencies have been removed from KINTSUGI. The project now uses pure Python implementations
> (CuPy/NumPy) for all processing including Extended Depth of Focus (EDF). Only libvips and VALIS
> remain as external dependencies.

## Executive Summary

This document analyzes the current dependency structure of the KINTSUGI project, focusing on how external dependencies (libvips, maven, java-jdk, VALIS, Stackview) are integrated and managed. It identifies gaps in the current setup and proposes improvements for proper package management.

---

## 1. Current Project Structure

```
KINTSUGI/
├── env.yml                    # Primary conda environment (580 deps)
├── env_new.yml                # Alternative environment (identical)
├── env_streamlined.yml        # Minimal environment (~150 deps)
├── .env                       # Windows-specific paths
├── notebooks/
│   ├── Kreg/                  # Registration module (16 files)
│   ├── Kview2/                # Visualization module (25 files)
│   ├── Kstitch/               # Stitching module (6 files)
│   ├── dl_refinement/         # Deep learning refinement (4 files)
│   └── *.ipynb                # Workflow notebooks
├── docs/
├── fiji_macros/
└── KDecon_v*.exe/m            # MATLAB deconvolution
```

---

## 2. External Dependencies Analysis

### 2.1 libvips (via pyvips)

| Attribute | Value |
|-----------|-------|
| **Package** | `pyvips==2.2.3` (pip) |
| **Native Library** | libvips 8.16 (from Zenodo) |
| **Purpose** | High-performance image I/O and processing |

**Current Integration:**
- Installed via pip as `pyvips==2.2.3`
- Native library distributed via Zenodo download
- No automated validation of native library presence

**Files Using pyvips:**
- `notebooks/Kreg/registration.py` (line 16)
- `notebooks/Kreg/slide_io.py` (lines 7, 38, 45)
- `notebooks/Kreg/valtils.py` (line 7)

**Configuration Pattern:**
```python
import pyvips
pyvips.cache_set_max(0)  # Disable caching for memory management
```

**Issues Identified:**
1. No system-level dependency check for libvips C library
2. Version mismatch risk between Python bindings and native library
3. Platform-specific installation not documented (Windows vs Linux)

---

### 2.2 Maven & Java-JDK

| Attribute | Maven | Java-JDK |
|-----------|-------|----------|
| **Package** | `maven=3.9.9` (conda) | `openjdk=11.0.25` (conda) |
| **Alternative** | Manual install | `java-jdk21` (Zenodo) |
| **Channel** | conda-forge | conda-forge |

**Current Integration:**
- Conda packages in env.yml
- Alternative: Zenodo pre-compiled binaries
- JPype bridge for Python-Java interop

**Files Using Java:**
- `notebooks/Kreg/slide_io.py`:
  - Line 24: `import jpype`
  - Line 27: `import scyjava`
  - JVM lifecycle management implemented

**JVM Management Pattern:**
```python
import jpype
import scyjava
# JVM started once, managed throughout session
```

**Issues Identified:**
1. Two different Java versions (OpenJDK 11 vs JDK 21)
2. No documentation on which Java version to use when
3. JVM startup overhead not optimized
4. Memory configuration for JVM not documented

---

### 2.3 VALIS (valis-wsi)

| Attribute | Value |
|-----------|-------|
| **Package** | `valis-wsi==1.2.0` (pip) |
| **Purpose** | Whole-slide image registration |
| **Status** | Standard pip package |

**Current Integration:**
- Installed via pip
- Used in registration workflow
- Well-integrated with Kreg module

**References:**
- `notebooks/IMPLEMENTATION_SUMMARY.md` (lines 54, 56)
- `notebooks/README_Registration.md` (lines 69-73)
- `notebooks/Kreg/registration.py` (docstrings, lines 1498-1537)

**Issues Identified:**
1. No indication if this is a modified version (per user query)
2. Version pinned but transitive dependencies not locked
3. No fallback if VALIS unavailable

---

### 2.4 Stackview

| Attribute | Value |
|-----------|-------|
| **Package** | `stackview==0.18.0` (pip) |
| **Purpose** | Multi-dimensional image visualization |
| **Status** | Standard pip package |

**Current Integration:**
- Installed via pip
- Integrated with Napari in Kview2 module

**Files Using Stackview:**
- `notebooks/Kview2/_animate.py`
- `notebooks/Kview2/_context.py`
- `notebooks/Kview2/_crop.py`

**Issues Identified:**
1. No indication if this is a modified version (per user query)
2. Napari plugin compatibility not verified

---

## 3. Dependency Management Assessment

### 3.1 Current State

| Mechanism | Description | Files |
|-----------|-------------|-------|
| Conda Environment | Primary dependency management | env.yml, env_new.yml, env_streamlined.yml |
| Pip (within Conda) | Python-only packages | Listed in env.yml pip section |
| Zenodo Archive | External binaries | README.md installation instructions |
| Manual Installation | Alternative setup | Documented in README.md |

### 3.2 What's Missing

| Component | Status |
|-----------|--------|
| `setup.py` | Not present |
| `pyproject.toml` | Not present |
| `requirements.txt` | Not present (embedded in env.yml) |
| Package metadata | Not defined |
| Entry points | Not configured |
| CI/CD testing | Not implemented |

---

## 4. Identified Issues

### 4.1 Critical Issues

1. **No Python Package Structure**
   - Cannot be installed via `pip install -e .`
   - Module imports require manual `sys.path` manipulation
   - No version tracking for the package itself

2. **External Dependencies Not Validated**
   - libvips native library presence not checked at import
   - Java/Maven availability not verified programmatically
   - No graceful degradation if optional dependencies missing

3. **Environment File Redundancy**
   - `env.yml` and `env_new.yml` are identical
   - Creates confusion about which to use
   - No clear documentation on differences

### 4.2 Moderate Issues

1. **Platform-Specific Paths**
   - `.env` file contains Windows paths
   - No Linux/macOS equivalent provided
   - Cross-platform setup not documented

2. **Version Pinning Inconsistency**
   - Some packages have exact versions
   - Others use flexible version specifiers
   - Risk of dependency conflicts on update

3. **No Dependency Lockfile**
   - `conda-lock` or equivalent not used
   - Reproducible builds not guaranteed
   - Environment drift possible over time

### 4.3 Minor Issues

1. **Documentation Gaps**
   - GPU memory requirements not documented
   - Performance tuning guidance missing
   - Troubleshooting guide incomplete

2. **No Automated Testing**
   - `check_setup.py` and `test_registration_setup.py` exist but are manual
   - No CI/CD pipeline
   - No unit tests for core modules

---

## 5. Proposed Improvements

### 5.1 Package Structure (Priority: High)

Create proper Python package with `pyproject.toml`:

```toml
[build-system]
requires = ["setuptools>=61.0", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "kintsugi"
version = "1.0.0"
description = "Multi-cycle immunofluorescence image registration and analysis"
requires-python = ">=3.10"
dependencies = [
    "numpy>=1.24.0",
    "scipy>=1.10.0",
    "scikit-image>=0.21.0",
    "opencv-python>=4.8.0",
    "pyvips>=2.2.0",
    "valis-wsi>=1.2.0",
    "stackview>=0.18.0",
    # ... additional dependencies
]

[project.optional-dependencies]
gpu = [
    "torch>=2.0.0",
    "torchvision>=0.15.0",
]
java = [
    "JPype1>=1.5.0",
    "scyjava>=1.0.0",
]
viz = [
    "napari>=0.4.19",
    "magicgui>=0.7.0",
]
dev = [
    "pytest>=7.0.0",
    "pytest-cov>=4.0.0",
    "black>=23.0.0",
    "ruff>=0.1.0",
]

[project.scripts]
kintsugi-register = "kintsugi.cli:main"

[tool.setuptools.packages.find]
where = ["notebooks"]
include = ["Kreg*", "Kview2*", "Kstitch*", "dl_refinement*"]
```

### 5.2 Dependency Validation (Priority: High)

Create runtime dependency checker:

```python
# kintsugi/dependency_check.py
def check_dependencies():
    """Validate all external dependencies at startup."""
    issues = []

    # Check pyvips/libvips
    try:
        import pyvips
        version = pyvips.version(0)
    except Exception as e:
        issues.append(f"libvips not found: {e}")

    # Check Java (optional)
    try:
        import jpype
        if not jpype.isJVMStarted():
            jpype.startJVM()
    except Exception as e:
        issues.append(f"Java/JVM not available: {e}")

    return issues
```

### 5.3 Environment Consolidation (Priority: Medium)

1. **Remove redundant env_new.yml**
2. **Create platform-specific variants:**
   - `env-linux.yml`
   - `env-windows.yml`
   - `env-macos.yml`
3. **Add conda-lock for reproducibility:**
   ```bash
   conda-lock -f env.yml -p linux-64 -p win-64
   ```

### 5.4 Installation Script (Priority: Medium)

Create cross-platform installation script:

```bash
#!/bin/bash
# install.sh - KINTSUGI installation script

# Detect platform
OS=$(uname -s)

# Install conda environment
conda env create -f env.yml

# Validate libvips
python -c "import pyvips; print(f'libvips version: {pyvips.version(0)}')"

# Validate Java (optional)
python -c "import jpype; print('Java available')" 2>/dev/null || echo "Java not available (optional)"

# Run setup validation
python notebooks/check_setup.py
```

### 5.5 Documentation Updates (Priority: Low)

1. Add platform-specific installation guides
2. Document GPU memory requirements
3. Create troubleshooting guide for common dependency issues
4. Add performance tuning documentation

---

## 6. Implementation Roadmap

### Phase 1: Package Structure (Week 1)
- [ ] Create `pyproject.toml`
- [ ] Reorganize source files
- [ ] Add `__version__` tracking
- [ ] Create entry points

### Phase 2: Dependency Management (Week 2)
- [ ] Implement dependency validation
- [ ] Consolidate environment files
- [ ] Create conda-lock files
- [ ] Add installation scripts

### Phase 3: Testing & CI (Week 3)
- [ ] Add unit tests for core modules
- [ ] Configure GitHub Actions CI
- [ ] Add dependency validation to CI
- [ ] Create release workflow

### Phase 4: Documentation (Week 4)
- [ ] Update README with new installation
- [ ] Add platform-specific guides
- [ ] Create troubleshooting guide
- [ ] Document GPU requirements

---

## 7. Summary

### Current State
- Dependencies managed via conda environment files
- External binaries distributed via Zenodo
- No proper Python package structure
- Manual setup validation required

### Recommended Actions
1. **Immediate**: Create `pyproject.toml` for package management
2. **Short-term**: Implement dependency validation and consolidate environments
3. **Medium-term**: Add CI/CD pipeline with automated testing
4. **Long-term**: Create cross-platform installation automation

### Risk Assessment
| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Dependency conflicts | Medium | High | Add conda-lock |
| Missing native libs | High | High | Add validation |
| Environment drift | Medium | Medium | Pin all versions |
| Cross-platform issues | High | Medium | Platform-specific docs |

---

*Document generated: 2025-12-03*
*Analysis based on: KINTSUGI repository exploration*
