# KINTSUGI Troubleshooting Guide

This guide covers common issues and their solutions when installing or running KINTSUGI.

## Table of Contents

- [Installation Issues](#installation-issues)
  - [Conda Environment Creation Fails](#conda-environment-creation-fails)
  - [Package Installation Errors](#package-installation-errors)
- [Dependency Issues](#dependency-issues)
  - [libvips Not Found](#libvips-not-found)
  - [Java/JVM Issues](#javajvm-issues)
  - [CUDA/GPU Issues](#cudagpu-issues)
- [Runtime Issues](#runtime-issues)
  - [Import Errors](#import-errors)
  - [Memory Issues](#memory-issues)
  - [Registration Failures](#registration-failures)
- [Platform-Specific Issues](#platform-specific-issues)
  - [Windows](#windows)
  - [Linux](#linux)
  - [macOS](#macos)

---

## Installation Issues

### Conda Environment Creation Fails

**Symptom**: `conda env create` command fails with dependency conflicts.

**Solutions**:

1. **Use libmamba solver** (faster and more reliable):
   ```bash
   conda install -n base conda-libmamba-solver
   conda config --set solver libmamba
   conda env create -f envs/env-linux.yml
   ```

2. **Update conda**:
   ```bash
   conda update -n base conda
   ```

3. **Try the streamlined environment**:
   ```bash
   conda env create -f env_streamlined.yml
   ```

4. **Create minimal environment and add packages**:
   ```bash
   conda create -n KINTSUGI python=3.10
   conda activate KINTSUGI
   pip install -e .
   # Install additional dependencies as needed
   ```

### Package Installation Errors

**Symptom**: `pip install -e .` fails.

**Solutions**:

1. **Ensure setuptools is updated**:
   ```bash
   pip install --upgrade pip setuptools wheel
   ```

2. **Check Python version**:
   ```bash
   python --version  # Should be 3.10+
   ```

3. **Install with verbose output**:
   ```bash
   pip install -e . -v
   ```

---

## Dependency Issues

### libvips Not Found

**Symptom**: Error message like `cannot find library libvips` or `pyvips.GObject.Error`.

#### Windows

1. **Download from Zenodo**:
   - Download PyVips-dev-8.16 from [Zenodo](https://zenodo.org/records/14969214)
   - Extract to KINTSUGI folder

2. **Set PATH manually**:
   ```powershell
   $env:PATH = "C:\Users\[username]\KINTSUGI\vips-dev-8.16\bin;$env:PATH"
   ```

3. **Add to system PATH permanently**:
   - Open System Properties > Environment Variables
   - Add `C:\Users\[username]\KINTSUGI\vips-dev-8.16\bin` to PATH

#### Linux

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y libvips-dev

# Fedora/RHEL
sudo dnf install vips-devel

# Verify installation
vips --version
```

#### macOS

```bash
brew install vips
```

### Verify libvips Installation

```python
import pyvips
print(f"libvips version: {pyvips.version(0)}.{pyvips.version(1)}.{pyvips.version(2)}")
```

---

### Java/JVM Issues

**Symptom**: `java.lang.Exception` or `JVMNotFoundException`.

#### Check Java Installation

```bash
java -version
echo $JAVA_HOME  # Linux/macOS
echo %JAVA_HOME%  # Windows
```

#### Fix Missing Java

1. **Install via conda** (recommended):
   ```bash
   conda install openjdk=11
   ```

2. **Download from Zenodo** (Windows):
   - Download java-jdk21 from Zenodo
   - Extract to KINTSUGI folder

3. **Set JAVA_HOME**:
   ```bash
   # Linux/macOS
   export JAVA_HOME=/path/to/java
   export PATH=$JAVA_HOME/bin:$PATH

   # Windows PowerShell
   $env:JAVA_HOME = "C:\Users\[username]\KINTSUGI\java-jdk21"
   $env:PATH = "$env:JAVA_HOME\bin;$env:PATH"
   ```

#### JVM Already Started Error

**Symptom**: `JVMAlreadyStarted` exception.

This occurs when trying to start JVM multiple times. Solutions:

1. **Restart Python kernel** (in Jupyter)
2. **Single JVM initialization**:
   ```python
   import jpype
   if not jpype.isJVMStarted():
       jpype.startJVM()
   ```

---

### CUDA/GPU Issues

**Symptom**: `torch.cuda.is_available()` returns `False`.

#### Verify GPU Detection

```python
import torch
print(f"CUDA available: {torch.cuda.is_available()}")
print(f"CUDA version: {torch.version.cuda}")
if torch.cuda.is_available():
    print(f"GPU: {torch.cuda.get_device_name(0)}")
```

#### Solutions

1. **Check NVIDIA driver**:
   ```bash
   nvidia-smi
   ```

2. **Reinstall PyTorch with CUDA**:
   ```bash
   pip uninstall torch torchvision
   pip install torch torchvision --index-url https://download.pytorch.org/whl/cu121
   ```

3. **Verify CUDA toolkit**:
   ```bash
   nvcc --version
   ```

4. **Use CPU fallback**:
   - Most KINTSUGI operations work without GPU
   - Set device manually: `device = "cpu"`

---

## Runtime Issues

### Import Errors

**Symptom**: `ModuleNotFoundError` or `ImportError`.

#### General Solutions

1. **Verify environment is activated**:
   ```bash
   conda activate KINTSUGI
   which python  # Should show conda environment path
   ```

2. **Run dependency check**:
   ```bash
   kintsugi check
   ```

3. **Reinstall package**:
   ```bash
   pip install -e . --force-reinstall
   ```

#### Specific Import Errors

**`ImportError: cannot import name 'Valis'`**:
```python
# Correct import path
from Kreg.registration import Valis
# or via kintsugi package
from kintsugi.kreg import Valis
```

**`ModuleNotFoundError: No module named 'Kreg'`**:
```python
import sys
sys.path.insert(0, '/path/to/KINTSUGI/notebooks')
```

---

### Memory Issues

**Symptom**: `MemoryError` or process killed during processing.

#### Solutions

1. **Reduce image dimensions**:
   ```python
   config = {
       "max_image_dim_px": 2048,  # Reduce from default
       "max_processed_image_dim_px": 1024,
   }
   ```

2. **Process images in tiles**:
   - Use tiled processing for large images
   - Reduce batch size

3. **Monitor memory usage**:
   ```python
   import psutil
   print(f"Memory: {psutil.virtual_memory().percent}%")
   ```

4. **Clear pyvips cache**:
   ```python
   import pyvips
   pyvips.cache_set_max(0)
   ```

5. **Use zarr for large datasets**:
   ```python
   import zarr
   # Process in chunks
   ```

---

### Registration Failures

**Symptom**: Registration produces poor results or fails.

#### Check Input Images

1. **Verify image format**:
   ```python
   import tifffile
   img = tifffile.imread('image.tif')
   print(f"Shape: {img.shape}, dtype: {img.dtype}")
   ```

2. **Check image quality**:
   - Ensure images are not corrupted
   - Verify adequate contrast and features

#### Adjust Parameters

1. **Try different registrars**:
   ```python
   config = {
       "micro_rigid_registrar_cls": "RigidRegistrar",  # or "AffineRegistrar"
   }
   ```

2. **Adjust image dimensions**:
   ```python
   config = {
       "max_image_dim_px": 4096,  # Increase for more detail
   }
   ```

3. **Disable non-rigid registration**:
   ```python
   config = {
       "compose_non_rigid": False,
   }
   ```

---

## Platform-Specific Issues

### Windows

#### Long Path Issues

**Symptom**: `FileNotFoundError` with long paths.

**Solution**:
1. Enable long paths in Windows:
   - Run `gpedit.msc`
   - Navigate to: Computer Configuration > Administrative Templates > System > Filesystem
   - Enable "Enable Win32 long paths"

2. Or move KINTSUGI to shorter path (e.g., `C:\KINTSUGI`)

#### DLL Load Failures

**Symptom**: `ImportError: DLL load failed`.

**Solutions**:
1. Install Visual C++ Redistributable
2. Verify all Zenodo dependencies are extracted
3. Add paths to system PATH

---

### Linux

#### Permission Issues

**Symptom**: Permission denied errors.

**Solutions**:
```bash
# Fix permissions
chmod +x scripts/install.sh
chmod -R u+rw KINTSUGI/

# Don't run as root
# If needed, use: sudo chown -R $USER:$USER KINTSUGI/
```

#### Display Issues (Napari)

**Symptom**: Napari doesn't display.

**Solutions**:
```bash
# Check display
echo $DISPLAY

# For headless servers, use Xvfb
Xvfb :99 -screen 0 1024x768x24 &
export DISPLAY=:99
```

---

### macOS

#### Apple Silicon (M1/M2)

**Symptom**: Package incompatibility on ARM.

**Solutions**:
1. Use Rosetta 2 for x86 compatibility:
   ```bash
   arch -x86_64 conda create -n KINTSUGI python=3.10
   ```

2. Use native ARM packages where available:
   ```bash
   conda config --add channels apple
   ```

#### Gatekeeper Blocks

**Symptom**: "Cannot be opened because the developer cannot be verified".

**Solution**:
```bash
xattr -d com.apple.quarantine /path/to/file
```

---

## Getting Help

If you're still experiencing issues:

1. **Run full diagnostics**:
   ```bash
   kintsugi check --verbose > diagnostics.txt
   python -c "import sys; print(sys.version)"
   conda list > packages.txt
   ```

2. **Check GitHub Issues**:
   [https://github.com/smith6jt-cop/KINTSUGI/issues](https://github.com/smith6jt-cop/KINTSUGI/issues)

3. **Create a new issue** with:
   - Operating system and version
   - Python version
   - Error message (full traceback)
   - Steps to reproduce
   - Output of `kintsugi check`

---

## Quick Reference

| Issue | Quick Fix |
|-------|-----------|
| libvips not found | Windows: Download from Zenodo; Linux: `apt install libvips-dev` |
| Java not found | `conda install openjdk=11` |
| GPU not detected | Check `nvidia-smi`; reinstall PyTorch with CUDA |
| Import errors | `conda activate KINTSUGI && pip install -e .` |
| Memory errors | Reduce `max_image_dim_px` in config |
| Environment conflicts | Use `env_streamlined.yml` |
