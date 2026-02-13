# API Reference

## Main Package

```python
import kintsugi
```

### Functions

#### kintsugi.check_dependencies

```python
kintsugi.check_dependencies(verbose: bool = True) -> dict
```

Check all KINTSUGI dependencies and return status.

**Parameters:**
- `verbose` (bool): If True, print status messages. Default: True

**Returns:**
- dict: Dictionary with dependency status information.

#### kintsugi.get_config_template

```python
kintsugi.get_config_template() -> dict
```

Get a template configuration dictionary for registration workflow.

**Returns:**
- dict: Template configuration with all available options.

## Kreg Module (Registration)

```python
from kintsugi import Kreg
# or
from kintsugi.kreg import Valis
```

### Valis Class

The main class for image registration.

```python
from kintsugi.kreg import Valis

registrar = Valis(
    src_dir="/path/to/images",
    dst_dir="/path/to/output",
    reference_img_f="reference.tif",
    **kwargs
)
```

**Parameters:**
- `src_dir` (str): Path to source image directory
- `dst_dir` (str): Path to output directory
- `reference_img_f` (str): Filename of reference image
- `align_to_reference` (bool): Align all images to reference. Default: True
- `max_image_dim_px` (int): Maximum image dimension. Default: 2048
- `compose_non_rigid` (bool): Enable non-rigid registration. Default: True
- `crop_to_overlap` (bool): Crop output to overlapping region. Default: True

**Methods:**
- `register()`: Perform initial registration
- `register_micro()`: Perform micro (fine) registration
- `warp_and_merge_slides()`: Warp images and save results

## Kview2 Module (Visualization)

```python
from kintsugi import Kview
# or
from kintsugi.kview2 import imshow, curtain, crop
```

### imshow

Display an image interactively.

```python
from kintsugi.kview2 import imshow

imshow(image, **kwargs)
```

### curtain

Compare two images with an interactive curtain view.

```python
from kintsugi.kview2 import curtain

curtain(image1, image2, **kwargs)
```

### crop

Interactively crop an image.

```python
from kintsugi.kview2 import crop

cropped = crop(image)
```

## Kstitch Module (Stitching)

```python
from kintsugi import Kstitch
```

Provides tile stitching functionality for large tiled images.

## signal Module

```python
from kintsugi.signal import (
    subtract_autofluorescence,
    analyze_for_subtraction,
    AutofluorescenceSubtractor,
)
```

Provides autofluorescence subtraction with intelligent parameter suggestion and learning capabilities.

## qc Module

```python
from kintsugi.qc import ImageQC, CellQC, MarkerQC, BatchQC
```

Quality control and assessment tools for image processing pipelines.

## denoise Module

```python
from kintsugi.denoise import (
    adaptive_denoise,
    denoise_median,
    denoise_nlm,
    denoise_bilateral,
)
```

Advanced denoising algorithms including adaptive, NLM, bilateral, and N2V methods.

## edf Module (Extended Depth of Focus)

```python
from kintsugi.edf import extended_depth_of_focus_variance
```

### extended_depth_of_focus_variance

Compute extended depth of focus from a z-stack using variance-based focus selection.

```python
from kintsugi.edf import extended_depth_of_focus_variance

edf = extended_depth_of_focus_variance(
    stack,               # (Z, Y, X) array
    radius_x=2,          # Variance window radius X
    radius_y=2,          # Variance window radius Y
    sigma=10.0,          # Smoothing sigma for variance map
    blend_depth=2,       # Adjacent z-slices to blend
    z_smooth_sigma=1.0,  # Gaussian smoothing for z-index map
)
```

Uses CuPy (GPU) when available, falls back to NumPy/SciPy (CPU). Device mode is controlled by `KINTSUGI_DEVICE_MODE` environment variable in SLURM jobs.

## kcorrect_gpu Module (Illumination Correction)

```python
from kintsugi.kcorrect_gpu import KCorrectGPU
```

GPU-accelerated BaSiC illumination correction using CuPy FFT.

```python
corrector = KCorrectGPU(device_id=0)
flatfield, darkfield = corrector.estimate(images)
corrected = corrector.transform(images, flatfield, darkfield)
```

## project Module

```python
from kintsugi.project import KintsugiProject, ExperimentConfig
```

### KintsugiProject Class

```python
project = KintsugiProject.load("/path/to/project")
print(project.name)
print(project.paths.raw)       # data/raw/
print(project.paths.processed) # data/processed/
print(project.config)          # ExperimentConfig
```

### ExperimentConfig Class

```python
config = ExperimentConfig.from_dict(json_data)
# Automatically translates CODEX field names to KINTSUGI equivalents
```

## gpu Module

```python
from kintsugi.gpu import get_gpu_manager
```

Multi-GPU detection and management.

```python
gpu = get_gpu_manager()
print(gpu.summary())         # List available GPUs
print(gpu.cupy_installed)    # True if CuPy package is importable
print(gpu.cupy_available)    # True if GPU hardware is accessible
```

## hpc Module

```python
from kintsugi.hpc import detect_multi_account_resources, detect_live_multi_account
```

HPC auto-detection for multi-account SLURM environments.

```python
# Static detection (from sacctmgr associations)
resources = detect_multi_account_resources()

# Live detection (accounts for running jobs)
live = detect_live_multi_account()
```

## parallel_io Module

```python
from kintsugi.parallel_io import load_tiles_parallel, save_stack_parallel
```

Parallel image loading and saving utilities for multi-threaded I/O.

## zarr_io Module

```python
from kintsugi.zarr_io import write_ome_zarr, read_ome_zarr
```

OME-Zarr format I/O with Dask lazy loading support.

## deps Module

```python
from kintsugi import deps
```

### DependencyChecker Class

```python
from kintsugi.deps import DependencyChecker

checker = DependencyChecker()
results = checker.check_all(verbose=True)
```

Check and report on all external dependencies.

## cli Module

```python
from kintsugi.cli import main
```

Click-based CLI with command groups: `main`, `mcp`, `slurm`, `workflow`. See [CLI Reference](cli.md) for full command documentation.
