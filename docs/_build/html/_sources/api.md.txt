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
from kintsugi import Kview2
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

## dl_refinement Module

```python
from kintsugi import dl_refinement
```

Provides deep learning-based refinement tools for image processing.

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
