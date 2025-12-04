# Processing Workflows

This guide describes the typical image processing workflows available in KINTSUGI.

## Overview

KINTSUGI provides a complete pipeline for multi-cycle immunofluorescence image analysis:

1. **Illumination Correction** - Fix uneven lighting across images
2. **Stitching** - Combine tiled images into single large images
3. **Deconvolution** - Remove optical blur and improve resolution
4. **Extended Depth of Focus (EDoF)** - Combine Z-stacks into 2D
5. **Registration** - Align images across multiple cycles
6. **Autofluorescence Removal** - Subtract background signal
7. **Segmentation** - Identify cells and structures
8. **Spatial Analysis** - Analyze spatial relationships

## Workflow 1: Single Channel Evaluation

Use this workflow to test and tune parameters on a single channel before batch processing.

**Notebook**: `notebooks/1_Single_Channel_Eval.ipynb`

### Steps:
1. Load a single channel image
2. Test illumination correction parameters
3. Test stitching parameters
4. Test deconvolution settings
5. Evaluate results visually

## Workflow 2: Cycle Processing

Batch process all channels across multiple cycles.

**Notebook**: `notebooks/2_Cycle_Processing.ipynb`

### Steps:
1. Configure input/output directories
2. Set processing parameters
3. Run illumination correction
4. Run stitching
5. Run deconvolution (optional)
6. Run EDoF (optional)
7. Run registration

## Workflow 3: Signal Isolation

Remove autofluorescence and isolate true signal.

**Notebook**: `notebooks/3_Signal_Isolation.ipynb`

### Steps:
1. Load registered images
2. Identify autofluorescence channels
3. Subtract autofluorescence
4. Apply filtering
5. Generate final signal images

## Workflow 4: Segmentation Analysis

Perform cell segmentation and spatial analysis.

**Notebook**: `notebooks/4_Segmentation_Analysis.ipynb`

### Steps:
1. Load processed images
2. Run InstanSeg segmentation
3. Extract cell features
4. Perform clustering
5. Analyze spatial relationships

## Registration Module (Kreg)

The Kreg module provides image registration functionality based on VALIS.

### Basic Usage

```python
from kintsugi.kreg import Valis

# Initialize registrar
registrar = Valis(
    src_dir="/path/to/images",
    dst_dir="/path/to/output",
    reference_img_f="cycle1.tif",
    align_to_reference=True,
)

# Run registration
registrar.register()
registrar.register_micro()  # Fine registration
registrar.warp_and_merge_slides()  # Output registered images
```

### Configuration Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `src_dir` | Source image directory | Required |
| `dst_dir` | Output directory | Required |
| `reference_img_f` | Reference image filename | Required |
| `max_image_dim_px` | Max image dimension for processing | 2048 |
| `compose_non_rigid` | Enable non-rigid registration | True |
| `crop_to_overlap` | Crop to overlapping region | True |

## Visualization Module (Kview2)

The Kview2 module provides interactive visualization tools.

### Basic Usage

```python
from kintsugi.kview2 import imshow, curtain, crop

# Display image
imshow(image)

# Compare two images with curtain view
curtain(image1, image2)

# Interactive crop
cropped = crop(image)
```

## Stitching Module (Kstitch)

The Kstitch module provides tile stitching functionality.

### Basic Usage

```python
from kintsugi.kstitch import stitch_tiles

# Stitch tiled images
stitched = stitch_tiles(
    tile_dir="/path/to/tiles",
    output_path="/path/to/output.tif",
    overlap=0.1,  # 10% overlap
)
```

## Tips for Large Datasets

1. **Start small**: Test parameters on a subset before full batch processing
2. **Monitor memory**: Use `kintsugi info` to check available resources
3. **Use tiled processing**: Enable tiled output for very large images
4. **GPU acceleration**: Enable CuPy for faster deconvolution
