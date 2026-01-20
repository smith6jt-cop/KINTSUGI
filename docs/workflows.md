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
5. **Quantitative QC analysis (with PDF export)**
6. Run deconvolution (optional)
7. Run EDoF (optional)
8. Run registration

### QC Plot Output

The notebook automatically generates QC plots and saves them as PDFs to `PROJECT_DIR/qc_plots/`:

| File | Description |
|------|-------------|
| `raw_summary_heatmaps.pdf` | SNR, CV, and intensity heatmaps for raw data |
| `stitched_summary_heatmaps.pdf` | Same metrics after stitching |
| `deconvolved_summary_heatmaps.pdf` | Same metrics after deconvolution |
| `edf_summary_heatmaps.pdf` | EDF projection statistics |
| `{stage}_zprofile_cyc{NN}_ch{N}.pdf` | Z-plane profiles for first cycle |

The plotting functions accept optional parameters:
- `stage_name`: Label for the processing stage (used in PDF filename)
- `save_pdf`: Whether to save PDF (default: True)

### QC Module

The QC statistics and plotting functions are available in `notebooks/Kprocess.py` for programmatic use:

```python
from Kprocess import (
    run_raw_qc,        # Complete raw data QC workflow
    run_stitched_qc,   # Complete stitched QC with raw comparison
    run_decon_qc,      # Complete decon QC with stitched comparison
    run_edf_qc,        # Complete EDF QC with decon comparison
    plot_summary_heatmaps,  # Individual heatmap plotting
    plot_zplane_profiles,   # Z-profile plotting
)

# Example: Run stitched QC programmatically
stitched_df = run_stitched_qc(
    stitch_dir=stitch_dir,
    cache_file=PROJECT_DIR / 'cache' / 'stitched_stats.pkl',
    start_cycle=1, end_cycle=13,
    start_channel=1, end_channel=4,
    n_zplanes=15,
    qc_output_dir=PROJECT_DIR / 'qc_plots',
    raw_stats_df=raw_df,  # Optional: for comparison
)
```

## Workflow 3: Signal Isolation

Remove autofluorescence and isolate true signal.

> **Note:** The notebook-based approach (`3_Signal_Isolation.ipynb`) is deprecated.
> Use the Claude-guided workflow or Python API instead.

### Option A: Claude-Guided Workflow (Recommended)

Use Claude Code with the KINTSUGI MCP server for interactive, AI-assisted signal isolation.

**Setup:**
```bash
pip install kintsugi[claude]
# If using kintsugi init, Claude config is automatic
# For existing projects: kintsugi mcp config /path/to/project
```

**Usage:**
Claude Code can:
- Load and analyze channels
- Suggest optimal parameters based on image characteristics
- Apply processing with real-time feedback
- Learn from successful parameters for future recommendations

### Option B: Python API

```python
from kintsugi.denoise import adaptive_denoise, denoise_nlm
from kintsugi.qc import ImageQC

# Load and assess image quality
qc = ImageQC()
result = qc.assess(image)
print(f"Quality: {result.quality_score}, Issues: {result.issues}")

# Adaptive denoising (auto-selects best method)
denoised = adaptive_denoise(image, strength="auto")

# Or use specific methods
denoised = denoise_nlm(image, patch_size=7, patch_distance=11)
```

### Option C: Legacy Notebook

**Notebook**: `notebooks/3_Signal_Isolation.ipynb` (DEPRECATED)

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

## Workflow 5: Quality Control

Assess and validate image quality at multiple levels.

### Image-Level QC

```python
from kintsugi.qc import ImageQC

qc = ImageQC()
result = qc.assess(image, marker="CD3", tissue="tonsil")

print(f"Passed: {result.passed}")
print(f"Quality Score: {result.quality_score}")
print(f"Issues: {result.issues}")
print(f"Recommendations: {result.recommendations}")
```

### Cell-Level QC

```python
from kintsugi.qc import CellQC

qc = CellQC()
result = qc.assess(
    cell_data,
    marker_columns=["CD3", "CD20", "DAPI"],
    morphology_columns=["area", "eccentricity"]
)

# Filter problematic cells
filtered_data = result.filtered_data
print(f"Removed {len(result.outliers)} outliers")
```

### Marker Validation

```python
from kintsugi.qc import MarkerQC

qc = MarkerQC()
result = qc.assess(intensities, marker_name="CD3", cell_types=cell_types)

if result.crosstalk_detected:
    print(f"Crosstalk with: {result.crosstalk_markers}")
```

### Batch Effects

```python
from kintsugi.qc import BatchQC

qc = BatchQC()
result = qc.assess(
    data=combined_data,
    batch_column="batch_id",
    marker_columns=["CD3", "CD20", "DAPI"]
)

if result.batch_effects_detected:
    normalized = qc.normalize_batches(data, method="quantile")
```

## Denoising Module

KINTSUGI provides multiple denoising algorithms:

### Traditional Filters

```python
from kintsugi.denoise import (
    denoise_median,
    denoise_gaussian,
    denoise_bilateral,
    denoise_nlm,
)

# Median filter
result = denoise_median(image, size=3)

# Non-local means
result = denoise_nlm(image, patch_size=7, patch_distance=11)
```

### Deep Learning Denoising

```python
from kintsugi.denoise import denoise_n2v, denoise_care

# Noise2Void (self-supervised, no clean targets needed)
denoiser = N2VDenoiser()
denoiser.train(noisy_images, n_epochs=50)
result = denoiser.predict(image)

# CARE (supervised, requires paired data)
denoiser = CAREDenoiser()
denoiser.train(noisy_images, clean_images)
result = denoiser.predict(image)
```

### Adaptive Denoising

```python
from kintsugi.denoise import adaptive_denoise

# Automatically selects best method and parameters
result = adaptive_denoise(image, strength="auto")
```

## Segmentation Module

### Classical Segmentation

```python
from kintsugi.segment import segment_nuclei_watershed, segment_cells_watershed

# Nuclei segmentation
nuclei = segment_nuclei_watershed(
    dapi_image,
    min_distance=10,
    threshold_method="otsu"
)

# Cell segmentation with membrane expansion
cells = segment_cells_watershed(
    membrane_image,
    nuclei_labels=nuclei,
    expansion_distance=5
)
```

### SAM Segmentation

```python
from kintsugi.segment import SAMSegmenter

segmenter = SAMSegmenter(model_type="vit_b")
masks = segmenter.segment(image)

# With box prompts
masks = segmenter.segment_boxes(image, boxes=[[x1, y1, x2, y2]])
```

### Post-Processing

```python
from kintsugi.segment import refine_masks, filter_masks_by_size

# Filter by size
filtered = filter_masks_by_size(labels, min_size=50, max_size=5000)

# Comprehensive refinement
refined = refine_masks(
    labels,
    min_size=50,
    fill_holes=True,
    smooth=True,
    split_touching=True
)
```

## Tips for Large Datasets

1. **Start small**: Test parameters on a subset before full batch processing
2. **Monitor memory**: Use `kintsugi info` to check available resources
3. **Use tiled processing**: Enable tiled output for very large images
4. **GPU acceleration**: Enable CuPy for faster deconvolution
5. **Parameter learning**: Use Claude Code to build a database of successful parameters
6. **Batch QC**: Run BatchQC to detect and correct batch effects before analysis
