# KINTSUGI Migration Guide

## Transitioning from Notebooks 3 and 5 to the New Unified Workflow

This guide helps you transition from the deprecated Notebooks 3 (Signal Isolation) and 5 (DL Channel Refinement) to the new unified `3_Signal_Isolation_QC.ipynb` notebook.

## What Changed

The old notebooks have been replaced with a single, more powerful notebook:

- **Old:** `3_Signal_Isolation.ipynb` + `5_DL_Channel_Refinement.ipynb`
- **New:** `3_Signal_Isolation_QC.ipynb`

## Overview of Changes

| Old Approach | New Approach |
|--------------|--------------|
| Notebook 3: Manual parameter tuning | Unified notebook with 3 workflow options |
| Notebook 5: DL quality assessment | Integrated QC at each processing step |
| Kutils.py functions | Same functions + new kintsugi.denoise/qc modules |
| Manual review | Automated assessment with human-in-the-loop |
| No parameter memory | Parameter learning by tissue/marker |

## Quick Start

### Option 0: New Unified Notebook (Simplest)

Open `notebooks/3_Signal_Isolation_QC.ipynb` and follow the step-by-step workflow. The notebook supports three approaches:

1. **Claude-Guided**: AI recommends parameters, learns from success
2. **Interactive Tuners**: Widget-based parameter adjustment
3. **Python API**: Direct programmatic control

### Option 1: Claude Code Integration (Recommended for AI-assisted workflow)

1. **Install dependencies:**
   ```bash
   pip install kintsugi[claude]
   ```

2. **Configure Claude Code:**
   - If you used `kintsugi init`, configuration is already set up
   - For existing projects: `kintsugi mcp config /path/to/your/project`

3. **Start using with Claude Code:**
   - Claude can now load channels, suggest parameters, and process images
   - Parameters are learned and improve over time

### Option 2: Jupyter Interactive Tuners

```python
from kintsugi.claude import (
    blank_subtraction_tuner,
    denoise_tuner,
    clahe_tuner,
    clean_tuner,
)

# Interactive blank subtraction
result = blank_subtraction_tuner(
    signal_channel="cyc001_CD3",
    blank_channel="cyc001_blank",
)

# Get optimized parameters for batch processing
params = result.get_params()
```

### Option 3: Pure Python API

```python
from kintsugi.denoise import denoise_n2v, denoise_nlm, adaptive_denoise
from kintsugi.qc import ImageQC, CellQC, MarkerQC

# Advanced denoising
denoised = adaptive_denoise(image, strength="auto")

# Quality control
qc = ImageQC()
result = qc.assess(image)
if not result.passed:
    print(f"Issues: {result.issues}")
```

---

## Detailed Migration

### Notebook 3: Signal Isolation

#### Old: ini_params function
```python
# OLD (Kutils.py)
from Kutils import ini_params
result = ini_params(
    signal, blank,
    blank_clip_factor=1000,
    blank_scale_factor=0.8,
    smooth_low=True,
    erosion=1
)
```

#### New: MCP Tool or Parameter Tuner
```python
# NEW Option 1: Jupyter tuner
from kintsugi.claude import blank_subtraction_tuner
result = blank_subtraction_tuner(signal_channel, blank_channel)
params = result.get_params()

# NEW Option 2: Direct API
from kintsugi.mcp.tools.signal_isolation import subtract_blank
result = await subtract_blank(
    signal_channel="cyc001_CD3",
    blank_channel="cyc001_blank",
    blank_clip_factor=1000,
    blank_scale_factor=0.8,
    smooth_low=True,
    erosion=1
)
```

#### Old: denoise function
```python
# OLD
from Kutils import denoise
result = denoise(image, medn=True, med_size=3)
```

#### New: Advanced Denoising
```python
# NEW: Traditional filters
from kintsugi.denoise import denoise_median, denoise_nlm

result = denoise_median(image, size=3)
# or
result = denoise_nlm(image, patch_size=7)

# NEW: Self-supervised (N2V)
from kintsugi.denoise import denoise_n2v
result = denoise_n2v(image, n_epochs=50)

# NEW: Adaptive (auto-select best method)
from kintsugi.denoise import adaptive_denoise
result = adaptive_denoise(image, strength="auto")
```

#### Old: CLAHE function
```python
# OLD
from Kutils import CLAHE
result = CLAHE(image, clip_limit=0.01, tileGridSize=70)
```

#### New: MCP Tool
```python
# NEW: Via MCP
result = await apply_clahe(
    channel="cyc001_CD3",
    clip_limit=0.01,
    tile_grid_size=70
)

# Or via Jupyter tuner
from kintsugi.claude import clahe_tuner
result = clahe_tuner(channel_name)
```

---

### Notebook 5: DL Channel Refinement

#### Old: ChannelAssessor
```python
# OLD
from dl_refinement import ChannelAssessor
assessor = ChannelAssessor(model_path="model.pth")
result = assessor.process_image(image, "CD3")
```

#### New: QC Module
```python
# NEW: Image-level QC
from kintsugi.qc import ImageQC
qc = ImageQC()
result = qc.assess(image, marker="CD3", tissue="tonsil")

print(f"Passed: {result.passed}")
print(f"Quality Score: {result.quality_score}")
print(f"Issues: {result.issues}")
print(f"Recommendations: {result.recommendations}")
```

#### Old: BatchChannelProcessor
```python
# OLD
from dl_refinement import BatchChannelProcessor
processor = BatchChannelProcessor(assessor)
results = processor.process_directory(input_dir, output_dir)
```

#### New: Batch QC
```python
# NEW: Batch effects detection
from kintsugi.qc import BatchQC, detect_batch_effects

qc = BatchQC()
result = qc.assess(
    data=cell_data,
    batch_column="batch_id",
    marker_columns=["CD3", "CD20", "DAPI"]
)

if result.batch_effects_detected:
    print(f"Affected markers: {result.affected_markers}")
```

---

## Parameter Learning

The new system learns from successful parameter choices:

```python
# Parameters are automatically recorded when you approve results
from kintsugi.claude import ParameterLearningEngine

engine = ParameterLearningEngine("/path/to/project")

# Get recommendations based on tissue/marker
recommendations = engine.recommend_parameters(
    tissue_type="tonsil",
    marker_name="CD3",
    operation="blank_subtraction"
)

print(f"Recommended params: {recommendations['recommended_parameters']}")
print(f"Confidence: {recommendations['confidence']}")
```

When using Claude Code, this happens automatically through the MCP tools.

---

## Available MCP Tools

| Tool | Description |
|------|-------------|
| `load_channel` | Load channel image from project |
| `subtract_blank` | Autofluorescence subtraction |
| `denoise` | Basic denoising filters |
| `denoise_advanced` | N2V, NLM, BM3D, etc. |
| `apply_clahe` | Contrast enhancement |
| `clean_background` | Background removal |
| `gaussian_subtract` | Structured background removal |
| `assess_quality` | Channel quality assessment |
| `compute_snr` | Signal-to-noise ratio |
| `suggest_parameters` | AI-powered suggestions |
| `suggest_with_learning` | Suggestions using learned history |
| `approve_and_learn` | Record successful params |

---

## New Modules

### kintsugi.denoise
- `denoise_median`, `denoise_gaussian`, `denoise_bilateral`
- `denoise_nlm` (Non-Local Means)
- `denoise_n2v` (Noise2Void)
- `denoise_bm3d_lite` (BM3D-inspired)
- `adaptive_denoise` (auto-select)

### kintsugi.qc
- `ImageQC` - Image-level quality control
- `CellQC` - Cell-level filtering
- `MarkerQC` - Marker validation
- `BatchQC` - Batch effects detection

### kintsugi.segment
- `SAMSegmenter` - SAM wrapper for microscopy
- `segment_nuclei_watershed` - Classical segmentation
- `refine_masks` - Post-processing utilities

---

## FAQ

**Q: What if I don't use Claude Code?**
A: You can use the Jupyter interactive tuners (Kview2 + Kutils) or direct Python API. The MCP server is optional.

**Q: Will my existing parameters work?**
A: Yes, the parameters are similar. See the mapping tables above.

**Q: How does parameter learning work?**
A: Successful parameters are stored in a SQLite database, indexed by tissue type and marker name. Future recommendations are weighted by past success.

**Q: Can I still use Kutils.py functions?**
A: Yes, `Kutils.py` is still available in the notebooks folder and is used by the interactive tuners in `3_Signal_Isolation_QC.ipynb`.

---

## Getting Help

- Documentation: See `docs/` directory
- Issues: https://github.com/smith6jt-cop/KINTSUGI/issues
- Claude Code: Use `/help` for available tools
