# Evaluation: PyImageJ/CLIJ2 vs Pure Python for Extended Depth of Field (EDF)

> **⚠️ DEPRECATION NOTICE (2025):** This evaluation document is historical. Based on this analysis,
> KINTSUGI has **adopted the pure Python approach** (CuPy/NumPy) as the default implementation.
> PyImageJ/CLIJ2 has been deprecated and removed as a required dependency. See `src/kintsugi/edf.py`
> for the current implementation.

## Executive Summary

This document evaluates the current PyImageJ/CLIJ2 implementation for 3D to 2D Extended Depth of Field (EDF) conversion and assesses pure Python alternatives for both notebook testing and batch processing workflows.

**UPDATE (2025)**: The pure Python implementation using CuPy (GPU) or NumPy (CPU) has been adopted
as the default. PyImageJ/CLIJ2 is now deprecated and no longer required.

~~**Recommendation**: A **hybrid approach** is recommended:~~
~~1. **Keep PyImageJ/CLIJ2** for production batch processing (GPU-accelerated, proven performance)~~
~~2. **Implement a pure Python fallback** for environments without Java/GPU or for testing~~

---

## 1. Current Implementation Analysis

### 1.1 Architecture Overview

The current EDF implementation uses:
- **PyImageJ** (`imagej >= 1.4.0`) as the Python-Java bridge
- **CLIJ2** plugin for GPU-accelerated processing
- **`extendedDepthOfFocusVarianceProjection`** algorithm

**Key Files:**
- `fiji_macros/EDF_Macro.ijm` - Standalone macro (72 lines)
- `notebooks/1_Single_Channel_Eval.ipynb` - Interactive testing (Section 6)
- `notebooks/2_Cycle_Processing.ipynb` - Batch processing (Section 5)

### 1.2 Algorithm Parameters

From `EDF_Macro.ijm`:
```
radius_x = 5.0      # Gaussian blur radius X
radius_y = 5.0      # Gaussian blur radius Y
sigma = 20.0        # Variance projection sigma
xTiles = 3          # X-dimension tile splits
yTiles = 3          # Y-dimension tile splits
z_start = 3         # Starting Z-slice
z_end = 15          # Ending Z-slice
```

### 1.3 Performance Benchmarks

From macro comments:
| Device | Time (2x2 tiling) |
|--------|-------------------|
| NVIDIA RTX A4500 | 42.1 seconds |
| Intel UHD Graphics 770 | 1m 59.4s |

Image size: (17, 9484, 9962) - approximately 17 Z-slices at ~9500x10000 pixels

### 1.4 Current Dependencies

```yaml
# From env-*.yml and pyproject.toml
pyimagej >= 1.4.0
jpype1 >= 1.5.0
scyjava >= 1.0.0
openjdk = 11 (or 21 on Windows)
maven = 3.9.9
```

Plus local Fiji installation with CLIJ2 plugin.

---

## 2. CLIJ2 Algorithm Details

### 2.1 What `extendedDepthOfFocusVarianceProjection` Does

The CLIJ2 variance projection algorithm:

1. **For each (x,y) position** in the image:
   - Calculates local variance in a window of size `(2*radius_x+1, 2*radius_y+1)` for each Z-slice
   - The variance is computed after Gaussian smoothing with `sigma`
   - Selects the pixel value from the Z-slice with **maximum local variance**

2. **Rationale**: High local variance indicates sharp edges/features (in-focus regions). Low variance indicates blur (out-of-focus regions).

3. **Mathematical formulation**:
   ```
   For each (x,y):
     variance[z] = Var(GaussianBlur(I[z], sigma)[x-rx:x+rx, y-ry:y+ry])
     z_best = argmax(variance)
     output[x,y] = I[z_best, x, y]
   ```

---

## 3. Pure Python Alternatives

### 3.1 Available Libraries

| Library | Method | GPU Support | Notes |
|---------|--------|-------------|-------|
| **NumPy/SciPy** | Variance-based | No (but vectorized) | Most flexible, requires custom implementation |
| **focus-stack** | Laplacian gradient | No | PyPI package, designed for microscopy |
| **OpenCV** | Laplacian pyramid | Yes (CUDA) | Well-tested, multi-platform |
| **Mahotas** | Sobel-based focus measure | No | Microscopy-focused, good documentation |
| **CuPy** | Custom variance | Yes (CUDA) | NumPy API with GPU acceleration |

### 3.2 Recommended Pure Python Implementation

A variance-based implementation matching CLIJ2's approach:

```python
import numpy as np
from scipy.ndimage import uniform_filter, gaussian_filter

def extended_depth_of_focus_variance(
    stack: np.ndarray,
    radius_x: int = 5,
    radius_y: int = 5,
    sigma: float = 20.0,
    use_gpu: bool = False
) -> np.ndarray:
    """
    Extended depth of focus using local variance projection.

    Matches CLIJ2's extendedDepthOfFocusVarianceProjection algorithm.

    Parameters
    ----------
    stack : np.ndarray
        3D image stack with shape (z, y, x)
    radius_x : int
        Half-width of variance calculation window in X
    radius_y : int
        Half-height of variance calculation window in Y
    sigma : float
        Gaussian smoothing sigma before variance calculation
    use_gpu : bool
        If True and CuPy available, use GPU acceleration

    Returns
    -------
    np.ndarray
        2D projected image with shape (y, x)
    """
    if use_gpu:
        try:
            import cupy as cp
            from cupyx.scipy.ndimage import uniform_filter as gpu_uniform_filter
            from cupyx.scipy.ndimage import gaussian_filter as gpu_gaussian_filter
            xp = cp
            stack = cp.asarray(stack)
            uf = gpu_uniform_filter
            gf = gpu_gaussian_filter
        except ImportError:
            xp = np
            uf = uniform_filter
            gf = gaussian_filter
    else:
        xp = np
        uf = uniform_filter
        gf = gaussian_filter

    n_slices, height, width = stack.shape
    variances = xp.zeros_like(stack, dtype=xp.float32)

    # Window size for variance calculation
    window_size = (2 * radius_y + 1, 2 * radius_x + 1)

    for z in range(n_slices):
        # Apply Gaussian smoothing
        smoothed = gf(stack[z].astype(xp.float32), sigma=sigma)

        # Calculate local variance: Var(X) = E[X^2] - E[X]^2
        local_mean = uf(smoothed, size=window_size)
        local_mean_sq = uf(smoothed ** 2, size=window_size)
        variances[z] = local_mean_sq - local_mean ** 2

    # Find Z-index with maximum variance at each (x,y)
    max_var_idx = xp.argmax(variances, axis=0)

    # Select pixels from best-focus slices
    y_idx, x_idx = xp.meshgrid(
        xp.arange(height),
        xp.arange(width),
        indexing='ij'
    )
    result = stack[max_var_idx, y_idx, x_idx]

    if use_gpu:
        result = cp.asnumpy(result)

    return result
```

### 3.3 Tiled Processing for Large Images

```python
def edf_tiled(
    stack: np.ndarray,
    tile_size: tuple = (3000, 3000),
    overlap: int = 50,
    **edf_kwargs
) -> np.ndarray:
    """
    Process large stacks in tiles to manage memory.

    Parameters
    ----------
    stack : np.ndarray
        3D image stack (z, y, x)
    tile_size : tuple
        (height, width) of each tile
    overlap : int
        Pixel overlap between tiles for seamless blending
    **edf_kwargs
        Arguments passed to extended_depth_of_focus_variance

    Returns
    -------
    np.ndarray
        2D EDF result
    """
    n_slices, height, width = stack.shape
    tile_h, tile_w = tile_size

    result = np.zeros((height, width), dtype=stack.dtype)
    weight = np.zeros((height, width), dtype=np.float32)

    # Create blending weights (raised cosine)
    blend_y = np.ones(tile_h)
    blend_x = np.ones(tile_w)
    if overlap > 0:
        ramp = np.linspace(0, 1, overlap)
        blend_y[:overlap] = ramp
        blend_y[-overlap:] = ramp[::-1]
        blend_x[:overlap] = ramp
        blend_x[-overlap:] = ramp[::-1]
    blend_weights = np.outer(blend_y, blend_x)

    # Process tiles
    for y_start in range(0, height, tile_h - overlap):
        for x_start in range(0, width, tile_w - overlap):
            y_end = min(y_start + tile_h, height)
            x_end = min(x_start + tile_w, width)

            tile = stack[:, y_start:y_end, x_start:x_end]
            tile_result = extended_depth_of_focus_variance(tile, **edf_kwargs)

            # Apply blending weights
            tile_blend = blend_weights[:y_end-y_start, :x_end-x_start]
            result[y_start:y_end, x_start:x_end] += tile_result * tile_blend
            weight[y_start:y_end, x_start:x_end] += tile_blend

    # Normalize by weights
    result = (result / np.maximum(weight, 1e-10)).astype(stack.dtype)

    return result
```

---

## 4. Comparison: PyImageJ/CLIJ2 vs Pure Python

### 4.1 Performance Comparison

| Aspect | PyImageJ/CLIJ2 | Pure Python (NumPy) | Pure Python (CuPy) |
|--------|----------------|---------------------|-------------------|
| **Speed (GPU)** | ~42s | N/A | ~45-60s estimated |
| **Speed (CPU)** | ~2 min | ~3-5 min | N/A |
| **Memory** | High (JVM + GPU) | Moderate | Moderate (GPU) |
| **Setup complexity** | High | Low | Medium |
| **Dependencies** | Java, Fiji, Maven | NumPy, SciPy | NumPy, CuPy, CUDA |

### 4.2 Pros and Cons

#### PyImageJ/CLIJ2
**Pros:**
- Proven, battle-tested CLIJ2 implementation
- Excellent GPU utilization via OpenCL
- Access to ImageJ's extensive plugin ecosystem
- Already integrated and working

**Cons:**
- Complex dependency chain (Java, Fiji, Maven, OpenCL)
- Long initialization time (30s - 3min on first run)
- JVM memory overhead
- Harder to debug/modify
- Platform-specific setup issues

#### Pure Python
**Pros:**
- Simple installation (`pip install numpy scipy`)
- Easy to understand, modify, and debug
- No JVM overhead
- Better integration with Python ML/analysis pipelines
- Portable across all platforms

**Cons:**
- CPU implementation ~2-3x slower
- Need to implement and validate algorithm
- No access to ImageJ plugins

---

## 5. Recommendations

### 5.1 Hybrid Architecture (Recommended)

```
┌─────────────────────────────────────────────────────────┐
│                    EDF Processing                        │
├─────────────────────────────────────────────────────────┤
│                                                          │
│   ┌─────────────┐         ┌─────────────────────────┐   │
│   │  Input:     │         │   Backend Selection:    │   │
│   │  3D Stack   │────────▶│   - GPU + CLIJ2 avail?  │   │
│   │  (Z,Y,X)    │         │   - CuPy available?     │   │
│   └─────────────┘         │   - CPU fallback        │   │
│                           └───────────┬─────────────┘   │
│                                       │                  │
│              ┌────────────────────────┼────────────┐    │
│              ▼                        ▼            ▼    │
│   ┌──────────────────┐  ┌─────────────────┐ ┌────────┐ │
│   │ PyImageJ/CLIJ2   │  │ CuPy (GPU)      │ │ NumPy  │ │
│   │ (Production)     │  │ (Alternative)   │ │ (CPU)  │ │
│   └────────┬─────────┘  └───────┬─────────┘ └───┬────┘ │
│            │                    │               │       │
│            └────────────────────┼───────────────┘       │
│                                 ▼                        │
│                        ┌─────────────┐                  │
│                        │  Output:    │                  │
│                        │  2D Image   │                  │
│                        └─────────────┘                  │
└─────────────────────────────────────────────────────────┘
```

### 5.2 Implementation Plan

#### Phase 1: Create Pure Python EDF Module
1. Implement `extended_depth_of_focus_variance()` function
2. Add tiled processing support
3. Add optional CuPy GPU acceleration
4. Validate output against CLIJ2 results

#### Phase 2: Create Unified Interface
```python
# src/kintsugi/edf.py
class EDFProcessor:
    def __init__(self, backend='auto'):
        """
        backend: 'clij2', 'cupy', 'numpy', or 'auto'
        """
        self.backend = self._detect_backend(backend)

    def process(self, stack, radius_x=5, radius_y=5, sigma=20.0,
                z_start=None, z_end=None, tiles=(3,3)):
        """Unified EDF processing interface."""
        pass
```

#### Phase 3: Integration
1. Update notebooks to use unified interface
2. Add backend selection in batch processing
3. Update documentation

### 5.3 Use Case Recommendations

| Use Case | Recommended Backend | Rationale |
|----------|---------------------|-----------|
| **Notebook testing** | Pure Python (NumPy) | Fast iteration, no setup needed |
| **Small batches (<10 images)** | Pure Python (CuPy if available) | Avoid JVM startup overhead |
| **Large batches (>10 images)** | PyImageJ/CLIJ2 | Best sustained throughput |
| **CI/CD testing** | Pure Python (NumPy) | No Java dependencies |
| **Cloud/HPC** | Pure Python (CuPy) | Better CUDA integration |
| **Local workstation** | PyImageJ/CLIJ2 | Leverage existing setup |

---

## 6. Validation Strategy

To ensure the pure Python implementation matches CLIJ2:

```python
def validate_edf_implementations():
    """Compare Python EDF output to CLIJ2 reference."""
    # Load test stack
    test_stack = load_test_data()

    # Run CLIJ2 (reference)
    clij2_result = run_clij2_edf(test_stack)

    # Run Python implementation
    python_result = extended_depth_of_focus_variance(test_stack)

    # Compare
    mse = np.mean((clij2_result - python_result) ** 2)
    ssim = structural_similarity(clij2_result, python_result)

    assert mse < 1e-3, f"MSE too high: {mse}"
    assert ssim > 0.99, f"SSIM too low: {ssim}"
```

---

## 7. Conclusion

A pure Python implementation of the EDF variance projection algorithm is **feasible and recommended** as a complement to the existing PyImageJ/CLIJ2 workflow. The implementation should:

1. **Match CLIJ2's algorithm** exactly for consistency
2. **Support GPU acceleration** via CuPy for competitive performance
3. **Provide a unified interface** that auto-selects the best backend
4. **Be validated** against CLIJ2 output to ensure correctness

This approach provides flexibility for different deployment scenarios while maintaining the proven performance of CLIJ2 for production batch processing.

---

## Appendix: Quick Reference

### Current EDF Call in Notebooks

```python
# From 1_Single_Channel_Eval.ipynb, cell 83
macro = """
run("CLIJ2 Macro Extensions", "cl_device=[" + device + "]");
Ext.CLIJ2_extendedDepthOfFocusVarianceProjection(input, output, radius_x, radius_y, sigma);
"""
ij.py.run_macro(macro, args)
```

### Proposed Python Equivalent

```python
# Proposed usage
from kintsugi.edf import EDFProcessor

processor = EDFProcessor(backend='auto')  # or 'clij2', 'cupy', 'numpy'
result = processor.process(
    stack,
    radius_x=5,
    radius_y=5,
    sigma=20.0,
    z_start=3,
    z_end=15,
    tiles=(3, 3)
)
```
