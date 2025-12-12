# Algorithm Optimization Plan: BM3D-lite & NLM Denoising

## Overview

This document outlines the plan for optimizing two critical image denoising algorithms:

1. **BM3D-lite**: O(n²) patch matching → O(n log n) with KD-tree spatial indexing
2. **NLM (Non-Local Means)**: O(n⁴) pixel-level → O(n²) with OpenCV's optimized implementation

## Current Implementation Analysis

### BM3D-lite (`denoise_bm3d_lite`)
- **Location**: `src/kintsugi/denoise/patch_based.py:97-256`
- **Current Complexity**: O(n²) where n = number of patches
- **Bottleneck**: Lines 167-187 - nested loop comparing all patches

```python
for i, (ref_patch, (ref_y, ref_x)) in enumerate(zip(patches, positions)):
    for j, (_other_patch, (other_y, other_x)) in enumerate(zip(patches, positions)):
        if abs(other_y - ref_y) <= half_window and abs(other_x - ref_x) <= half_window:
            diff = np.sum((patch_dcts[j] - ref_dct) ** 2)
```

### NLM (`denoise_patch_similarity`)
- **Location**: `src/kintsugi/denoise/patch_based.py:259-395`
- **Current Complexity**: O(n⁴) for non-fast mode, O(n²) for fast mode
- **Bottleneck**: Lines 369-388 - four nested loops for pixel-by-pixel processing

```python
for y in range(h):
    for x in range(w):
        for dy in range(-search_radius, search_radius + 1):
            for dx in range(-search_radius, search_radius + 1):
```

---

## Optimization Strategy

### 1. BM3D-lite: KD-Tree Spatial Indexing

**Approach**: Use `scipy.spatial.cKDTree` to find patches within the search window in O(log n) instead of O(n).

**Implementation**:
```python
from scipy.spatial import cKDTree

# Build spatial index on patch positions (one-time cost: O(n log n))
positions_array = np.array(positions)
tree = cKDTree(positions_array)

for i, (ref_patch, (ref_y, ref_x)) in enumerate(zip(patches, positions)):
    # Query patches within search window: O(log n + k)
    candidate_indices = tree.query_ball_point([ref_y, ref_x], r=half_window)

    # Filter by DCT similarity (only on candidates)
    for j in candidate_indices:
        if i == j:
            continue
        diff = np.sum((patch_dcts[j] - ref_dct) ** 2)
        ...
```

**Expected Speedup**: 10-50x for typical image sizes (1000+ patches)

**Risks**:
- Numerical differences from different patch ordering
- Edge cases with very few or very many patches in window

### 2. NLM: OpenCV Replacement

**Approach**: Replace custom implementation with `cv2.fastNlMeansDenoising()` which uses optimized C++ with integral images.

**Implementation**:
```python
import cv2

def denoise_nlm_optimized(image, h_param=None, patch_size=7, search_radius=21):
    img = _to_numpy(image)

    # Convert to uint8 for OpenCV (with proper scaling)
    if img.dtype != np.uint8:
        img_min, img_max = img.min(), img.max()
        img_scaled = ((img - img_min) / (img_max - img_min) * 255).astype(np.uint8)
    else:
        img_scaled = img
        img_min, img_max = 0, 255

    # OpenCV NLM
    template_window = patch_size
    search_window = search_radius * 2 + 1
    denoised = cv2.fastNlMeansDenoising(img_scaled, None, h_param, template_window, search_window)

    # Scale back to original range
    result = denoised.astype(np.float64) / 255 * (img_max - img_min) + img_min
    return result.astype(image.dtype)
```

**Expected Speedup**: 100-1000x

**Risks**:
- OpenCV only supports uint8 input (requires scaling)
- Parameter mapping differences (h_param interpretation)
- Slight numerical differences due to implementation details

---

## Validation Framework

### Test Categories

1. **Synthetic Images**: Known ground truth for quantitative metrics
2. **Real Images**: Representative microscopy images from KINTSUGI workflows
3. **Edge Cases**: Very noisy, very clean, small, large images

### Quality Metrics

| Metric | Description | Acceptance Criteria |
|--------|-------------|---------------------|
| PSNR | Peak Signal-to-Noise Ratio | Δ < 0.5 dB from original |
| SSIM | Structural Similarity Index | Δ < 0.01 from original |
| MSE | Mean Squared Error | Δ < 5% from original |
| Visual | Side-by-side comparison | No visible artifacts |

### Performance Metrics

| Metric | Description | Target |
|--------|-------------|--------|
| Wall Time | Execution time | >5x improvement |
| Memory Peak | Maximum memory usage | No increase |
| Scaling | Time vs image size | Linear or better |

### Test Images

| Image | Size | Noise Level | Purpose |
|-------|------|-------------|---------|
| Synthetic gradient | 256×256 | σ=15, 25, 50 | Quantitative PSNR |
| Synthetic shapes | 512×512 | σ=25 | Edge preservation |
| Real DAPI channel | 1024×1024 | Natural | Practical validation |
| Real marker channel | 2048×2048 | Natural | Large image scaling |
| Stress test | 4096×4096 | σ=25 | Performance limits |

---

## Test Implementation

### Phase 1: Baseline Measurements
1. Run current implementations on all test images
2. Record timing, PSNR, SSIM for each
3. Save denoised outputs for visual comparison

### Phase 2: Optimized Implementation
1. Implement KD-tree version of BM3D-lite
2. Implement OpenCV wrapper for NLM
3. Ensure API compatibility (same function signatures)

### Phase 3: Comparative Testing
1. Run optimized versions on identical test images
2. Compare quality metrics to baseline
3. Compare performance metrics to baseline
4. Generate visual difference maps

### Phase 4: Validation Criteria

**PASS conditions**:
- PSNR difference < 0.5 dB on all synthetic images
- SSIM difference < 0.01 on all images
- No visible artifacts in visual inspection
- Performance improvement > 5x

**FAIL conditions**:
- PSNR drops by > 1 dB on any image
- Visible artifacts (blocking, ringing, etc.)
- Edge blurring significantly worse than original

---

## Rollback Plan

If optimizations fail validation:
1. Keep original implementations as default
2. Add optimized versions as `_fast` variants with appropriate warnings
3. Document known limitations in docstrings

---

## Timeline

| Phase | Tasks | Est. Time |
|-------|-------|-----------|
| 1 | Create test framework, generate test images | 1-2 hours |
| 2 | Implement BM3D-lite optimization | 1-2 hours |
| 3 | Implement NLM replacement | 1 hour |
| 4 | Run validation tests | 1 hour |
| 5 | Document results, create PR | 1 hour |

---

## Files to Create/Modify

### New Files
- `tests/test_denoise_optimization.py` - Validation test suite
- `src/kintsugi/denoise/patch_based_optimized.py` - Optimized implementations (temporary for testing)

### Modified Files (after validation)
- `src/kintsugi/denoise/patch_based.py` - Replace with optimized versions

---

## Success Criteria Summary

| Algorithm | Quality | Performance | Decision |
|-----------|---------|-------------|----------|
| BM3D-lite KD-tree | PSNR ±0.5dB, SSIM ±0.01 | >10x faster | Merge if both pass |
| NLM OpenCV | PSNR ±0.5dB, SSIM ±0.01 | >50x faster | Merge if both pass |
