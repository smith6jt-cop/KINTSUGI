# KINTSUGI Processing Speed Optimization - Session Summary

**Date:** December 9, 2025
**Context:** Addressing processing time regression from ~15 min to ~60-175 min after recent changes

---

## Background

### Previous Work Completed
1. Optimized GPU BaSiC SVD with power iteration (10x faster than full SVD)
2. Restored ThreadPoolExecutor parallelism with GPU BaSiC (8 workers)
3. Added timing per stage and progress logging (ProgressCounter class)
4. Changed output to TIFF for fair speed comparison (excluding zarr overhead)

### Current Implementation Files
- `src/kintsugi/kcorrect_gpu.py` - GPU BaSiC illumination correction
- `notebooks/Kstitch/stitching.py` - Image stitching with GPU FFT
- `notebooks/2_Cycle_Processing.ipynb` - Main processing pipeline
- `notebooks/1_Single_Channel_Eval.ipynb` - Parameter tuning

---

## Repository Review Findings

### Repositories Analyzed
1. **cecelia** (smith6jt-cop) - R/napari image analysis framework
2. **PyBaSiCCellprofilerPlugin** (smith6jt-cop) - BaSiC for CellProfiler
3. **m2stitch** (smith6jt-cop) - MIST-inspired stitching
4. **mcmicro** (smith6jt-cop) - Nextflow multiplex pipeline
5. **ashlar** (smith6jt76) - Stitching/registration tool
6. **cylinter** (smith6jt) - QC for multiplex imaging
7. **RAPID** (smith6jt-cop) - MATLAB processing pipeline
8. **BaSiCPy** (peng-lab) - JAX-based BaSiC implementation

### Key Optimization Techniques Found

#### From BaSiCPy (peng-lab)
- JAX provides ~6x speedup over CuPy/Numba through JIT compilation
- Device-agnostic arrays with `jnp` (JAX NumPy)
- 3D DCT transforms with `JaxDCT.dct3d()`
- Multi-worker parallelization parameter
- **User Decision:** Stick with CuPy optimization only (no JAX dependency)

#### From MIST (NIST)
- Hybrid CPU-GPU pipelining: **24x speedup** (59x42 tiles in 26 seconds)
- cuFFT with CUDA kernels for NCC computation
- Pipelining overlaps CPU data loading with GPU computation (11.2x speedup)
- Source: [NIST MIST Paper](https://www.nature.com/articles/s41598-017-04567-y)

#### From PyBaSiCCellprofilerPlugin
- **Caching mechanism**: Compute flatfield/darkfield once, reuse across cycles
- Memory efficiency: Store baseline drift as per-image array

#### From RAPID
- Parallel Computing with configurable worker counts
- GPU acceleration via CUDA for deconvolution
- FFT-based convolution (ConvFFT3_S) for 3D operations
- Batch processing of tiles, cycles, regions

---

## Proposed Implementation Plan

### Phase 1: BaSiC Illumination Correction
- Implement caching for flatfield/darkfield
- Optimize DCT/IDCT with CuPy
- Keep power iteration SVD (already implemented)

### Phase 2: Stitching Optimization
- GPU pipelining (overlap CPU/GPU work)
- Batch FFT with async memory transfers
- Pre-allocated GPU memory buffers

### Phase 3: Parallel Processing
- Optimized ThreadPoolExecutor
- GPU for heavy compute, CPU for I/O
- Memory-mapped file I/O

### Phase 4: Output Format
- TIFF as default for speed
- OME-Zarr as optional post-processing

### Expected Gains
| Optimization | Expected Speedup |
|--------------|------------------|
| BaSiC caching | 2-3x |
| GPU pipelining | 2-4x |
| Batch FFT | 1.5-2x |
| Memory pre-allocation | 1.2-1.5x |
| **Combined** | **5-10x potential** |

---

## IMPORTANT: BaSiC Caching Validation Required

### User Concern
The current per-z-plane BaSiC processing may produce better results than caching. Need objective validation before implementing caching.

### Test Design

#### Test Cases
1. **DAPI (CH1)** - High signal, present in all cycles
2. **Blank channels (Cycle 1 & 13, CH2-4)** - Background/noise only
3. **Sparse marker (Cycle 2 CH3)** - Fewer positive cells
4. **Dense marker (Cycle 3 CH3)** - More positive cells

#### Processing Modes to Compare
- **Mode A (Current)**: Compute BaSiC per z-plane individually
- **Mode B (Cached)**: Compute BaSiC once from reference plane, apply to all z-planes

#### Metrics to Evaluate
1. **Intensity Statistics** - Mean, std, min, max, CV across tiles
2. **Flatfield Quality** - Uniformity (std/mean), center-to-edge ratio
3. **Darkfield Quality** - Magnitude and pattern
4. **Corrected Image Quality** - Tile boundary artifacts, inter-tile variation, SNR
5. **Biological Signal** - Segmentation consistency, positive cell detection, blank residuals

#### Success Criteria
- **Pass**: Cached mode within 5% of individual mode on all metrics
- **Fail**: Any metric differs >10% or visible artifacts

---

## Files to Modify

| File | Changes |
|------|---------|
| `src/kintsugi/kcorrect_gpu.py` | Add caching, optimize DCT with CuPy |
| `notebooks/Kstitch/stitching.py` | GPU pipelining, batch FFT, memory pre-allocation |
| `notebooks/Kstitch/_translation_computation.py` | GPU-optimized NCC |
| `notebooks/1_Single_Channel_Eval.ipynb` | Use optimized BaSiC |
| `notebooks/2_Cycle_Processing.ipynb` | Parallel processing, caching, memory mapping |

---

## Constraints
- Must remain pure Python (CuPy for CUDA)
- Must maintain transparent, tunable parameters
- Must preserve current API and notebook structure
- GPU acceleration optional (CPU fallback required)

---

## Next Steps
1. Create and run BaSiC caching validation test
2. Based on results, decide whether to implement caching
3. Implement other optimizations (GPU pipelining, batch FFT, etc.)
4. Benchmark against ~15 min baseline

---

## References
- [MIST: Accurate and Scalable Microscopy Image Stitching Tool](https://www.nature.com/articles/s41598-017-04567-y)
- [BaSiCPy GitHub](https://github.com/peng-lab/BaSiCPy)
- [CuPy Performance Best Practices](https://docs.cupy.dev/en/stable/user_guide/performance.html)
- [JAX GPU Performance Tips](https://docs.jax.dev/en/latest/gpu_performance_tips.html)
