"""
Multi-GPU utilities for KINTSUGI.

Provides centralized GPU detection, device selection, and multi-GPU coordination
for both PyTorch and CuPy workloads.

Features:
- Auto-detect all available non-integrated NVIDIA GPUs
- Distribute work across multiple GPUs
- DataParallel wrapper for PyTorch models
- Per-GPU CuPy accelerator instances
- Detailed timing and progress logging for hardware detection

Usage:
    from kintsugi.gpu import GPUManager, get_gpu_manager

    # Get singleton manager (with progress logging)
    gpu = get_gpu_manager(verbose=True)

    # Check available GPUs
    print(f"Available GPUs: {gpu.device_count}")
    print(f"Device IDs: {gpu.device_ids}")

    # Wrap PyTorch model for multi-GPU inference
    model = gpu.wrap_model(model)

    # Get CuPy device for specific GPU
    with gpu.cupy_device(device_id=0):
        # CuPy operations on GPU 0
        pass
"""

from __future__ import annotations

import time
import warnings
from collections.abc import Callable
from contextlib import contextmanager
from dataclasses import dataclass, field
from functools import lru_cache
from typing import TYPE_CHECKING, TypeVar

if TYPE_CHECKING:
    import torch
    import torch.nn as nn

# Type variable for model wrapping
T = TypeVar("T")


def _timed_step(
    step_name: str,
    func: Callable[[], T],
    verbose: bool = False,
    threshold_warn: float = 5.0,
) -> tuple[T, float]:
    """
    Execute a function with timing, optionally logging progress.

    Parameters
    ----------
    step_name : str
        Name of the step for logging
    func : Callable
        Function to execute
    verbose : bool
        If True, print progress messages
    threshold_warn : float
        Warn if step takes longer than this (seconds)

    Returns
    -------
    tuple
        (result, elapsed_time)
    """
    if verbose:
        print(f"  → {step_name}...", end="", flush=True)

    t0 = time.perf_counter()
    try:
        result = func()
        elapsed = time.perf_counter() - t0

        if verbose:
            status = "✓" if elapsed < threshold_warn else "⚠"
            print(f" {status} ({elapsed:.1f}s)")
        elif elapsed > threshold_warn:
            warnings.warn(
                f"GPU detection step '{step_name}' took {elapsed:.1f}s", stacklevel=3
            )

        return result, elapsed

    except Exception as e:
        elapsed = time.perf_counter() - t0
        if verbose:
            print(f" ✗ ({elapsed:.1f}s) - {type(e).__name__}")
        raise

# Global singleton
_gpu_manager: GPUManager | None = None


@dataclass
class GPUInfo:
    """Information about a single GPU device."""

    device_id: int
    name: str
    total_memory: int  # bytes
    compute_capability: tuple[int, int]
    is_integrated: bool = False
    multi_processor_count: int = 0


@dataclass
class GPUManager:
    """
    Centralized GPU management for multi-GPU workloads.

    Automatically detects all available NVIDIA GPUs and provides utilities
    for distributing work across them.

    Attributes
    ----------
    device_ids : list[int]
        List of available GPU device IDs (non-integrated only)
    device_count : int
        Number of available GPUs
    gpu_info : list[GPUInfo]
        Detailed information about each GPU
    pytorch_available : bool
        Whether PyTorch CUDA is available
    cupy_available : bool
        Whether CuPy is available
    detection_times : dict[str, float]
        Timing breakdown for each detection step (for diagnostics)
    """

    device_ids: list[int] = field(default_factory=list)
    gpu_info: list[GPUInfo] = field(default_factory=list)
    pytorch_available: bool = False
    cupy_available: bool = False
    detection_times: dict[str, float] = field(default_factory=dict)
    _initialized: bool = False
    _verbose: bool = False

    def __post_init__(self):
        if not self._initialized:
            self._detect_gpus(verbose=self._verbose)
            self._initialized = True

    @property
    def device_count(self) -> int:
        """Number of available non-integrated GPUs."""
        return len(self.device_ids)

    @property
    def has_gpu(self) -> bool:
        """Whether any GPU is available."""
        return self.device_count > 0

    @property
    def primary_device(self) -> int:
        """Primary GPU device ID (first available)."""
        return self.device_ids[0] if self.device_ids else 0

    def _detect_gpus(self, verbose: bool = False) -> None:
        """
        Detect all available NVIDIA GPUs, excluding integrated graphics.

        Parameters
        ----------
        verbose : bool
            If True, print detailed progress with timing for each step.
            Useful for diagnosing slow hardware detection.
        """
        total_start = time.perf_counter()
        self.device_ids = []
        self.gpu_info = []
        self.detection_times = {}

        if verbose:
            print("Hardware Detection (detailed timing):")

        # =====================================================================
        # Step 1: Try PyTorch CUDA detection (typically faster, more common)
        # =====================================================================
        try:
            # Import PyTorch (can be slow on networked filesystems)
            _, import_time = _timed_step(
                "Import PyTorch",
                lambda: __import__("torch"),
                verbose=verbose,
                threshold_warn=10.0,
            )
            self.detection_times["pytorch_import"] = import_time

            import torch

            # Check CUDA availability (initializes CUDA context)
            def check_cuda():
                return torch.cuda.is_available()

            cuda_available, cuda_time = _timed_step(
                "Check CUDA availability",
                check_cuda,
                verbose=verbose,
                threshold_warn=5.0,
            )
            self.detection_times["cuda_check"] = cuda_time

            if cuda_available:
                self.pytorch_available = True

                # Get device count
                count, count_time = _timed_step(
                    "Get device count",
                    torch.cuda.device_count,
                    verbose=verbose,
                    threshold_warn=2.0,
                )
                self.detection_times["device_count"] = count_time

                # Query each device (can be slow for many GPUs)
                for i in range(count):

                    def get_props(dev_id=i):
                        return torch.cuda.get_device_properties(dev_id)

                    props, prop_time = _timed_step(
                        f"Query GPU {i} properties",
                        get_props,
                        verbose=verbose,
                        threshold_warn=3.0,
                    )
                    self.detection_times[f"gpu_{i}_props"] = prop_time

                    name = props.name.lower()

                    # Skip integrated GPUs (Intel, AMD APU, etc.)
                    is_integrated = any(
                        x in name
                        for x in ["intel", "integrated", "amd radeon graphics", "vega"]
                    )

                    if is_integrated:
                        if verbose:
                            print(f"    [Skip] GPU {i}: integrated graphics")
                        continue

                    info = GPUInfo(
                        device_id=i,
                        name=props.name,
                        total_memory=props.total_memory,
                        compute_capability=(props.major, props.minor),
                        is_integrated=False,
                        multi_processor_count=props.multi_processor_count,
                    )
                    self.device_ids.append(i)
                    self.gpu_info.append(info)

        except ImportError:
            if verbose:
                print("  → PyTorch not installed, skipping")
            self.detection_times["pytorch_import"] = 0.0

        # =====================================================================
        # Step 2: Check CuPy availability (needed for GPU-accelerated operations)
        # =====================================================================
        try:
            # Import CuPy (can trigger JIT compilation on first use)
            _, cupy_import_time = _timed_step(
                "Import CuPy",
                lambda: __import__("cupy"),
                verbose=verbose,
                threshold_warn=30.0,  # CuPy first import can be very slow
            )
            self.detection_times["cupy_import"] = cupy_import_time

            import cupy

            self.cupy_available = True

            # If PyTorch didn't find GPUs, try CuPy detection
            if not self.device_ids:
                if verbose:
                    print("  → PyTorch found no GPUs, trying CuPy detection...")

                count, count_time = _timed_step(
                    "CuPy get device count",
                    cupy.cuda.runtime.getDeviceCount,
                    verbose=verbose,
                    threshold_warn=5.0,
                )
                self.detection_times["cupy_device_count"] = count_time

                for i in range(count):

                    def get_cupy_props(dev_id=i):
                        with cupy.cuda.Device(dev_id):
                            return cupy.cuda.runtime.getDeviceProperties(dev_id)

                    props, prop_time = _timed_step(
                        f"CuPy query GPU {i}",
                        get_cupy_props,
                        verbose=verbose,
                        threshold_warn=5.0,
                    )
                    self.detection_times[f"cupy_gpu_{i}_props"] = prop_time

                    name = (
                        props["name"].decode()
                        if isinstance(props["name"], bytes)
                        else props["name"]
                    )

                    # Skip integrated
                    if any(x in name.lower() for x in ["intel", "integrated"]):
                        if verbose:
                            print(f"    [Skip] GPU {i}: integrated graphics")
                        continue

                    info = GPUInfo(
                        device_id=i,
                        name=name,
                        total_memory=props["totalGlobalMem"],
                        compute_capability=(props["major"], props["minor"]),
                        is_integrated=False,
                        multi_processor_count=props["multiProcessorCount"],
                    )
                    self.device_ids.append(i)
                    self.gpu_info.append(info)

        except ImportError:
            if verbose:
                print("  → CuPy not installed")
            self.detection_times["cupy_import"] = 0.0

        except Exception as e:
            warnings.warn(f"Error detecting CuPy GPUs: {e}", stacklevel=2)

        # =====================================================================
        # Summary
        # =====================================================================
        total_time = time.perf_counter() - total_start
        self.detection_times["total"] = total_time

        if verbose:
            print(f"\n  Total detection time: {total_time:.1f}s")
            if total_time > 30:
                print("  ⚠ Detection took >30s. Consider:")
                print("    - Pre-importing torch/cupy in earlier cells")
                print("    - Checking network filesystem latency")
                print("    - Using 'verbose=True' to identify slow steps")

    def wrap_model(
        self,
        model: nn.Module,
        device_ids: list[int] | None = None,
    ) -> nn.Module:
        """
        Wrap a PyTorch model for multi-GPU inference using DataParallel.

        Parameters
        ----------
        model : nn.Module
            PyTorch model to wrap
        device_ids : list[int], optional
            Specific GPU IDs to use. If None, uses all available GPUs.

        Returns
        -------
        nn.Module
            Wrapped model (DataParallel if multiple GPUs, original if single/no GPU)
        """
        if not self.pytorch_available:
            return model

        import torch
        import torch.nn as nn

        if device_ids is None:
            device_ids = self.device_ids

        if len(device_ids) == 0:
            # No GPU - return model as-is (CPU)
            return model
        elif len(device_ids) == 1:
            # Single GPU - just move to device
            device = torch.device(f"cuda:{device_ids[0]}")
            return model.to(device)
        else:
            # Multiple GPUs - use DataParallel
            device = torch.device(f"cuda:{device_ids[0]}")
            model = model.to(device)
            model = nn.DataParallel(model, device_ids=device_ids)
            return model

    def get_torch_device(self, device_id: int | None = None) -> torch.device:
        """
        Get a PyTorch device object.

        Parameters
        ----------
        device_id : int, optional
            Specific GPU ID. If None, uses primary GPU or CPU.

        Returns
        -------
        torch.device
            PyTorch device object
        """
        import torch

        if not self.has_gpu:
            return torch.device("cpu")

        if device_id is None:
            device_id = self.primary_device

        return torch.device(f"cuda:{device_id}")

    @contextmanager
    def cupy_device(self, device_id: int | None = None):
        """
        Context manager for CuPy operations on a specific GPU.

        Parameters
        ----------
        device_id : int, optional
            GPU device ID. If None, uses primary GPU.

        Yields
        ------
        cupy.cuda.Device
            CuPy device context
        """
        if not self.cupy_available:
            yield None
            return

        import cupy

        if device_id is None:
            device_id = self.primary_device

        with cupy.cuda.Device(device_id):
            yield cupy.cuda.Device(device_id)

    def distribute_work(
        self,
        items: list[T],
        device_ids: list[int] | None = None,
    ) -> list[tuple[int, list[T]]]:
        """
        Distribute work items across GPUs.

        Parameters
        ----------
        items : list
            Items to distribute
        device_ids : list[int], optional
            GPUs to use. If None, uses all available.

        Returns
        -------
        list[tuple[int, list]]
            List of (device_id, items) tuples for each GPU
        """
        if device_ids is None:
            device_ids = self.device_ids

        if not device_ids:
            # No GPUs - return all items for CPU
            return [(-1, items)]

        n_devices = len(device_ids)
        n_items = len(items)

        # Distribute items evenly
        items_per_device = n_items // n_devices
        remainder = n_items % n_devices

        result = []
        start = 0
        for i, device_id in enumerate(device_ids):
            # Give extra items to first GPUs if not evenly divisible
            count = items_per_device + (1 if i < remainder else 0)
            end = start + count
            if count > 0:
                result.append((device_id, items[start:end]))
            start = end

        return result

    def get_memory_info(self, device_id: int | None = None) -> dict[str, int]:
        """
        Get GPU memory information.

        Parameters
        ----------
        device_id : int, optional
            GPU device ID. If None, uses primary GPU.

        Returns
        -------
        dict
            Dictionary with 'total', 'used', 'free' memory in bytes
        """
        if device_id is None:
            device_id = self.primary_device

        if self.pytorch_available:
            import torch

            if torch.cuda.is_available():
                with torch.cuda.device(device_id):
                    total = torch.cuda.get_device_properties(device_id).total_memory
                    allocated = torch.cuda.memory_allocated(device_id)
                    cached = torch.cuda.memory_reserved(device_id)
                    return {
                        "total": total,
                        "used": allocated,
                        "cached": cached,
                        "free": total - cached,
                    }

        if self.cupy_available:
            import cupy

            with cupy.cuda.Device(device_id):
                mempool = cupy.get_default_memory_pool()
                total = cupy.cuda.runtime.getDeviceProperties(device_id)["totalGlobalMem"]
                used = mempool.used_bytes()
                return {
                    "total": total,
                    "used": used,
                    "cached": mempool.total_bytes(),
                    "free": total - used,
                }

        return {"total": 0, "used": 0, "cached": 0, "free": 0}

    def summary(self) -> str:
        """Get a summary string of available GPUs."""
        if not self.has_gpu:
            return "No GPUs available (CPU mode)"

        lines = [f"Found {self.device_count} GPU(s):"]
        for info in self.gpu_info:
            mem_gb = info.total_memory / (1024**3)
            lines.append(
                f"  [{info.device_id}] {info.name} "
                f"({mem_gb:.1f} GB, CC {info.compute_capability[0]}.{info.compute_capability[1]})"
            )
        return "\n".join(lines)

    def timing_report(self) -> str:
        """
        Get a formatted timing report from the last GPU detection.

        Useful for diagnosing slow hardware detection on HPC systems.

        Returns
        -------
        str
            Formatted timing breakdown string
        """
        if not self.detection_times:
            return "No timing data available (detection not run with verbose=True)"

        lines = ["GPU Detection Timing Report:"]
        lines.append("-" * 40)

        # Sort by time (slowest first)
        sorted_times = sorted(
            [(k, v) for k, v in self.detection_times.items() if k != "total"],
            key=lambda x: x[1],
            reverse=True,
        )

        for step, elapsed in sorted_times:
            if elapsed > 0.01:  # Only show meaningful times
                status = "⚠ SLOW" if elapsed > 5.0 else ""
                lines.append(f"  {step}: {elapsed:.2f}s {status}")

        total = self.detection_times.get("total", 0)
        lines.append("-" * 40)
        lines.append(f"  TOTAL: {total:.2f}s")

        if total > 30:
            lines.append("")
            lines.append("⚠ Detection took >30s. Possible causes:")
            lines.append("  - First PyTorch/CuPy import on networked filesystem")
            lines.append("  - CuPy JIT compilation (first use)")
            lines.append("  - CUDA driver initialization")
            lines.append("  Tip: Pre-import torch/cupy in an earlier cell")

        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"GPUManager(devices={self.device_ids}, "
            f"pytorch={self.pytorch_available}, "
            f"cupy={self.cupy_available})"
        )


def get_gpu_manager(verbose: bool = False) -> GPUManager:
    """
    Get the global GPU manager singleton.

    Parameters
    ----------
    verbose : bool
        If True and this is the first call (creating the manager),
        print detailed timing information for each detection step.
        Useful for diagnosing slow hardware detection.

    Returns
    -------
    GPUManager
        The global GPU manager instance

    Examples
    --------
    >>> # First call with verbose to see timing breakdown
    >>> gpu = get_gpu_manager(verbose=True)
    Hardware Detection (detailed timing):
      → Import PyTorch... ✓ (0.8s)
      → Check CUDA availability... ✓ (0.3s)
      ...

    >>> # Subsequent calls return cached singleton instantly
    >>> gpu = get_gpu_manager()  # No detection, returns cached
    """
    global _gpu_manager
    if _gpu_manager is None:
        # Create with verbose flag - only affects first creation
        _gpu_manager = GPUManager(_verbose=verbose)
    return _gpu_manager


def reset_gpu_manager() -> None:
    """Reset the global GPU manager (for testing or re-detection)."""
    global _gpu_manager
    _gpu_manager = None


def get_detection_times() -> dict[str, float]:
    """
    Get timing breakdown from the last GPU detection.

    Returns
    -------
    dict[str, float]
        Dictionary mapping step names to elapsed seconds.
        Useful for diagnosing slow hardware detection.

    Examples
    --------
    >>> times = get_detection_times()
    >>> print(f"PyTorch import: {times.get('pytorch_import', 0):.1f}s")
    >>> print(f"CuPy import: {times.get('cupy_import', 0):.1f}s")
    >>> print(f"Total: {times.get('total', 0):.1f}s")
    """
    return get_gpu_manager().detection_times.copy()


# Convenience functions with caching


@lru_cache(maxsize=1)
def get_device_count() -> int:
    """Get number of available GPUs (cached after first call)."""
    return get_gpu_manager().device_count


@lru_cache(maxsize=1)
def get_device_ids() -> tuple[int, ...]:
    """
    Get tuple of available GPU device IDs (cached after first call).

    Note: Returns tuple instead of list for hashability/caching.
    """
    return tuple(get_gpu_manager().device_ids)


@lru_cache(maxsize=1)
def has_multi_gpu() -> bool:
    """Check if multiple GPUs are available (cached after first call)."""
    return get_gpu_manager().device_count > 1


def wrap_model_multi_gpu(model: nn.Module) -> nn.Module:
    """Wrap a PyTorch model for multi-GPU inference."""
    return get_gpu_manager().wrap_model(model)


def clear_detection_cache() -> None:
    """
    Clear cached detection results.

    Call this if GPU configuration changes during runtime
    (rare, but possible in HPC environments with dynamic allocation).
    """
    get_device_count.cache_clear()
    get_device_ids.cache_clear()
    has_multi_gpu.cache_clear()


def prewarm_gpu_imports(verbose: bool = False) -> dict[str, float]:
    """
    Pre-import PyTorch and CuPy to avoid slow imports during hardware detection.

    On HPC systems with networked filesystems, the first import of torch/cupy
    can take 30-60+ seconds. Calling this function early (e.g., in the imports
    cell) amortizes this cost so the hardware detection cell runs quickly.

    Parameters
    ----------
    verbose : bool
        If True, print timing for each import

    Returns
    -------
    dict[str, float]
        Dictionary with import times: 'torch', 'cupy', 'total'

    Examples
    --------
    >>> # Call this early in your notebook imports cell
    >>> from kintsugi.gpu import prewarm_gpu_imports
    >>> times = prewarm_gpu_imports(verbose=True)
    Pre-warming GPU imports...
      → PyTorch... ✓ (1.2s)
      → CuPy... ✓ (0.8s)
    Done (2.0s total)
    """
    times = {}
    total_start = time.perf_counter()

    if verbose:
        print("Pre-warming GPU imports...")

    # Import PyTorch
    try:
        _, torch_time = _timed_step(
            "PyTorch",
            lambda: __import__("torch"),
            verbose=verbose,
            threshold_warn=10.0,
        )
        times["torch"] = torch_time
    except ImportError:
        if verbose:
            print("  → PyTorch not installed")
        times["torch"] = 0.0

    # Import CuPy
    try:
        _, cupy_time = _timed_step(
            "CuPy",
            lambda: __import__("cupy"),
            verbose=verbose,
            threshold_warn=30.0,
        )
        times["cupy"] = cupy_time
    except ImportError:
        if verbose:
            print("  → CuPy not installed")
        times["cupy"] = 0.0

    total = time.perf_counter() - total_start
    times["total"] = total

    if verbose:
        print(f"Done ({total:.1f}s total)")

    return times


def cleanup_gpu_memory(device_ids: list[int] | None = None, verbose: bool = False) -> None:
    """
    Explicitly release GPU memory pools to prevent Jupyter cell hangs.

    Call this at the end of GPU-intensive processing cells to force
    memory cleanup and prevent the cell from hanging.

    Parameters
    ----------
    device_ids : list[int], optional
        GPU device IDs to clean up. If None, cleans up all available GPUs.
    verbose : bool
        If True, print memory freed

    Examples
    --------
    >>> # At the end of a processing cell:
    >>> from kintsugi.gpu import cleanup_gpu_memory
    >>> cleanup_gpu_memory(verbose=True)
    GPU memory freed: 2.45 GB on device 0
    """
    gpu = get_gpu_manager()

    if device_ids is None:
        device_ids = gpu.device_ids

    if not device_ids:
        return

    if not gpu.cupy_available:
        return

    try:
        import gc

        import cupy as cp

        for dev_id in device_ids:
            try:
                with cp.cuda.Device(dev_id):
                    mempool = cp.get_default_memory_pool()
                    pinned_mempool = cp.get_default_pinned_memory_pool()

                    used_before = mempool.used_bytes()

                    mempool.free_all_blocks()
                    pinned_mempool.free_all_blocks()

                    if verbose and used_before > 0:
                        print(f"GPU memory freed: {used_before / 1e9:.2f} GB on device {dev_id}")

            except Exception as e:
                if verbose:
                    print(f"Warning: Could not clean GPU {dev_id}: {e}")

        gc.collect()

    except ImportError:
        pass
    except Exception as e:
        if verbose:
            print(f"Warning: GPU cleanup error: {e}")
