"""
Multi-GPU utilities for KINTSUGI.

Provides centralized GPU detection, device selection, and multi-GPU coordination
for both PyTorch and CuPy workloads.

Features:
- Auto-detect all available non-integrated NVIDIA GPUs
- Distribute work across multiple GPUs
- DataParallel wrapper for PyTorch models
- Per-GPU CuPy accelerator instances

Usage:
    from kintsugi.gpu import GPUManager, get_gpu_manager

    # Get singleton manager
    gpu = get_gpu_manager()

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

import warnings
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, TypeVar

if TYPE_CHECKING:
    import torch
    import torch.nn as nn

# Type variable for model wrapping
T = TypeVar("T")

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
    """

    device_ids: list[int] = field(default_factory=list)
    gpu_info: list[GPUInfo] = field(default_factory=list)
    pytorch_available: bool = False
    cupy_available: bool = False
    _initialized: bool = False

    def __post_init__(self):
        if not self._initialized:
            self._detect_gpus()
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

    def _detect_gpus(self) -> None:
        """Detect all available NVIDIA GPUs, excluding integrated graphics."""
        self.device_ids = []
        self.gpu_info = []

        # Try PyTorch first (more common)
        try:
            import torch

            if torch.cuda.is_available():
                self.pytorch_available = True
                count = torch.cuda.device_count()

                for i in range(count):
                    props = torch.cuda.get_device_properties(i)
                    name = props.name.lower()

                    # Skip integrated GPUs (Intel, AMD APU, etc.)
                    is_integrated = any(
                        x in name for x in ["intel", "integrated", "amd radeon graphics", "vega"]
                    )

                    if is_integrated:
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
            pass

        # Check CuPy availability
        try:
            import cupy

            self.cupy_available = True

            # If PyTorch didn't find GPUs, try CuPy
            if not self.device_ids:
                count = cupy.cuda.runtime.getDeviceCount()
                for i in range(count):
                    with cupy.cuda.Device(i):
                        props = cupy.cuda.runtime.getDeviceProperties(i)
                        name = (
                            props["name"].decode()
                            if isinstance(props["name"], bytes)
                            else props["name"]
                        )

                        # Skip integrated
                        if any(x in name.lower() for x in ["intel", "integrated"]):
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
            pass
        except Exception as e:
            warnings.warn(f"Error detecting CuPy GPUs: {e}", stacklevel=2)

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

    def __repr__(self) -> str:
        return (
            f"GPUManager(devices={self.device_ids}, "
            f"pytorch={self.pytorch_available}, "
            f"cupy={self.cupy_available})"
        )


def get_gpu_manager() -> GPUManager:
    """
    Get the global GPU manager singleton.

    Returns
    -------
    GPUManager
        The global GPU manager instance
    """
    global _gpu_manager
    if _gpu_manager is None:
        _gpu_manager = GPUManager()
    return _gpu_manager


def reset_gpu_manager() -> None:
    """Reset the global GPU manager (for testing)."""
    global _gpu_manager
    _gpu_manager = None


# Convenience functions


def get_device_count() -> int:
    """Get number of available GPUs."""
    return get_gpu_manager().device_count


def get_device_ids() -> list[int]:
    """Get list of available GPU device IDs."""
    return get_gpu_manager().device_ids


def has_multi_gpu() -> bool:
    """Check if multiple GPUs are available."""
    return get_gpu_manager().device_count > 1


def wrap_model_multi_gpu(model: nn.Module) -> nn.Module:
    """Wrap a PyTorch model for multi-GPU inference."""
    return get_gpu_manager().wrap_model(model)
