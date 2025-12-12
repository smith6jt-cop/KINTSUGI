"""
RAPIDS GPU-Accelerated Data Processing for KINTSUGI.

Provides transparent GPU acceleration using NVIDIA RAPIDS libraries:
- cuDF: GPU DataFrame library (drop-in pandas replacement)
- cuML: GPU machine learning (drop-in scikit-learn replacement)
- cuGraph: GPU graph analytics

The module provides smart fallbacks to CPU when RAPIDS is not available
or when data size doesn't warrant GPU transfer overhead.

Usage:
    from kintsugi.rapids import get_dataframe_backend, RAPIDSManager

    # Get DataFrame module (cudf if available, pandas otherwise)
    pd = get_dataframe_backend()
    df = pd.read_csv("data.csv")

    # Or use RAPIDSManager for more control
    rapids = RAPIDSManager()
    df = rapids.to_gpu(pandas_df)  # Convert to cuDF if available

    # cuDF-pandas integration (automatic acceleration)
    rapids.enable_cudf_pandas()
    import pandas as pd  # Now GPU-accelerated transparently

References:
    - https://rapids.ai/cudf-pandas/
    - https://docs.rapids.ai/
"""

from __future__ import annotations

import logging
from contextlib import contextmanager
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, TypeVar

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import cudf
    import cuml

logger = logging.getLogger("kintsugi.rapids")

# Type variable for generic DataFrame operations
DF = TypeVar("DF")

# Track RAPIDS availability
_RAPIDS_AVAILABLE = {
    "cudf": False,
    "cuml": False,
    "cugraph": False,
    "cudf_pandas": False,
}

# Try importing RAPIDS components
try:
    import cudf

    _RAPIDS_AVAILABLE["cudf"] = True
except ImportError:
    cudf = None

try:
    import cuml

    _RAPIDS_AVAILABLE["cuml"] = True
except ImportError:
    cuml = None

try:
    import cugraph

    _RAPIDS_AVAILABLE["cugraph"] = True
except ImportError:
    cugraph = None

# Check for cudf.pandas (RAPIDS 24.02+)
try:
    import cudf.pandas

    _RAPIDS_AVAILABLE["cudf_pandas"] = True
except (ImportError, AttributeError):
    pass


@dataclass
class RAPIDSConfig:
    """Configuration for RAPIDS acceleration."""

    # Enable/disable RAPIDS
    enabled: bool = True

    # Minimum data size to use GPU (smaller data not worth transfer overhead)
    min_rows_for_gpu: int = 10000
    min_cols_for_gpu: int = 10

    # Memory management
    spill_to_host: bool = True  # Spill GPU memory to host when full
    pool_allocator: bool = True  # Use memory pool for faster allocations

    # cuDF-pandas settings
    cudf_pandas_mode: str = "auto"  # "auto", "always", "never"

    # Device selection
    device_id: int | None = None  # None = use primary from GPUManager


@dataclass
class RAPIDSInfo:
    """Information about RAPIDS availability and configuration."""

    cudf_available: bool = False
    cuml_available: bool = False
    cugraph_available: bool = False
    cudf_pandas_available: bool = False
    cudf_version: str = ""
    cuml_version: str = ""
    gpu_memory_total: int = 0
    gpu_memory_free: int = 0


class RAPIDSManager:
    """
    Centralized RAPIDS management for GPU-accelerated data processing.

    Provides:
    - Smart switching between cuDF and pandas based on data size
    - cuML-accelerated ML operations
    - Memory management utilities
    - cuDF-pandas transparent acceleration

    Example
    -------
    >>> rapids = RAPIDSManager()
    >>> print(rapids.info())

    >>> # Convert pandas DataFrame to cuDF if beneficial
    >>> gpu_df = rapids.to_gpu(pandas_df)

    >>> # Use cuML for clustering
    >>> if rapids.cuml_available:
    ...     from cuml.cluster import KMeans
    ...     kmeans = KMeans(n_clusters=10)
    ...     labels = kmeans.fit_predict(gpu_df)
    """

    def __init__(self, config: RAPIDSConfig | None = None):
        """
        Initialize RAPIDS manager.

        Parameters
        ----------
        config : RAPIDSConfig, optional
            Configuration options
        """
        self.config = config or RAPIDSConfig()
        self._cudf_pandas_enabled = False

        # Set device if specified
        if self.config.device_id is not None and _RAPIDS_AVAILABLE["cudf"]:
            import cupy as cp

            cp.cuda.Device(self.config.device_id).use()

        # Configure memory pool if requested
        if self.config.pool_allocator and _RAPIDS_AVAILABLE["cudf"]:
            self._setup_memory_pool()

    @property
    def cudf_available(self) -> bool:
        """Check if cuDF is available."""
        return _RAPIDS_AVAILABLE["cudf"] and self.config.enabled

    @property
    def cuml_available(self) -> bool:
        """Check if cuML is available."""
        return _RAPIDS_AVAILABLE["cuml"] and self.config.enabled

    @property
    def cugraph_available(self) -> bool:
        """Check if cuGraph is available."""
        return _RAPIDS_AVAILABLE["cugraph"] and self.config.enabled

    @property
    def cudf_pandas_available(self) -> bool:
        """Check if cuDF-pandas integration is available."""
        return _RAPIDS_AVAILABLE["cudf_pandas"] and self.config.enabled

    def _setup_memory_pool(self):
        """Configure RAPIDS memory pool."""
        try:
            import rmm

            # Use pool allocator for faster allocations
            rmm.reinitialize(pool_allocator=True)

            if self.config.spill_to_host:
                # Enable managed memory that can spill to host
                rmm.reinitialize(managed_memory=True)

            logger.info("RAPIDS memory pool configured")
        except Exception as e:
            logger.warning(f"Could not configure RAPIDS memory pool: {e}")

    def info(self) -> RAPIDSInfo:
        """
        Get RAPIDS availability and configuration info.

        Returns
        -------
        RAPIDSInfo
            Information about RAPIDS components
        """
        info = RAPIDSInfo(
            cudf_available=self.cudf_available,
            cuml_available=self.cuml_available,
            cugraph_available=self.cugraph_available,
            cudf_pandas_available=self.cudf_pandas_available,
        )

        if self.cudf_available:
            import cudf as _cudf

            info.cudf_version = _cudf.__version__

            try:
                import cupy as cp

                mem = cp.cuda.Device().mem_info
                info.gpu_memory_free = mem[0]
                info.gpu_memory_total = mem[1]
            except Exception:
                pass

        if self.cuml_available:
            import cuml as _cuml

            info.cuml_version = _cuml.__version__

        return info

    def summary(self) -> str:
        """Get a summary string of RAPIDS availability."""
        info = self.info()
        lines = ["RAPIDS Status:"]

        if info.cudf_available:
            mem_gb = info.gpu_memory_total / (1024**3)
            free_gb = info.gpu_memory_free / (1024**3)
            lines.append(f"  cuDF: v{info.cudf_version} (GPU: {free_gb:.1f}/{mem_gb:.1f} GB free)")
        else:
            lines.append("  cuDF: Not available")

        if info.cuml_available:
            lines.append(f"  cuML: v{info.cuml_version}")
        else:
            lines.append("  cuML: Not available")

        if info.cugraph_available:
            lines.append("  cuGraph: Available")
        else:
            lines.append("  cuGraph: Not available")

        if info.cudf_pandas_available:
            status = "Enabled" if self._cudf_pandas_enabled else "Available (not enabled)"
            lines.append(f"  cuDF-pandas: {status}")
        else:
            lines.append("  cuDF-pandas: Not available")

        return "\n".join(lines)

    def should_use_gpu(self, data: Any) -> bool:
        """
        Determine if GPU should be used for given data.

        Parameters
        ----------
        data : array-like or DataFrame
            Input data

        Returns
        -------
        bool
            True if GPU should be used
        """
        if not self.cudf_available:
            return False

        # Get data shape
        if hasattr(data, "shape"):
            shape = data.shape
            n_rows = shape[0]
            n_cols = shape[1] if len(shape) > 1 else 1
        elif hasattr(data, "__len__"):
            n_rows = len(data)
            n_cols = 1
        else:
            return False

        return n_rows >= self.config.min_rows_for_gpu or n_cols >= self.config.min_cols_for_gpu

    def to_gpu(self, data: pd.DataFrame) -> cudf.DataFrame | pd.DataFrame:
        """
        Convert pandas DataFrame to cuDF if beneficial.

        Parameters
        ----------
        data : pd.DataFrame
            Input pandas DataFrame

        Returns
        -------
        cudf.DataFrame or pd.DataFrame
            GPU DataFrame if RAPIDS available and data large enough,
            otherwise returns input unchanged
        """
        if not self.should_use_gpu(data):
            return data

        try:
            import cudf

            return cudf.DataFrame.from_pandas(data)
        except Exception as e:
            logger.warning(f"Could not convert to cuDF: {e}")
            return data

    def to_cpu(self, data: Any) -> pd.DataFrame:
        """
        Convert cuDF DataFrame to pandas.

        Parameters
        ----------
        data : cudf.DataFrame or pd.DataFrame
            Input DataFrame

        Returns
        -------
        pd.DataFrame
            Pandas DataFrame
        """
        if hasattr(data, "to_pandas"):
            return data.to_pandas()
        return data

    def get_dataframe_module(self) -> Any:
        """
        Get the appropriate DataFrame module (cudf or pandas).

        Returns
        -------
        module
            cudf if available and enabled, pandas otherwise
        """
        if self.cudf_available:
            import cudf

            return cudf
        return pd

    def enable_cudf_pandas(self) -> bool:
        """
        Enable cuDF-pandas transparent acceleration.

        After calling this, importing pandas will use cuDF acceleration
        transparently. This must be called BEFORE importing pandas in
        user code for full effect.

        Returns
        -------
        bool
            True if enabled successfully
        """
        if not self.cudf_pandas_available:
            logger.warning("cuDF-pandas not available")
            return False

        try:
            import cudf.pandas

            cudf.pandas.install()
            self._cudf_pandas_enabled = True
            logger.info("cuDF-pandas acceleration enabled")
            return True
        except Exception as e:
            logger.warning(f"Could not enable cuDF-pandas: {e}")
            return False

    @contextmanager
    def gpu_context(self, device_id: int | None = None):
        """
        Context manager for GPU operations.

        Parameters
        ----------
        device_id : int, optional
            GPU device ID

        Yields
        ------
        RAPIDSManager
            Self for chaining
        """
        if device_id is not None and self.cudf_available:
            import cupy as cp

            with cp.cuda.Device(device_id):
                yield self
        else:
            yield self

    # cuML convenience methods

    def pairwise_distances(
        self,
        X: np.ndarray,
        Y: np.ndarray | None = None,
        metric: str = "euclidean",
    ) -> np.ndarray:
        """
        Compute pairwise distances using cuML if available.

        Parameters
        ----------
        X : np.ndarray
            First set of points (n_samples_X, n_features)
        Y : np.ndarray, optional
            Second set of points (n_samples_Y, n_features)
        metric : str
            Distance metric

        Returns
        -------
        np.ndarray
            Distance matrix
        """
        if self.cuml_available and self.should_use_gpu(X):
            try:
                import cupy as cp
                from cuml.metrics import pairwise_distances as cuml_pdist

                X_gpu = cp.asarray(X)
                Y_gpu = cp.asarray(Y) if Y is not None else None

                result = cuml_pdist(X_gpu, Y_gpu, metric=metric)
                return cp.asnumpy(result)
            except Exception as e:
                logger.warning(f"cuML pairwise_distances failed, using scipy: {e}")

        # CPU fallback
        from scipy.spatial.distance import cdist

        return cdist(X, Y if Y is not None else X, metric=metric)

    def nearest_neighbors(
        self,
        X: np.ndarray,
        n_neighbors: int = 5,
        metric: str = "euclidean",
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Find nearest neighbors using cuML if available.

        Parameters
        ----------
        X : np.ndarray
            Points to find neighbors for (n_samples, n_features)
        n_neighbors : int
            Number of neighbors
        metric : str
            Distance metric

        Returns
        -------
        distances : np.ndarray
            Distances to neighbors (n_samples, n_neighbors)
        indices : np.ndarray
            Indices of neighbors (n_samples, n_neighbors)
        """
        if self.cuml_available and self.should_use_gpu(X):
            try:
                import cupy as cp
                from cuml.neighbors import NearestNeighbors

                X_gpu = cp.asarray(X)
                nn = NearestNeighbors(n_neighbors=n_neighbors, metric=metric)
                nn.fit(X_gpu)
                distances, indices = nn.kneighbors(X_gpu)

                return cp.asnumpy(distances), cp.asnumpy(indices)
            except Exception as e:
                logger.warning(f"cuML NearestNeighbors failed, using sklearn: {e}")

        # CPU fallback
        from sklearn.neighbors import NearestNeighbors

        nn = NearestNeighbors(n_neighbors=n_neighbors, metric=metric)
        nn.fit(X)
        return nn.kneighbors(X)

    def kmeans(
        self,
        X: np.ndarray,
        n_clusters: int,
        random_state: int = 42,
        max_iter: int = 300,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        K-means clustering using cuML if available.

        Parameters
        ----------
        X : np.ndarray
            Data to cluster (n_samples, n_features)
        n_clusters : int
            Number of clusters
        random_state : int
            Random seed
        max_iter : int
            Maximum iterations

        Returns
        -------
        labels : np.ndarray
            Cluster labels
        centers : np.ndarray
            Cluster centers
        """
        if self.cuml_available and self.should_use_gpu(X):
            try:
                import cupy as cp
                from cuml.cluster import KMeans

                X_gpu = cp.asarray(X)
                kmeans = KMeans(
                    n_clusters=n_clusters,
                    random_state=random_state,
                    max_iter=max_iter,
                )
                labels = kmeans.fit_predict(X_gpu)
                centers = kmeans.cluster_centers_

                return cp.asnumpy(labels), cp.asnumpy(centers)
            except Exception as e:
                logger.warning(f"cuML KMeans failed, using sklearn: {e}")

        # CPU fallback
        from sklearn.cluster import KMeans

        kmeans = KMeans(n_clusters=n_clusters, random_state=random_state, max_iter=max_iter)
        labels = kmeans.fit_predict(X)
        return labels, kmeans.cluster_centers_

    def pca(
        self,
        X: np.ndarray,
        n_components: int = 2,
        random_state: int = 42,
    ) -> np.ndarray:
        """
        PCA dimensionality reduction using cuML if available.

        Parameters
        ----------
        X : np.ndarray
            Data to transform (n_samples, n_features)
        n_components : int
            Number of components

        Returns
        -------
        np.ndarray
            Transformed data (n_samples, n_components)
        """
        if self.cuml_available and self.should_use_gpu(X):
            try:
                import cupy as cp
                from cuml.decomposition import PCA

                X_gpu = cp.asarray(X)
                pca = PCA(n_components=n_components, random_state=random_state)
                result = pca.fit_transform(X_gpu)

                return cp.asnumpy(result)
            except Exception as e:
                logger.warning(f"cuML PCA failed, using sklearn: {e}")

        # CPU fallback
        from sklearn.decomposition import PCA

        pca = PCA(n_components=n_components, random_state=random_state)
        return pca.fit_transform(X)

    def umap(
        self,
        X: np.ndarray,
        n_components: int = 2,
        n_neighbors: int = 15,
        min_dist: float = 0.1,
        random_state: int = 42,
    ) -> np.ndarray:
        """
        UMAP dimensionality reduction using cuML if available.

        Parameters
        ----------
        X : np.ndarray
            Data to transform (n_samples, n_features)
        n_components : int
            Number of components
        n_neighbors : int
            Number of neighbors for manifold construction
        min_dist : float
            Minimum distance in embedding space

        Returns
        -------
        np.ndarray
            Transformed data (n_samples, n_components)
        """
        if self.cuml_available and self.should_use_gpu(X):
            try:
                import cupy as cp
                from cuml.manifold import UMAP

                X_gpu = cp.asarray(X)
                umap = UMAP(
                    n_components=n_components,
                    n_neighbors=n_neighbors,
                    min_dist=min_dist,
                    random_state=random_state,
                )
                result = umap.fit_transform(X_gpu)

                return cp.asnumpy(result)
            except Exception as e:
                logger.warning(f"cuML UMAP failed, trying CPU: {e}")

        # CPU fallback - requires umap-learn package
        try:
            import umap

            reducer = umap.UMAP(
                n_components=n_components,
                n_neighbors=n_neighbors,
                min_dist=min_dist,
                random_state=random_state,
            )
            return reducer.fit_transform(X)
        except ImportError:
            raise ImportError("UMAP requires either cuML or umap-learn package")


# Global singleton
_rapids_manager: RAPIDSManager | None = None


def get_rapids_manager() -> RAPIDSManager:
    """
    Get the global RAPIDS manager singleton.

    Returns
    -------
    RAPIDSManager
        The global RAPIDS manager instance
    """
    global _rapids_manager
    if _rapids_manager is None:
        _rapids_manager = RAPIDSManager()
    return _rapids_manager


def reset_rapids_manager() -> None:
    """Reset the global RAPIDS manager (for testing)."""
    global _rapids_manager
    _rapids_manager = None


# Convenience functions


def get_dataframe_backend() -> Any:
    """
    Get the best available DataFrame backend.

    Returns cudf if available and enabled, pandas otherwise.

    Returns
    -------
    module
        DataFrame module (cudf or pandas)
    """
    return get_rapids_manager().get_dataframe_module()


def is_rapids_available() -> bool:
    """Check if any RAPIDS component is available."""
    return any(_RAPIDS_AVAILABLE.values())


def enable_cudf_pandas() -> bool:
    """
    Enable cuDF-pandas transparent acceleration.

    Should be called early, before importing pandas.

    Returns
    -------
    bool
        True if enabled successfully
    """
    return get_rapids_manager().enable_cudf_pandas()


def rapids_info() -> dict[str, Any]:
    """
    Get RAPIDS availability information.

    Returns
    -------
    dict
        Dictionary with RAPIDS component status
    """
    manager = get_rapids_manager()
    info = manager.info()
    return {
        "cudf": {
            "available": info.cudf_available,
            "version": info.cudf_version,
        },
        "cuml": {
            "available": info.cuml_available,
            "version": info.cuml_version,
        },
        "cugraph": {
            "available": info.cugraph_available,
        },
        "cudf_pandas": {
            "available": info.cudf_pandas_available,
        },
        "gpu_memory": {
            "total_gb": info.gpu_memory_total / (1024**3) if info.gpu_memory_total else 0,
            "free_gb": info.gpu_memory_free / (1024**3) if info.gpu_memory_free else 0,
        },
    }
