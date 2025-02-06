import zarr
import os
from skimage import io


def convert_tiff_to_zarr(tiff_path, zarr_path, chunk_size=(1000, 1000)):
    """Convert TIFF to Zarr format with optimal chunking"""
    # Load TIFF
    tiff_data = io.imread(tiff_path)
    
    # Create zarr array with compression
    z = zarr.open(zarr_path, mode='w', shape=tiff_data.shape,
                  chunks=chunk_size, dtype=tiff_data.dtype,
                  compressor=zarr.Blosc(cname='zstd', clevel=3))
    
    # Write data
    z[:] = tiff_data
    return z

def load_images_as_zarr(path, out_dir, chunk_size=(1000, 1000)):
    """Load both signal and blank images as zarr arrays"""
    name = os.path.basename(path).split('.')[0]
    os.makedirs(os.path.join(out_dir, "zarr"),exist_ok=True)
    zarr_dir = os.path.join(out_dir, "zarr", f'{name}.zarr')
    zarr_out = convert_tiff_to_zarr(path, zarr_dir, chunk_size)

    return zarr_out

def convert_array_to_zarr(array, zarr_path, chunk_size=(1000, 1000)):
    """Convert array to Zarr format with optimal chunking"""

    # Create zarr array with compression
    z = zarr.open(zarr_path, mode='w', shape=array.shape,
                  chunks=chunk_size, dtype=array.dtype,
                  compressor=zarr.Blosc(cname='zstd', clevel=3))
    
    # Write data
    z[:] = array
    return z