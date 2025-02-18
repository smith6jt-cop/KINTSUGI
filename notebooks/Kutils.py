import zarr
import os
from skimage import io
from skimage import morphology, filters
from skimage.exposure import rescale_intensity
import numpy as np
from scipy.ndimage import median_filter, uniform_filter
import dask.array as da


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
    zarr_dir = os.path.join(out_dir, "zarr", f'{name}')
    zarr_out = convert_tiff_to_zarr(path, zarr_dir, chunk_size)

    return zarr_out

def convert_array_to_zarr(array, out_dir, name, chunk_size=(1000, 1000)):
    """Convert array to Zarr format with optimal chunking"""
    os.makedirs(os.path.join(out_dir, "zarr"),exist_ok=True)
    zarr_dir = os.path.join(out_dir, "zarr", f'{name}')
    # Create zarr array with compression
    z = zarr.open(zarr_dir, mode='w', shape=array.shape,
                  chunks=chunk_size, dtype=array.dtype,
                  compressor=zarr.Blosc(cname='zstd', clevel=3))
    
    # Write data
    z[:] = array
    return z

def clean(im_cl, background_threshold:int=1000, smooth:bool=False, smooth_threshold:int = 1000, footprint:int=1, remove_small:bool=False, View_original:bool=False):

    background_threshold = max(1, background_threshold)
    smooth_threshold = max(1, smooth_threshold)
    footprint = max(1, footprint)
    
    processed = im_cl.copy()
    
    background = processed <= background_threshold
    processed[background] = 0

    if smooth:
        transition_mask = (processed > 0) & (processed < 2*smooth_threshold)
        processed[transition_mask] = median_filter(processed, size=footprint)[transition_mask]
    
    if remove_small:
        signal_mask = processed > 0
        signal_mask = morphology.remove_small_objects(signal_mask, min_size=30, connectivity=2)
        processed[~signal_mask] = 0
        processed = morphology.closing(processed, morphology.disk(1))
    
    result = im_cl if View_original else processed

    print(f"Original range: {im_cl.max()}/{im_cl.min()}")
    print(f"Output range: {result.max()}/{result.min()}")
    
    return result

def process_chunk(image, subtracted, smooth_low, low_size, smooth_high, high_size, low_percentile_value, high_percentile_value):
    
    if smooth_low:
        low_mask = image < low_percentile_value
        low_mask = morphology.binary_dilation(low_mask, morphology.disk(1))
        subtracted = da.where(low_mask, (uniform_filter(subtracted, size=low_size)), subtracted)

    if smooth_high:
        high_mask = image > high_percentile_value
        high_mask = morphology.binary_dilation(high_mask, morphology.disk(high_size))
        subtracted = da.where(high_mask, 0, subtracted)

    return subtracted