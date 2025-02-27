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

def clean(im_cl, background_threshold:int=1000, smooth:bool=False, smooth_threshold:int = 1000, footprint:int=1, remove_small:bool=False, small_size:int=30, View_original:bool=False):

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
        signal_mask = morphology.remove_small_objects(signal_mask, min_size=small_size, connectivity=2)
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

def ini_params(im1, im2, blank_clip_factor:int=0, blank_scale_factor:float=0.0, smooth_low:bool=False, low_size:int=1, smooth_high:bool=False, high_size:int=1, erosion:int=0, View_original:bool=False):
    im2_clip = np.clip(im2, blank_clip_factor, im2.max())
    im2_clip[im2_clip <= blank_clip_factor] = 0
    im3 = im1 - (np.minimum(im1, im2_clip * blank_scale_factor))
    if smooth_low:
        low_values_mask = im1 < np.percentile(im1.ravel(), 60)
        low_values_mask = morphology.binary_dilation(low_values_mask, morphology.disk(1))
        # im3[low_values_mask] = 0
        im3[low_values_mask] = uniform_filter(im3, size=low_size)[low_values_mask]
    
    if smooth_high:
        high_values_mask = im1 > np.percentile(im1.ravel(), 90)
        high_values_mask = morphology.binary_dilation(high_values_mask, morphology.disk(high_size))
        im3[high_values_mask] = 0
        # im3[high_values_mask] = uniform_filter(im3, size=smoothing_size)[high_values_mask]
    if View_original:
        im4 =im1.copy()
    else:
        im4 = morphology.erosion(im3, morphology.disk(erosion))  
    return im4


def ini_params_full(im1, im2, blank_clip_factor:int=0, blank_scale_factor:float=0.0, smooth_low:bool=False, smooth_high:bool=False, smoothing_size:int=1):
    
    if smoothing_size <= 0:
        smoothing_size = 1

    im3 = da.empty(
        shape=im2.shape,
        dtype=im2.dtype
    )

    def process_chunk(blank_chunk, signal_chunk):
        blank_chunk = da.clip(blank_chunk, blank_clip_factor, blank_chunk.max())
        blank_chunk = da.where(blank_chunk <= blank_clip_factor, 0, blank_chunk)

        processed_chunk = signal_chunk - da.minimum(signal_chunk, blank_chunk * blank_scale_factor)

        if smooth_low:
            low_mask = processed_chunk < da.percentile(processed_chunk.ravel(), 60)
            low_mask = morphology.binary_opening(low_mask, morphology.disk(2))
            low_smoothed = median_filter(processed_chunk, size=smoothing_size)
            processed_chunk = da.where(low_mask, low_smoothed, processed_chunk)

        if smooth_high:
            high_mask = processed_chunk > da.percentile(processed_chunk.ravel(), 90)
            high_mask = morphology.binary_opening(high_mask, morphology.disk(1))
            high_smoothed = median_filter(processed_chunk, size=smoothing_size)
            processed_chunk = da.where(high_mask, high_smoothed, processed_chunk)

        return processed_chunk

    im3 = da.map_blocks(process_chunk, im2, im1, dtype=im2.dtype)
    im4 = im3.copy()
    return im4
