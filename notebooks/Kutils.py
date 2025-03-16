import zarr
import os
from skimage import io
from skimage import morphology, filters
from skimage.exposure import rescale_intensity
import numpy as np
from scipy.ndimage import median_filter, uniform_filter
import dask_image.ndfilters
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

def clean(im_cl, backgrnd_thresh:int=100, smooth:bool=False, smooth_thresh:int = 100, footprint:int=1, remove_small:bool=False, small_size:int=30, View_original:bool=False):

    backgrnd_thresh = max(1, backgrnd_thresh)
    smooth_thresh = max(1, smooth_thresh)
    footprint = max(1, footprint)
    
    processed = da.Array.copy(im_cl)
    
    processed = da.where(processed <= backgrnd_thresh, 0, processed)

    if smooth:
        transition_mask = da.where((processed < smooth_thresh), processed, 0)
        processed = da.where(transition_mask, dask_image.ndfilters.median_filter(processed, size=footprint), processed)
    
    if remove_small:
        # Create a boolean mask
        binary_mask = processed > 0
        
        # Process each chunk with remove_small_objects
        def remove_small_chunk(chunk):
            return morphology.remove_small_objects(chunk, min_size=small_size, connectivity=2)
        
        # Apply to each block and maintain boolean dtype
        filtered_mask = binary_mask.map_blocks(remove_small_chunk, dtype=bool)
        
        # Apply mask to processed data
        processed = da.where(filtered_mask, processed, 0)
        
        # Apply closing operation through map_blocks
        def apply_closing(chunk):
            return morphology.closing(chunk, morphology.disk(1))
        
        processed = processed.map_blocks(apply_closing, dtype=processed.dtype)
    
    result = im_cl if View_original else processed

    print(f"Original range: {im_cl.compute().max()}/{im_cl.compute().min()}")
    print(f"Output range: {result.compute().max()}/{result.compute().min()}")
    
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
from skimage.filters import threshold_otsu


def ini_params_SNR(im1, im2, blank_clip_factor:int=10000, blank_scale_factor:float=0.8, smooth_low:bool=True, smooth_high:bool=True, show_labels:bool=True, smoothing_size:int=3):
    from skimage.measure import label, regionprops
    import pandas as pd
    im2_clip = np.clip(im2, blank_clip_factor, im2.max())
    im2_clip[im2_clip <= blank_clip_factor] = 0

    im3 = im1 - (np.minimum(im1, im2_clip * blank_scale_factor))
    if smooth_low:
        low_values_mask = im3 < np.percentile(im3, 25)
        im3[low_values_mask] = median_filter(im3, size=smoothing_size)[low_values_mask]
    
    if smooth_high:
        high_values_mask = im3 > np.percentile(im3, 95)
        im3[high_values_mask] = median_filter(im3, size=smoothing_size)[high_values_mask]
    mask = im3 == 0
    im4 = im3.copy()
    
    binary_image = im4  <= threshold_otsu(im4 )
    # edge_image = sobel(binary_image)
    seg_mask = label(binary_image, background=0)
    properties = regionprops(seg_mask, intensity_image=im4 )
    statistics = {
   
    'area':       [p.area               for p in properties if p.area<800],
    'mean':       [p.mean_intensity     for p in properties if p.area<800]
    }
    df = pd.DataFrame(statistics)
    
    MFI = np.asarray(df['mean'])
    aX = MFI.flatten()
    # compute 20 largest values in aX
    top20 = np.sort(aX)[-20:]
    # compute the mean of bottom 10th percentile of aX
    btm10 = np.sort(aX)[:int(len(aX)*0.1)]
    top20btm10 = np.mean(top20)/np.mean(btm10)
    

    print(f'SNR is {top20btm10}') 
    print(df.describe())

    if show_labels:
        return label(binary_image)
    else:
        return im4
