import math

import numpy as np
import xarray as xr


def stitch_images(data_xr, num_cols):
    """Stitch together a stack of different channels from different FOVs into a single 2D image
    for each channel

    Args:
        data_xr (xarray.DataArray):
            xarray containing image data from multiple fovs and channels
        num_cols (int):
            number of images stitched together horizontally

    Returns:
        xarray.DataArray:
            the stitched image data
    """

    num_imgs = data_xr.shape[0]
    num_rows = math.ceil(num_imgs / num_cols)
    row_len = data_xr.shape[1]
    col_len = data_xr.shape[2]

    total_row_len = num_rows * row_len
    total_col_len = num_cols * col_len

    stitched_data = np.zeros(
        (1, total_row_len, total_col_len, data_xr.shape[3]), dtype=data_xr.dtype
    )

    img_idx = 0
    for row in range(num_rows):
        for col in range(num_cols):
            stitched_data[
                0, row * row_len : (row + 1) * row_len, col * col_len : (col + 1) * col_len, :
            ] = data_xr[img_idx, ...]
            img_idx += 1
            if img_idx == num_imgs:
                break

    stitched_xr = xr.DataArray(
        stitched_data,
        coords=[["stitched_image"], range(total_row_len), range(total_col_len), data_xr.channels],
        dims=["fovs", "rows", "cols", "channels"],
    )
    return stitched_xr
