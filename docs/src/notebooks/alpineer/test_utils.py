import os

import numpy as np
import xarray as xr

from alpineer import image_utils, io_utils, tiff_utils


def _make_blank_file(folder, name) -> None:
    with open(os.path.join(folder, name), "w"):
        pass


def _make_small_file(folder: str, name: str) -> None:
    """Creates small file.  Creating a blank file will cause a stall for 0-size checking

    Args:
        folder (str):
            Folder to store file in
        name (str):
            Name of small file
    """
    with open(os.path.join(folder, name), "w") as f:
        f.write("a")


def gen_fov_chan_names(num_fovs, num_chans, return_imgs=False, use_delimiter=False):
    """Generate fov and channel names

    Names have the format 'fov0', 'fov1', ..., 'fovN' for fovs and 'chan0', 'chan1', ...,
    'chanM' for channels.

    Args:
        num_fovs (int):
            Number of fov names to create
        num_chans (int):
            Number of channel names to create
        return_imgs (bool):
            Return 'chanK.tiff' as well if True.  Default is False
        use_delimiter (bool):
            Appends '_otherinfo' to the first fov. Useful for testing fov id extraction from
            filenames.  Default is False

    Returns:
        tuple (list, list) or (list, list, list):
            If return_imgs is False, only fov and channel names are returned
            If return_imgs is True, image names will also be returned
    """
    fovs = [f"fov{i}" for i in range(num_fovs)]
    if use_delimiter:
        fovs[0] = f"{fovs[0]}_otherinfo"
    chans = [f"chan{i}" for i in range(num_chans)]

    if return_imgs:
        imgs = [f"{chan}.tiff" for chan in chans]
        return fovs, chans, imgs
    else:
        return fovs, chans


# required metadata for mibitiff writing
MIBITIFF_METADATA = {
    "run": "20180703_1234_test",
    "date": "2017-09-16T15:26:00",
    "coordinates": (12345, -67890),
    "size": 500.0,
    "slide": "857",
    "fov_id": "fov1",
    "fov_name": "R1C3_Tonsil",
    "folder": "fov1/RowNumber0/Depth_Profile0",
    "dwell": 4,
    "scans": "0,5",
    "aperture": "B",
    "instrument": "MIBIscope1",
    "tissue": "Tonsil",
    "panel": "20170916_1x",
    "mass_offset": 0.1,
    "mass_gain": 0.2,
    "time_resolution": 0.5,
    "miscalibrated": False,
    "check_reg": False,
    "filename": "20180703_1234_test",
    "description": "test image",
    "version": "alpha",
}


def _gen_tif_data(fov_number, chan_number, img_shape, fills, dtype):
    """Generates random or set-filled image data

    Args:
        fov_number (int):
            Number of fovs required
        chan_number (int):
            Number of channels required
        img_shape (tuple):
            Single image dimensions (x pixels, y pixels)
        fills (bool):
            If False, data is randomized.  If True, each single image will be filled with a value
            one less than that of the next channel.  If said image is the last channel, then the
            value is one less than that of the first channel in the next fov.
        dtype (type):
            Data type for generated data

    Returns:
        numpy.ndarray:
            Image data with shape (fov_number, img_shape[0], img_shape[1], chan_number)

    """
    if not fills:
        tif_data = np.random.randint(0, 100, size=(fov_number, *img_shape, chan_number)).astype(
            dtype
        )
    else:
        tif_data = np.full(
            (*img_shape, fov_number, chan_number),
            (np.arange(fov_number * chan_number) % 256).reshape(fov_number, chan_number),
            dtype=dtype,
        )
        tif_data = np.moveaxis(tif_data, 2, 0)

    return tif_data


def _gen_label_data(fov_number, comp_number, img_shape, dtype):
    """Generates quadrant-based label data

    Args:
        fov_number (int):
            Number of fovs required
        comp_number (int):
            Number of components
        img_shape (tuple):
            Single image dimensions (x pixels, y pixesl)
        dtype (type):
            Data type for generated labels

    Returns:
        numpy.ndarray:
            Label data with shape (fov_number, img_shape[0], img_shape[1], comp_number)
    """
    label_data = np.zeros((fov_number, *img_shape, comp_number), dtype=dtype)

    right = (img_shape[1] - 1) // 2
    left = (img_shape[1] + 2) // 2
    up = (img_shape[0] - 1) // 2
    down = (img_shape[0] + 2) // 2

    counter = 1
    for fov in range(fov_number):
        for comp in range(comp_number):
            label_data[fov, :up, :right, comp] = counter
            counter = (counter % 255) + 1
            label_data[fov, :up, left:, comp] = counter
            counter = (counter % 255) + 1
            label_data[fov, down:, :right, comp] = counter
            counter = (counter % 255) + 1
            label_data[fov, down:, left:, comp] = counter
            counter = (counter % 255) + 1

    return label_data


def _write_tifs(base_dir, fov_names, img_names, shape, sub_dir, fills, dtype, single_dir=False):
    """Generates and writes single tifs to into base_dir/fov_name/sub_dir

    Args:
        base_dir (str):
            Path to base directory
        fov_names (list):
            List of fov folders to create/fill
        img_names (list):
            Channel names
        shape (tuple):
            Single image shape (x pixels, y pixels)
        sub_dir (str):
            Subdirectory to write images into
        fills (bool):
            If False, data is randomized.  If True, each single image will be filled with a value
            one less than that of the next channel.  If said image is the last channel, then the
            value is one less than that of the first channel in the next fov.
        dtype (type):
            Data type for generated images
        single_dir (bool):
            whether to create single image dir with fov names prepended to image file

    Returns:
        tuple (dict, numpy.ndarray):

        - File locations, indexable by fov names
        - Image data as an array with shape (num_fovs, shape[0], shape[1], num_channels)
    """
    tif_data = _gen_tif_data(len(fov_names), len(img_names), shape, fills, dtype)

    if sub_dir is None:
        sub_dir = ""

    filelocs = {}

    for i, fov in enumerate(fov_names):
        filelocs[fov] = []
        # write to individual fov folders
        if not single_dir:
            fov_path = os.path.join(base_dir, fov, sub_dir)
            os.makedirs(fov_path)
        for j, name in enumerate(img_names):
            # prepend fov name to single directory images
            if single_dir:
                img_path = os.path.join(base_dir, f"{fov}_{name}")
            else:
                img_path = os.path.join(fov_path, name)
            image_utils.save_image(img_path + ".tiff", tif_data[i, :, :, j])
            filelocs[fov].append(img_path)

    return filelocs, tif_data


def _write_multitiff(
    base_dir, fov_names, channel_names, shape, sub_dir, fills, dtype, channels_first=False
):
    """Generates and writes multitifs to into base_dir

    Args:
        base_dir (str):
            Path to base directory
        fov_names (list):
            List of fov files to write
        channel_names (list):
            Channel names
        shape (tuple):
            Single image shape (x pixels, y pixels)
        sub_dir (str):
            Ignored.
        fills (bool):
            If False, data is randomized.  If True, each single image will be filled with a value
            one less than that of the next channel.  If said image is the last channel, then the
            value is one less than that of the first channel in the next fov.
        dtype (type):
            Data type for generated images
        channels_first(bool):
            Indicates whether the data should be saved in channels_first format. Default: False

    Returns:
        tuple (dict, numpy.ndarray):

        - File locations, indexable by fov names
        - Image data as an array with shape (num_fovs, shape[0], shape[1], num_channels)
    """
    tif_data = _gen_tif_data(len(fov_names), len(channel_names), shape, fills, dtype)

    filelocs = {}

    for i, fov in enumerate(fov_names):
        tiffpath = os.path.join(base_dir, f"{fov}.tiff")
        v = tif_data[i, :, :, :]
        if channels_first:
            v = np.moveaxis(v, -1, 0)
        image_utils.save_image(tiffpath, v)
        filelocs[fov] = tiffpath

    return filelocs, tif_data


def _write_mibitiff(base_dir, fov_names, channel_names, shape, sub_dir, fills, dtype):
    """Generates and writes mibitiffs to into base_dir

    Args:
        base_dir (str):
            Path to base directory
        fov_names (list):
            List of fov files to write
        channel_names (list):
            Channel names
        shape (tuple):
            Single image shape (x pixels, y pixels)
        sub_dir (str):
            Ignored.
        fills (bool):
            If False, data is randomized.  If True, each single image will be filled with a value
            one less than that of the next channel.  If said image is the last channel, then the
            value is one less than that of the first channel in the next fov.
        dtype (type):
            Data type for generated images

    Returns:
        tuple (dict, numpy.ndarray):

        - File locations, indexable by fov names
        - Image data as an array with shape (num_fovs, shape[0], shape[1], num_channels)
    """
    tif_data = _gen_tif_data(len(fov_names), len(channel_names), shape, fills, dtype)

    filelocs = {}

    mass_map = tuple(enumerate(channel_names, 1))

    for i, fov in enumerate(fov_names):
        tiffpath = os.path.join(base_dir, f"{fov}.tiff")
        tiff_utils.write_mibitiff(tiffpath, tif_data[i, :, :, :], mass_map, MIBITIFF_METADATA)
        filelocs[fov] = tiffpath

    return filelocs, tif_data


def _write_reverse_multitiff(base_dir, fov_names, channel_names, shape, sub_dir, fills, dtype):
    """Generates and writes 'reversed' multitifs to into base_dir

    Saved images have shape (num_channels, shape[0], shape[1]).  This is mostly useful for
    testing deepcell-input loading.

    Args:
        base_dir (str):
            Path to base directory
        fov_names (list):
            List of fov files to write
        channel_names (list):
            Channel names
        shape (tuple):
            Single image shape (x pixels, y pixels)
        sub_dir (str):
            Ignored.
        fills (bool):
            If False, data is randomized.  If True, each single image will be filled with a value
            one less than that of the next channel.  If said image is the last channel, then the
            value is one less than that of the first channel in the next fov.
        dtype (type):
            Data type for generated images

    Returns:
        tuple (dict, numpy.ndarray):

        - File locations, indexable by fov names
        - Image data as an array with shape (num_fovs, shape[0], shape[1], num_channels)
    """
    tif_data = _gen_tif_data(len(channel_names), len(fov_names), shape, fills, dtype)

    filelocs = {}

    for i, fov in enumerate(fov_names):
        tiffpath = os.path.join(base_dir, f"{fov}.tiff")
        image_utils.save_image(tiffpath, tif_data[:, :, :, i])
        filelocs[fov] = tiffpath

    tif_data = np.swapaxes(tif_data, 0, -1)

    return filelocs, tif_data


def _write_labels(base_dir, fov_names, comp_names, shape, sub_dir, fills, dtype, suffix=""):
    """Generates and writes label maps to into base_dir

    Args:
        base_dir (str):
            Path to base directory
        fov_names (list):
            List of fov files to write
        comp_names (list):
            Component names
        shape (tuple):
            Single image shape (x pixels, y pixels)
        sub_dir (str):
            Ignored.
        fills (bool):
            Ignored.
        dtype (type):
            Data type for generated labels
        suffix (str):
            Suffix for label datafiles

    Returns:
        tuple (dict, numpy.ndarray):

        - File locations, indexable by fov names
        - Label data as an array with shape (num_fovs, shape[0], shape[1], num_components)
    """
    label_data = _gen_label_data(len(fov_names), len(comp_names), shape, dtype)

    filelocs = {}

    for i, fov in enumerate(fov_names):
        tiffpath = os.path.join(base_dir, f"{fov}{suffix}.tiff")
        image_utils.save_image(tiffpath, label_data[i, :, :, 0])
        filelocs[fov] = tiffpath

    return filelocs, label_data


TIFFMAKERS = {
    "tiff": _write_tifs,
    "multitiff": _write_multitiff,
    "reverse_multitiff": _write_reverse_multitiff,
    "mibitiff": _write_mibitiff,
    "labels": _write_labels,
}


def create_paired_xarray_fovs(
    base_dir,
    fov_names,
    channel_names,
    img_shape=(10, 10),
    mode="tiff",
    delimiter=None,
    sub_dir=None,
    fills=False,
    dtype="int8",
    channels_first=False,
    single_dir=False,
):
    """Writes data to file system (images or labels) and creates expected xarray for reloading
    data from said file system.

    Args:
        base_dir (str):
            Path to base directory.  All data will be written into this folder.
        fov_names (list):
            List of fovs
        channel_names (list):
            List of channels/components
        img_shape (tuple):
            Single image shape (x pixels, y pixels)
        mode (str):
            The type of data to generate. Current options are:

            - 'tiff'
            - 'multitiff'
            - 'reverse_multitiff'
            - 'mibitiff'
            - 'labels'
        delimiter (str or None):
            Delimiting character or string separating fov_id from rest of file/folder name.
            Default is None.
        sub_dir (str):
            Only active for 'tiff' mode.  Creates another sub directory in which tiffs are stored
            within the parent fov folder.  Default is None.
        fills (bool):
            Only active for image data (not 'labels'). If False, data is randomized.  If True,
            each single image will be filled with a value one less than that of the next channel.
            If said image is the last channel, then the value is one less than that of the first
            channel in the next fov.
        dtype (type):
            Data type for generated images/labels.  Default is int16
        channels_first (bool):
            Indicates whether the data should be saved in channels_first format when
            mode is 'multitiff'. Default: False
        single_dir (bool):
            whether to create single image dir with fov names prepended to image file

    Returns:
        tuple (dict, xarray.DataArray):

        - File locations, indexable by fov names
        - Image/label data as an xarray with shape
          (num_fovs, im_shape[0], shape[1], num_channels)
    """

    # validation checks
    io_utils.validate_paths(base_dir)

    if fov_names is None or fov_names is []:
        raise ValueError("No fov names were given...")

    if channel_names is None or channel_names is []:
        raise ValueError("No image names were given...")

    if not isinstance(fov_names, list):
        fov_names = [fov_names]

    if not isinstance(channel_names, list):
        channel_names = [channel_names]

    if mode == "multitiff":
        filelocs, tif_data = TIFFMAKERS[mode](
            base_dir,
            fov_names,
            channel_names,
            img_shape,
            sub_dir,
            fills,
            dtype,
            channels_first=channels_first,
        )
    elif mode == "tiff":
        filelocs, tif_data = TIFFMAKERS[mode](
            base_dir, fov_names, channel_names, img_shape, sub_dir, fills, dtype, single_dir
        )
    else:
        filelocs, tif_data = TIFFMAKERS[mode](
            base_dir, fov_names, channel_names, img_shape, sub_dir, fills, dtype
        )

    if delimiter is not None:
        fov_ids = [fov.split(delimiter)[0] for fov in fov_names]
    else:
        fov_ids = fov_names

    if "multitiff" in mode:
        channel_names = range(len(channel_names))

    if mode == "labels":
        data_xr = make_labels_xarray(tif_data, fov_ids, channel_names, *img_shape)
    else:
        data_xr = make_images_xarray(tif_data, fov_ids, channel_names, *img_shape)

    return filelocs, data_xr


def make_images_xarray(
    tif_data, fov_ids=None, channel_names=None, row_size=10, col_size=10, dtype="int16"
):
    """Generate a correctly formatted image data xarray

    Args:
        tif_data (numpy.ndarray or None):
            Image data to embed within the xarray.  If None, randomly generated image data is used,
            but fov_ids and channel_names must not be None.
        fov_ids (list or None):
            List of fov names.  If None, fov id's will be generated based on the shape of tif_data
            following the scheme 'fov0', 'fov1', ... , 'fovN'. Default is None.
        channel_names (list or None):
            List of channel names.  If None, channel names will be generated based on the shape of
            tif_data following the scheme 'chan0', 'chan1', ... , 'chanM'.  Default is None.
        row_size (int):
            Horizontal size of individual image.  Default is 10.
        col_size (int):
            Vertical size of individual image. Default is 10.
        dtype (type):
            Data type for generated images.  Default is int16.

    Returns:
        xarray.DataArray:
            Image data with standard formatting
    """
    if tif_data is None:
        tif_data = _gen_tif_data(
            len(fov_ids), len(channel_names), (row_size, col_size), False, dtype=dtype
        )
    else:
        row_size, col_size = tif_data.shape[1:3]

        buf_fov_ids, buf_chan_names = gen_fov_chan_names(tif_data.shape[0], tif_data.shape[-1])
        if fov_ids is None:
            fov_ids = buf_fov_ids
        if channel_names is None:
            channel_names = buf_chan_names

    coords = [fov_ids, range(row_size), range(col_size), channel_names]
    dims = ["fovs", "rows", "cols", "channels"]
    return xr.DataArray(tif_data, coords=coords, dims=dims)


def make_labels_xarray(
    label_data, fov_ids=None, compartment_names=None, row_size=10, col_size=10, dtype="int16"
):
    """Generate a correctly formatted label data xarray

    Args:
        label_data (numpy.ndarray or None):
            Label data to embed within the xarray.  If None, automatically generated label data is
            used, but fov_ids and compartment_names must not be None.
        fov_ids (list or None):
            List of fov names.  If None, fov id's will be generated based on the shape of tif_data
            following the scheme 'fov0', 'fov1', ... , 'fovN'. Default is None.
        compartment_names (list or None):
            List of compartment names.  If None, compartment names will be ['whole_cell'] or
            ['whole_cell', 'nuclear'] if label_data.shape[-1] is 1 or 2 respecticely. Default is
            None.
        row_size (int):
            Horizontal size of individual image.  Default is 10.
        col_size (int):
            Vertical size of individual image. Default is 10.
        dtype (type):
            Data type for generated labels.  Default is int16.

    Returns:
        xarray.DataArray:
            Label data with standard formatting
    """
    if label_data is None:
        label_data = _gen_label_data(
            len(fov_ids), len(compartment_names), (row_size, col_size), dtype=dtype
        )
    else:
        row_size, col_size = label_data.shape[1:3]

        buf_fov_ids, _ = gen_fov_chan_names(label_data.shape[0], 0)
        if fov_ids is None:
            fov_ids = buf_fov_ids
        if compartment_names is None:
            comp_dict = {1: ["whole_cell"], 2: ["whole_cell", "nuclear"]}
            compartment_names = comp_dict[label_data.shape[-1]]

    coords = [fov_ids, range(row_size), range(col_size), compartment_names]
    dims = ["fovs", "rows", "cols", "compartments"]
    return xr.DataArray(label_data, coords=coords, dims=dims)
