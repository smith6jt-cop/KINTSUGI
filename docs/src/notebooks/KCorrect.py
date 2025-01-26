import time
import numpy as np
from typing import List
from skimage.transform import resize as skresize
from scipy.fftpack import dct, idct

RESIZE_ORDER = 1
RESIZE_MODE = "symmetric"
PRESERVE_RANGE = True
OUTPUT_IMAGE = "OutputImage"

def _resize_images_list(images_list: List, side_size: float = None, x_side_size: float = None, y_side_size: float = None):
    if side_size is not None:
        y_side_size = x_side_size = side_size
    resized_images_list = []
    for im in images_list:
        if im.shape[0] != x_side_size or im.shape[1] != y_side_size:
            resized_images_list.append(skresize(
                im, 
                (x_side_size, y_side_size), 
                order = RESIZE_ORDER, 
                mode = RESIZE_MODE,
                preserve_range = PRESERVE_RANGE
                )
            )
        else:
            resized_images_list.append(im)
    return resized_images_list

def _dct2d(mtrx: np.array):
    """
    Calculates 2D discrete cosine transform.
    
    Parameters
    ----------
    mtrx
        Input matrix.  
        
    Returns
    -------    
    Discrete cosine transform of the input matrix.
    """
     
    # Check if input object is 2D.
    if mtrx.ndim != 2:
        raise ValueError("Passed object should be a matrix or a numpy array with dimension of two.")

    return dct(dct(mtrx.T, norm='ortho').T, norm='ortho')

def _idct2d(mtrx: np.array):
    """
    Calculates 2D inverse discrete cosine transform.
    
    Parameters
    ----------
    mtrx
        Input matrix.  
        
    Returns
    -------    
    Inverse of discrete cosine transform of the input matrix.
    """
     
    # Check if input object is 2D.
    if mtrx.ndim != 2:
        raise ValueError("Passed object should be a matrix or a numpy array with dimension of two.")
 
    return idct(idct(mtrx.T, norm='ortho').T, norm='ortho')

def _shrinkageOperator(matrix, epsilon):
    temp1 = matrix - epsilon
    temp1[temp1 < 0] = 0
    temp2 = matrix + epsilon
    temp2[temp2 > 0] = 0
    res = temp1 + temp2
    return res

def _inexact_alm_rspca_l1(
    images, 
    lambda_flatfield, 
    if_darkfield, 
    lambda_darkfield, 
    optimization_tolerance, 
    max_iterations,
    weight=None, 
    ):

    if weight is not None and weight.size != images.size:
            raise IOError('weight matrix has different size than input sequence')

    # Shape variables
    p = images.shape[0]           # Image height
    q = images.shape[1]           # Image width
    m = p*q                       # Total pixels per image
    n = images.shape[2]          # Number of images
    images = np.reshape(images, (m, n), order='F')

    if weight is not None:
        weight = np.reshape(weight, (m, n), order='F')
    else:
        weight = np.ones_like(images)

    # SVD and Norm variables    
    _, svd, _ = np.linalg.svd(images, full_matrices=False) # SVD decomposition
    norm_two = svd[0] # Largest singular value
    d_norm = np.linalg.norm(images, ord='fro') # Frobenius norm of images

    # Optimization parameters
    dual_var_lowrank = 0 # Dual variable for low-rank component (Y1)
    lagrange_mult1 = 1 # Lagrange multiplier for first constraint (ent1)
    lagrange_mult2 = 10 # Lagrange multiplier for second constraint (ent2)
    penalty_factor = 12.5 / norm_two # (mu)
    penalty_factor_bar = penalty_factor * 1e7 # Upper bound for penalty (mu_bar)
    scale_ratio = 1.5 # Scale factor for penalty updates

    # Component matrices
    A1_hat = np.zeros_like(images) # Estimated flat-field
    A1_coeff = np.ones((1, images.shape[1])) # Flat-field coefficients
    E1_hat = np.zeros_like(images) # Estimated error/noise
    W_hat = _dct2d(np.zeros((p, q)).T) # DCT coefficients

    # Offset and mask variables
    A_offset = np.zeros((m, 1)) # Offset for flatfield
    B1_uplimit = np.min(images) # Upper limit for darkfield
    B1_offset = 0 # Darkfield offset
    A_inmask = np.zeros((p, q)) # Inner mask for optimization
    # Mask covers central region (5/6 x 5/6)
    A_inmask[int(np.round(p / 6) - 1): int(np.round(p*5 / 6)), int(np.round(q / 6) - 1): int(np.round(q * 5 / 6))] = 1

    # main iteration loop starts
    iter = 0
    total_svd = 0
    converged = False

    while not converged:
 
        iter += 1

        if len(A1_coeff.shape) == 1:
            A1_coeff = np.expand_dims(A1_coeff, 0)
        if len(A_offset.shape) == 1:
            A_offset = np.expand_dims(A_offset, 1)

        W_idct_hat = _idct2d(W_hat.T)
        A1_hat = np.dot(np.reshape(W_idct_hat, (-1,1), order='F'), A1_coeff) + A_offset

        temp_W = (images - A1_hat - E1_hat + (1 / penalty_factor) * dual_var_lowrank) / lagrange_mult1
        temp_W = np.reshape(temp_W, (p, q, n), order='F')
        temp_W = np.mean(temp_W, axis=2)

        W_hat = W_hat + _dct2d(temp_W.T)
        W_hat = np.maximum(W_hat - lambda_flatfield / (lagrange_mult1 * penalty_factor), 0) + np.minimum(W_hat + lambda_flatfield / (lagrange_mult1 * penalty_factor), 0)
        W_idct_hat = _idct2d(W_hat.T)

        if len(A1_coeff.shape) == 1:
            A1_coeff = np.expand_dims(A1_coeff, 0)
        if len(A_offset.shape) == 1:
            A_offset = np.expand_dims(A_offset, 1)

        A1_hat = np.dot(np.reshape(W_idct_hat, (-1,1), order='F'), A1_coeff) + A_offset
        E1_hat = images - A1_hat + (1 / penalty_factor) * dual_var_lowrank / lagrange_mult1
        E1_hat = _shrinkageOperator(E1_hat, weight / (lagrange_mult1 * penalty_factor))
        R1 = images - E1_hat
        A1_coeff = np.mean(R1, 0) / np.mean(R1)
        A1_coeff[A1_coeff < 0] = 0

        if if_darkfield:
            validA1coeff_idx = np.where(A1_coeff < 1)

            B1_coeff = (np.mean(R1[np.reshape(W_idct_hat, -1, order='F') > np.mean(W_idct_hat) - 1e-6][:, validA1coeff_idx[0]], 0) - \
            np.mean(R1[np.reshape(W_idct_hat, -1, order='F') < np.mean(W_idct_hat) + 1e-6][:, validA1coeff_idx[0]], 0)) / np.mean(R1)
            k = np.array(validA1coeff_idx).shape[1]
            temp1 = np.sum(A1_coeff[validA1coeff_idx[0]]**2)
            temp2 = np.sum(A1_coeff[validA1coeff_idx[0]])
            temp3 = np.sum(B1_coeff)
            temp4 = np.sum(A1_coeff[validA1coeff_idx[0]] * B1_coeff)
            temp5 = temp2 * temp3 - temp4 * k
            if temp5 == 0:
                B1_offset = 0
            else:
                B1_offset = (temp1 * temp3 - temp2 * temp4) / temp5

            B1_offset = np.maximum(B1_offset, 0)
            B1_offset = np.minimum(B1_offset, B1_uplimit / np.mean(W_idct_hat))
            B_offset = B1_offset * np.reshape(W_idct_hat, -1, order='F') * (-1)
            B_offset = B_offset + np.ones_like(B_offset) * B1_offset * np.mean(W_idct_hat)

            A1_offset = np.mean(R1[:, validA1coeff_idx[0]], axis=1) - np.mean(A1_coeff[validA1coeff_idx[0]]) * np.reshape(W_idct_hat, -1, order='F')
            A1_offset = A1_offset - np.mean(A1_offset)
            A_offset = A1_offset - np.mean(A1_offset) - B_offset

            # smooth A_offset
            W_offset = _dct2d(np.reshape(A_offset, (p,q), order='F').T)
            W_offset = np.maximum(W_offset - lambda_darkfield / (lagrange_mult2 * penalty_factor), 0) + \
                np.minimum(W_offset + lambda_darkfield / (lagrange_mult2 * penalty_factor), 0)
            A_offset = _idct2d(W_offset.T)
            A_offset = np.reshape(A_offset, -1, order='F')

            # encourage sparse A_offset
            A_offset = np.maximum(A_offset - lambda_darkfield / (lagrange_mult2 * penalty_factor), 0) + \
                np.minimum(A_offset + lambda_darkfield / (lagrange_mult2 * penalty_factor), 0)
            A_offset = A_offset + B_offset


        Z1 = images - A1_hat - E1_hat # Constraint violation
        dual_var_lowrank = dual_var_lowrank + penalty_factor * Z1 # Dual variable updates
        penalty_factor = np.minimum(penalty_factor * scale_ratio, penalty_factor_bar)
        
        # Stop Criterion
        stopCriterion = np.linalg.norm(Z1, ord='fro') / d_norm
        # print(f'Iteration {iter}, stopCriterion: {stopCriterion}')
        if stopCriterion < optimization_tolerance:
            converged = True

        if not converged and iter >= max_iterations:
            converged = True
    A_offset = np.squeeze(A_offset)
    A_offset = A_offset + B1_offset * np.reshape(W_idct_hat, -1, order='F')
    return A1_hat, E1_hat, A_offset

def _resize_image(image: np.ndarray, side_size: float  = None, x_side_size: float = None, y_side_size: float = None):
    if side_size is not None:
        y_side_size = x_side_size = side_size
    if image.shape[0] != x_side_size or image.shape[1] != y_side_size:
        return skresize(
            image,
            (x_side_size, y_side_size), 
            order = RESIZE_ORDER, 
            mode = RESIZE_MODE,
            preserve_range = PRESERVE_RANGE
        )
    else:
        return image
    
def validate_correction(original, corrected):
    # Check relative intensity preservation
    orig_ratio = np.max(original) / np.mean(original)
    corr_ratio = np.max(corrected) / np.mean(corrected)
    if abs(orig_ratio - corr_ratio) > 0.1 * orig_ratio:  # 10% threshold
        print("Significant change in relative intensities detected") 

def KCorrect(
        images_list, 
        if_darkfield,
        max_iterations,
        optimization_tolerance,
        max_reweight_iterations,
        reweight_tolerance 
    ):
    _saved_size = images_list[0].shape
    working_size = 128
    nrows = ncols = working_size
    downsized_image = np.dstack(_resize_images_list(images_list, side_size=working_size))
    mean_image = np.mean(downsized_image, axis=2)
    mean_div = mean_image/np.mean(mean_image)
    Dis_Cos_Trans_mean = _dct2d(mean_div.T)
    sorted_images = np.sort(downsized_image, axis=2)

    lambda_flatfield = np.sum(np.abs(Dis_Cos_Trans_mean)) / 400 * 0.5
    lambda_darkfield = lambda_flatfield * 0.2

    XAoffset = np.zeros((nrows, ncols))
    weight = np.ones(sorted_images.shape)
    eplson = 0.1

    reweighting_iter = 0
    flag_reweighting = True
    flatfield_last = np.ones((nrows, ncols))
    darkfield_last = np.random.randn(nrows, ncols)

    while flag_reweighting:
        reweighting_iter += 1

        initial_flatfield = False
        if initial_flatfield:
            raise IOError('Initial flatfield option not implemented yet!')
        else:
            X_k_A, X_k_E, X_k_Aoffset = _inexact_alm_rspca_l1(
                images = sorted_images, 
                lambda_flatfield = lambda_flatfield,
                if_darkfield = if_darkfield, 
                lambda_darkfield = lambda_darkfield, 
                optimization_tolerance = optimization_tolerance, 
                max_iterations = max_iterations,
                weight=weight
            )

        XA = np.reshape(X_k_A, [nrows, ncols, -1], order='F')
        XE = np.reshape(X_k_E, [nrows, ncols, -1], order='F')
        XAoffset = np.reshape(X_k_Aoffset, [nrows, ncols], order='F')
        XE_norm = XE / np.mean(XA, axis=(0, 1))

        weight = np.ones_like(XE_norm) / (np.abs(XE_norm) + eplson)

        weight = weight * weight.size / np.sum(weight)

        temp = np.mean(XA, axis=2) - XAoffset
        flatfield_current = temp / np.mean(temp)
        darkfield_current = XAoffset
        mad_flatfield = np.sum(np.abs(flatfield_current - flatfield_last)) / np.sum(np.abs(flatfield_last))
        temp_diff = np.sum(np.abs(darkfield_current - darkfield_last))
        
        if temp_diff < 1e-7:
            mad_darkfield = 0
        else:
            mad_darkfield = temp_diff / np.maximum(np.sum(np.abs(darkfield_last)), 1e-6)

        # print(f"Re-weighting iteration {reweighting_iter}: MAD flatfield = {mad_flatfield}; MAD darkfield = {mad_darkfield}")
        flatfield_last = flatfield_current
        darkfield_last = darkfield_current
        if np.maximum(mad_flatfield,
                        mad_darkfield) <= reweight_tolerance or \
                reweighting_iter >= max_reweight_iterations:
            flag_reweighting = False

    shading = np.mean(XA, 2) - XAoffset
    flatfield = _resize_image(
        image = shading, 
        x_side_size = _saved_size[0], 
        y_side_size = _saved_size[1]
    )

    flatfield = flatfield / np.mean(flatfield)

    if if_darkfield:
        darkfield = _resize_image(
            image = XAoffset, 
            x_side_size = _saved_size[0], 
            y_side_size = _saved_size[1]
        )
    else:
        darkfield = np.zeros_like(flatfield)
    return flatfield, darkfield

# def save_model(model_dir: PathLike, overwrite: bool = False) -> None:
#     """Save current model to folder.

#     Args:
#         model_dir: path to model directory

#     Raises:
#         FileExistsError: if model directory already exists
#     """
#     path = Path(model_dir)

#     try:
#         path.mkdir()
#     except FileExistsError:
#         if not overwrite:
#             raise FileExistsError("Model folder already exists.")

#     # save settings
#     with open(path / "settings.json", "w") as fp:
#         # see pydantic docs for output options
#         fp.write(self.json())

#     # NOTE emit warning if profiles are all zeros? fit probably not run
#     # save profiles
#     np.savez(
#         path / _PROFILES_FNAME,
#         flatfield=np.array(self.flatfield),
#         darkfield=np.array(self.darkfield),
#         baseline=np.array(self.baseline),
#     )

# def load_model(cls, model_dir: PathLike) -> BaSiC:
#     """Create a new instance from a model folder."""
#     path = Path(model_dir)

#     if not path.exists():
#         raise FileNotFoundError("Model directory not found.")

#     with open(path / "settings.json") as fp:
#         model = json.load(fp)

#     profiles = np.load(path / "profiles.npz")
#     model["flatfield"] = profiles["flatfield"]
#     model["darkfield"] = profiles["darkfield"]
#     model["baseline"] = profiles["baseline"]

#     return BaSiC(**model)