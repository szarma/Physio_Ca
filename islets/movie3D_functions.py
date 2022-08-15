from typing import Union, Tuple, List
import numpy as np
import cv2
from .numeric import multi_map
from skimage.transform import warp, AffineTransform
from skimage.feature import match_template
import logging
from matplotlib.pyplot import plot

def extract_frameshift_3d(
    frame: np.typing.NDArray[np.float_],
    template: np.typing.NDArray[np.float_],
    max_shifts: Union[ np.typing.ArrayLike, List[int]],
    method: str = 'skimage',
    log: bool = False,
) -> Tuple[np.typing.NDArray[np.float_], np.float_]:
    """
    Performs 3D motion correction for a single frame at a subpixel level by comparing it to a reference template.

    Args:
        max_shifts: maximum pixel shifts allowed when correcting in the z-, height, and width direction

        template: reference template (Z,H,W)

        method: 'skimage' (default) or 'opencv' (not implemented yet).

        log: whether to log the arrays when determining the correlation

    Returns:
        shifts : tuple, contains shifts in z,h, and w dimension and correlation with template

        xcorrs: cross correlation of the movies with the template (estimate)

    Raises:
        Exception 'Unknown motion correction method!'

    """
    assert template.shape == frame.shape
    z_i, h_i, w_i = frame.shape
    
    template = template[
               max_shifts[0]: z_i - max_shifts[0],
               max_shifts[1]: h_i - max_shifts[1],
               max_shifts[2]: w_i - max_shifts[2]].copy()
    if log:
        mint = np.percentile(template[template>0],1)
        template = np.maximum(template, mint)
        template = np.log(template)

        minf = np.percentile(frame[frame>0],1)
        frame = np.maximum(frame, minf)
        frame = np.log(frame)


    if method == 'opencv':
        raise NotImplementedError("not yet implemented in opencv")
        # res = cv2.matchTemplate(movie, template, cv2.TM_CCORR_NORMED)
        # idx_max = cv2.minMaxLoc(res)[3]
    elif method == 'skimage':
        res = match_template(frame, template)
        if res.max()<=0:
            logging.warning("No reasonable shift found. Defaulting to 0.")
            return np.ones(3, dtype="float"), 1
        idx_max = np.unravel_index(np.argmax(res), res.shape)
        # idx_max = idx_max[::-1]
    else:
        raise Exception('Unknown motion correction method!')


    shifts = np.array([ np.clip(idx_max[j], 1, 2 * max_shifts[j] - 1) if max_shifts[j]>0 else 0 for j in range(3) ])
    shift_subpx = shifts.astype("float")
    log_cor = np.log(1e-10 + np.maximum(res,0))
    # now we fit quadratic for the logged correlation coeficients to check for subpixel shifts
    # y = a*x^2 + bx + c in each of the dimensions
    for dim, sh in enumerate(shifts):
        if max_shifts[dim]==0:
            continue
        slice_ = list(shifts)
        # n = shifts[dim]
        slice_[dim] = slice(sh - 1, sh + 2)
        y = log_cor[tuple(slice_)]
        a = (y[2] + y[0] - 2 * y[1]) / 2
        b = (y[2] - y[0]) / 2 - 2 * a * shifts[dim]
        # position of the minimum in dimension dim is redefined,
        # but only if the new value is close to the original
        if np.abs(a)<1e-10:
            shift_subpx[dim] = 0
        else:
            if np.abs(-b / 2 / a - shifts[dim]) < 1:
                shift_subpx[dim] = -b / 2 / a
    corr_est = res.max()
    for j in range(3):
        shift_subpx[j] = np.clip( - max_shifts[j] + shift_subpx[j], -max_shifts[j]-.5, max_shifts[j]+.5)
    return shift_subpx, corr_est

def extract_shifts_3d(
    movie: np.typing.NDArray[np.float_],
    max_shifts: Union[np.typing.ArrayLike, List[int], None] = None,
    template: Union[np.typing.NDArray[np.float_], None] = None,
    n_processes: int = 1,
    method: str = 'skimage',
    mode: str = "consecutive+template"
    ) -> Tuple[np.typing.ArrayLike, np.typing.ArrayLike]:
    """
    Performs motion correction using the opencv matchtemplate function. At every iteration a template is built by taking the median of all frames and then used to align the other frames.

    Args:
        max_shifts[2],max_shifts[1]: maximum pixel shifts allowed when correcting in the width and height direction

        template: if a good template for frame by frame correlation is available it can be passed. If None it is automatically computed

        method: depends on what is installed 'opencv' or 'skimage'. 'skimage' is an order of magnitude slower

    Returns:
        shifts : tuple, contains shifts in x and y and correlation with template

        xcorrs: cross correlation of the movies with the template

    Raises:
        Exception 'Unknown motion correction method!'

    """

    global iterf_
    valid_modes = ["consecutive","template","consecutive+template"]
    if mode not in valid_modes:
        raise ValueError(f"Mode can only assume one of the following values: {valid_modes}")
    if max_shifts is None:
        max_shifts = np.array([2,5,5])
    else:
        max_shifts = np.asanyarray(max_shifts)
    infer_template = template is None
    # _, z_i, h_i, w_i = movie.shape
    shifts = np.zeros(shape = (len(movie),3), dtype = float)
    if "consecutive" in mode:
        def iterf_(pair_of_frames):
            return extract_frameshift_3d(pair_of_frames[0], pair_of_frames[1], max_shifts, method = method)[0]
        # for i in range(1):
        #     print (f"step {i} of a consecutive mode")
        dshifts = np.array([[0] * 3] + list(multi_map(iterf_, zip(movie, movie[1:]), processes = n_processes)))
        shifts_ = np.cumsum(dshifts, axis = 0)
        shifts_ = shifts_ - np.median(shifts_, axis = 0)
        movie = apply_shifts_3d(movie, shifts_, )
        shifts += shifts_
        # c = plot([])[0].get_color()
        # plot(shifts,c=c)
        # if all(np.abs(dshifts).max()<max_shifts):
        #     break
    else:
        movie = movie.copy()
    corrs = None

    if "template" in mode:
        def iterf_(frame):
            return extract_frameshift_3d(frame, template, max_shifts, method = method)[0]
        for i in range(10): # only because I am afraid of while True
            # print (f"step {i} of a template mode")
            if infer_template:
                median_movie = np.median(movie, axis = 0)
                L1_dist_to_median = np.sum(np.abs(movie - median_movie), axis = (1, 2, 3))
                ichoose = np.argmin(L1_dist_to_median)
                logging.warning(f"Template not provided. Choosing frame {ichoose} as a template.")
                template = movie[ichoose].copy()
            # shifts_and_corrs = multi_map(iterf_, movie, processes = n_processes)
            # corrs = np.array([el[1] for el in shifts_and_corrs])
            # shifts_and_corrs = list(shifts_and_corrs)
            # dshifts = np.array([el[0] for el in shifts_and_corrs])
            dshifts = np.array(list(multi_map(iterf_, movie, processes = n_processes)))
            shifts += dshifts
            # c = plot([])[0].get_color()
            # plot(shifts,c=c)
            maxDshifts = np.abs(dshifts).max(0)
            maxDshifts = maxDshifts[max_shifts>0]
            if all(maxDshifts<max_shifts[max_shifts>0]):
                break
            else:
                movie = apply_shifts_3d(movie, dshifts)
        corrs = np.array([np.corrcoef(frame.flat,template.flat)[0,1] for frame in movie])

    return (shifts, corrs)


def apply_shifts_3d(
        movie: np.typing.NDArray,
        shifts: np.typing.NDArray,
        interpolation: str = 'linear',
        method: str = 'opencv'
) -> np.typing.NDArray[np.float_]:
    """
    Apply precomputed shifts to a movie, using subpixels adjustment (cv2.INTER_CUBIC function)

    Args:
        shifts: array of tuples representing x and y shifts for each frame

        interpolation: 'linear', 'cubic', 'nearest' or cvs.INTER_XXX

        method: (undocumented)

        remove_blanks: (undocumented)

    Returns:
        self

    Raise:
        Exception 'Interpolation method not available'

        Exception 'Method not defined'
    """

    # if type(self.flat[0]) is not np.float32:
    #     logging.warning('Casting the array to float32')
    #     self = np.asanyarray(self, dtype = np.float32)

    if interpolation == 'cubic':
        if method == 'opencv':
            interpolation = cv2.INTER_CUBIC
        else:
            interpolation = 3
        logging.debug('cubic interpolation')

    elif interpolation == 'nearest':
        if method == 'opencv':
            interpolation = cv2.INTER_NEAREST
        else:
            interpolation = 0
        logging.debug('nearest interpolation')

    elif interpolation == 'linear':
        if method == 'opencv':
            interpolation = cv2.INTER_LINEAR
        else:
            interpolation = 1
        logging.debug('linear interpolation')
    elif interpolation == 'area':
        if method == 'opencv':
            interpolation = cv2.INTER_AREA
        else:
            raise Exception('Method not defined')
        logging.debug('area interpolation')
    elif interpolation == 'lanczos4':
        if method == 'opencv':
            interpolation = cv2.INTER_LANCZOS4
        else:
            interpolation = 4
        logging.debug('lanczos/biquartic interpolation')
    else:
        raise Exception('Interpolation method not available')

    t, z, h, w = movie.shape

    output = movie.astype(np.float32)

    for i in range(t):
        sh_z, sh_h, sh_w = shifts[i]

        if method == 'opencv':
            # the usual 2D:
            Mwh = np.float32([
                [1, 0, sh_w],
                [0, 1, sh_h]
            ])
            min_, max_ = np.min(output[i]), np.max(output[i])
            if sh_z != 0: # 3D part
                Mwz = np.float32([
                    [1, 0, 0],
                    [0, 1, sh_z]
                ])
                for ih in range(h):
                    # frame_shifted = cv2.warpAffine(frame[:, ih].astype(np.float32), Mwz, (w, z),
                    frame_shifted = cv2.warpAffine(output[i,:, ih], Mwz, (w, z),
                                                   flags = interpolation, borderMode = cv2.BORDER_CONSTANT, borderValue=0)
                    output[i, :, ih, :] = np.clip(frame_shifted, min_, max_)
            for iz in range(z):
                frame_shifted = cv2.warpAffine(output[i,iz], Mwh, (w, h), flags = interpolation,
                                               borderMode = cv2.BORDER_CONSTANT, borderValue=0)
                output[i, iz, :] = np.clip(frame_shifted, min_, max_)

        elif method == 'skimage':
            raise NotImplementedError("Not yet implemented using scikit-image. Try using opencv instead.")
            # tform = AffineTransform( translation = (-sh_z, -sh_w, -sh_h), dimensionality = 3 )
            # output[i] = warp(frame, tform, preserve_range = True, order = interpolation)

        else:
            raise Exception('Unknown shift application method')
    return output