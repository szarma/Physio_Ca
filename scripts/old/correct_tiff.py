#!/opt/tljh/user/envs/physio/bin/python

import argparse
import logging
from pathlib import Path
from sys import stdout

import numpy as np
from skimage.exposure import rescale_intensity
from skimage.io import imread, imsave


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog='correct_tiff',
                                     description='Utility to swap channels in a TIFF-file')
    parser.add_argument('files',
                        metavar='N',
                        help='File(s) to be corrected.',
                        nargs='+')
    parser.add_argument('--swap-red-green', '-rg',
                        help='Specifies whether to swap the red and green channel.',
                        action='store_true',
                        dest='swap_red_green')
    parser.add_argument('--lut',
                        dest='LUT',
                        nargs=3,
                        default=None,
                        help='LUT values for R, G, B channels of the image.')
    parser.add_argument('--recursive', '-R',
                        action='store_true',
                        dest='is_recursive',
                        help='Specifies if paths are given from directories and every TIFF per directory shall be '
                             'considered.')
    parser.add_argument('--verbose',
                        action='store_true',
                        dest='verbose',
                        help='Trigger to display more information while executing the script.')
    return parser.parse_args()


def prepare_logger(_args: argparse.Namespace) -> logging.Logger:
    log_level = logging.INFO
    if _args.verbose:
        log_level = logging.DEBUG
    logger = logging.getLogger(__name__)
    logger.setLevel(log_level)
    handler = logging.StreamHandler(stdout)
    handler.setLevel(log_level)
    logger.addHandler(handler)
    return logger


if __name__ == '__main__':
    args = parse_args()
    log = prepare_logger(args)
    files = args.files
    if args.is_recursive:
        files_tmp = []
        log.debug(f'Searching files in directories:')
        for entry in files:
            r_files = []
            cmd_dir = Path(entry)
            if cmd_dir.is_dir():
                log.debug(f'Searching .tif files in "{entry}"')
                res = cmd_dir.glob('*.tif')
                for file in res:
                    if file.is_file():
                        r_files.append(file)
                        log.debug(f'Added "{file}" to list.')
            [files_tmp.append(file) for file in r_files]
        files = files_tmp

    if args.swap_red_green:
        log.info(f'Swapping channels (R <-> G)')
        for file in files:
            path = Path(file)
            if path.is_file():
                image: np.ndarray = imread(path.as_posix())
                r = image[:, :, 0].copy()
                g = image[:, :, 1].copy()
                b = image[:, :, 2].copy()
                if args.LUT is not None:
                    log.info(f'LUT was set for adjustments.')
                    log.info(f'R: {args.LUT[0]}, G: {args.LUT[1]}, B: {args.LUT[2]}')
                    r = rescale_intensity(image[:, :, 0], in_range=(r.min(), int(args.LUT[0])))
                    g = rescale_intensity(image[:, :, 1], in_range=(g.min(), int(args.LUT[1])))
                    b = rescale_intensity(image[:, :, 2], in_range=(b.min(), int(args.LUT[2])))

                converted_image = np.dstack((g, r, b))
                save_path = path.as_posix().replace(path.suffix, '_converted.tif')
                imsave(save_path, converted_image, check_contrast=False)
                log.info(f'Converted "{path.as_posix()}" and saved it to "{save_path}".')
            else:
                log.debug(f'"{path.as_posix()}" does not exist or is no file! Skipped.')
