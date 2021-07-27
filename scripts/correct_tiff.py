#!/opt/tljh/user/envs/physio/bin/python

import argparse
import logging
from pathlib import Path
from sys import stdout

from PIL import Image

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
parser.add_argument('--verbose',
                    action='store_true',
                    dest='verbose',
                    help='Trigger to display more information while executing the script.')


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
    args = parser.parse_args()
    log = prepare_logger(args)

    if args.swap_red_green:
        log.info(f'Swapping channels (R <-> G)')
        for file in args.files:
            path = Path(file)
            if path.is_file():
                image = Image.open(path.as_posix())
                r, g, b = image.split()
                converted_image = Image.merge('RGB', (g, r, b))
                save_path = path.as_posix().replace(path.suffix, '_converted.tif')
                converted_image.save(save_path)
                log.info(f'Converted "{path.as_posix()}" and saved it to "{save_path}".')
            else:
                log.debug(f'"{path.as_posix()}" does not exist or is no file! Skipped.')
