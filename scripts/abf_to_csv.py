#!/opt/tljh/user/envs/physio/bin/python

import argparse
import logging
import sys
from pathlib import Path

import dask.dataframe
from dask.diagnostics import ProgressBar
import numpy as np
import pyabf


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Utility to convert files using Axon Bioformat (ABF) to CSV format.')
    parser.add_argument('--input-file', '-i',
                        required=True,
                        nargs=1,
                        dest='input_file',
                        help='Path to the file, which shall be converted.')
    parser.add_argument('--verbose',
                        dest='verbose',
                        action='store_true',
                        help='Triggers output of more information while running.')
    return parser.parse_args()


def prepare_logger(_args: argparse.Namespace) -> logging.Logger:
    log_level = logging.INFO
    if _args.verbose:
        log_level = logging.DEBUG
    logger = logging.getLogger(__name__)
    logger.setLevel(log_level)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(log_level)
    logger.addHandler(handler)
    return logger


if __name__ == '__main__':
    cmd_args = parse_args()
    log = prepare_logger(cmd_args)
    file_path = Path(cmd_args.input_file[0])
    if not file_path.is_file():
        raise FileNotFoundError

    abf = pyabf.ABF(file_path.as_posix())
    log.info(abf.__str__())

    for i in range(0, abf.sweepCount):
        log.info(f'Converting sweep ({i + 1}/{abf.sweepCount})...')
        abf.setSweep(sweepNumber=i)

        dask_df = dask.dataframe.from_array(np.vstack((abf.sweepX, abf.sweepY, abf.sweepC)).T,
                                            columns=[abf.sweepLabelX, abf.sweepLabelY, abf.sweepLabelC])
        log.info(f'Saving sweep #{i}...')
        with ProgressBar():
            dask.dataframe.to_csv(df=dask_df,
                                  filename=file_path.as_posix().replace(file_path.suffix, f'_sweep{i}.csv'),
                                  single_file=True)
    log.info('Process finished')

    log.info('Used protocols:')
    log.info('-' * 20)
    [log.info(f'{t} sec: {s}') for t, s in zip(abf.tagTimesSec, abf.tagComments)]
