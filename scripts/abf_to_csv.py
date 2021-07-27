#!/opt/tljh/user/envs/physio/bin/python

import argparse
import logging
import os
import sys
import uuid
from pathlib import Path

import numpy as np
import pandas as pd
import pyabf
from pandas.core.internals import BlockManager


class BlockManagerUnconsolidated(BlockManager):
    def __init__(self, *args, **kwargs):
        BlockManager.__init__(self, *args, **kwargs)
        self._is_consolidated = False
        self._known_consolidated = False

    def _consolidate_inplace(self): pass

    def _consolidate(self): return self.blocks


def df_from_arrays(arrays, columns, index) -> pd.DataFrame:
    from pandas.core.internals import make_block

    def gen():
        _len = None
        p = 0
        for a in arrays:
            if _len is None:
                _len = len(a)
                assert len(index) == _len
            assert _len == len(a)
            yield make_block(values=a.reshape((1, _len)), placement=(p,))
            p += 1

    blocks = tuple(gen())
    mgr = BlockManagerUnconsolidated(blocks=blocks, axes=[columns, index])
    return pd.DataFrame(mgr, copy=False)


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
    cmd_args = parser.parse_args()
    log = prepare_logger(cmd_args)
    file_path = Path(cmd_args.input_file[0])
    if not file_path.is_file():
        raise FileNotFoundError

    abf = pyabf.ABF(file_path.as_posix())
    log.info(abf.__str__())

    for i in range(0, abf.sweepCount):
        log.info(f'Converting sweep ({i + 1}/{abf.sweepCount})...')
        abf.setSweep(sweepNumber=i)

        # Create memory map for control:
        tmp_c_filename = '/tmp/' + str(uuid.uuid4()) + '.tmp'
        log.debug(f'Creating temporary file "{tmp_c_filename}" serving as memory map...')
        mem_c = np.memmap(tmp_c_filename, mode='w+', shape=abf.sweepC.shape, dtype=abf.sweepC.dtype)
        mem_c[:] = abf.sweepC[:]
        mem_c.flush()
        mem_c.flags['WRITEABLE'] = False

        # Create memory map for time:
        tmp_x_filename = '/tmp/' + str(uuid.uuid4()) + '.tmp'
        log.debug(f'Creating temporary file "{tmp_x_filename}" serving as memory map...')
        mem_x = np.memmap(tmp_x_filename, mode='w+', shape=abf.sweepX.shape, dtype=abf.sweepX.dtype)
        mem_x[:] = abf.sweepX[:]
        del abf.sweepX
        mem_x.flush()
        mem_x.flags['WRITEABLE'] = False

        # Create memory map for measured data:
        tmp_y_filename = '/tmp/' + str(uuid.uuid4()) + '.tmp'
        log.debug(f'Creating temporary file "{tmp_y_filename}" serving as memory map...')
        mem_y = np.memmap(tmp_y_filename, mode='w+', shape=abf.sweepY.shape, dtype=abf.sweepY.dtype)
        mem_y[:] = abf.sweepY[:]
        del abf.sweepY
        mem_y.flush()
        mem_y.flags['WRITEABLE'] = False

        log.debug('Creating memory mapped dataframe...')
        df = df_from_arrays([mem_x, mem_y, mem_c],
                            columns=[abf.sweepLabelX, abf.sweepLabelY, abf.sweepLabelC],
                            index=range(len(mem_x)))
        df = df.T
        log.info(f'Saving sweep #{i}...')
        df.to_csv(file_path.as_posix().replace(file_path.suffix, f'_sweep{i}.csv'), index=False)
        del df

        log.info('Deleting temporary files...')
        os.remove(tmp_c_filename)
        log.debug(f'Removed {tmp_c_filename}')
        os.remove(tmp_x_filename)
        log.debug(f'Removed {tmp_x_filename}')
        os.remove(tmp_y_filename)
        log.debug(f'Removed {tmp_y_filename}')
    log.info('Process finished')

    log.info('Used protocols:')
    log.info('-' * 20)
    [log.info(f'{t} sec: {s}') for t, s in zip(abf.tagTimesSec, abf.tagComments)]
