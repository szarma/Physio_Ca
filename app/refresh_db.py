import argparse
import logging
import pandas as pd
from pathlib import Path
import re
import sqlite3
from sys import stdout
import tqdm
from typing import List, NamedTuple, Optional
from DatabaseManager import DatabaseManager


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Utility to swap channels in a TIFF-file')
    parser.add_argument('--root-dir', '-rd',
                        default='/data/',
                        dest='root_dir',
                        help='Specifies the root directory, where the algorithm should look for updates.')
    parser.add_argument('--verbose',
                        action='store_true',
                        dest='verbose',
                        help='Trigger to display more information while executing the script.')
    return parser.parse_args()


def prepare_logger(is_verbose: bool) -> logging.Logger:
    log_level = logging.DEBUG if is_verbose else logging.INFO
    logger = logging.getLogger(__name__)
    logger.setLevel(log_level)
    handler = logging.StreamHandler(stdout)
    handler.setLevel(log_level)
    logger.addHandler(handler)
    return logger


def get_file_paths(root_dir: Path, log_out: logging.Logger) -> List[Path]:
    file_list: List[Path] = []
    for path_object in root_dir.rglob('*'):
        try:
            if path_object.suffix in ['.lif', '.nd2']:
                if path_object.is_file():
                    log_out.debug(f'Found file while indexing: "{path_object.as_posix()}"')
                    file_list.append(path_object)
        except PermissionError as e:
            log_out.error(f'Permission denied accessing file: "{e.filename}"')
    log_out.info(f'{len(file_list)} files (.lif or .nd2) found under "{root_dir.as_posix()}"')
    return file_list


class SqlEntry(NamedTuple):
    date: str
    experiment_directory: str
    experiment: str
    series: str
    frequency: float
    protocol: str
    protocol_path: str
    pickles: List[str]
    pixel_size: float
    size_x: int
    size_y: int
    size_z: int
    size_t: int
    # comments: str
    # comment_file: str
    duration: float


if __name__ == '__main__':
    args = parse_args()
    log = prepare_logger(args.verbose)
    log.info('Establishing connection to SQLite database...')
    dbm = DatabaseManager(Path('/home/jupyter-johannes/test_db.db'), verbose=True)

    file_paths = get_file_paths(Path(args.root_dir), log)
    file_paths = sorted(file_paths)
    series_regex = r'\w+\d+-?\d+'
    protocol_regex = r'.*protocol[.]txt'
    pkl_regex = r'(?P<pickle_param>\d(?:\.\d)?)_rois[.]pkl'
    compiled_ser_regex = re.compile(series_regex)
    compiled_pro_regex = re.compile(protocol_regex)
    compiled_pkl_regex = re.compile(pkl_regex)

    for path in tqdm.tqdm(file_paths):
        recording_type = ''
        if path.suffix == '.lif':
            recording_type = 'Leica'
        elif path.suffix == '.nd2':
            recording_type = 'Nikon'
        else:
            raise IOError('Module is not able to read the specified file format!')

        metadata: Optional[pd.DataFrame] = None
        metadata_exists = path.parent.joinpath(f'.{path.name}.meta').is_file()
        if metadata_exists:
            meta_path = path.parent.joinpath(f'.{path.name}.meta')
            try:
                log.debug(f'Reading metadata of file "{meta_path.as_posix()}"')
                metadata = pd.read_csv(meta_path, index_col=0)
            except IOError:
                log.error(f'Could not parse metadata of file "{meta_path.as_posix()}" !')
                log.debug(f'Setting metadata to None.')

        analysis_directory = path.parent.joinpath(path.name + '_analysis/')
        if analysis_directory.is_dir():
            log.debug(f'Found analysis directory at "{analysis_directory.as_posix()}"')
            existing_series = []
            for child in analysis_directory.iterdir():
                if child.is_dir():
                    if compiled_ser_regex.match(child.name) is not None:
                        log.debug(f'Matched {child.name} as a valid series of the analysis directory.')
                        existing_series.append(child.name)
            log.debug(f'Series found: {existing_series.__str__()}')
            for series in existing_series:
                pickles: List[str] = []
                protocol_file: Optional[Path] = None
                protocol: Optional[pd.DataFrame] = None
                series_directory = analysis_directory.joinpath(series)
                item_list = series_directory.glob('**/*')
                time: Optional[float] = None
                frequency: Optional[float] = None

                files = [f for f in item_list if f.is_file()]
                for file in files:
                    protocol_match = compiled_pro_regex.match(file.name)
                    if protocol_match:
                        protocol_file = file
                        log.debug(f'Found protocol file "{file.as_posix()}".')
                        try:
                            protocol = pd.read_csv(protocol_file.as_posix(), dtype=str)
                            protocol.dropna(how='all', inplace=True)
                            protocol['t_begin'] = pd.to_timedelta(['00:' + el if type(el) == str else '00:00:00'
                                                                   for el in protocol['begin']]).total_seconds()
                            if metadata is not None:
                                time = pd.to_timedelta(metadata['Duration']).sum().total_seconds()
                                frequency = metadata.loc[0, 'Frequency']
                                protocol['t_end'] = pd.to_timedelta(['00:' + el if type(el) == str else time /
                                                                    frequency * 1e9 for el in
                                                                    protocol['end']]).total_seconds()
                            else:
                                # TODO: Add functionality to parse metadata again to be able to add protocol information
                                protocol = None

                        except FileNotFoundError as fe:
                            log.error(f'File "{protocol_file.as_posix()}" could not be found! ({fe}')
                        except ValueError as ve:
                            log.error(f'Protocol file "{protocol_file.as_posix()}" holds values, which could not be '
                                      f'parsed! ({ve})')
                        except PermissionError as pe:
                            log.error(f'No permissions to read file "{protocol_file.as_posix()}"! ({pe})')
                    else:
                        pickle = compiled_pkl_regex.match(file.name)
                        if pickle:
                            pickles.append(str(pickle.group('pickle_param')))

                # TODO: Add functionality that series get their corresponding metadata
                size_x: Optional[int] = metadata['SizeX'].unique()[0] if metadata is not None else None
                size_y: Optional[int] = metadata['SizeY'].unique()[0] if metadata is not None else None
                size_z: Optional[int] = None
                if metadata is not None:
                    size_z = metadata['SizeZ'].unique()[0] if 'SizeZ' in metadata.columns else None
                size_t: Optional[int] = metadata['SizeT'].sum() if metadata is not None else None
                px_size: Optional[float] = metadata.loc[0, 'pxSize'] if metadata is not None else None
                entry = SqlEntry(
                    date=metadata.loc[0, 'Start time'] if metadata is not None else pd.NaT,
                    experiment_directory=series_directory.as_posix(),
                    protocol_path=None if protocol_file is None else protocol_file.as_posix(),
                    pickles=pickles,
                    series=series,
                    experiment=path.name,
                    frequency=frequency,
                    protocol=protocol.to_string() if protocol is not None else '',
                    size_x=size_x,
                    size_y=size_y,
                    size_z=size_z,
                    size_t=size_t,
                    pixel_size=px_size,
                    duration=time
                )
                try:
                    cursor.execute(
                        '''insert or replace into experiments(date, experiment_directory, protocol_path, pickles, series,
                           experiment, frequency, protocol, pxSize, SizeX, SizeY, SizeZ, SizeT, duration)
                           values(?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',
                        (entry.date, entry.experiment_directory, entry.protocol_path, str(entry.pickles), entry.series,
                         entry.experiment, entry.frequency, entry.protocol, entry.pixel_size, entry.size_x, entry.size_y,
                         entry.size_z, entry.size_t, entry.duration)
                    )
                except sqlite3.InterfaceError as sql_error:
                    log.error(f'Could not add experiment to database!.')
                    log.error(f'Experiment path: {path.as_posix()}')
                    log.error(f'Experiment: {path.name}')
                    log.error(f'Series: {series}')
                log.debug(f'Available pickles for series "{series}": {pickles}')

    log.info('Serializing database to file...')
    connection.commit()
    connection.close()
    log.info('Database refreshed.')
