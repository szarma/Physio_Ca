import bioformats
import datetime
import javabridge
import logging
from nd2reader import ND2Reader
import pandas as pd
from pathlib import Path
import re
import sqlite3
from typing import List, NamedTuple, Optional, Tuple, Union, Pattern
from sys import stdout

SUPPORTED_FORMATS = ['.nd2', '.lif']


class DatabaseManager:
    _series_regex: str
    _protocol_regex: str
    _pkl_regex: str
    _compiled_ser_regex: Pattern[str]
    _compiled_pro_regex: Pattern[str]
    _compiled_pkl_regex: Pattern[str]
    _logger: logging.Logger
    connection: sqlite3.Connection
    _cursor: sqlite3.Cursor

    def __init__(self, database_path: Path, verbose: bool):
        """
        Standard constructor of a database manager for CTN experiments/recordings.

        :param database_path: Path of the database, which contains the CTN data.
        :param verbose: Sets whether debug log is printed.
        """

        # Initialize logging functions for output
        log_level = logging.DEBUG if verbose else logging.INFO
        self._logger = logging.getLogger(__name__)
        self._logger.setLevel(log_level)
        handler = logging.StreamHandler(stdout)
        handler.setLevel(log_level)
        self._logger.addHandler(handler)

        # Initialize connection to database
        if not database_path.is_file():
            self._logger.error(f'Could not find file "{database_path.as_posix()}". New database will be created!')
        try:
            self.connection = sqlite3.connect(database_path)
            self._cursor = self.connection.cursor()
        except sqlite3.DatabaseError as e:
            self._logger.error(f'"Unable to open {database_path.as_posix()}"!')
            self._logger.error(f'Error: {e}')
            raise IOError

        # Check if tables already exist and create them if they are not there.
        self._logger.info('Checking if tables exist (create if not existing)...')
        self._cursor.execute('''create table if not exists recording
        (id integer primary key autoincrement, path text not null , metadata integer not null)''')
        self._cursor.execute('''create table if not exists metadata
        (id integer primary key autoincrement, date date, frequency real, pixel_size real, size_x integer,
        size_y integer, size_z int, size_t int)''')
        self._cursor.execute('''create table if not exists experiment
        (id integer primary key autoincrement, pickles text, series text, path_to_movie text, path_to_pickles text,
        path_to_additional_info text, recording integer, duration real)''')
        self.connection.commit()

        # Create regex entries and compile them to boost performance:
        self._series_regex = r'\w+\d+-?\d+'
        self._protocol_regex = r'.*protocol[.]txt'
        self._pkl_regex = r'(?P<pickle_param>\d(?:\.\d)?)_rois[.]pkl'
        self._compiled_ser_regex = re.compile(self._series_regex)
        self._compiled_pro_regex = re.compile(self._protocol_regex)
        self._compiled_pkl_regex = re.compile(self._pkl_regex)

    def __del__(self):
        self.connection.commit()
        self._cursor.close()
        self.connection.close()

    def refresh(self, root_dir: Path, wipe: bool = False):
        """
        Function refreshes the whole database structure under the given file tree.

        :param root_dir: Directory, which shall be parsed recursively for all experiments.
        :param wipe: Indicates whether existing data should be deleted before refreshing.
        """
        recordings: List[Path]

        recordings = self.__get_recordings(root_dir=root_dir)
        self.__set_recordings(recording_list=recordings, wipe=wipe)
        self.__set_experiments(recording_list=recordings, wipe=wipe)
        raise NotImplementedError

    def __get_experiments(self, recording_path: Path) -> List[Tuple[str, Path]]:
        experiments: List[Tuple[str, Path]] = []

        analysis_directory = recording_path.parent.joinpath(recording_path.name + '_analysis/')
        if analysis_directory.is_dir():
            self._logger.debug(f'Analysis directory: {analysis_directory.as_posix()}')
            existing_series = []
            for child in analysis_directory.iterdir():
                if child.is_dir():
                    if self._compiled_ser_regex.match(child.name) is not None:
                        self._logger.debug(f'Matched {child.name} as a valid series of the analysis directory.')
                        experiments.append((child.name, recording_path.joinpath(child.name)))
                self._logger.debug(f'Series found: {existing_series.__str__()}')

        return experiments

    def __set_experiments(self, recording_list: List[Path], wipe: bool = False) -> None:
        if wipe:
            self._logger.info(f'Deleting existing data in order to refresh...')
            self._cursor.execute('''delete from experiment''')
        for recording in recording_list:
            experiments = self.__get_experiments(recording)
            for experiment_name, experiment_path in experiments:
                pass
        raise NotImplementedError

    def __get_recordings(self, root_dir: Path) -> List[Path]:
        file_list: List[Path] = []
        for path_object in root_dir.rglob('*'):
            try:
                if path_object.suffix in SUPPORTED_FORMATS:
                    self._logger.debug(f'Found file while indexing: "{path_object.as_posix()}"')
                    file_list.append(path_object)
            except PermissionError as e:
                self._logger.error(f'Permission denied accessing file: "{e.filename}"')
        self._logger.debug(f'Found {len(file_list)} recordings.')
        return file_list

    def __set_recordings(self, recording_list: List[Path], wipe: bool = False) -> None:
        if wipe:
            self._logger.info(f'Deleting existing data in order to refresh...')
            self._cursor.execute('''delete from recording''')
        for recording in recording_list:
            metadata = self.__get_metadata(recording_path=recording, identifier=None)
            metadata_entry = SqlMetadataEntry(
                id=0,
                date=metadata['date'],
                frequency=.0,
                pixel_size=.0,
                size_x=0,
                size_y=0,
                size_z=0,
                size_t=0
            )
            self._cursor.execute('''insert or replace into metadata(date) values (?)''',
                                 (metadata_entry.date,))
            metadata_id = self._cursor.lastrowid

            recording_entry = SqlRecordingEntry(
                id=0,
                path=recording,
                metadata=metadata_id
            )
            self._cursor.execute('''insert or replace into recording(path, metadata) values (?,?)''',
                                 (recording_entry.path.as_posix(), recording_entry.metadata))
        self.connection.commit()

    def __get_metadata(self, recording_path: Path, identifier: Optional[Union[str, int]] = None) -> pd.DataFrame:
        metadata = pd.DataFrame(data={'Name': None,
                                      'SizeT': 0,
                                      'SizeX': 0,
                                      'SizeY': 0,
                                      'SizeZ': 0,
                                      'pxSize': 0.0,
                                      'bit depth': None,
                                      'Frequency': 0.0,
                                      'Start time': datetime.datetime.now(),
                                      'End time': datetime.datetime.now(),
                                      'Duration': None,
                                      'line scan': 0,
                                      'date': datetime.datetime.now()})
        metadata_exists = recording_path.parent.joinpath(f'.{recording_path.name}.meta').is_file()
        if metadata_exists:
            meta_path = recording_path.parent.joinpath(f'.{recording_path.name}.meta')
            try:
                self._logger.debug(f'Reading metadata of file "{meta_path.as_posix()}"')
                metadata = pd.read_csv(meta_path, index_col=0)
            except IOError:
                self._logger.error(f'Could not parse metadata of file "{meta_path.as_posix()}" !')
                self._logger.debug(f'Setting metadata to None.')
        else:
            try:
                javabridge.start_vm(class_path=bioformats.JARS, run_headless=True)
            except RuntimeError:
                self._logger.error(f'Error occurred while starting JVM! Please restart the kernel and try again.')
                exit()
            xml = bioformats.OMEXML(bioformats.get_omexml_metadata(recording_path))
            num_images = xml.get_image_count()
            if num_images > 0:
                if recording_path.suffix == '.nd2':
                    with ND2Reader(recording_path.as_posix()) as reader:
                        metadata['date'] = reader.metadata['date']
                else:
                    metadata['date'] = xml.image().AcquisitionDate
            else:
                self._logger.info(f'No recording data found in "{recording_path.as_posix()}". Skipped file.')

        return metadata


class SqlRecordingEntry(NamedTuple):
    id: int
    path: Path
    metadata: int


class SqlExperimentEntry(NamedTuple):
    id: int
    pickles: List[str]
    series: Union[List[str], str]
    movie_path: Path
    protocol_path: Path
    additional_info_path: Path
    recording: int
    duration: float


class SqlMetadataEntry(NamedTuple):
    id: int
    date: datetime.datetime
    frequency: float
    pixel_size: float
    size_x: int
    size_y: int
    size_z: int
    size_t: int
