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

import islets.Recording
from . import SUPPORTED_FORMATS
from ..Recording import Recording, parse_leica
from ..utils import get_series_dir

class ExperimentDatabase:
    _logger: logging.Logger
    connection: sqlite3.Connection
    _cursor: sqlite3.Cursor

    def __init__(self, database_path: Path, verbose: bool = False):
        """Standard constructor of a database manager for CTN experiments/recordings.

        Constructor of an instance of the ExperimentDatabase. The database uses sqlite as an engine.

        Raises
        ------
        IOError
            Gets raised if the database cannot be opened due to failures during IO.

        Parameters
        ----------
        database_path : Path
            Path of the database, which contains the CTN data.
        verbose : bool, default: False
            Sets whether debug log is printed.
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

    def __del__(self):
        """Standard destructor.

        Standard destructor of the database manager. To prevent losing unsaved changes it commits all actions of the
        current connection and closes the cursor for this instance.
        """
        self.connection.commit()
        self._cursor.close()
        self.connection.close()
        self._logger.info('Connection closed.')

    def __get_recording(self, root_dir: Path) -> List[Path]:
        """Function getting all recordings.

        Parameters
        ----------
        root_dir : Path
            Directory, which is considered as root to search for experiment files.

        Returns
        -------
        List[Path]
            List of paths to all found recordings.
        """
        file_list: List[Path] = []
        for path_object in root_dir.rglob('*'):
            try:
                if path_object.suffix in SUPPORTED_FORMATS:
                    self._logger.debug(f'Found file while indexing: "{path_object.as_posix()}"')
                    file_list.append(path_object)
            except PermissionError as e:
                self._logger.error(f'Permission denied accessing file: "{e.filename}"')

        return file_list

    def refresh(self, root_dir: Path, wipe: bool = False) -> None:
        """Refresh the database.

        Function, which can be called to refresh the database. It can either be rewritten by wiping all previous data or
        if can just add newly found entries to the existing database.
        Experiments will be searched in each subfolder of the root directory.

        Parameters
        ----------
        root_dir : Path
            Path of the directory, which will be used as root directory to refresh the database.
        wipe : bool, default: False
            Indicator whether the database should be rewritten completely. Default: False

        Returns
        -------
        None
        """
        recordings: List[Path]

        recordings = self.__get_recording(root_dir=root_dir)
        for recording in recordings:
            rec_data = Recording(recording.as_posix())
            analysis_folder = recording.parent / (recording.name + '_analysis')
            if not analysis_folder.is_dir():
                self._logger.debug(f'No analysis folder could be found for "{recording.as_posix()}".')
            else:
                pass

        raise NotImplementedError
