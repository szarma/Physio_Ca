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

class ExperimentDatabase:
    _logger: logging.Logger
    connection: sqlite3.Connection
    _cursor: sqlite3.Cursor

    def __init__(self, database_path: Path, verbose: bool = False):
        """Standard constructor of a database manager for CTN experiments/recordings.

        Constructor of an instance of the ExperimentDatabase. The database uses sqlite as an engine.

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
