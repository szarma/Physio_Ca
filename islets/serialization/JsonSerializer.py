from copy import deepcopy
from datetime import date
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Optional, Union

from .IDataSerializer import IDataSerializer
from .. import Regions


class JsonSerializer(IDataSerializer):
    __file_path: Path

    def __init__(self, file_path: Union[Path, str]):
        """Standard constructor of the Hdf5Serializer.

        Parameters
        ----------
        file_path : Path
            Path of the file to be loaded/saved.
        """
        self.file_path = Path(file_path) if type(file_path) == str else file_path

    @property
    def filepath(self) -> Path:
        return self.__file_path

    @filepath.setter
    def filepath(self, val: Path):
        val = Path(val.as_posix().replace(' ', '_'))
        self.__file_path = val

    def load(self) -> Regions:
        """Saves the regions into a JSON file at the specified location.

        Parameters
        """
        pass

    def save(self,
             regions: Regions,
             movie: Optional[np.ndarray] = None,
             columns_to_save: Optional[List[str]] = None,
             add_date: bool = True,
             use_compression: bool = True) -> None:
        """Saves the regions into a JSON file at the specified location.

        Parameters
        ----------
        """
        pass
