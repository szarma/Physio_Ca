from copy import deepcopy
from datetime import date
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Optional, Union

from .IDataSerializer import IDataSerializer
from .. import Regions


class Hdf5Serializer(IDataSerializer):
    """Class for serialization of regions created by the islets module into HDF5-Format.

    The dataframe will be stored into "/df" of the HDF5-file.
    Every other value of the object will be stored in "/__dict__".

    Documentation of the HDF5 file format: https://www.hdfgroup.org/solutions/hdf5
    """

    __file_path: Path

    def __init__(self, file_path: Union[Path, str]):
        """Standard constructor of the Hdf5Serializer.

        Parameters
        ----------
        file_path : Path
            Path of the file to be loaded/saved.
        """
        self.filepath = Path(file_path) if type(file_path) == str else file_path

    @property
    def filepath(self) -> Path:
        return self.__file_path

    @filepath.setter
    def filepath(self, val: Path) -> None:
        val = Path(val.as_posix().replace(' ', '_'))
        self.__file_path = val

    def load(self) -> Regions:
        """Load the specified file from HDF5 into regions.

        Returns
        -------
        Regions
            Regions object loaded from the HDF5-file.

        Raises
        ------
        FileNotFoundError
            If the file cannot be found or does not exist.
        """
        if not self.__file_path.is_file():
            raise FileNotFoundError

        regions = Regions({(0, 0): [0, 0]}, mode='dummy')
        store = pd.HDFStore(self.__file_path.as_posix())
        series = store['__dict__']
        regions.__dict__ = series.to_dict()
        regions.df = store['df']
        return regions

    def save(self,
             regions: Regions,
             movie: Optional[np.ndarray] = None,
             columns_to_save: Optional[List[str]] = None,
             add_date: bool = True) -> None:
        """Saves the regions into a HDF5 file at the specified location.

        Parameters
        ----------
        regions : Regions
            Regions object, which shall be saved into the specified file.
        movie : numpy.ndarray, optional
            Movie, which shall be applied as an update before saving regions.
        columns_to_save : List[str], optional
            Columns of the dataframe to save except the necessary ones.
        add_date : bool, default: True
            Determines whether the current date should get added to the name of the saved file.
        """
        if columns_to_save is None:
            columns_to_save = ['trace']
        if movie is not None:
            regions.update(movie)

        if not self.__file_path.parent.is_dir():
            self.__file_path.parent.mkdir(parents=True)
        if add_date:
            today = date.today()
            self.__file_path = Path(self.__file_path.parent,
                                    '_'.join([today.strftime('%Y_%m_%d'), self.__file_path.stem]),
                                    self.__file_path.suffix)

        store = pd.HDFStore(self.__file_path.as_posix())
        saving = ['statImages', 'mode', 'image', 'filterSize', 'trange', 'FrameRange',
                  'analysisFolder', 'time', 'Freq', 'metadata']

        copy_movie = hasattr(self, 'movie')
        tmp_movie = None
        if copy_movie:
            tmp_movie = regions.movie
            del regions.movie
        regions_to_save = deepcopy(regions)
        all_attributes = list(regions_to_save.__dict__.keys())

        # Serialize dataframe:
        for k in regions.df.columns:
            if k not in ['peak', 'pixels', 'peakValue', 'tag', 'interest'] + columns_to_save:
                del regions_to_save.df[k]
        store['df'] = regions_to_save.df

        # Serialize other object values
        for k in all_attributes:
            if k not in saving:
                del regions_to_save.__dict__[k]
        save_dict = pd.Series(regions_to_save.__dict__)
        store['__dict__'] = save_dict
        store.close()

        if copy_movie:
            regions.movie = tmp_movie
            del tmp_movie
