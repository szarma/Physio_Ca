from copy import deepcopy
from datetime import date
import json

import mgzip
import numpy as np
import orjson
import pandas as pd
from pathlib import Path
from typing import List, Optional, Union

from .IDataSerializer import IDataSerializer
from islets import Regions


class JsonSerializer(IDataSerializer):
    __file_path: Path

    def __init__(self, file_path: Union[Path, str]):
        """Standard constructor of the JsonSerializer.

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
    def filepath(self, val: Path):
        val = Path(val.as_posix().replace(' ', '_'))
        self.__file_path = val

    def load(self) -> Regions:
        """Saves the regions into a JSON file at the specified location.

        Returns
        -------
        Deserialized regions object.

        Raises
        ------
        FileNotFoundError
            Gets raised when the specified file has not been found.
        """
        if not self.__file_path.is_file():
            raise FileNotFoundError

        use_compression = True if self.__file_path.suffix == '.gz' else False
        file_content : str
        if use_compression:
            with mgzip.open(self.__file_path.as_posix(), 'rt') as file:
                file_content = file.read()
        else:
            with open(self.__file_path.as_posix(), 'rt') as file:
                file_content = file.read()

        json_obj = orjson.loads(file_content.encode())

        df = pd.DataFrame.from_dict(json_obj['df'])
        df["trace"] = [np.array(tr) for tr in df["trace"]]
        df['peak'] = df['peak'].apply(lambda x: tuple(x))
        df['pixels'] = df['pixels'].apply(lambda x: [tuple(y) for y in x])
        try:
            df.index = df.index.astype(int)
        except:
            pass
        del json_obj['df']

        json_obj['metadata'] = pd.Series(json_obj['metadata'])
        json_obj['Freq'] = np.float64(json_obj['Freq'])
        json_obj['image'] = np.array(json_obj['image'])

        for key, value in json_obj['statImages'].items():
            json_obj['statImages'][key] = np.array(value)
        
        regions = Regions(dict(zip(df['peak'], df['pixels'])))
        for key, value in json_obj.items():
            regions.__setattr__(key, value)
        regions.df = df
        regions.update()
        regions.time = np.array(regions.time)

        regions = Regions(regions)
        protocol_files = list(self.__file_path.parent.glob('*protocol*.*'))
        if len(protocol_files) > 0:
            regions.import_protocol(protocol_files[0].as_posix())
        regions.pathToPickle = self.__file_path.as_posix()

        return regions

    def save(self,
             regions: Regions,
             movie: Optional[np.ndarray] = None,
             columns_to_save: Optional[List[str]] = None,
             add_date: bool = True,
             use_compression: bool = True) -> None:
        """Saves the regions into a JSON file at the specified location.

        Parameters
        ----------
        regions : Regions
            Regions object, which shall be saved.
        movie : numpy.ndarray, optional
            Movie, which is used to update regions before saving. If movie is None no update will be performed before.
        columns_to_save : List[str], optional
            List of columns indicating columns to save. If None only standard columns get saved.
        add_date : bool, default: True
            Indicates whether the date shall be added to the name of the file.
        use_compression: bool, default: True
            Indicates whether compression should be used for serialization of the data.
        """
        if columns_to_save is None:
            columns_to_save = ['trace']
        if movie is not None:
            regions.update(movie)
        if add_date:
            today = date.today()
            self.__file_path = Path(self.__file_path.parent,
                                    '_'.join([today.strftime('%Y_%m_%d'), self.__file_path.stem]),
                                    self.__file_path.suffix)
        if not self.__file_path.parent.exists():
            self.__file_path.parent.mkdir(parents=True)

        saving = ['statImages', 'mode', 'image', 'filterSize', 'df', 'trange', 'FrameRange', 'analysisFolder', 'time',
                  'Freq', 'metadata']

        tmp_movie = hasattr(regions, 'movie')
        if tmp_movie:
            movie = regions.movie
            del regions.movie
        all_keys = list(regions.__dict__.keys())
        tmp_regions = deepcopy(regions)
        if tmp_movie:
            regions.movie = movie
            del movie

        for k in all_keys:
            if k not in saving:
                del tmp_regions.__dict__[k]
        for k in regions.df.columns:
            if k not in ['peak', 'pixels', 'peakValue', 'tag', 'interest'] + columns_to_save:
                del tmp_regions.df[k]

        json_dict = {}
        for k, v in tmp_regions.__dict__.items():
            value = v
            if isinstance(v, (pd.DataFrame, pd.Series)):
                value = json.JSONDecoder().decode(v.to_json(double_precision=5))
            if isinstance(v, (np.float64, np.int64)):
                value = str(v)
            if isinstance(v, np.ndarray):
                value = v.tolist()
            if isinstance(v, dict):
                for k_, v_ in v.items():
                    if isinstance(v_, np.ndarray):
                        v[k_] = v_.tolist()
                value = v

            json_dict[k] = value

        json_string = orjson.dumps(json_dict).decode()

        if use_compression:
            if not self.__file_path.suffix == '.gz':
                self.__file_path = self.__file_path.parent / (self.__file_path.name + '.gz')
            with mgzip.open(self.__file_path.as_posix(), 'wt') as file:
                file.write(json_string)
        else:
            with open(self.__file_path.as_posix(), 'wt') as file:
                file.write(json_string)

        regions.pathToPickle = self.__file_path.as_posix()
