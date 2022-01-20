from abc import ABC, abstractmethod
from pathlib import Path
from typing import Union
from .. import Regions


class IDataSerializer(ABC):
    __file_path: Path

    @abstractmethod
    def __init__(self, file_path: Union[Path, str]):
        pass

    @property
    @abstractmethod
    def filepath(self) -> Path:
        pass

    @filepath.setter
    @abstractmethod
    def filepath(self, val: Path) -> None:
        pass

    @abstractmethod
    def load(self, **kwargs) -> Regions:
        pass

    @abstractmethod
    def save(self, regions: Regions, **kwargs) -> None:
        pass
