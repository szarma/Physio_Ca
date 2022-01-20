from pathlib import Path
from typing import Union

from .IDataSerializer import IDataSerializer
from .Hdf5Serializer import Hdf5Serializer
from .JsonSerializer import JsonSerializer

__supported_file_types__ = {
    'hdf5': ['.hdf5', '.h5'],
    'json': ['.json', '.gz'],
    'pickle': ['.pkl']
}


def get_serializer(file_path: Union[str, Path], file_format: str = 'infer') -> IDataSerializer:
    """Function to get a serializer for loading or saving islet data.

    Parameters
    ----------
    file_path : Union[str, Path]
        Path of the file, which shall be deserialized.
    file_format : {'hdf5', 'json', 'infer'}, default: 'infer'
        Format, which shall be used to load the file. If 'infer' is passed, the function will automatically try to
        use the right format by guessing based on the file suffix.

    Raises
    ------
    KeyError
        Gets raised when the file type cannot be inferred or is not supported.

    Returns
    -------
    Serializer implementing the IDataSerializer interface for serialization and deserialization of data files.
    """

    serializer: IDataSerializer

    if type(file_path) == 'str':
        file_path = Path(file_path)

    if file_format == 'infer':
        for key, value in __supported_file_types__.items():
            if file_path.suffix in value:
                file_format = key

    if file_format == 'hdf5':
        serializer = Hdf5Serializer(file_path=file_path)
    elif file_format == 'json':
        serializer = JsonSerializer(file_path=file_path)
    elif file_format == 'pickle':
        raise NotImplementedError
    else:
        raise KeyError

    return serializer
