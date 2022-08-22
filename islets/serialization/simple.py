import warnings
import logging
import h5py
import io
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict

def df2hdf(df: pd.DataFrame,
            buffer,
            key: str = 'dataframe',
            verbose: bool = False,
            overwrite: bool = False) -> None:
    assert isinstance(df, pd.DataFrame)
    if key in buffer.keys():
        if overwrite:
            del buffer[key]
        else:
            raise ValueError(
                f"Key {key} exists. if you want to overwrite it, you will "
                "need to specify `overwrite=True` in the function call."
            )
    df_handle = buffer.create_group(key)
    if verbose:
        print(f"'{key}' created in '{buffer.name}'")
    df_handle.create_dataset(name = "index", data = df.index.values)
    if verbose:
        print(f"\t'{key}/index' saved.")
    for col in df.columns:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error")
                v = df[col].values
                if hasattr(v[0], "shape") and len(v[0].shape):
                    try:
                        v = np.vstack(v)
                    except:
                        if verbose:
                            print(
                                f"Attempted to save {col} as a multidimensional "
                                "dataset, but it was unsuccessful."
                            )
                df_handle.create_dataset(data = v, name = col, compression="gzip")
                if verbose:
                    print(f"\t'{key}/{col}' saved")
        except (TypeError, pd.errors.PerformanceWarning):
            df_handle.create_dataset(data = df[col].apply(str).values, name = col, compression="gzip")
            if verbose:
                print(f"\t'{key}/{col}' saved (as strings)", end = ' ')
            coltypes = df[col].apply(type).apply(lambda xi: xi.__name__)
            if len(coltypes.value_counts()) == 1:
                k = f"_{col}_type"
                df_handle.attrs[k] = coltypes.iloc[0]
                if verbose:
                    # noinspection SqlNoDataSourceInspection
                    print(f"with type saved as an attribute '{k}' ")
            else:
                k = f"_{col}_types"
                df_handle.create_dataset(data = coltypes.values, name = k, compression="gzip")
                if verbose:
                    print(f"with types saved in '{k}'.")

def dict2hdf(dict_: Dict,
              buffer,
              key: str,
              in_attrs: bool = True,
              verbose: bool = False,
              overwrite: bool = False,
              ) -> None:

    if key in buffer.keys():
        if overwrite:
            del buffer[key]
        else:
            raise ValueError(
                f"Key {key} exists. if you want to overwrite it, you will "
                "need to specify `overwrite=True` in the function call."
            )
    df_handle = buffer.create_group(key)
    for k in dict_:
        if in_attrs:
            df_handle.attrs[k] = dict_[k]
        else:
            df_handle.create_dataset(data = dict_[k], name=k, compression="gzip")
    if verbose:
        print(f"{list(dict_.keys())} saved in {buffer.name}/{key}")

def fig2hdf(fig, buffer, key, overwrite = False):
    if key in buffer.keys():
        if overwrite:
            del buffer[key]
        else:
            raise ValueError(
                f"Key {key} exists. if you want to overwrite it, you will "
                "need to specify `overwrite=True` in the function call."
            )
    from PIL import Image
    import numpy as np
    my_stringIObytes = io.BytesIO()
    fig.savefig(my_stringIObytes, format = 'tif')
    my_stringIObytes.seek(0)
    im = Image.open(my_stringIObytes)
    imdata = np.asarray(im)
    write_to_hdf(imdata, buffer, key, stringify = False, overwrite = overwrite)
    Xset = buffer[key]
    Xset.attrs.create('CLASS', 'IMAGE', dtype = 'S6')
    # Xset.attrs.create('IMAGE_MINMAXRANGE', [0, 255], dtype = np.uint8)
    Xset.attrs.create('IMAGE_SUBCLASS', 'IMAGE_TRUECOLOR', dtype = 'S16')
    # Xset.attrs.create('IMAGE_VERSION', '1.2', dtype = 'S4')
    # Xset.attrs.create('INTERLACE_MODE', 'INTERLACE_PIXEL', dtype = 'S16')



def write_to_hdf(value, buff, key: str, stringify: bool = False, overwrite = False, verbose = False) -> None:
    if not overwrite:
        if key in buff.attrs or key in buff.keys():
            raise FileExistsError(
                f"Key '{key}' already exists. If you wish to overwrite it,"
                "you need to explicitly include overwrite=True when calling the function."
            )
    if stringify:
        value = str(value)
    try:
        buff.attrs[key] = value
        if verbose:
            print(f"'{key}' saved in '{buff.name}' as an attribute")
    except:
        buff.create_dataset(data = value, name=key, compression="gzip")
        if verbose:
            print(f"'{key}' saved in '{buff.name}' as a dataset")

def hdf2df(filename: [str, Path], key='_dataframe_', verbose=False, keep_open=False) -> pd.DataFrame:
    filename = Path(filename)
    hf = h5py.File(filename, "r")
    df_handle = hf[key]
    index = df_handle['index']
    df = pd.DataFrame(index=index)
    # if verbose:
    #     print("DataFrame instantiated")
    columns = [col for col in df_handle.keys() if col!="index" and not col.endswith("_types")]
    attrs = dict(df_handle.attrs)
    if verbose:
        print(f"Attributes are {attrs}")
    # allowed = "[],.()" for let in punctuation if let not in []
    for col in columns:
        if len(df_handle[col].shape)==1:
            df[col] = df_handle[col][:]
        else:
            df.loc[:,col] = [df_handle[col][j] for j in range(df.shape[0])]
        if df_handle[col].dtype == np.dtype("O"):
            df[col] = df[col].apply(lambda xi: xi.decode("UTF8"))
        if verbose:
            print(f"Column '{col}' successfully imported.")

        if f"_{col}_types" in df_handle.keys():
            logging.warning(f"Column {col} should consist of multiple types. This is an experimental",
                            "feature. You should check the result to assure it's correct.")
        if f"_{col}_type" in attrs or f"_{col}_types" in df_handle.keys():
            #coltype = df_handle.attrs[f"_{col}_type"]
            try:
                df[col] = [el if el.isalpha() else eval(el) for el in df[col]]
                ## making sure there are no alphabetical characters in the string
                ## this prevents us from storing dictionaries in this simple way.
                ## might be worth reconsidering...
                coltypes = df[col].apply(type)
                # if future, we might opt for automatic import of necessary module to recast
                # the values into the appropriate type (but seems to be too risky)
                # coltype = eval(coltype)
                # df[col] = df[col].apply(coltype)
                vc = coltypes.value_counts()
                if len(vc)==1:
                    coltype = vc.index[0]
                    if verbose:
                        print(f"Column '{col}' successfully coerced into '{coltype}'.")
                else:
                    assert False
            except:
                logging.warning(f"Column '{col}' could not be (fully) coerced into appropriate type(s).")
    if keep_open:
        df.hdf_file_handle = hf
    else:
        hf.close()

    return df