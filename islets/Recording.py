import os
from sys import exc_info
from warnings import warn
import errno
import bioformats as bf
import numpy as np
import pandas as pd
from . import cmovie
from nd2reader import ND2Reader
from tifffile import memmap as tiffmemmap
from tifffile import TiffFile
from .utils import autocorr2d
from IPython.display import display

def parse_leica(rec,
                merge=True,
                skipTags=None,
                index=False,
                verbose=False,
               ):
    if not hasattr(rec,"metadata"):
        if verbose>0:
            print("Recording has no metadata.")
        return []
    if skipTags is None:
        skipTags = ["crop", "split", "frame", "every", "half", "snapshot", "-", "proj", "max", "resize"]

    toDrop = [i for i,row in rec.metadata.iterrows() if "Series" not in row["Name"] 
#               or row["SizeY"]*row["SizeT"]<2000 
              #or row["SizeT"]/row["Frequency"]<2
             ]
    for tag in skipTags:
        toDrop += [i for i,row in rec.metadata.iterrows() if tag in row["Name"].lower()]
    rec.metadata.drop(index=np.unique(toDrop), inplace=True)
    if len(rec.metadata)==0:
        warn("Can not parse any series from %s"%rec.path)
        return []
    if not merge:
        if index:
            return  list(zip(rec.metadata.index,rec.metadata.Name.values))
        else:
            return  list(rec.metadata.Name.values)
    if len(rec.metadata)>1:
        rec.calc_gaps()
        ff = np.any(
                [rec.metadata["gap"]>5]+
                [rec.metadata["gap"]<-5]+
                [~np.isfinite(rec.metadata["gap"])]+
                [rec.metadata[c].diff().abs()>0 for c in ["pxSize", "SizeX", "SizeY"]]+
                [rec.metadata["Frequency"].diff().abs()/rec.metadata["Frequency"]>.02],
            axis=0)
        sers = np.split(rec.metadata.Name.values, np.where(ff)[0])[1:]
        idxs = np.split(rec.metadata.index.values, np.where(ff)[0])[1:]
    else:
        sers = [rec.metadata.Name.values]
        idxs = [list(rec.metadata.index)]
    if verbose:
        print ("series grouped as following:", sers)
    outSer = []
    for serlist in sers:
#         if len(serlist)==len(rec.metadata):
#             ser="all"
#         else:
        serrange = [int(el.replace("Series","")) for el in serlist]
        if len(serrange)>1:
            ser = "Series%03i-%i"%(serrange[0],serrange[-1])
        else:
            ser = "Series%03i"%(serrange[0])
        outSer += [ser]
    if index:
        return list(zip(idxs,outSer))
    return outSer


class Recording:
    def __init__(self, pathToExperiment):
        if not os.path.isfile(pathToExperiment):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), pathToExperiment)
        self.folder, self.Experiment = os.path.split(pathToExperiment)
        self.path = pathToExperiment
        self.metafile = os.path.join(self.folder,  f".{self.Experiment}.meta")
        self.is_nd2 = os.path.splitext(self.path)[1] == '.nd2'
        if os.path.isfile(self.metafile):
            try:
                self.metadata = pd.read_csv(self.metafile, index_col=0)
                for c in ["Start time","End time"]:
                    self.metadata[c] = pd.to_datetime(self.metadata[c])
                for c in ["Duration"]:
                    self.metadata[c] = pd.to_timedelta(self.metadata[c])
                from .general_functions import td_nanfloor
                self.metadata["End time"] = self.metadata["End time"].apply(td_nanfloor)
            except:
                warn("Metadata not read correctly. WIll try to parse them now.")
        else:
            try:
                print (f"Recording {pathToExperiment} not yet preprocessed. Preprocessing takes a few seconds and will speed up the usage later...", end=" ")
                from sys import stdout
                stdout.flush()
                self.parse_metadata()
                self.save_metadata()
                print (f"Finished.", end="")
                stdout.flush()
            except:
                warn(f"Could not parse metadata. {exc_info()}")
            finally:
                print("")
        if hasattr(self, "metadata"):
            self.nSeries = len(self.metadata)
            if "Name" in self.metadata:
                self.allSeries = self.metadata.Name.values
                
    def __del__(self):
        if hasattr(self, "tempdir") and not self.tempdir is None:
            from os import system
            system(f"rm -rf {self.tempdir}")
        
    def parse_metadata(self,verbose=False):
        metadata = pd.DataFrame(columns=[
            "Name","SizeT","SizeX","SizeY","SizeZ","Nchannels",
            "pxSize","pxUnit","bit depth", "Frequency",
            "Start time", "End time", "Duration"
        ])
        if self.path.endswith(".tif") or self.path.endswith(".tiff"):
            self.metadata = metadata
            self.metadata.loc[0,"Name"] = "all"
            try:
                self.metadata.loc[0,"Frequency"] = float(os.path.split(self.path)[1].split("Hz")[0].strip("_").strip().split( "_")[-1].split()[-1])
            except:
                print(exc_info())
                pass
            return None
        if self.path.endswith(".stk"):
            stkFile = TiffFile(self.path)
            metadata.loc[0,"SizeT"] = stkFile.stk_metadata['NumberPlanes']
            metadata.loc[0,"Start time"] = stkFile.stk_metadata['CreateTime']
            metadata.loc[0,"End time"] = stkFile.stk_metadata['LastSavedTime']
            metadata.loc[0,"Duration"] = (stkFile.stk_metadata['LastSavedTime'] - stkFile.stk_metadata['CreateTime']).total_seconds()
            page = stkFile.pages[0].asarray()
            metadata.loc[0,"SizeX"] = page.shape[1]
            metadata.loc[0,"SizeY"] = page.shape[2]
            metadata.loc[0,"bit depth"] = page.dtype.name
            metadata.loc[0,"pxSize"] = "_"
            metadata.loc[0,"pxUnit"] = "_"
            metadata.loc[0,"Name"] = "all"
            # dts = np.diff(stkFile.stk_metadata["DatetimeCreated"])
            # dt = np.median(dts) / np.timedelta64(1, 's')
            # metadata.loc[0,"Frequency"] = 1/dt
            metadata["Frequency"] = metadata['SizeT']/metadata['Duration']
            self.metadata = metadata
            return None

        else:
            md = bf.get_omexml_metadata(self.path)
            self.xml = bf.OMEXML(md)
            self.nSeries = self.xml.get_image_count()
            self.allSeries = [self.xml.image(i).Name for i in range(self.nSeries)]


            for i in range(self.nSeries):
                im = self.xml.image(i)
                for c in metadata.columns:
                    try:
                        metadata.loc[i,c] = getattr(im.Pixels, c)
                        metadata.loc[i,c] = int(metadata.loc[i,c])
                    except:
                        pass
                metadata.loc[i,"Name"] = im.Name
                metadata.loc[i,"pxSize"] = im.Pixels.get_PhysicalSizeX()
                metadata.loc[i,"pxUnit"] = im.Pixels.get_PhysicalSizeXUnit()
                metadata.loc[i,"bit depth"] = im.Pixels.get_PixelType()
                metadata.loc[i,"Nchannels"] = im.Pixels.get_channel_count()
                if self.is_nd2:
                    try:
                        with ND2Reader(self.path) as nd2data:
                            date = nd2data.metadata['date']
                    except:
                        date = pd.to_datetime("")
                        print (f"Error reading date from {self.path}: ", exc_info())
                else:
                    date = im.get_AcquisitionDate()
                metadata.loc[i,"Start time"] = date
                if metadata.loc[i,"SizeZ"]>1:
                    metadata.loc[i,"Z-stack height"] = im.Pixels.get_PhysicalSizeZ()
                if metadata.loc[i,"SizeT"]>1:
                    lastT = im.Pixels.Plane(int(metadata.loc[i,"Nchannels"] * metadata.loc[i,"SizeT"]-1)).DeltaT
                    if lastT!=0 and lastT is not None:
                        if verbose:
                            print (im.Name, metadata.loc[i,"SizeT"], type(metadata.loc[i,"SizeT"]), lastT)
                        metadata.loc[i,"Frequency"] = (metadata.loc[i,"SizeT"])/lastT

            for c,t in list(zip(metadata.columns,["str"]+["int"]*5+["float","str","str","float"])):
                metadata[c] = metadata[c].astype(t)
            metadata["Start time"] = pd.to_datetime(metadata["Start time"])
            metadata["Duration"] = pd.to_timedelta(metadata["SizeT"]/metadata["Frequency"], unit="s")
            metadata["End time"] = metadata["Start time"]+metadata["Duration"]
            self.metadata = metadata
            try:
                self.tag_linescans()
            except:
                self.metadata["line scan"]="none"
        
    def calc_gaps(self):
        metadata = self.metadata
        x = [(metadata["Start time"].iloc[i+1]-metadata["End time"].iloc[i]).total_seconds() for i in range(len(metadata)-1)]
#         x = [el if el>=pd.Timedelta(0) else pd.Timedelta(0) for el in x]
        metadata["gap"] = [np.nan]+x
    
    def tag_linescans(rec, min_size=6,verbose=False,indices=None):
        metadata = rec.metadata
        if "line scan" not in metadata.index:
            metadata["line scan"]="none"
        if indices is None:
            indices = metadata.index
        for i in indices:
            if metadata.loc[i,"SizeT"]==1: continue
            freq = rec.metadata.loc[i,"Frequency"]
            if not np.isfinite(freq): continue
            scale = max(1,int(2*min_size/rec.metadata.loc[i,"pxSize"]))
            if metadata.loc[i,"SizeY"]<=scale:
                rec.metadata.loc[i,"line scan"]="multi"
            else:
                data = rec.import_series(metadata.loc[i,"Name"], restrict=(0,1.5/freq),save=False)[1]
                assert len(data)==1
                data = data[0]
#                 if verbose: print (data.shape)
                try:
                    stat = autocorr2d(data, dxrange=[scale], dyrange=[0])[0, 0] / np.abs(
                        autocorr2d(data, dyrange=[scale], dxrange=[0])[0, 0])
                    if stat>2:
                        rec.metadata.loc[i,"line scan"]="single"
                except:
                    if verbose: print("autocorr2d failed")
        
    def save_metadata(self):
        from os import umask
        umask(0o002)
        self.metadata.to_csv(self.metafile, float_format = "%#.3g")

    def import_series(self, Series, onlyMeta=False, isLineScan=False, restrict=None, mode="auto", save=True,
                      pathToTif=None, channel=0):
        # print ("+"*20,Series, self.metadata.Name)

        if Series == "all":
            SeriesList = self.allSeries
        elif Series in self.metadata.Name.values:
            SeriesList = [Series]
        elif "Series" in Series[0]:
            SeriesList = list(Series)
            serrange = [int(el.replace("Series", "")) for el in SeriesList]
            if len(serrange) > 1:
                Series = "Series%03i-%i" % (serrange[0], serrange[-1])
            else:
                Series = "Series%03i" % (serrange[0])
        else:
            try:
                serrange = Series.split("Series")[1].split("-")
                serrange = [int(el) for el in serrange]
                singleFile = len(serrange) == 1
                if singleFile:
                    SeriesList = [Series]
                else:
                    serrange = range(serrange[0], serrange[-1] + 1)
                    SeriesList = ["Series%03i" % i for i in serrange]
            except:
                warn(f"Cannot import {Series} from {self.path}. Probably a modified name.")
                return None
        try:
            self.Series
        except:
            self.Series = {}
        self.Series[Series] = {}

        metadata = self.metadata.copy()
        toDrop = [i for i, r in metadata.iterrows() if r.Name not in SeriesList]
        metadata.drop(index=toDrop, inplace=True)
        if len(metadata) == 0:
            return None
        if len(metadata) > 1:
            assert (np.abs(metadata["Frequency"] / metadata["Frequency"].mean() - 1) < 1e-1).all()

        for c in ["SizeX", "SizeY", "pxSize", "pxUnit", "bit depth"]:
            assert all([x == metadata[c].iloc[0] for x in metadata[c]])
        tsum = metadata[["Name", "SizeT", "Start time", "Duration", "End time"]].copy()
        metadata1 = metadata.iloc[0].copy()
        metadata1['SizeT'] = tsum["SizeT"].sum()
        if restrict is not None:
            t_begin, t_end = restrict
            metadata1["time_range"] = restrict
            time = np.arange(metadata1['SizeT']) / metadata1["Frequency"]
            if t_end < 0:
                t_end = time.max() + t_end
            if t_end < t_begin:
                t_end = time.max()
            frame_begin = np.where(time >= t_begin)[0][0]
            frame_end = np.where(time <= t_end)[0][-1]
            metadata1["frame_range"] = frame_begin, frame_end
            metadata1['SizeT'] = frame_end - frame_begin
        else:
            frame_begin, frame_end = 0, metadata1['SizeT']
            metadata1["frame_range"] = frame_begin, frame_end
        metadata1["individual Series"] = tsum
        metadata1["Name"] = Series

        if isLineScan:
            from copy import deepcopy
            metadata2 = deepcopy(metadata1)
            metadata2["Frequency"] = metadata2["Frequency"] * metadata2["SizeY"]
            metadata2["SizeT"] = metadata2["SizeT"] * metadata2["SizeY"]
            metadata2["SizeY"] = 1
            self.Series[Series]["metadata"] = metadata2
        else:
            metadata2 = metadata1

        if save:
            self.Series[Series]["metadata"] = metadata2

        if onlyMeta:
            if save:
                return None
            else:
                return metadata2
        ###################### if not onlyMeta ##################

        if pathToTif is not None:
            if mode != "auto":
                warn("If you supply a tiff file, the mode will be ignored.")
                mode = "auto"
            print("loading")
            # load motion corrected
            if pathToTif.lower().endswith("tif") or pathToTif.lower().endswith("tiff"):
                # from caiman import load as cload
                try:
                    data = cmovie(tiffmemmap(pathToTif))
                except:
                    from .utils import  load_tif
                    data = cmovie(load_tif(pathToTif))
            elif pathToTif.lower().endswith("npy"):
                data = cmovie(np.load(pathToTif))
            else:
                raise ValueError("Unknown file extension")
            data.fr = None
            FrameRange = metadata2.frame_range
            data = data[FrameRange[0]:FrameRange[1]]
        else:
            if mode == "auto":
                Nbytes = np.zeros(1, dtype="float32").nbytes
                Nbytes = Nbytes * metadata1.SizeT * metadata1.SizeY * metadata1.SizeX
                if Nbytes > 1e9:
                    mode = "tif"
                else:
                    mode = "ram"
            if mode in ["memmap","tif"]:
                if not hasattr(self, "tempdir"):
                    user = os.environ["USER"].split("jupyter-")[-1]
                    from datetime import date
                    today = date.today().strftime("%Y_%m_%d")
                    rndnumber = np.random.randint(int(1e10))
                    today_dir = f"/data/.tmp/{today}"
                    if not os.path.isdir(today_dir):
                        os.makedirs(today_dir)
                        #os.chown(today_dir, -1, "jupyterhub-users")
                        os.chmod(today_dir, 0o777)
                    tempdir = os.path.join(today_dir, f"{user}/{rndnumber}")
                    if not os.path.isdir(tempdir):
                        os.makedirs(tempdir)
                    self.tempdir = tempdir
                    # tempdir = "/data/.tmp"
                    # self.tempdir = os.path.join(tempdir, f"{np.random.randint(int(1e10))}")
                    # os.makedirs(self.tempdir)
                filename = os.path.join(self.tempdir,
                                        f"{Series.replace('/','_')}_d1_{metadata1.SizeX}_d2_{metadata1.SizeY}_d3_{metadata1.SizeZ}_order_C_frames_{metadata1.SizeT}_.{mode}")
                if mode=="memmap":
                    data = np.memmap(
                        filename,
                        dtype = "float32",
                        mode = "w+",
                        shape = (metadata1.SizeT, metadata1.SizeY, metadata1.SizeX)
                    )
                else: # it is tiff
                    data = tiffmemmap(
                        filename,
                        dtype = "float32",
                        mode = "w+",
                        photometric = "minisblack",
                        shape = (int(metadata1.SizeT), int(metadata1.SizeY), int(metadata1.SizeX))
                    )
            elif mode == "ram":
                data = np.zeros(
                    shape=(metadata1.SizeT, metadata1.SizeY, metadata1.SizeX),
                    dtype="float32",
                )
            else:
                raise ValueError("mode can only take the following values: 'ram', 'memmap', 'tif', and 'auto'.")
            if self.path.endswith(".stk"):
                stkFile = TiffFile(self.path)
                data = stkFile.pages[0].asarray()[frame_begin: frame_end]
            else:
                try:
                    self.rdr
                except:
                    self.rdr = bf.ImageReader(self.path, perform_init=True)
                assert metadata.SizeT.iloc[0] > frame_begin
                for i in metadata.index:
                    if metadata.index[0] == i:
                        dt = frame_begin
                    else:
                        dt = 0
                    offset = metadata.loc[:i - 1, "SizeT"].sum() - frame_begin
                    for t in range(dt, metadata.loc[i, "SizeT"]):
                        if t + offset >= len(data):
                            break
                        data[t + offset] = self.rdr.read(series=i, rescale=False, t=t, c=channel)
                if isLineScan:
                    data = data.reshape((np.prod(data.shape[:2]), 1, data.shape[-1]))
        if save:
            self.Series[Series]["data"] = data
        else:
            return metadata2, data
    def __repr__(self):
        display(self.metadata)
        return ""

    def raw_import(self, ser, mode="auto", save=True):
        try:
            self.Series
        except:
            self.Series = {}
        self.Series[ser] = {}
        if isinstance(ser,str):
            j = self.metadata.index[self.metadata.Name==ser][0]
        else:
            j = ser
        metadata1 = self.metadata.loc[j]
        if "Nchannels" not in metadata1.index:
            metadata1.Nchannels = 1

        if mode=="auto":
            Nbytes = np.zeros(1,dtype=metadata1["bit depth"]).nbytes
            Nbytes = Nbytes * metadata1.SizeT * metadata1.SizeY * metadata1.SizeX * metadata1.SizeZ * metadata1.Nchannels
            if Nbytes>1e9:
                mode="tif"
            else:
                mode="ram"
        shape = (metadata1.SizeX, metadata1.SizeY)
        if metadata1.SizeZ > 1:
            shape = (metadata1.SizeZ, ) + shape
        if metadata1.SizeT > 1:
            shape = (metadata1.SizeT, ) + shape
        if metadata1.Nchannels > 1:
            shape = shape + (metadata1.Nchannels,)

        if mode == "ram":
            data = np.zeros(
                shape=shape,
                dtype=metadata1["bit depth"],
            )
        elif mode in ["memmap","tif"]:
            if not hasattr(self, "tempdir"):
                user = os.environ["USER"].split("jupyter-")[-1]
                from datetime import date
                today = date.today().strftime("%Y_%m_%d")
                rndnumber = np.random.randint(int(1e10))
                today_dir = f"/data/.tmp/{today}"
                if not os.path.isdir(today_dir):
                    os.makedirs(today_dir)
                    #os.chown(today_dir, -1, "jupyterhub-users")
                    os.chmod(today_dir, "777")
                tempdir = os.path.join(today_dir, f"{user}/{rndnumber}")
                if not os.path.isdir(tempdir):
                    os.makedirs(tempdir)
                self.tempdir = tempdir

            if mode=="memmap":
                filename = os.path.join(self.tempdir, f"{ser}_d1_{metadata1.SizeX}_d2_{metadata1.SizeY}_d3_{metadata1.SizeZ}_order_C_frames_{metadata1.SizeT}_.mmap")
                data =np.memmap(filename,
                                dtype=metadata1["bit depth"],
                                mode="w+",
                                shape=shape
                                )
            if mode=="tif":
                filename = os.path.join(self.tempdir, "tmp.tif")
                shape = tuple(int(x) for x in shape)
                data =tiffmemmap(filename,
                                 shape=shape,
                                 dtype=metadata1["bit depth"],
                                 photometric="minisblack"
                                 )
        else:
            raise ValueError("mode can only take the following values: 'ram', 'memmap', 'tif', and 'auto'.")

        try:
            self.rdr
        except:
            self.rdr = bf.ImageReader(self.path, perform_init=True)

        from itertools import product
        for c, t, z in product(range(metadata1.Nchannels), range(metadata1.SizeT), range(metadata1.SizeZ)):
            slicer = slice(None), slice(None)
            if metadata1.SizeZ > 1:
                slicer = (z, ) + slicer
            if metadata1.SizeT > 1:
                slicer = (t, ) + slicer
            if metadata1.Nchannels > 1:
                slicer = slicer + (c,)
            data[slicer] = self.rdr.read(series=j, rescale=False, t=t, c=c,z=z)

        if save:
            self.Series[ser]["data"] = data
        else:
            return data

