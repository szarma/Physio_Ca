import numpy as np
import os
import bioformats as bf
import pandas as pd

def parse_leica(rec,
                merge=True,
                skipTags = ["crop","split","frame", "every","half", "snapshot","-","proj","max"],
                index=False,
                verbose=False,
               ):
    from pandas import Timedelta
    toDrop = [i for i,row in rec.metadata.iterrows() if "Series" not in row["Name"] or row["SizeY"]*row["SizeT"]<2000 or row["SizeT"]/row["Frequency"]<2]
    for tag in skipTags:
        toDrop += [i for i,row in rec.metadata.iterrows() if tag in row["Name"].lower()]
    rec.metadata.drop(index=np.unique(toDrop), inplace=True)
    if len(rec.metadata)==0:
        print ("Can not parse any series.")
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
                [rec.metadata[c].diff().abs()>0 for c in ["pxSize", "SizeX", "SizeY"]],
            axis=0)
        sers = np.split(rec.metadata.Name.values, np.where(ff)[0])[1:]
        idxs = np.split(list(rec.metadata.index), np.where(ff)[0])[1:]
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

def saveMovie(movie, filename, maxFreq=2, frameRate=60,dpi=100):
    from .utils import show_movie
    from .numeric import rebin
    if maxFreq<movie.fr:
        nrebin = int(np.ceil((movie.fr/maxFreq)))
        if nrebin>1:
            showMovie = rebin(movie, nrebin)
            from caiman import movie as cmovie
            showMovie = cmovie(showMovie, fr=movie.fr/nrebin)
        else:
            showMovie = movie+1
    else:
        showMovie = movie+1
    if filename=="embed":
        return show_movie(showMovie, 
               fps=frameRate,
               NTimeFrames=len(showMovie),
               # out="save",
               # saveName=filename,
               tmax=len(showMovie)/showMovie.fr,
               dpi=dpi
              )
    else:
        show_movie(showMovie, 
               fps=frameRate,
               NTimeFrames=len(showMovie),
               out="save",
               saveName=filename,
               tmax=len(showMovie)/showMovie.fr,
               dpi=dpi
              )
        return 0

class Recording:
    def __init__(self, pathToExperiment):
        self.folder, self.Experiment = os.path.split(pathToExperiment)
        self.path = pathToExperiment
        self.metafile = os.path.join(self.folder,  f".{self.Experiment}.meta")
        try:
            self.metadata = pd.read_csv(self.metafile, index_col=0)
            for c in ["Start time","End time"]:
                self.metadata[c] = pd.to_datetime(self.metadata[c])
            for c in ["Duration"]:
                self.metadata[c] = pd.to_timedelta(self.metadata[c])
            from general_functions import td_nanfloor
            self.metadata["End time"] = self.metadata["End time"].apply(td_nanfloor)
            # print (f"metadata imported from {self.metafile}.")
        except:
            print (f"Recording {pathToExperiment} not yet preprocessed. Preprocessing takes a few seconds and will speed up the usage later...", end=" ")
            from sys import stdout
            stdout.flush()
            self.parse_metadata()
            self.save_metadata()
            print (f"Finished.")
            stdout.flush()
        self.nSeries = len(self.metadata)
        self.allSeries = self.metadata.Name.values
                
#     def __del__(self):
#         import javabridge
#         javabridge.kill_vm()
        
    def parse_metadata(self,verbose=False):
        md = bf.get_omexml_metadata(self.path)
        self.xml = bf.OMEXML(md)
        self.nSeries = self.xml.get_image_count()
        self.allSeries = [self.xml.image(i).Name for i in range(self.nSeries)]
#         if Series is None:
#             SeriesList = self.allSeries
#             Series = "all"
#         else:
#             serrange = Series.split("Series")[1].split("-")
#             serrange = [int(el) for el in serrange]
#             singleFile = len(serrange)==1
#             if singleFile:
#                 SeriesList = [Series]
#             else:
#                 serrange = range(serrange[0],serrange[-1]+1)
#                 SeriesList = ["Series%03i"%i for i in serrange]
#         self.Series = Series
        metadata = pd.DataFrame(columns=["Name","SizeT","SizeX","SizeY","SizeZ","pxSize","pxUnit",
                                         "bit depth", "Frequency", "Start time", "End time", "Duration"])
        for i in range(self.nSeries):
#             try:
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
                metadata.loc[i,"Start time"] = im.get_AcquisitionDate()
                if metadata.loc[i,"SizeT"]>10:
                    if verbose:
                        print (im.Name, metadata.loc[i,"SizeT"], type(metadata.loc[i,"SizeT"]))
                    lastT = im.Pixels.Plane(int(metadata.loc[i,"SizeT"]-1)).DeltaT
                    metadata.loc[i,"Frequency"] = (metadata.loc[i,"SizeT"]-1)/lastT
#                     metadata.loc[i,"Frequency"] = 1./np.diff([im.Pixels.Plane(i).DeltaT for i in range(min(metadata.loc[i,"SizeT"],100))]).mean()
                if metadata.loc[i,"SizeZ"]>1:
                    metadata.loc[i,"Z-stack height"] = im.Pixels.get_PhysicalSizeZ()
#             except:
#                 print (str(exc_info()))
#                 pass
        for c,t in list(zip(metadata.columns,["str"]+["int"]*4+["float","str","str","float"])):
            metadata[c] = metadata[c].astype(t)
        metadata["Start time"] = pd.to_datetime(metadata["Start time"])
        metadata["Duration"] = pd.to_timedelta(metadata["SizeT"]/metadata["Frequency"], unit="s")
        metadata["End time"] = metadata["Start time"]+metadata["Duration"]
        self.metadata = metadata
        
    def calc_gaps(self):
        metadata = self.metadata
        x = [(metadata["Start time"].iloc[i+1]-metadata["End time"].iloc[i]).total_seconds() for i in range(len(metadata)-1)]
#         x = [el if el>=pd.Timedelta(0) else pd.Timedelta(0) for el in x]
        metadata["gap"] = [np.nan]+x
        
    def save_metadata(self):
        self.metadata.to_csv(self.metafile, float_format = "%#.3g")
        
    def import_series(self,Series, onlyMeta=False,isLineScan=False,restrict=None,memmap=False):
        if Series=="all":
            SeriesList = self.allSeries
        elif Series in self.metadata.Name.values:
            SeriesList = [Series]
        elif "Series" in Series[0]:
            SeriesList = list(Series)
            serrange = [int(el.replace("Series","")) for el in SeriesList]
            if len(serrange)>1:
                Series = "Series%03i-%i"%(serrange[0],serrange[-1])
            else:
                Series = "Series%03i"%(serrange[0])
        else:
            serrange = Series.split("Series")[1].split("-")
            serrange = [int(el) for el in serrange]
            singleFile = len(serrange)==1
            if singleFile:
                SeriesList = [Series]
            else:
                serrange = range(serrange[0],serrange[-1]+1)
                SeriesList = ["Series%03i"%i for i in serrange]
        try:
            self.Series
        except:
            self.Series = {}
        self.Series[Series] = {}
        
        metadata = self.metadata.copy()
        toDrop = [i for i,r in metadata.iterrows() if r.Name not in SeriesList]
        metadata.drop(index=toDrop, inplace=True)
        if len(metadata)>1:
            assert (np.abs(metadata["Frequency"]/metadata["Frequency"].mean()-1)<1e-1).all()

        for c in ["SizeX","SizeY","pxSize","pxUnit","bit depth"]:
            assert all([x==metadata[c].iloc[0] for x in metadata[c]])
        
        tsum = metadata[["Name","SizeT","Start time","Duration","End time"]].copy()
        metadata1 = metadata.iloc[0].copy()
        metadata1['SizeT'] = tsum["SizeT"].sum()
        if restrict is not None:
            t_begin, t_end = restrict
            metadata1["time_range"] = restrict
            time = np.arange(metadata1['SizeT'])/metadata1["Frequency"]
            if t_end<0:
                t_end = time.max()+t_end
            if t_end<t_begin:
                t_end = time.max()
            
            frame_begin = np.where(time>=t_begin)[0][0]
            frame_end   = np.where(time<=t_end)[0][-1]
            metadata1["frame_range"] = frame_begin, frame_end
            metadata1['SizeT'] = frame_end-frame_begin
        else:
            frame_begin, frame_end = 0, metadata1['SizeT']
            metadata1["frame_range"] = frame_begin, frame_end
        metadata1["individual Series"] = tsum
        metadata1["Name"] = Series
        self.Series[Series]["metadata"] = metadata1
        
        if isLineScan:
            from copy import deepcopy
            metadata2 = deepcopy(metadata1)
            metadata2["Frequency"] = metadata2["Frequency"]*metadata2["SizeY"]
            metadata2["SizeT"] = metadata2["SizeT"]*metadata2["SizeY"]
            metadata2["SizeY"] = 1
            self.Series[Series]["metadata"] = metadata2
        
        if onlyMeta:
            return None
#         data = np.zeros((tsum["SizeT"].sum(), metadata1.SizeY, metadata1.SizeX), dtype=metadata1["bit depth"])
        data = np.zeros((metadata1.SizeT, metadata1.SizeY, metadata1.SizeX), dtype=metadata1["bit depth"])
        
        try: self.rdr
        except: self.rdr = bf.ImageReader(self.path, perform_init=True)
        
        assert metadata.SizeT.iloc[0]>frame_begin
        
        for i in metadata.index:
#             firstFrame = self.rdr.read(series=i, rescale=False, t=0)
#             if len(firstFrame.shape)==3:
            if metadata.index[0]==i:
                dt = frame_begin
            else:
                dt = 0
                
            offset = metadata.loc[:i-1,"SizeT"].sum()-frame_begin
            for t in range(dt,metadata.loc[i,"SizeT"]):
#                 print (t, t+offset)
                if t+offset>=len(data):break
                data[t+offset] = self.rdr.read(series=i, rescale=False, t=t, c=0)
        
        if isLineScan:
            data = data.reshape((np.prod(data.shape[:2]),1,data.shape[-1]))
#         if restrict is not None:
#             data = data[frame_begin:frame_end]
        self.Series[Series]["data"] = data
