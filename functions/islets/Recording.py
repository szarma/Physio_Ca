import numpy as np
import os
import bioformats as bf
import pandas as pd

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
            # print (f"metadata imported from {self.metafile}.")
        except:
            print (f"Recording not yet preprocessed. Preprocessing takes a few seconds and will speed up the usage later...", end=" ")
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
        
    def parse_metadata(self,Series=None):
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
        metadata = pd.DataFrame(columns=["Name","SizeT","SizeX","SizeY","pxSize","pxUnit",
                                         "bit depth", "Frequency", "Start time", "End time", "Duration"])
        for i in range(self.nSeries):
            im = self.xml.image(i)
#             if im.Name not in SeriesList: continue
            for c in metadata.columns:
                try:
                    metadata.loc[i,c] = getattr(im.Pixels, c)
                except:
                    pass
            metadata.loc[i,"Name"] = im.Name
            metadata.loc[i,"pxSize"] = im.Pixels.get_PhysicalSizeX()
            metadata.loc[i,"pxUnit"] = im.Pixels.get_PhysicalSizeXUnit()
            metadata.loc[i,"bit depth"] = im.Pixels.get_PixelType()
            metadata.loc[i,"Start time"] = im.get_AcquisitionDate()
            if metadata.loc[i,"SizeT"]>1:
#             metadata.loc[i,"Frequency"] = 1/np.diff([im.Pixels.Plane(i).DeltaT for i in range(max(metadata.loc[i,"SizeT"],100))]).mean()
                metadata.loc[i,"Frequency"] = (metadata.loc[i,"SizeT"]-1)/im.Pixels.Plane(metadata.loc[i,"SizeT"]-1).DeltaT

        for c,t in list(zip(metadata.columns,["str"]+["int"]*3+["float","str","str","float"])):
            metadata[c] = metadata[c].astype(t)
        metadata["Start time"] = pd.to_datetime(metadata["Start time"])
        metadata["Duration"] = pd.to_timedelta(metadata["SizeT"]/metadata["Frequency"], unit="s")
        metadata["End time"] = metadata["Start time"]+metadata["Duration"]
        self.metadata = metadata
        
    def calc_gaps(self):
        metadata = self.metadata
        x = [metadata["Start time"].iloc[i+1]-metadata["End time"].iloc[i] for i in range(len(metadata)-1)]
        x = [el if el>=pd.Timedelta(0) else pd.Timedelta(0) for el in x]
        metadata["gap"] = [pd.NaT]+x
        
    def save_metadata(self):
        self.metadata.to_csv(self.metafile, float_format = "%#.3g")
        
    def import_series(self,Series):
        if Series=="all":
            SeriesList = self.allSeries
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
        
        assert (np.abs(metadata["Frequency"]/metadata["Frequency"].mean()-1)<1e-1).all()

        for c in ["SizeX","SizeY","pxSize","pxUnit","bit depth"]:
            assert all([x==metadata[c].iloc[0] for x in metadata[c]])
        
        tsum = metadata[["Name","SizeT","Start time","Duration","End time"]].copy()
        metadata1 = metadata.iloc[0].copy()
        metadata1["individual Series"] = tsum
        metadata1['SizeT'] = tsum["SizeT"].sum()
        metadata1["Name"] = Series
        self.Series[Series]["metadata"] = metadata1
        
        
        data = np.zeros((metadata1.SizeT, metadata1.SizeX, metadata1.SizeY), dtype=metadata1["bit depth"])
        
        try: self.rdr
        except: self.rdr = bf.ImageReader(self.path, perform_init=True)
        
        for i in metadata.index:
            offset = metadata.loc[:i-1,"SizeT"].sum()
            for t in range(metadata.loc[i,"SizeT"]):
                data[t+offset] = self.rdr.read(series=i, rescale=False, t=t)
        self.Series[Series]["data"] = data
