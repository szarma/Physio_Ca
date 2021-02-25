#!/opt/tljh/user/envs/physio/bin/python
## importing stuff
import os
import pickle
from sys import path as syspath
fundir = os.path.expanduser("~/srdjan_functs/")
if os.path.isdir(fundir):
    syspath.append(fundir)
else:
    syspath.append("./functions/")
from general_functions import suppress_stdout, suppress_stderr

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from islets.Recording import Recording, saveMovie
from islets.Regions import Regions
from islets.numeric import power_spectrum, rebin
from islets.utils import saveRois
from islets.Regions import getPeak2BoundaryDF, getGraph_of_ROIs_to_Merge, mergeBasedOnGraph
from matplotlib.colors import LogNorm

def get_filterSizes(px, physSize=None):
    if physSize is None:
        try:
            physSize = args.__dict__["filter_size"]
        except:
            physSize = 5.5
    base = int(np.ceil(physSize/px))
    wider = int(np.ceil(base*1.1))
    if wider==base: wider += 1
    toComb = int(np.ceil(base*1.3))
    if toComb <= wider: toComb += 1
    return [(base,), (wider,), (base,wider), (base,toComb)]

def parse_leica(rec, merge=True, skipTags = ["crop","split","frame", "every","half", "snapshot","-"], index=False):
    from pandas import Timedelta
    toDrop = [i for i,row in rec.metadata.iterrows() if "Series" not in row["Name"] or row["SizeY"]*row["SizeT"]<2000]
    for tag in skipTags:
        toDrop += [i for i,row in rec.metadata.iterrows() if tag in row["Name"].lower()]
    rec.metadata.drop(index=np.unique(toDrop), inplace=True)
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
                [rec.metadata[c].diff().abs()>0 for c in ["pxSize", "SizeX", "SizeY"]],
            axis=0)
        sers = np.split(rec.metadata.Name.values, np.where(ff)[0])
        idxs = np.split(list(rec.metadata.index), np.where(ff)[0])
    else:
        sers = [rec.metadata.Name.values]
        idxs = [list(rec.metadata.index)]
    
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

fracSaturTh = .05
cmap = plt.cm.Greys
cmap.set_bad("lime")

if __name__=="__main__":
    import warnings
    warnings.simplefilter("ignore")
    import argparse
    parser = argparse.ArgumentParser(description='Extract series and process from a recording')
    parser.add_argument('--recording', '-rec', type=str,
                        help='path to the recording')
    parser.add_argument('--series', '-ser', type=str,
                        help='name of the series', default="")
    parser.add_argument('--filter-size', '-fs', type=float,
                        help='physical size of desired rois (in µm)', default=5.5)
    parser.add_argument('--interactive', '-i', const=True, default=False,action="store_const",
                        help='toggle interactive mode [the series argument is then ignored]')
    parser.add_argument('--leave-movie', const=True, default=False,action="store_const",
                        help='if the movie exists, do not attempt to overwrite it')
    parser.add_argument('--leave-pickles', const=True, default=False,action="store_const",
                        help='if the pickles exist, do not attempt to overwrite them')
    parser.add_argument('--test', const=True, default=False,action="store_const",
                        help='toggle test mode')

    args = parser.parse_args()

    # for k in args.__dict__.keys():
    #     print("%20s"%k, args.__dict__[k])

    if not args.test:
        import javabridge
        from bioformats import JARS as bfJARS
    #     log_config = os.path.join(os.path.split(__file__)[0], "resources", "log4j.properties")
        javabridge.start_vm(
    #         args=[
    #             "-Dlog4j.configuration=file:{}".format(log_config),
    #         ],
            class_path=bfJARS,
            max_heap_size="20G"
        )
        from islets import cmovie
    rec = Recording(args.recording)
    rec.calc_gaps()
    restrict=(0,-2)
    timeRestriction = ""
#     print (rec.metadata.gap)
    metadata = rec.metadata.copy()
    for c in metadata:
        if "time" in c:
            metadata[c] = metadata[c].dt.strftime('%H:%M:%S')
    metadata["Duration"] = pd.to_datetime(metadata["Duration"]).dt.strftime('%H:%M:%S')
#     metadata["gap"] = pd.to_timedelta(metadata["gap"],unit="s")#.strftime('%H:%M:%S')
    metadata["gap [s]"] = ["%.3g"%t for t in metadata["gap"]]
    del metadata["gap"]
    print ("The recording contains the following series/images:\n")
    from IPython.display import display
    display(metadata)

    if args.interactive:
        ser_idxs = input("Please input the index (e.q. 0) corresponding to the series you wish to be processed, or a range of indices if you wish the series to be stitched and processed together (e.g. 0-3).\t")
        if "-" in ser_idxs:
            indices = ser_idxs.split("-")
            assert len(indices)==2
            indices = metadata.index[slice(int(indices[0]), int(indices[1])+1)]
            serlist = metadata.Name.loc[indices[0]], metadata.Name.loc[indices[-1]]
        else:
            serlist = metadata.Name.loc[int(ser_idxs)]
        timeRestriction = input("Please input the time range (in seconds!) you wish to restrict your analysis to (e.g. 10-1800). In case you wish to analyse the whole series, you can also leave the fiels empty. \t")
        if len(timeRestriction):
            restrict = timeRestriction.split("-")
            restrict = tuple( float(t) for t in restrict)
        else:
            restrict = None
        rec.import_series(serlist, onlyMeta=True, restrict=restrict)
        for ser in rec.Series:
            break
    else:
        ser = args.series
    if ser=="all":
        if rec.Experiment.endswith("lif"):
            sers = parse_leica(rec, merge=True)
        else:
            sers = ["all"]
    else:
        sers = [ser]
    print (sers, "will be attempted to process")
    
    for ser in sers:
        print (f"Attempting to process {ser}...")

        try:
        #     with suppress_stdout():
        #     with suppress_stderr():
            j = rec.metadata[rec.metadata.Name==ser.split("-")[0]].index[0]
            rec.import_series(ser, isLineScan=rec.metadata.loc[j,"SizeY"]>550,restrict=restrict)
#             quitnow=False
        except:
            print (f"could not import {ser} from {rec.path}. exiting...")
            continue
#             quitnow=True
#         finally:
#             javabridge.kill_vm()
#         if quitnow:
#             quit()
        saveDir = os.path.join(rec.folder, rec.Experiment+"_analysis", ser)
        if len(timeRestriction):
            saveDir += "_"+timeRestriction+"s"
        if not os.path.isdir(saveDir):
            print ("creating directory", saveDir)
            os.makedirs(saveDir)
        else:
            print (f"{saveDir} exists already.")
        metadata = rec.Series[ser]["metadata"]
        movie = cmovie(rec.Series[ser]["data"],fr=metadata["Frequency"])
        if len(rec.metadata)==1:
            movieFilename = os.path.join(saveDir, ".".join(rec.Experiment.split(".")[:-1]+["mp4"]))
        else:
            movieFilename = os.path.join(saveDir, rec.Experiment+"_"+ser+".mp4")
        if len(rec.metadata)==1:
            protocolFilename = os.path.join(saveDir, ".".join(rec.Experiment.split(".")[:-1]+["_protocol","txt"]))
        else:
            protocolFilename = os.path.join(saveDir, rec.Experiment+"_"+ser+"_protocol.txt")

        if (not os.path.isfile(movieFilename)) or (not args.leave_movie):
            print("Writing the movie...")
            saveMovie(movie,movieFilename)
        if not os.path.isfile(protocolFilename):
            pd.DataFrame([[None]*4],columns=["compound","concentration","begin","end"]).to_csv(protocolFilename,index=False)
        # anull saturated above threshold
        Nsatur = (movie==movie.max()).sum(0)
        toAnull = np.where(Nsatur>len(movie)*fracSaturTh)
        movie[(slice(None), )+toAnull] = 0
        
        # go thorugh filter sizes
        regions = Regions(movie, full=False, diag=True, debleach=False, processes=5)
        filtSizes = get_filterSizes(metadata.pxSize)

        for spFilt in filtSizes:
            print ("\t"*2,"#"*5,spFilt)
            pickleFile = os.path.join(saveDir, ".".join(map(str,spFilt))+"_rois.pkl")
            if os.path.isfile(pickleFile) and args.leave_pickles:
                print ("already exists, skipping.")
                continue
            else:
                print ("processing with filter size of ", spFilt)
            regions.constructRois(gSig_filt=spFilt,img_th=0, processes=5)
            regions.update(movie)
            while True:
                peak2bnd = getPeak2BoundaryDF(regions.df)
                df = peak2bnd.query("dist<1.")[["i","j"]]
                if len(df)==0: break
                gRois = getGraph_of_ROIs_to_Merge(df,regions,plot=False)
                dropped = mergeBasedOnGraph(gRois,regions)
                if dropped == 0: break

            regions.purge_lones((min(spFilt)*.4)**2)
            regions.sortFromCenter()
            saveRois(regions, saveDir, filename= ".".join(map(str,spFilt)), add_date=False, formats=["vienna"])

            fig = plt.figure(figsize=(5,5*np.divide(*movie.shape[1:])))
            ax = fig.add_axes([0.01,0.01,.98,.98])
            regions.plotEdges(imkw_args={"cmap":cmap},color="darkred", ax = ax)
            fig.savefig(os.path.join(saveDir, ".image_"+".".join(map(str,spFilt))+'.png'),dpi=150)
            plt.close(fig)

    if not args.test:
        javabridge.kill_vm()
    quit()
