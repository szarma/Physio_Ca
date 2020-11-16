#!/opt/tljh/user/envs/physio/bin/python
import argparse

from contextlib import contextmanager

@contextmanager
def suppress_stdout():
    import os, sys
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout
            
@contextmanager
def suppress_stderr():
    import os, sys
    with open(os.devnull, "w") as devnull:
        old_stderr = sys.stderr
        sys.stderr = devnull
        try:  
            yield
        finally:
            sys.stderr = old_stderr
            
parser = argparse.ArgumentParser(description='Extract series and process from a recording')
parser.add_argument('--recording', '-rec', type=str,
                    help='path to the recording')
parser.add_argument('--series', '-ser', type=str,
                    help='name of the series', default="")
parser.add_argument('--restrict', type=str,
                    help='restrict analysis to the time interval (in seconds!), e.g. "0-100" will only process first 100 seconds of the movie', default="")
parser.add_argument('--corrected', type=str,
                    help='path to a motion-corrected movie, bypasses the original', default="")
parser.add_argument('--leave-movie', const=True, default=False,action="store_const",
                    help='if the movie exists, do not attempt to overwrite it')
parser.add_argument('--verbose', const=True, default=False,action="store_const",
                    help='toggle verbose output')
parser.add_argument('--leave-pickles', const=True, default=False,action="store_const",
                    help='if the pickles exist, do not attempt to overwrite them')
parser.add_argument('--test', const=True, default=False,action="store_const",
                    help='toggle test mode on')
parser.add_argument('--only-movie', const=True, default=False, action="store_const",
                    help='only do movie')
parser.add_argument('--line-scan', default="none", type=str,
                    help='indicate if it is a line scan, and if yes, what kind ("single" or "multi")')
parser.add_argument('--debug', const=True, default=False, action="store_const",
                    help='toggle debug mode')

args = parser.parse_args()


if args.test:
    for k in args.__dict__.keys():
        print("%20s"%k, args.__dict__[k])
    
if args.verbose:
    print("importing modules...")
    
import os
import warnings
import numpy as np
import pickle
from sys import exc_info
    
from pandas import DataFrame
from sys import path as syspath
syspath.append(os.path.expanduser("~/srdjan_functs/"))

from islets.Recording import Recording, saveMovie, parse_leica
from islets.LineScan import LineScan
from islets.numeric import rebin
from islets.utils import saveRois, get_filterSizes
from islets.Regions1 import Regions, getPeak2BoundaryDF, getGraph_of_ROIs_to_Merge, mergeBasedOnGraph

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category= FutureWarning,)
    from caiman import movie as cmovie
    from caiman import load as cload

import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from copy import copy
fracSaturTh = .05
cmap = copy(plt.cm.Greys)
cmap.set_bad("lime")

pathToCorrected = args.corrected
recFile = args.recording
ser = args.series




if len(pathToCorrected)==0:
    import javabridge
    import bioformats as bf
#     log_config = os.path.join(os.path.split(__file__)[0], "resources", "log4j.properties")

    javabridge.start_vm(
        class_path=bf.JARS,
        max_heap_size="20G",
#         args=["-Dlog4j.configuration=file:{}".format(log_config),],
    )
#     bf.init_logger()


rec = Recording(recFile)
isNikon = recFile.endswith(".nd2")

restrict = tuple([int(t) for t in args.restrict.split("_")]) if len(args.restrict) else None
if args.verbose:
    print("importing series...")


if isNikon:
    if len(ser):
        warnings.warn("Nikon files (.nd2) can apparently contain only one series, so the series argument is ignored")
    serToImport = rec.metadata.loc[0,"Name"]
    ser = os.path.split(recFile)[1].split(".nd2")[0]
else:
    serToImport = ser


saveDir = os.path.join(rec.folder, rec.Experiment+"_analysis", ser)    

if args.line_scan!="none": 
    restrict = None
else:
    if restrict is not None:
        saveDir += f"_{restrict[0]}-{restrict[1]}s"
    else:
        restrict = (0,-2)


if not os.path.isdir(saveDir):
    if args.verbose:
        print ("creating directory", saveDir)
    os.makedirs(saveDir)
else:
    if args.verbose:
        print (f"{saveDir} exists already.")

if args.debug:
    assert False ########################## debug stop ###########

if args.line_scan!="none":
    nameDict = dict([(name, list(ii)) for ii, name in parse_leica(rec, index=True)])
    indices = nameDict[serToImport]
    serNames = rec.metadata.loc[indices,"Name"]
    for ix, name in zip(indices,serNames):
        lsname = "%s: %s"%(rec.Experiment[:-4], name)
        rec.import_series(name, isLineScan=(args.line_scan=="single") )
        data = rec.Series[name]["data"].astype("float32")
        if args.line_scan=="multi":
            data = data.sum(1)
        else:
            assert data.shape[1]==1
            data = data[:,0]
        linescan = LineScan(
            data = data.T,
            metadata = rec.Series[name]["metadata"],
            name = lsname
            )
        linescan.plot(save=os.path.join(saveDir,lsname.replace(": ","_")+".png"), Npoints=2000)   
        
        
        
else:

    rec.import_series(serToImport, restrict=restrict, 
                      onlyMeta= bool(len(pathToCorrected))
                     )
    metadata = rec.Series[serToImport]['metadata']

    try: javabridge.kill_vm()
    except: pass
    
    if args.debug: assert False
    ########################## debug stop ###########

    t0, te = restrict
    if len(pathToCorrected):
        movie = cload(
            pathToCorrected,
            fr=metadata.Frequency,
            outtype=metadata['bit depth']
        )
        time = np.arange(len(movie))/movie.fr
        if te<0:
            te = time[-1]+te
        FrameRange = (np.searchsorted(time, t0), np.searchsorted(time, te))
        movie = movie[FrameRange[0]:FrameRange[1]]
    else:
        movie = cmovie(
            rec.Series[serToImport]['data'],
            fr=metadata.Frequency
        )
    if args.test:
        movie = movie[:,:20,:20]
    #### movie saving (or not)
    if len(rec.metadata)==1:
        movieFilename = os.path.join(saveDir, ".".join(rec.Experiment.split(".")[:-1]+["mp4"]))
    else:
        movieFilename = os.path.join(saveDir, rec.Experiment+"_"+ser+".mp4")


    if metadata.pxSize<.8:
        if args.verbose: print ("Resizing the movie resolution by 2...")
        movie = movie.resize(1/2,1/2,1)
        metadata.pxSize *= 2
        metadata.SizeX /= 2
        metadata.SizeY /= 2

    writeMovie = True
    if os.path.isfile(movieFilename):
        if args.verbose: print("Movie already exists, ", end="")
        if args.leave_movie:
            writeMovie = False
            if args.verbose: print("and I leave it be.")
        else:
            if args.verbose: print("and I'll rewrite it.")
    if writeMovie:
        if args.verbose: print("Writing the movie...")
        if not args.test: saveMovie(movie,movieFilename)

    if not args.only_movie: 


        #### protocol filename
        protocolFilename = movieFilename.replace(".mp4", "_protocol.txt")
        if not os.path.isfile(protocolFilename):
            if args.verbose: print("placed dummy protocol file at", protocolFilename)
            if not args.test:
                DataFrame([[None]*4],columns=["compound","concentration","begin","end"]).to_csv(protocolFilename,index=False)

        filtSizes = get_filterSizes(metadata.pxSize)

        # anull saturated above threshold
        Nsatur = (movie==movie.max()).sum(0)
        toAnull = np.where(Nsatur>len(movie)*fracSaturTh)
        movie[(slice(None), )+toAnull] = 0


        for spFilt in filtSizes:
            if args.verbose: print ("\t"*2,"#"*5,spFilt)

            pickleFile = os.path.join(saveDir, ".".join(map(str,spFilt))+"_rois.pkl")
            if os.path.isfile(pickleFile) and args.leave_pickles:
                if args.verbose: print ("already exists, skipping.")
                continue
            else:
                if args.verbose: print ("processing with filter size of ", spFilt)

            regions = Regions(movie,gSig_filt=spFilt,diag=True)
            regions.time += t0
            if args.verbose:
                print (f"initiallized with {len(regions.df)} rois.")

            regions.purge_lones((min(spFilt)*.4)**2, verbose=args.verbose)
            regions.sortFromCenter()
            regions.calcTraces()
            regions.metadata = metadata
            if not args.test: 
                saveRois(regions, saveDir, filename= ".".join(map(str,spFilt)), add_date=False, formats=["vienna"])

            # preview image
            fig = plt.figure(figsize=(5,5*np.divide(*movie.shape[1:])))
            ax = fig.add_axes([0.01,0.01,.98,.98])
            regions.plotEdges(imkw_args={"cmap":cmap},color="darkred", ax=ax)
            text = ax.text(.98,.02,len(regions.df),size=16,transform = ax.transAxes, ha="right",color="goldenrod")
            text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'),
                                   path_effects.Normal()])
            if not args.test:
                fig.savefig(os.path.join(saveDir, ".image_"+".".join(map(str,spFilt))+'.png'),dpi=150)
            plt.close(fig)
            del regions



# except:
#     print (exc_info())
# finally:
try: javabridge.kill_vm()
except: pass
exit()