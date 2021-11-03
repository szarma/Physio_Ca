#!/opt/tljh/user/envs/physio/bin/python
import argparse
parser = argparse.ArgumentParser(description='Extract series and process from a recording')
parser.add_argument('--recording', '-rec', type=str,
                    help='path to the recording')
parser.add_argument('--series', '-ser', type=str,
                    help='name of the series', default="")
parser.add_argument('--restrict', type=str,
                    help='restrict analysis to the time interval (in seconds!), e.g. "0-100" will only process first 100 seconds of the movie', default="")
parser.add_argument('--shifted', type=str,
                    help='path to a motion-corrected movie, bypasses the original', default="")
parser.add_argument('--leave-movie', const=True, default=False,action="store_const",
                    help='if the movie exists, do not attempt to overwrite it')
parser.add_argument('--verbose', const=True, default=False,action="store_const",
                    help='toggle verbose output')
parser.add_argument('--leave-pickles', const=True, default=False,action="store_const",
                    help='if the pickles exist, do not attempt to overwrite them')
parser.add_argument('--test', const=True, default=False,action="store_const",
                    help='toggle test mode on')
parser.add_argument('--only-movie', const=True, default=False,action="store_const",
                    help='only do movie')

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

    
from pandas import DataFrame
from sys import path as syspath
syspath.append(os.path.expanduser("~/srdjan_functs/"))

from islets.Recording import Recording
from islets import saveMovie
from islets.numeric import rebin
from islets.utils import show_movie, saveRois, get_filterSizes

from islets.Regions1 import Regions, getPeak2BoundaryDF, getGraph_of_ROIs_to_Merge, mergeBasedOnGraph

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category= FutureWarning,)
    from caiman import movie as cmovie
    from caiman import load as cload

import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects

fracSaturTh = .05
cmap = plt.cm.Greys
cmap.set_bad("lime")

recFile = args.recording#"/data/Sandra/2020/2020_08_05/Experiment63d.lif"
ser = args.series#"Series004"
restrict = tuple([int(t) for t in args.restrict.split("-")]) if len(args.restrict) else None
pathToShifted = args.shifted


if len(pathToShifted)==0:
    import javabridge
    from bioformats import JARS as bfJARS
    log_config = os.path.join(os.path.split(__file__)[0], "resources", "log4j.properties")

    javabridge.start_vm(
        class_path=bfJARS,
        max_heap_size="20G",
#         args=["-Dlog4j.configuration=file:{}".format(log_config),],
                       )


if args.verbose:
    print("importing series...")
    
rec = Recording(recFile)
isNikon = recFile.endswith(".nd2")

if isNikon:
    if len(ser):
        warnings.warn("Nikon files (.nd2) can apparently contain only one series, so the series argument is ignored")
    serToImport = rec.metadata.loc[0,"Name"]
    ser = os.path.split(recFile)[1].split(".nd2")[0]
else:
    serToImport = ser

    
saveDir = os.path.join(rec.folder, rec.Experiment+"_analysis", ser)    

if restrict is not None:
    saveDir += f"_{restrict[0]}-{restrict[1]}s"
else:
    restrict = (0,-2)

rec.import_series(serToImport, restrict=restrict, onlyMeta=len(pathToShifted))
metadata = rec.Series[serToImport]['metadata']


if len(pathToShifted)==0: javabridge.kill_vm()

if not os.path.isdir(saveDir):
    if args.verbose:
        print ("creating directory", saveDir)
    os.makedirs(saveDir)
else:
    if args.verbose:
        print (f"{saveDir} exists already.")
        

t0, te = restrict
if len(pathToShifted):
    movie = cload(
        pathToShifted,
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

# if args.leave_movie:
#     if args.verbose: print("Not even checking whether the movie exists")
# else:
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

if args.only_movie: 
    quit()
    
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


quit()
