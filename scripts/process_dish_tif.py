#!/opt/tljh/user/envs/physio/bin/python
import argparse
parser = argparse.ArgumentParser(description='Extract series and process from a recording')
parser.add_argument('--recording', '-rec', type=str,
                    help='path to the recording')
parser.add_argument('--skip', '-skip', type=str,
                    help='skip every kth frame if int else infer', default="auto")
# parser.add_argument('--restrict', type=str,
#                     help='restrict analysis to the time interval (in seconds!), e.g. "0-100" will only process first 100 seconds of the movie', default="")
parser.add_argument('--leave-movie', const=True, default=False, action="store_const",
                    help='if the movie exists, do not attempt to overwrite it')
parser.add_argument('--verbose', const=True, default=False, action="store_const",
                    help='toggle verbose output')
parser.add_argument('--leave-pickles', const=True, default=False, action="store_const",
                    help='if the pickles exist, do not attempt to overwrite them')
parser.add_argument('--test', const=True, default=False, action="store_const",
                    help='toggle test mode on')
parser.add_argument('--only-movie', const=True, default=False, action="store_const",
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
from PIL import Image
from sys import exc_info
    
from pandas import DataFrame
from sys import path as syspath
syspath.append(os.path.expanduser("~/srdjan_functs/"))

from islets.Recording import Recording, saveMovie, autocorr2d
from islets.utils import saveRois, get_filterSizes
from islets.Regions import Regions, getStatImages

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category= FutureWarning,)
    from islets import cload

import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from copy import copy
fracSaturTh = .05
cmap = copy(plt.cm.Greys)
cmap.set_bad("lime")

recFile = args.recording
# restrict = tuple([int(t) for t in args.restrict.split("_")]) if len(args.restrict) else None
if args.verbose:
    print("importing series...")

saveDir = recFile+"_analysis/all"


if not os.path.isdir(saveDir):
    if args.verbose:
        print ("creating directory", saveDir)
    os.makedirs(saveDir)
else:
    if args.verbose:
        print (f"{saveDir} exists already.")

Frequency = float(os.path.split(recFile)[1].split("Hz")[0].strip("_").strip().split( "_")[-1].split()[-1])
# with Image.open(recFile) as im:
#     im.seek(0)
#     outtype = np.array(im).dtype
outtype = "float32"

movie = cload(
    recFile,
    outtype=outtype
)
if args.skip == "auto":
    allstd = movie.std(0)
    for skip in [3,2,1]:
        meanSecond = np.std([movie[j::skip].mean(0) for j in range(skip)],0)
        if np.mean(meanSecond/allstd)>.9:
            if args.verbose:
                print("inferred skip =",skip)
            break
else:
    skip = int(args.skip)
movie = movie[::skip]
movie.fr = Frequency
# time = np.arange(len(movie))/movie.fr
# if te<0:
#     te = time[-1]+te
# FrameRange = (np.searchsorted(time, t0), np.searchsorted(time, te))
# movie = movie[FrameRange[0]:FrameRange[1]]

statImages = getStatImages(movie[::])

tpcx = autocorr2d(statImages["std"],np.arange(0,20),[0]).flatten()
tpcy = autocorr2d(statImages["std"],[0],np.arange(0,20)).flatten()
x = (tpcx+tpcy)/2
cellSizePx = np.searchsorted(-x,-.1,)
if args.verbose:
    print (f"Assuming the cells are approx {cellSizePx} pixels in radius")

if args.test:
    movie = movie[:,:5*cellSizePx,:5*cellSizePx]
    
    
movieFilename = os.path.join(saveDir,os.path.splitext(os.path.split(recFile)[1])[0]+".mp4")

nRebSpace = cellSizePx//6
if nRebSpace>1:
    if args.verbose:
        print ("Resizing the movie resolution by", nRebSpace)
    movie = movie.resize(1/nRebSpace,1/nRebSpace,1)
    cellSizePx /= nRebSpace
        

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

filtSizes = get_filterSizes(8./cellSizePx)

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

    regions = Regions(movie,gSig_filt=spFilt,diag=True,use_restricted=True)
    regions.time #+= t0
    if args.verbose:
        print (f"initiallized with {len(regions.df)} rois.")

    regions.purge_lones((min(spFilt)*.4)**2, verbose=args.verbose)
    regions.sortFromCenter()
    regions.calcTraces()
    if "metadata" in globals():
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