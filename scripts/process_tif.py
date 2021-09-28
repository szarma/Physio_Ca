#!/opt/tljh/user/envs/physio/bin/python
import argparse
parser = argparse.ArgumentParser(description='Extract series and process from a recording')
parser.add_argument('--recording'     , '-rec' , type=str,   help='path to the recording')
parser.add_argument('--frequency'     , '-freq', type=float, help='recording frequency',)
parser.add_argument('--spatial-filter', '-sp'  , type=str, 
                    help='''produce roi pickles with exactly these filter sizes,
                    e.g. -sp="5" or -sp="5+6" to produce simple rois with indicated sizes,
                    or sp="5,5+6" to produce both 5 and 5+6. 
                    Filter size determines the size of final rois.
                    If you specify rebinning, the numbers here refer to the movie upon rebinning.''',
                    default="infer"
                    )

parser.add_argument('--rebin','-rb',type=int,
                    help= '''rebin recording with large resolution by a number (e.g. 2,3,5).
                    Default is not to rebin.''', 
                    default=1
                   )
# parser.add_argument('--skip', '-skip', type=str,
#                     help='skip every kth frame if int else infer', default="auto")
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
# import pickle
# from PIL import Image
from sys import exit
    
from pandas import DataFrame

from islets.Recording import saveMovie, autocorr2d
from islets.utils import saveRois, getStatImages
from islets.Regions import Regions

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category= FutureWarning,)
    from islets import cload

import matplotlib.pyplot as plt
from copy import copy
fracSaturTh = .05
cmap = copy(plt.cm.Greys)
cmap.set_bad("lime")

if args.verbose:
    print("importing recording...")
tiffile = args.recording
if not os.path.isfile(tiffile):
    raise ValueError(f"{tiffile} cannot be found")
freq = float(args.frequency)

# restrict = tuple([int(t) for t in args.restrict.split("_")]) if len(args.restrict) else None
movie = cload(tiffile, fr=freq)

if args.test:
    movie = movie[:,200:300,200:400]
    print ("reducing the movie for testing purposes.")
    
if args.rebin>1:
    movie = movie.resize(1/args.rebin, 1/args.rebin, 1)
    print (f"resized the movie as requested. The new resolution is {movie.shape[1:]}.")
    
saveDir = tiffile+"_analysis/all"

if not os.path.isdir(saveDir):
    if args.verbose:
        print ("creating directory", saveDir)
    os.makedirs(saveDir)
else:
    if args.verbose:
        print (f"{saveDir} exists already.")

# with Image.open(recFile) as im:
#     im.seek(0)
#     outtype = np.array(im).dtype
# outtype = "float32"

# movie = cload(
#     recFile,
#     outtype=outtype
# )

# if args.skip == "auto":
#     allstd = movie.std(0)
#     for skip in [3,2,1]:
#         meanSecond = np.std([movie[j::skip].mean(0) for j in range(skip)],0)
#         if np.mean(meanSecond/allstd)>.9:
#             if args.verbose:
#                 print("inferred skip =",skip)
    #             break
# else:
#     skip = int(args.skip)
# movie = movie[::skip]
# movie.fr = Frequency
# time = np.arange(len(movie))/movie.fr
# if te<0:
#     te = time[-1]+te
# FrameRange = (np.searchsorted(time, t0), np.searchsorted(time, te))
# movie = movie[FrameRange[0]:FrameRange[1]]




    
movieFilename = os.path.join(saveDir,os.path.splitext(os.path.split(tiffile)[1])[0]+".mp4")

        

writeMovie = True
if os.path.isfile(movieFilename):
    if args.verbose: print("Movie already exists, ", end="")
    if args.leave_movie:
        writeMovie = False
        if args.verbose: print("and I leave it be.")
    else:
        if args.verbose: print("and I'll rewrite it.")
else:
    if args.verbose: print("Movie does not exist. ", end="")
if writeMovie:
    if args.verbose: print("Writing the movie...")
    if not args.test: saveMovie(movie,movieFilename)

if args.only_movie: 
    exit()
    
if args.verbose: print("Calculating movie statistics...")
statImages = getStatImages(movie[::])

if args.spatial_filter=="infer":
    tpcx = autocorr2d(statImages["highperc"],np.arange(0,20),[0]).flatten()
    tpcy = autocorr2d(statImages["highperc"],[0],np.arange(0,20)).flatten()
    x = (tpcx+tpcy)/2
    cellSizePx = np.searchsorted(-x,-.1,)
    cellSizePx = max(cellSizePx,3)
    if args.verbose:
        print (f"Assuming the cells are approx {cellSizePx} pixels in radius, and will use this as a spatial filter.")
        
    # nRebSpace = cellSizePx//6
    # if nRebSpace>1:
    #     if args.verbose:
    #         print ("Resizing the movie resolution by", nRebSpace)
    #     movie = movie.resize(1/nRebSpace,1/nRebSpace,1)
    #     cellSizePx /= nRebSpace
    filtSizes = [[cellSizePx]]
else:
    filtSizes = args.spatial_filter.split(",")
    filtSizes = [eval(el.replace("+",",")) if "+" in el else (int(el),) for el in filtSizes]
    

    
#### protocol filename
protocolFilename = movieFilename.replace(".mp4", "_protocol.txt")
if not os.path.isfile(protocolFilename):
    if args.verbose: print("placed dummy protocol file at", protocolFilename)
    if not args.test:
        DataFrame([[None]*4],columns=["compound","concentration","begin","end"]).to_csv(protocolFilename,index=False)

# if args.test:
#     raise ValueError("stopped because of testing")

for spFilt in filtSizes:
    if args.verbose: print ("\t"*2,"#"*5,spFilt)

    pickleFile = os.path.join(saveDir, ".".join(map(str,spFilt))+"_rois.pkl")
    if os.path.isfile(pickleFile) and args.leave_pickles:
        if args.verbose: print ("already exists, skipping.")
        continue
    else:
        if args.verbose: print ("processing with filter size of ", spFilt)

    regions = Regions(statImages,gSig_filt=spFilt, verbose=args.verbose)
    if args.verbose:
        print (f"initiallized with {len(regions.df)} rois.")
    regions.merge_closest(verbose=args.verbose)
    regions.sortInOrder()
    regions.calcTraces(movie)
    #regions.time += metadata.frame_range[0]/metadata.Frequency
    regions.infer_gain()
    regions.calc_interest()
    if not args.test: 
        saveRois(regions, saveDir, filename= ".".join(map(str,spFilt)), add_date=False, formats=["vienna"])

exit()