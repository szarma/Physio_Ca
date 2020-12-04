#!/opt/tljh/user/envs/physio/bin/python

import argparse
            
parser = argparse.ArgumentParser(description='Extract series and process from a recording')
parser.add_argument('--recording', '-rec', type=str,
                    help='path to the recording')
parser.add_argument('--series', '-ser', type=str,
                    help='name of the series', default="")
parser.add_argument('--restrict', type=str, default="",
                    help='restrict analysis to the time interval (in seconds!), e.g. "0-100" will only process first 100 seconds of the movie')
parser.add_argument('--use-tif', type=str, default=None, 
                    help='path to a tif file to use instad of the original, useful when passing motion-corrected movie')
parser.add_argument('--leave-movie', const=True, default=False,action="store_const",
                    help='if the movie exists, do not attempt to overwrite it')
parser.add_argument('--verbose', const=True, default=False,action="store_const",
                    help='toggle verbose output')
# parser.add_argument('--mostly-blank', const=True, default=False,action="store_const",
#                     help='use if the FOV contains large areas without cells')
parser.add_argument('--leave-pickles', const=True, default=False,action="store_const",
                    help='if the pickles exist, do not attempt to overwrite them')
parser.add_argument('--only-movie', const=True, default=False, action="store_const",
                    help='only do movie')
parser.add_argument('--spatial-filter',"-sp", default=None, 
                    help='''produce roi pickles with exactly these filter sizes,
                    e.g. -sp="5" or -sp="5+6" to produce simple rois with indicated sizes,
                    or sp="5,5+6" to produce both 5 and 5+6. Default (None) will guess four
                    sizes based on pxSize in the metadata if there are.''')
parser.add_argument('--line-scan', default="none", type=str,
                    help='indicate if it is a line scan, and if yes, what kind ("single" or "multi")')
parser.add_argument('--debug', const=True, default=False, action="store_const",
                    help='toggle debug mode')

args = parser.parse_args()

# def process_as_linescans():
#     nameDict = dict([(name, list(ii)) for ii, name in parse_leica(rec, index=True)])
#     if serToImport in nameDict:
#         indices = nameDict[serToImport]
#     else:
#         serBegin, serEnd = [int(part.strip("Series"))  for part in serToImport.split("-")]
#         serEnd += 1
#         possibleNames = ["Series%03i"%jSer for jSer in range(serBegin, serEnd)]
#         indices = rec.metadata.Name.isin(possibleNames)
#         indices = indices[indices].index
#     serNames = rec.metadata.loc[indices,"Name"]
#     for ix, name in zip(indices,serNames):
#         lsname = "%s: %s"%(rec.Experiment[:-4], name)
#         rec.import_series(name, isLineScan=(args.line_scan=="single") )
#         data = rec.Series[name]["data"].astype("float32")
#         if args.line_scan=="multi":
#             data = data.sum(1)
#         else:
#             assert data.shape[1]==1
#             data = data[:,0]
#         linescan = LineScan(
#             data = data.T,
#             metadata = rec.Series[name]["metadata"],
#             name = lsname
#             )
#         linescan.plot(save=os.path.join(saveDir,lsname.replace(": ","_")+".png"), Npoints=2000)   

if args.debug:
    for k in args.__dict__.keys():
        print("%20s"%k, args.__dict__[k])
    
if args.verbose:
    print("importing modules...")
            
import os
import warnings
import numpy as np
np.corrcoef(*np.random.randn(2,3))
import pickle
from sys import exc_info, exit

from pandas import DataFrame
from sys import path as syspath
syspath.append(os.path.expanduser("~/srdjan_functs/"))

from islets.Recording1 import Recording, saveMovie
from islets.LineScan import LineScan
from islets.numeric import rebin
from islets.utils import saveRois, get_filterSizes
from islets.Regions1 import Regions, getStatImages

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category= FutureWarning,)
    from caiman import movie as cmovie
    from caiman import load as cload


pathToTif = args.use_tif
if pathToTif is not None:
    args.leave_movie=True
recFile = args.recording
ser = args.series

start_vm = (pathToTif is None) or args.line_scan!="none"

try:
    rec = Recording(recFile)
except:
    # probably means the file is not yet preprocessed
    # and we need VM
    start_vm = True

if start_vm:
    import bioformats as bf
    #log_config = os.path.join(
    #    os.path.split(__file__)[0], 
    #    "resources", 
    #    "log4j.properties"
    #)
    bf.javabridge.start_vm(
        class_path=bf.JARS,
        max_heap_size="20G",
        # args=["-Dlog4j.configuration=file:{}".format(log_config),],
    )
    # bf.init_logger()


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

if not os.path.isdir(saveDir):
    if args.verbose:
        print ("creating directory", saveDir)
    os.makedirs(saveDir)
else:
    if args.verbose:
        print (f"{saveDir} exists already.")

if args.line_scan!="none":
    process_as_linescans()
    bf.javabridge.kill_vm()
    if 'rec' in globals():
        del rec
    exit()
################### LINESCANS STOP HERE ################    
    

# if restrict is None:
#     restrict = (0,-2)
rec.import_series(serToImport,
                  restrict=restrict, 
                  pathToTif=pathToTif
                 )
metadata = rec.Series[serToImport]['metadata']
if start_vm:
    bf.javabridge.kill_vm()

movie = cmovie(
    rec.Series[serToImport]['data'][:-1],
)
movie.fr=metadata.Frequency

# if args.debug:
#     movie = movie[:,:50,:50]

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
    if not args.debug: saveMovie(movie,movieFilename)

if args.only_movie:
    if 'rec' in globals():
        del rec
    exit()


#### protocol filename
protocolFilename = movieFilename.replace(".mp4", "_protocol.txt")
if not os.path.isfile(protocolFilename):
    if args.verbose: print("placed dummy protocol file at", protocolFilename)
    if not args.debug:
        DataFrame([[None]*4],columns=["compound","concentration","begin","end"]).to_csv(protocolFilename,index=False)

if args.spatial_filter is None:
    filtSizes = get_filterSizes(metadata.pxSize)
else:
    filtSizes = args.spatial_filter.split(",")
    filtSizes = [eval(el.replace("+",",")) if "+" in el else (int(el),) for el in filtSizes]

if args.debug: assert False  #### debug stop ###

statistics = getStatImages(movie)

for spFilt in filtSizes:
    if args.verbose: print ("\t"*2,"#"*5,spFilt)

    pickleFile = os.path.join(saveDir, ".".join(map(str,spFilt))+"_rois.pkl")
    if os.path.isfile(pickleFile) and args.leave_pickles:
        if args.verbose: print ("already exists, skipping.")
        continue
    else:
        if args.verbose: print ("processing with filter size of ", spFilt)

    regions = Regions(statistics,gSig_filt=spFilt)
    if args.verbose:
        print (f"initiallized with {len(regions.df)} rois.")
    regions.merge_closest(verbose=args.verbose)
    regions.sortInOrder()
    regions.metadata = metadata
    regions.calcTraces(movie)
    regions.time += metadata.frame_range[0]/metadata.Frequency
    regions.infer_gain()
    regions.calc_interest()
    if not args.debug: 
        saveRois(regions, saveDir, filename= ".".join(map(str,spFilt)), add_date=False, formats=["vienna"])


if 'rec' in globals():
    del rec
exit()
