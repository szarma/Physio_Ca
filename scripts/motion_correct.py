#!/opt/tljh/user/envs/physio/bin/python

import argparse

parser = argparse.ArgumentParser(description='Correct series for translational motion')
parser.add_argument('--recording', '-rec', type=str, default="", help='path to the recording')
parser.add_argument(   '--series', '-ser', type=str, default="", help='name of the series')
parser.add_argument( '--restrict', type=str, default="", help='restrict analysis to the time interval (in seconds!), e.g. "0_100" will only correct first 100 seconds of the movie')
parser.add_argument(  '--verbose', const=True, default=False, action="store_const",
                    help='toggle verbose output')
parser.add_argument('--debug', const=True, default=False, action="store_const",
                    help='toggle debug mode')
parser.add_argument('--correct-phase', const=True, default=False, action="store_const",
                    help='whether to correct the phase also, default is False. In case phase is to be corrected, note that it will take first 1s of the recording as a template onto which to match. If that is not desirable, more input is needed.')


args = parser.parse_args()


if args.debug:
    for k in args.__dict__.keys():
        print("%20s"%k, args.__dict__[k])

if args.verbose:
    print("importing modules...")

import numpy as np
np.corrcoef([0,1,2],[0,2,1])

import warnings
import os
from islets import cmovie, saveMovie
from islets.Recording import Recording
from tifffile import memmap as tifmemmap
from islets.utils import gentle_motion_correct
import bioformats as bf
bf.javabridge.start_vm(class_path=bf.JARS, max_heap_size="20G")

recFile = args.recording
rec = Recording(recFile)
ser = args.series
isNikon = recFile.endswith(".nd2")

restrict = tuple([int(t) for t in args.restrict.split("_")]) if len(args.restrict) else None

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

if not os.path.isdir(saveDir):
    if args.verbose:
        print ("creating directory", saveDir)
    os.makedirs(saveDir)
else:
    if args.verbose:
        print (f"{saveDir} exists already.")

if args.verbose:
    print("importing series...")

rec.import_series(serToImport, restrict=restrict, )
metadata = rec.Series[serToImport]['metadata']
freq = metadata.Frequency

if len(rec.metadata)==1:
    movieFilename = os.path.join(saveDir, ".".join(rec.Experiment.split(".")[:-1]+["mp4"]))
else:
    movieFilename = os.path.join(saveDir, rec.Experiment+"_"+ser+".mp4")

outMovieFilename =  movieFilename.replace(".mp4", "_rigshifted.mp4")
outTifFilename = outMovieFilename.replace(".mp4", ".tif")

orig_movie = cmovie(rec.Series[ser]["data"][:-1], fr=freq)
bf.javabridge.kill_vm()

if args.debug:
    assert False ########################## debug stop ###########

m_rshifted = tifmemmap(
    outTifFilename,
    shape=orig_movie.shape,
    dtype=orig_movie.dtype,
    photometric="minisblack"
)
cell_half_width_in_px = int(np.ceil(6 / metadata.pxSize))
max_dev_rig = cell_half_width_in_px // 2 + 1

if args.correct_phase:
    ####################
    even_template = orig_movie[:int(orig_movie.fr)+1,::2].mean(0)
    gentle_motion_correct(
        orig_movie[:,::2],
        m_rshifted[:,::2],
        plot_name=movieFilename.split(".mp4")[0]+"_even_shifts.png",
        template=even_template,
        max_dev_rig = (max_dev_rig, max_dev_rig)
    )
    odd_template = odd_template = (even_template[:-1]+even_template[1:])/2
    gentle_motion_correct(
        orig_movie[:,1::2],
        m_rshifted[:,1::2],
        plot_name=movieFilename.split(".mp4")[0]+"_odd_shifts.png",
        template=odd_template,
        max_dev_rig = (max_dev_rig, max_dev_rig)
    )
else:
    gentle_motion_correct(
        orig_movie,
        m_rshifted,
        plot_name=movieFilename.split(".mp4")[0]+"_shifts.png",
        max_dev_rig = (max_dev_rig, max_dev_rig)
    )

m_rshifted.flush()

del orig_movie, rec, m_rshifted

saveMovie(m_rshifted, outMovieFilename)

quit()
