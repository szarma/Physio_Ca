#!/opt/tljh/user/envs/physio/bin/python

import argparse
            
parser = argparse.ArgumentParser(description='Correct series for translational motion')
parser.add_argument('--recording', '-rec', type=str, default="", help='path to the recording')
parser.add_argument(   '--series', '-ser', type=str, default="", help='name of the series')
parser.add_argument( '--restrict', type=str, default="", help='restrict analysis to the time interval (in seconds!), e.g. "0_100" will only process first 100 seconds of the movie')
parser.add_argument(  '--verbose', const=True, default=False, action="store_const",
                    help='toggle verbose output')
# parser.add_argument('--test', const=True, default=False, action="store_const",
#                     help='toggle test mode on')
parser.add_argument('--debug', const=True, default=False, action="store_const",
                    help='toggle debug mode')

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
import matplotlib.pyplot as plt
import PIL

from islets import cmovie
from islets.numeric import rebin
from islets.Recording import Recording, saveMovie

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
# metadata   = rec.Series[ser]["metadata"]
freq = metadata.Frequency
orig_movie = cmovie(rec.Series[ser]["data"][:-1], fr=freq)

if args.debug:
    assert False ########################## debug stop ###########
    
freqMC = 1
n_rebin = int(np.round(freq/freqMC))
if n_rebin<1:
    n_rebin=1
reb_movie = rebin(orig_movie,n_rebin,dtype='float32')
reb_movie = cmovie(reb_movie, fr=freq/n_rebin)
cell_half_width_in_px = int(np.ceil(6/metadata.pxSize))
max_dev_rig = cell_half_width_in_px//2+1

figs, axs = plt.subplots(2, 1, figsize=(8, 8), sharex=True, sharey=True)
shifts = np.zeros((len(reb_movie), 2))
template = reb_movie[len(reb_movie) // 2]
while True:
    tmp = reb_movie.motion_correct(
        max_shift_w=max_dev_rig,
        max_shift_h=max_dev_rig,
        template=template
    )
    dshifts = np.array(tmp[1])
    template = template
    shifts += dshifts
    maxshift = np.abs(dshifts).max()

    axs[0].plot(np.arange(len(shifts)) / reb_movie.fr, shifts[:, 0])
    axs[1].plot(np.arange(len(shifts)) / reb_movie.fr, shifts[:, 1])
    if maxshift < max_dev_rig:
        break
plt.xlabel("time index")
axs[0].set_title("shifts (in pixels)")
axs[0].set_ylabel("vertical")
axs[1].set_ylabel("horizontal")
plt.savefig(os.path.join(saveDir, "shifts.png"))

if n_rebin>1:
    shifts = cmovie(shifts.reshape((1,)+shifts.shape)).resize(n_rebin,1,1)[0]

m_rshifted = np.zeros_like(orig_movie[:len(shifts)])

di = len(shifts)//20

mtype = metadata["bit depth"]
maxv = 2**int(str(mtype).strip("uint"))-1

for i in range(0,len(shifts),di):
    sl = slice(i,min(len(shifts),i+di))
    print(sl.start, sl.stop)
    tmpm = orig_movie[sl].apply_shifts(shifts[sl])
    tmpm[tmpm<0] = 0
    tmpm[tmpm>maxv] = maxv
    tmpm = np.round(tmpm)
    m_rshifted[sl] = tmpm

del orig_movie, rec

if len(rec.metadata)==1:
    movieFilename = os.path.join(saveDir, ".".join(rec.Experiment.split(".")[:-1]+["mp4"]))
else:
    movieFilename = os.path.join(saveDir, rec.Experiment+"_"+ser+".mp4")

movieFilename = movieFilename.replace(".mp4", "_rigshifted.tif")

im = PIL.Image.fromarray(m_rshifted[0])

im.save(movieFilename,
        save_all=True,
        append_images=[PIL.Image.fromarray(m_rshifted[i]) for i in range(1,len(m_rshifted))],
        compression="tiff_deflate"
       )

saveMovie(m_rshifted, movieFilename.replace(".tif",".mp4"), )
bf.javabridge.kill_vm()
quit()
