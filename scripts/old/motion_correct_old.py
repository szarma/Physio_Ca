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
            
parser = argparse.ArgumentParser(description='Correct series for translational motion')
parser.add_argument('--recording', '-rec', type=str, default="", help='path to the recording')
parser.add_argument(   '--series', '-ser', type=str, default="", help='name of the series')
parser.add_argument( '--restrict', type=str, default="", help='restrict analysis to the time interval (in seconds!), e.g. "0-100" will only process first 100 seconds of the movie')
parser.add_argument('--verbose', const=True, default=False, action="store_const",
                    help='toggle verbose output')
parser.add_argument('--test', const=True, default=False, action="store_const",
                    help='toggle test mode on')
parser.add_argument('--debug', const=True, default=False, action="store_const",
                    help='toggle debug mode')

args = parser.parse_args()


if args.test:
    for k in args.__dict__.keys():
        print("%20s"%k, args.__dict__[k])
    
if args.verbose:
    print("importing modules...")

    
import warnings
import cv2
try:
    cv2.setNumThreads(0)
    print ("NumThreads for cv2 ... set!")
except:
    pass

import numpy as np
import os
import PIL
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category= FutureWarning,)
    import caiman as cm
    from caiman.motion_correction import MotionCorrect
    from caiman.source_extraction.cnmf import params as params

from sys import path as syspath
syspath.append(os.path.expanduser("~/srdjan_functs/"))

from islets.numeric import rebin
from islets.Recording import Recording
from islets import saveMovie

import javabridge
from bioformats import JARS as bfJARS
javabridge.start_vm(class_path=bfJARS, max_heap_size="20G")

sessID = str(np.random.randint(1e9,1e10))
tmpfolder = f"/data/.tmp/{sessID}/"
os.makedirs(tmpfolder)

recFile = args.recording
rec = Recording(recFile)
ser = args.series

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
if restrict is not None:
    saveDir += f"_{restrict[0]}-{restrict[1]}s"


if not os.path.isdir(saveDir):
    if args.verbose:
        print ("creating directory", saveDir)
    os.makedirs(saveDir)
else:
    if args.verbose:
        print (f"{saveDir} exists already.")

# if args.debug:
#     assert False ########################## debug stop ###########
    
    

rec.import_series(serToImport, restrict=restrict, )
metadata = rec.Series[serToImport]['metadata']

#if restrist is not None:
#    t0, te = restrict


orig_images = rec.Series[ser]["data"]
metadata    = rec.Series[ser]["metadata"]

freq = metadata.Frequency
freqMC = 2
n_rebin = int(freq/freqMC)
if n_rebin>1:
    rebinned_images = rebin(orig_images,n_rebin,dtype='float32')
else:
    rebinned_images = orig_images.astype("float32")

#%% start a cluster for parallel processing (if a cluster already exists it will be closed and a new session will be opened)
if 'dview' in locals():
    cm.stop_server(dview=dview)
c, dview, n_processes = cm.cluster.setup_cluster(
    backend='local', n_processes=None, single_thread=False)

fname_rebinned = cm.save_memmap([rebinned_images], base_name=tmpfolder+"tmp",  order='C', border_to_0=0, dview=dview)    
fname_original = cm.save_memmap([orig_images],     base_name=tmpfolder+"orig", order='C', border_to_0=0, dview=dview)



# to free up memory:
del orig_images, rec.Series[ser]["data"]



cell_half_width_in_px = int(np.ceil(6/metadata.pxSize))
# motion correction parameters
opts = params.CNMFParams(params_dict={
    'fnames'              : fname_rebinned,
    "max_deviation_rigid" : cell_half_width_in_px//2+1,
    'border_nan'          : "copy",
    'pw_rigid'            : False,
    'gSig_filt'           : (cell_half_width_in_px, cell_half_width_in_px),  # size of high pass spatial filtering, used in 1p data,
    'nonneg_movie'        : True
}) 
mc = MotionCorrect(fname_rebinned, dview=dview, **opts.get_group('motion'))


# ### Shifts


while True:
    print ("constraining mc.max_deviation_rigid at ", mc.max_deviation_rigid)
    mc.motion_correct(save_movie=True);
    maxshift = np.abs(mc.shifts_rig).max()
    if maxshift<mc.max_deviation_rigid:
        break
    else:
        mc.max_deviation_rigid += 1


mc.shifts_rig = np.array(mc.shifts_rig)

if n_rebin>1:
    mc.shifts_rig = cm.movie(mc.shifts_rig.reshape((1,)+mc.shifts_rig.shape)).resize(n_rebin,1,1)[0]
m_rshifted_name = mc.apply_shifts_movie(fname_original, save_memmap=True, save_base_name=tmpfolder+"MC")

cm.stop_server(dview=dview)
c, dview, n_processes = cm.cluster.setup_cluster(
    backend='local', n_processes=None, single_thread=True)


Y, xy_res, time_res = cm.load_memmap(m_rshifted_name, mode="r+")
m_rshifted = cm.movie(Y.reshape(xy_res+(time_res,)).T, fr=freq)


mtype = metadata["bit depth"]
maxv = 2**int(str(mtype).strip("uint"))-1
m_rshifted[m_rshifted<0] = 0
if m_rshifted.max()>maxv:
    print (f"some pixels in the shifted movie are larger than {maxv}, that's strange. I'll put them to the max value.")
    m_rshifted[m_rshifted>maxv] = maxv
m_rshifted = np.round(m_rshifted)
m_rshifted = m_rshifted.astype(mtype)
cm.stop_server(dview=dview)


# ### Save
# exp = os.path.split(rec.path)[1].split(".")[0]
# if rec.path.endswith("nd2"):
#     shiftFilename = f"{rec.path}_analysis/{exp}_rigidshifted.tif"
# else:
#     shiftFilename = f"{rec.path}_analysis/{ser}/{exp}_{ser}_rigidshifted.tif"


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


os.system(f"rm -rf {tmpfolder}")
javabridge.kill_vm()
quit()
