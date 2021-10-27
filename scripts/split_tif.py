#!/opt/tljh/user/envs/physio/bin/python

import argparse
parser = argparse.ArgumentParser(description='Split large tif files into individual pages for easier loading. NB: If you do not need original files afterwards, you need to delete them manually, this script will not do it for you.')
parser.add_argument('path', type=str, help='full path to file; in case of multi-tiffs, a path to the first one.', metavar="--path")
parser.add_argument('--silent', '-s', default=False, help='do not display the progress bar',action="store_const",const=True)
parser.add_argument('--frames', type=str, default="None", help='which frames to take out. For example: range(100), range(100,200,10)')

args = parser.parse_args()

from glob import glob
import os
import tifffile
from tqdm import tqdm

pathToFirst = args.path#"/data/Cosmin/2021_10_21_Human_slice/Slice_1_Troubleshoot/Untitled_1_MMStack.ome.tif"

folder = os.path.split(pathToFirst)[0]

def convert_multitif(fn, output_template, ser=0, frames=eval(args.frames)):
    tif = tifffile.TiffFile(fn)
    ser = tif.series[0]
    n_images=len(ser)
    if frames is None:
        frames = range(n_images)
    if not args.silent:
        frames = tqdm(frames)
    for i in frames:
        pg = ser[i]
        tifffile.imsave(output_template%(i),pg.asarray())

# for f in glob(folder+"/splitted_*.tif"):
#     os.remove(f)
convert_multitif(pathToFirst, folder+"/splitted_%06i.tif")