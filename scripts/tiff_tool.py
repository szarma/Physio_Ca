#!/opt/tljh/user/envs/physio/bin/python

import argparse
parser = argparse.ArgumentParser(description='Split or merge large tif(s) for easier loading. NB: If you do not need original files afterwards, you need to delete them manually, this script will not do it for you.')
parser.add_argument('action', type=str, help='merge or split', metavar="--action")
parser.add_argument('path', type=str, help='full path to file; in case of multi-tiffs, a path to the first one.', metavar="--path")
parser.add_argument('--silent', '-s', default=False, help='do not display the progress bar',action="store_const",const=True)
parser.add_argument('--frames', type=str, default=None, help='which frames to include. For example: range(100), range(100,200,10)')
parser.add_argument('--drop-even', default=False, help='whether to drop even frames',action="store_const",const=True)
parser.add_argument('--rebin', type=int, default=2, help='how much to rebin frames, default is 2')
parser.add_argument('--outtype', type=str, default="uint8", help='output type, default is uint8')
#parser.add_argument('--output-name',type=str, default="merged.tif", help="output name")
parser.add_argument('--rescale',type=float, default=None, help="rescale each frame by")


args = parser.parse_args()

if args.action not in ["merge","split"]:
    raise ValueError("action can be either merge or split.")

from glob import glob
import os
import tifffile
from tqdm import tqdm
from islets.numeric import rebin

pathToFirst = args.path
folder = os.path.split(pathToFirst)[0]

# Figure out the shape of the output:
allfiles = sorted(glob(args.path))
if len(allfiles)>1:
    # action needs to be merge
    if args.action=="split":
        raise ValueError("You gave me a series of files, and asked to split them. Point me only to the first one in the series.")
    x = tifffile.imread(allfiles[0])
    n_images = len(allfiles)
    mode="merge many"
else:
    # action cab be either
    tif = tifffile.TiffFile(allfiles[0])
    ser = tif.series[0]
    n_images=len(ser)
    pg = ser[0]
    x = pg.asarray()
    mode = args.action

nrebin = args.rebin

if args.drop_even:
    ser_length = n_images//2
else:
    ser_length = n_images

if args.frames is None:
    frames = range(ser_length)
else:
    frames = eval(args.frames)


shape = (len(frames), x.shape[0]//nrebin, x.shape[1]//nrebin, )

if not args.silent:
    frames = tqdm(frames)

if mode=="split":
    output_template = os.path.splitext(args.path)[0] + "_splitted_%06i.tif"
else: # "merge" and "merge many"
    output_filename = os.path.splitext(allfiles[0])[0] + "_merged.tif"
    data = tifffile.memmap(
        output_filename,
        shape=shape,
        dtype=args.outtype,
        photometric="minisblack"
    )

maxValid = 2**int(args.outtype.split("int")[-1])-1
if args.rescale is None:
    maxv = x.max()*2
    rescale = maxValid/maxv
else:
    rescale=args.rescale

print ("processing...")
for i in frames:
    if args.drop_even: ichoose = i*2
    else: ichoose = i
    if mode=="merge many":
        x = tifffile.imread(allfiles[ichoose])
    else:
        x = ser[ichoose].asarray()
    y = rebin(x*rescale, (nrebin, nrebin), axis=(0,1), dtype="uint8")
    if mode=="split":
        tifffile.imsave(output_template % (i), y)
    else:
        data[i,:] = y

if "merge" in mode:
    data.flush()

quit()