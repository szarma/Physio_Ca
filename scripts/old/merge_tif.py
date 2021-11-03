#!/opt/tljh/user/envs/physio/bin/python

import argparse
parser = argparse.ArgumentParser(description='Merge and compress large tif file(s) into a single file. NB: If you do not need original files afterwards, you need to delete them manually, this script will not do it for you.')
parser.add_argument('path', type=str, help='full path to file; in case of multi-tiffs, a path to the first one.', metavar="--path")
parser.add_argument('--silent', '-s', default=False, help='do not display the progress bar',action="store_const",const=True)
parser.add_argument('--drop-even', default=False, help='whether to drop even frames',action="store_const",const=True)
# parser.add_argument('--frames', type=str, default="None", help='which frames to take out. For example: range(100), range(100,200,10)')
parser.add_argument('--rebin', type=int, default=2, help='how much to rebin frames')
parser.add_argument('--outtype', type=str, default="uint8", help='output type')
parser.add_argument('--output-name',type=str, default="merged.tif", help="output name")
parser.add_argument('--rescale',type=float, default=None, help="rescale each frame by")
parser.add_argument('--frames', default=None, help="frames to include")

args = parser.parse_args()

from glob import glob
import os
import tifffile
from tqdm import tqdm
from islets.numeric import rebin
from islets.utils import save_tiff
from numpy import zeros
pathToFirst = args.path
folder = os.path.split(pathToFirst)[0]

tif = tifffile.TiffFile(args.path)
ser = tif.series[0]
n_images=len(ser)
pg = ser[0]
x = pg.asarray()
nrebin = args.rebin

if args.drop_even:
    ser_length = len(ser)//2
else:
    ser_length = len(ser)

if args.frames is None:
    frames = range(ser_length)
else:
    frames = eval(args.frames)

shape = (len(frames), x.shape[0]//nrebin, x.shape[1]//nrebin, )
data = zeros(shape=shape,dtype=args.outtype)

if args.rescale is None:
    maxv = x.max()*2
    rescale = 255/maxv
else:
    rescale=args.rescale

outputFilename = os.path.join(folder, args.output_name)
if os.path.isfile(outputFilename):
    raise ValueError(f"{outputFilename} already exists. please remove it, or pass a different output name.")



if not args.silent:
    frames = tqdm(frames)
print ("loading data...")
for i in frames:
    if args.drop_even:
        ichoose = i*2
    else:
        ichoose = i
    x = ser[ichoose].asarray()
    data[i,:] = rebin(x*rescale, (nrebin, nrebin), axis=(0,1), dtype="uint8")

print ("writing data... (this will probably take longer than loading, just saying) ")
save_tiff(data,outputFilename)