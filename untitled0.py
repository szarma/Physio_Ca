### Preface
import javabridge
import bioformats as bf
javabridge.start_vm(class_path=bf.JARS)
import numpy as np
# from os.path import isdir
# from os import makedirs

from sys import path as syspath
syspath.append("./functions/")
from physio_def_1 import getApparentFreq, importFrames#, getTimes
pathToFile = './data/256x32_image_slice3_0la_2Ca001.nd2'
md = bf.get_omexml_metadata(pathToFile)
xml = bf.OMEXML(md)
Nimages = xml.get_image_count()
rdr = bf.ImageReader(pathToFile, perform_init=True)
extension = "."+pathToFile.split(".")[-1]

iSeries = 0

im = xml.image(iSeries)
Name = im.Name
dimensions = dict(zip("TXY",(getattr(im.Pixels, "Size"+dim) for dim in "TXY")))

if dimensions["T"]>1:
    dimensions['freq'] = getApparentFreq(iSeries,xml)

howManyFirstFrames = 1
firstFrames = importFrames(rdr,idx=iSeries, which=(howManyFirstFrames,))

extraDims = len(firstFrames.shape)
if extraDims>3:
    print("There are more dimensions here than expected.\nI'll remedy this for now, but be shure to keep this in mind.")
    retainDims = np.where(firstFrames.std(axis=tuple(range(len(firstFrames.shape)-1)))>0)[0]
    assert len(retainDims)==1
    retainDims = retainDims[0]
    firstFrames = firstFrames.T[retainDims].T
else:
    retainDims = None


from pandas import to_numeric
mintype = to_numeric(firstFrames.flatten(),downcast="unsigned").dtype

firstFrames = firstFrames[:,:2,:2].astype(mintype)

bf.write_image(
    "/Users/srdjan/tmp.tiff",
    firstFrames,
    pixel_type=mintype.name,
#     c=0,
#     z=0,
#     t=np.arange(len(firstFrames))*.003,
#     size_c=1,
#     size_z=1,
#     size_t=len(firstFrames),
#     channel_names="None",
)
#########################
javabridge.kill_vm()
raise SystemExit
#########################

meanFirstFrames = firstFrames.mean(axis=0)
stdFirstFrames = firstFrames.std(axis=0)
pxSize = im.Pixels.get_PhysicalSizeX()
pxUnit = im.Pixels.get_PhysicalSizeXUnit()
text = "\n".join([" %s:%i"%(c,dimensions[c])  for c in "XYT"])
if dimensions["T"]>1:
    text += "\n f:%.1f Hz"%dimensions['freq']






allimage = importFrames(
    rdr,
    idx=iSeries,
    dtype=mintype,
    which=(2,(10,11),(10,11))+tuple([retainDims][:int(extraDims)])
)


bf.write_image(
    "/Users/srdjan/tmp.tiff",
    allimage,
    pixel_type=mintype.name,
#     c=0,
#     z=0,
#     t=np.arange(len(firstFrames))*.003,
#     size_c=1,
#     size_z=1,
#     size_t=len(firstFrames),
#     channel_names="None",
)
