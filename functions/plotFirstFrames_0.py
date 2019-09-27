import javabridge
import bioformats as bf
javabridge.start_vm(class_path=bf.JARS)
# from collections import OrderedDict
# import pandas as pd
# from scipy.stats import distributions as dst
import numpy as np

from sys import argv
filename = argv[1]

try:
    iSeries = argv[2]
except:
    iSeries = "all"

def plotImage(IM,pxSize=1,pxUnit="µm",rescale=1./30,log=True,savePath = None, addInfo = ""):
    from matplotlib.colors import LogNorm
    figsize = np.array([IM.shape[1]*pxSize,IM.shape[0]*pxSize])*rescale
#     figsize[1] += 1
    fig = plt.figure(figsize=figsize)
    plt.subplot(111)
    plt.imshow(IM,
               norm=LogNorm() if log else None,#vmin=1, vmax=np.log10(IM.max())),
               extent=(0,IM.shape[1]*pxSize,IM.shape[0]*pxSize,0))
    plt.colorbar(shrink = .8)
    plt.xlabel("x [%s]"%pxUnit)
    plt.ylabel("y [%s]"%pxUnit)
    yt = np.arange(0,np.ceil(IM.shape[0]*pxSize)+1e-10,10)
    xt = np.arange(0,np.ceil(IM.shape[1]*pxSize)+1e-10,20)
    plt.yticks(yt,fontsize = 14)
    plt.xticks(xt,fontsize = 14)
    ax = plt.gca()
    plt.text(0,1,"\n"+addInfo,color="white",va="top",fontfamily="Consolas",transform = ax.transAxes)
    if savePath is not None:
        fig.tight_layout()
        fig.savefig(savePath,dpi=120)
        plt.close(fig)

def importFirstFrames(rdr,idx,n=100):
    from warnings import warn
    image = []
    for i in range(n):
        try:
            image += [rdr.read(series=idx, rescale=False, t=t)]
        except:
            warn(f"You asked for {n} first frames, but there was only {len(images)} so I give you that")
            break
    image = np.stack(image)
    return image


md = bf.get_omexml_metadata(filename)
xml = bf.OMEXML(md)
Nimages = xml.get_image_count()

if iSeries is not "all":
    iSeries = int(iSeries)
    if iSeries>=Nimages:
        raise ValueError(f"""
                         iSeries ({iSeries}) needs to be smaller than the
                         total number of total series {Nimages}""")

rdr = bf.ImageReader(dataDir+filename, perform_init=True)

im = xml.image(iSeries)

    
    
pngFile = filename.split(".")[0]+'_%i.png'%iSeries


image = np.mean(,axis=0)


if iSeries == "all":
    for 