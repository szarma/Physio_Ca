from sys import argv
import javabridge
import bioformats as bf
javabridge.start_vm(class_path=bf.JARS)
import numpy as np
import matplotlib.pyplot as plt
from os.path import isdir
from os import makedirs

def importFirstFrames(rdr,idx,n=100):
    from warnings import warn
    image = []
    for i in range(n):
        try:
            image += [rdr.read(series=idx, rescale=False, t=i)]
        except:
            warn(f"For {idx} series, you asked for {n} first frames, but there was only {len(image)} so I give you that")
            break
    image = np.stack(image)
    return image

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

def getApparentFreq(idx_,xml_=None):
    if xml_ is None:
        # assert "xml" in globals()
        xml_ = xml
    Pixels = xml_.image(index=idx_).Pixels
    nPlanes = Pixels.get_plane_count()
    if nPlanes==1:
        return 0
    else:
        return (nPlanes-1)/Pixels.Plane(nPlanes-1).DeltaT

pathToFile = argv[1]
try:
    iSeries = argv[2]
except:
    iSeries = "all"
separateFolder=True

md = bf.get_omexml_metadata(pathToFile)
xml = bf.OMEXML(md)
Nimages = xml.get_image_count()

if iSeries is not "all":
    iSeries = int(iSeries)
    if iSeries>=Nimages:
        raise ValueError(f"""
                         iSeries ({iSeries}) needs to be smaller than the
                         total number of total series {Nimages}""")

rdr = bf.ImageReader(pathToFile, perform_init=True)
extension = "."+pathToFile.split(".")[-1]

if iSeries!="all":
    nSeries = [iSeries]
else:
    nSeries = range(Nimages)

fn = pathToFile.split("/")[-1]

if separateFolder:
    saveDir = pathToFile.replace(extension,"/")
else:
    saveDir = pathToFile.replace(extension,"/")

if not isdir(saveDir):
    makedirs(saveDir)

for iSeries in nSeries:
    im = xml.image(iSeries)
    Name = im.Name
    if fn not in Name:
        Name = fn.rstrip(extension) + "_"+ Name
    else:
        Name = Name.replace(extension,"").strip().replace("(","").replace(")","").replace("  "," ").replace(" ","_")

    outName = ("%02i_"%iSeries)+Name+".png"
    dimensions = dict(zip("TXY",(getattr(im.Pixels, "Size"+dim) for dim in "TXY")))
    image = importFirstFrames(rdr,iSeries)
    image = np.mean(image,axis=0)
    ## only for Nikon
    if len(image.shape)!=2:
        varImage = image.var(axis=(0,1))
        assert varImage[0]==0 and varImage[2]==0 and len(varImage)==3
        image = image.T[1]

    pxSize = im.Pixels.get_PhysicalSizeX()
    pxUnit = im.Pixels.get_PhysicalSizeXUnit()
    
    text = "\n".join([" %s:%i"%(c,dimensions[c])  for c in "XYT"])
    if dimensions["T"]>1:
        text += "\n f:%.1f Hz"%getApparentFreq(iSeries,xml)

    if extension!=".nd2":
        image = image.T
    plotImage(image.T,pxSize=pxSize,pxUnit=pxUnit,
              savePath=saveDir+outName,
              addInfo=text)
    