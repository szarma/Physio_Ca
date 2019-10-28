import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import distributions as dst
from scipy.optimize import curve_fit#,minimize,basinhopping

omeTag = "http://www.openmicroscopy.org/Schemas/OME/2016-06"

def rebin(a,n,axis=0,norm=True):
    ashape = a.shape
    newShape = ashape[:axis]+(ashape[axis]//n,n)+ashape[axis+1:]
    idx = tuple([slice(None)] * axis + [slice(ashape[axis]//n*n)] + [slice(None)]*(len(ashape)-axis-1))
    out = a[idx].reshape(newShape)
    out = out.sum(axis=axis+1)
    if norm:
        out = out/n
    return out

def getDimensions(idx_,xml_=None):
    if xml_ is None:
        # assert "xml" in globals()
        xml_ = xml
    Pixels = xml_.image(index=idx_).Pixels
    return sum([[
        (coord       ,getattr(Pixels,"PhysicalSize"+coord)),
        ("[%s]"%coord,getattr(Pixels,"PhysicalSize%sUnit"%coord))]
        for coord in "XYZ"],[])

def getTimes(idx_,xml_=None):
    if xml_ is None:
        # assert "xml" in globals()
        xml_ = xml
    Pixels = xml_.image(index=idx_).Pixels
    nPlanes = Pixels.get_plane_count()
    return np.array([Pixels.Plane(i).DeltaT for i in range(nPlanes)])

def getApparentFreq(idx_,xml_=None):
    if xml_ is None:
        # assert "xml" in globals()
        xml_ = xml
    Pixels = xml_.image(index=idx_).Pixels
    nPlanes = Pixels.get_plane_count()
    if nPlanes==1:
        return 0
    for j in range(1,10):
        frq = (nPlanes-j)/Pixels.Plane(nPlanes-j).DeltaT
        if frq>0:
            return frq

def importFrames(rdr,idx=0,which=None,dtype=None):
    from warnings import warn
    firstImage = rdr.read(series=idx, rescale=False, t=0)
    dims = firstImage.shape
    if which is None:
        which = [None]*len(dims)
    if len(which)>len(dims)+1:
        raise ValueError(f"the number of dimensions to import ({len(which)}) is larger than the number of available dimensions {len(dims)}")
    if len(which)<len(dims)+1:
        which = list(which)+[None]*(len(dims)-len(which)+1)
    ix = tuple()
    for ix_ in which[1:]:
        try:
            if type(ix_) is int:
                ix += (ix_,)
            else:
                try:
                    ix += (slice(*ix_),)
                except:
                    ix += (slice(ix_),)
        except:
            raise ValueError("`which` needs to NaN(s) or iterable")
    t=0
    if which[0] is None:
        Ts = None
        image = []
    else:
        if type(which[0]) is int:
            Ts = range(which[0]) 
        else:
            Ts = which[0]
            t=Ts[0]
        image = np.zeros( (len(Ts),) + firstImage[ix].shape)
        if dtype is not None:
            image = image.astype(dtype)
    
    while True:       ######## dangerous!
        try:
            nextImage = rdr.read(series=idx, rescale=False, t=t)
            if dtype is not None:
                nextImage = nextImage.astype(dtype)
            if Ts is None:
                image += [nextImage[ix]]
            else:
                image[t-Ts[0]] = nextImage[ix]
        except:
            warn(f"Could not give all required time points. I advise you double check the output")
            break
        if Ts is not None:
            if t not in Ts:
                break
        t += 1
    if Ts is None:
        image = np.stack(image)
        if dtype is not None:
            image = image.astype(dtype)
    return image
# def importFrames(rdr,idx=0,which=None,dtype=None):
#     from warnings import warn
#     firstImage = rdr.read(series=idx, rescale=False, t=0)
#     dims = firstImage.shape
#     if which is None:
#         ix = tuple(slice(None) for j in dims)
#         image = []
#         t=0
#         while True:       ######## dangerous!
#             try:
#                 nextImage = rdr.read(series=idx, rescale=False, t=t)
#                 if dtype is not None:
#                     nextImage = nextImage.astype(dtype)
#                 image += [nextImage[ix]]
#             except:
#                 break
#             t += 1
#         image = np.stack(image)
#     else:
#         try:
#             if type(which[0]) is int:
#                 Ts = range(which[0]) 
#             else:
#                 Ts = which[0]
#         except:
#             raise ValueError("`which` needs to NaN(s) or iterable")
#         ix = tuple()
#         for ix_ in which[1:]:
#             if hasattr(ix_,"__index__"):
#                 ix += (ix_,)
#             else:
#                 try:
#                     ix += (slice(*ix_),)
#                 except:
#                     ix += (slice(ix_),)
                
#         image = np.zeros( (len(Ts),) + firstImage[ix].shape)
#         if dtype is not None:
#             image = image.astype(dtype)
#         for i,t in enumerate(Ts):
#             try:
#                 nextImage = rdr.read(series=idx, rescale=False, t=t)        
#                 image[i] = nextImage[ix]
#             except:
#                 warn(f"Could not give all required time points. I advise you double check the output")
#                 break
#     return image

    
def plotImageWithRois(
    pxShows=[],
    pxWin=1,
    image_=None,
    stdDev=False,
    imgHeight = None,
    axs = None,
    label=True
):
    from collections import OrderedDict
    if imgHeight is None:
        imgHeight = image_.shape[-1]/20
    if axs is None:
        fig,axs = plt.subplots(int(stdDev)+1,1,figsize=(15,(1+int(stdDev))*imgHeight))
    try: axs[0]
    except: axs = [axs]
    axs[0].imshow(np.log10(10+image_.mean(axis=(0)).T))
    if stdDev:
        axs[1].imshow(image_.std(axis=(0)).T)
    for ax in axs:
        addRoisToImage(
            pxShows=pxShows,
            pxWin=pxWin,
            ax = ax,
            label=label
        )
#     if isinstance(pxShows,dict):
#         roiLabels,roiCoords = pxShows.keys(), pxShows.values()
#     else:
#         roiLabels = range(len(pxShows))
#         roiCoords = pxShows
#     roiProfiles = OrderedDict()
#     for roiLabel,pxShow in zip(roiLabels,roiCoords):
#         xx = pxShow[0], min(pxShow[0]+pxWin,image_.shape[1])
#         yy = pxShow[1], min(pxShow[1]+pxWin,image_.shape[2])
#         roiProfiles[roiLabel] = image_[:,slice(*xx),slice(*yy)].mean(axis=(1,2))
#         for ax in axs:
#             roi = Rectangle(
#                         (xx[0]-.5+dd/2,yy[0]-.5+dd/2),
#                         width=xx[1]-xx[0]-dd,
#                         height=yy[1]-yy[0]-dd,
#                         fill=False,
#                         edgecolor="C1"
#                     )
#             ax.add_patch(roi)
#             if label:
#                 ax.text(xx[0],yy[0],roiLabel,va="top",color="C1")
#     # fig.tight_layout()
    return axs

def addRoisToImage(
    pxShows,
    pxWin=1,
    stdDev=False,
    ax = None,
    label=True,
    color="C1",
    rectangle_kw=None,
    dd = .1
):
    from matplotlib.patches import Rectangle
    if isinstance(pxShows,dict):
        roiLabels,roiCoords = pxShows.keys(), pxShows.values()
    else:
        roiLabels = range(len(pxShows))
        roiCoords = pxShows
    if rectangle_kw is None:
        rectangle_kw = dict(fill=False,edgecolor=color)
    for roiLabel,pxShow in zip(roiLabels,roiCoords):
        xx = pxShow[0], pxShow[0]+pxWin
        yy = pxShow[1], pxShow[1]+pxWin
        roi = Rectangle(
                    (xx[0]-.5+dd/2,yy[0]-.5+dd/2),
                    width=xx[1]-xx[0]-dd,
                    height=yy[1]-yy[0]-dd,
                    **rectangle_kw
                )
        ax.add_patch(roi)
        if label:
            ax.text(xx[0],yy[0],roiLabel,va="center",color=color)

def getRoiProfiles(
    pxShows,
    pxWin=1,
    image_=None,
    tWin_ = 1
):
    from collections import OrderedDict
    if isinstance(pxShows,dict):
        roiLabels,roiCoords = pxShows.keys(), pxShows.values()
    else:
        roiLabels = range(len(pxShows))
        roiCoords = pxShows
    roiProfiles = OrderedDict()
    for roiLabel,pxShow in zip(roiLabels,roiCoords):
        xx = pxShow[0], min(pxShow[0]+pxWin,image_.shape[1])
        yy = pxShow[1], min(pxShow[1]+pxWin,image_.shape[2])
        roiProfiles[roiLabel] = rebin(image_[:,slice(*xx),slice(*yy)].mean(axis=(1,2)),tWin_)
    return roiProfiles
