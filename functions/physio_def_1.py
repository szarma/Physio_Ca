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

# def getApparentFreq(idx_,xml_=None):
#     if xml_ is None:
#         # assert "xml" in globals()
#         xml_ = xml
#     Pixels = xml_.image(index=idx_).Pixels
#     nPlanes = Pixels.get_plane_count()
#     if nPlanes==1:
#         return 0
#     for j in range(1,10):
#         frq = (nPlanes-j)/Pixels.Plane(nPlanes-j).DeltaT
#         if frq>0:
#             return frq

def importFrames(rdr,idx=0,which=None,dtype=None):
    from warnings import warn
    firstImage = rdr.read(series=idx, rescale=False, t=0)
    dims = firstImage.shape
    if which is None:
        which = [None]*len(dims)
    if len(which)>len(dims)+1:
        raise ValueError("the number of dimensions to import (%i) is larger than the number of available dimensions %i"%(len(which),len(dims)))
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
        if dtype is None:
            dtype = "float32"
        image = np.zeros( (len(Ts),) + firstImage[ix].shape, dtype=dtype)
    
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
            warn("Could not give all required time points. I advise you double check the output")
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

def showMovie(m_show, figsize = (6,6), out="jshtml",fps = 30, saveName=None, NTimeFrames=100,log=True,additionalPlot=None):
    import matplotlib.pyplot as plt
    from matplotlib import animation
    if NTimeFrames is not None:
        n_rebin = len(m_show)//NTimeFrames
        if n_rebin>1:
            m_show = rebin(m_show, n_rebin)
    if log:
#         while True:
        for p in range(1,5):
            baseline = np.percentile(m_show,p)
            m_show = np.maximum(m_show, baseline)
            if np.all(m_show>0): break
        m_show = np.log(m_show)
    fig, ax = plt.subplots(figsize=figsize,dpi=150)
    im = ax.imshow(m_show[0].T, cmap="Greys", vmin=0, vmax=m_show.max())
    if additionalPlot is not None:
        additionalPlot(ax)
    plt.close(fig)
    def init():
        im.set_data(m_show[0].T)
        return (im,)
    def animate(i):
        im.set_data(m_show[i].T)
        return (im,)
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(m_show),
                                   interval=1000/fps,
                                   blit=True)
    if out=="html5":
        from IPython.display import HTML
        return HTML(anim.to_html5_video())
    if out=="jshtml":
        from IPython.display import HTML
        return HTML(anim.to_jshtml())
    if out=="save" or saveName is not None:
        try:
            anim.save(saveName)
        except:
            saveName = input("please enter a valid filename. Otherwise, I'll save it as 'video.mp4'.")
            try: anim.save(saveName)
            except:
                saveName = "video.mp4"
                anim.save(saveName)
        return None
    

def getBaseName(npzFiles):
    from os.path import split
    dfName = npzFiles[0]
    dfName = dfName.replace(".npz","")
    if len(npzFiles)>1:
        tmp = [split(fn)[1].split("_")[-1].split(".")[0].replace("Series","") for fn in npzFiles]
        dfName = dfName.replace(tmp[0],tmp[0]+"-"+tmp[-1])
    return dfName

def import_npz_files(npzFiles):
    orig_images = []
    for npzFile in npzFiles:
        npzData = np.load(npzFile)
        orig_images += [npzData["data"]]
    orig_images = np.concatenate(orig_images)
    return orig_images

# plotting definition
def plotInteresting(investigateIndex_, onlyImage=False, figsize = (10,10)):
    fig, axs = plt.subplots(1,2,figsize=[figsize[0],figsize[1]/2], sharex=True, sharey=True)
    for ax,im in zip(axs,["std_dres","template"]):
        showImage = images[im].astype("float")
        if im=="std_dres":
            ax.imshow(.01+showImage,cmap="Greys",norm=LogNorm(),vmax=showImage.max()*10)
        else:
            x = np.arctan(30*showImage)
            v = np.percentile(np.abs(x),99)*4
            im = ax.imshow(x,cmap="bwr", vmin=-v, vmax=v)
        regions.plotEdges(ix=investigateIndex_,separate=True,image=False, ax=ax)
        regions.plotPeaks(ix=investigateIndex_,image=False, ax=ax,labels=True)
    fig.show()
    
    if onlyImage: return None

    fig, axs = plt.subplots(2,1,figsize=figsize, sharex=True)
    ia = 0
    n = 1
    ns = 2

    # inax = inset_axes(axs[0],
    #                     width="30%", # width = 30% of parent_bbox
    #                     height=1.2, # height : 1 inch
    #                     loc=2)
    # showImage = images["std_dres"].astype("float")
    # inax.imshow(.1+showImage,cmap="Greys",norm=LogNorm(),vmax=showImage.max()*10)
    # regions.plotEdges(ix=investigateIndex_,separate=True,ax=inax)
    t = rebin(time,n)
    for roi in investigateIndex_:
        x  = rebin(C.loc[roi,"trace"],n)
        xf = rebin(C.loc[roi,"faster_%g"%ironScale],n)
        xsd= rebin(C.loc[roi,"faster_%g_std"%ironScale],n)/n**.5
        xs = rebin(C.loc[roi,"slower_%g"%ironScale],n)
        eventFilter = posRender[np.where(C.index==roi)[0][0]]

        xsd= xsd/xf.std()
        xf = xf/xf.std()

        c  = "C%i"%(roi%10)
        axs[0].plot(t,x [::],c=c,lw=.4,alpha = .3)
        axs[0].plot(t,xs[::],c=c,lw=.7,alpha = 1)
        axs[0].plot(t[eventFilter],xs[eventFilter],".",c=c,ms=2,alpha = .2)
        yoffset = 13*ia
        axs[1].plot(t,xf[::]+yoffset,c=c,lw=.3)
        axs[1].plot(t,+ns*xsd[::]+yoffset,c=c,lw=.7)
        axs[1].plot(t,-ns*xsd[::]+yoffset,c=c,lw=.7)
        axs[1].axhline(yoffset,color=c,label=roi)
        axs[1].text(0,yoffset,str(roi)+" ",fontdict={"color":c},ha="right",va="center")
        axs[0].text(0,xs[0],str(roi)+" ",fontdict={"color":c},ha="right",va="center")
    # axs[1].legend(loc=(1.01,.01))
    # #     x,y = C.loc[roi,"peak"]#np.mean(C.loc[i,"pixels"],axis=0)
    # #     inax.plot(x,y,"o",mfc="none",ms=C.loc[i,"size"]**.5/2,c=c,mew=.3)
    #     bb = list(C.loc[roi,"boundary"])
    #     bb += [bb[0]]
    #     x,y = np.array(bb).T
    #     inax.plot(x,y,c=c,lw=.5)
        ia += 1
    # yl = axs[0].get_ylim()[1]
    # axs[0].set_ylim(None,yl*1.2)

    # plt.xticks([])
    # plt.yticks([])
    fig.tight_layout()
    fig.show()