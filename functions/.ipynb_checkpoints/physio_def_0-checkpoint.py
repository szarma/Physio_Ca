import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import distributions as dst
from scipy.optimize import curve_fit#,minimize,basinhopping

omeTag = "http://www.openmicroscopy.org/Schemas/OME/2016-06"

def myLogNorm(t,Ampl,offset,loc,scale=.1,s=1):
    yd = dst.lognorm.pdf(t,loc=loc,scale=scale,s=s)
    yd = Ampl*yd/yd.max()+offset
    return yd

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
    else:
        return (nPlanes-1)/Pixels.Plane(nPlanes-1).DeltaT

def guessPars(t,x):
    x0 = np.percentile(x,5)
    crit = x>x0+4
    if np.any(crit):
        iloc = np.where(crit)[0][0]
    else:
        iloc = np.argmax(x)
#     if iloc<0:iloc=0
#     x0 = np.mean(x[:iloc])
    loc = t[iloc]
    imax = iloc+np.argmax(x[iloc:min(len(x),iloc+10)])
    xmaxAvg = x[max(imax-1,0):min(len(x),imax+2)].mean()
    if xmaxAvg<np.percentile(x,10):
        return False
    Ampl = xmaxAvg-x0
    scale = .07
    s = 1.
    return Ampl,x0,loc,scale,s

def multiFun(t,ps):
    if len(ps)==0:
        return np.zeros_like(t)
    else:
        return np.sum([myLogNorm(t,*p) for p in ps], axis=0)

def firstOccurConseq(a_,conseq=1):
    a_=a_.copy()
    icc = np.where(a_)[0]
    if conseq==1:
        return icc[0]
    else:
        for i in range(conseq-1):
            try:
                icc = icc[:-1][np.diff(icc)==1]
            except:
                return -1
        return icc[0]

def EvalModel(ps,t,x):
    xf = multiFun(t,ps)
    return np.sum((x-xf)**2)

def Fit(t,y_,ax=None,nPeaks=30,verbose=False, showFailed=True):
    yoffset = np.percentile(y_,5)
    y = y_-yoffset
    yscale = y.std()*2
    y = y/yscale
    if ax:
        ax.plot(t,y,lw=.5)
    iBegin = 0
    dt = np.diff(t).mean()
#     print (len(t),len(y))
    ps = []
    for ip in range(nPeaks):
        if iBegin>len(y)-10:
            break
        y0 = multiFun(t,ps)
#         if y0 = 0: y0 = y.mean()
        tr,yr = t[iBegin:],(y-y0)[iBegin:]
        p0 = guessPars(tr,yr)
        if not p0:
            if verbose: print ("could not guess, p0=",p0)
            break
        p0 = p0[:-1]
        crit = yr-myLogNorm(tr,*p0)>3
#         iStop = firstOccurConseq(crit,3)
        try:
            iStop = firstOccurConseq(crit,3)
            if iStop<0:
                iStop = 10
            if iStop>len(yr):
                iStop = len(yr)
        except: iStop = len(yr)
        iStop += iBegin
#         if iStop<0: 
#             iStop = min(len(yr),iBegin+1./dt)
#         else:
#         if iStop==iBegin:
#             iStop = iBegin+int(.3/dt)
#         iStop = min(len(y),iStop)
        if verbose: print (ip,iBegin,iStop,p0,t[iBegin],t[iStop-1])
        tr = tr[:iStop-iBegin]
        yr = yr[:iStop-iBegin]
        if ax: 
            ax.plot(tr,multiFun(tr,ps+[p0]),"lightgrey",lw=.8)
        try:
            t0 = p0[2]
            p = curve_fit(myLogNorm,tr,yr,p0=p0,
                          bounds=np.array([
                              (1,np.inf),
                              (-np.inf,np.inf),
                              (t0-.15,t0+.05),
                              (.05,.2),
                              (0.8,2),
                          ])[:len(p0)].T 
                         )[0]
        except:
            if verbose: print ("something was off with fitting")
            iBegin += int(.1/dt)
            continue
#         print (p)
        score    = EvalModel([ ],t[:iStop], y[:iStop]-multiFun(t[:iStop],ps)-y.mean()*int(len(ps)==0))
        newScore = EvalModel([p],t[:iStop], y[:iStop]-multiFun(t[:iStop],ps))
#         score    = EvalModel([ ],tr,yr-yr.mean())
#         newScore = EvalModel([p],tr,yr)

        if newScore>score*.99:
            if verbose: print ("refused since score0 = %.1f and suggested score = %.1f"%(score,newScore))
            if ax and showFailed: ax.plot(tr,multiFun(tr,ps+[p]),color="grey")
            iBegin += int(.1/dt)
            continue
        if verbose: print ("accepted since score0 = %.1f and suggested score = %.1f"%(score,newScore))
        ps += [p]
        score=newScore
        if ax: ax.plot(t,y0+myLogNorm(t,*p))
        iBegin = max(0,iStop - int(.2/dt))
#         iBegin += int(.2/dt)
    return ps

def manyFit(t,x,verbose=False,toll=.1):
    var0 = x.var()
    xf = x.copy()
    PS_ = []
    for i in range(10):
        ps = Fit(t,xf,nPeaks=3)
        if not len(ps):
            print ("stopped because there was nothing else to find")
            break
        xf = xf-multiFun(t,ps)
        var1 = xf.var()
        dr2 = 1-var1/var0
        if dr2<toll:
            print ("stopped because the fve difference was too small: %f --> %f (%f)"%(var0,var1,dr2))
            break
        PS_ += ps
        print ("accepted because the fve difference was large enough: %f --> %f (%f)"%(var0,var1,dr2))
        var0 = var1
#         print (PR_)
    return PS_

def multiProfilePlotWithRoi(pxShows,pxWin,
                          timeWindow=None,
                          timeAverageWindow=1,
                          times_=None,
                          image_=None,
                          firstDeriv=False,
                          stdDev=False,
                          showRaw=False
                          ):
    global dimensions
    from matplotlib.patches import Rectangle
    if times_ is None: times_ = times
    if image_ is None: image_ = image
    if timeWindow is None:
        timeWindow = (times_.min(), times_.max())
    twSlice = slice(*(np.where(times_>=tl)[0][0] for tl in timeWindow))
    times_ = times_[twSlice]
    rebTime    = rebin(times_ , timeAverageWindow,0)

    dd = .2
    fig,axs = plt.subplots(int(stdDev)+1,1,figsize=(15,(1+int(stdDev))*dimensions["Y"]/20))
    try: axs[0]
    except: axs = [axs]
    axs[0].imshow(np.log10(10+image_.mean(axis=(0)).T))
    if stdDev:
        axs[1].imshow(image_.std(axis=(0)).T)
    for j,pxShow in enumerate(pxShows):
        for ax in axs:
            xx = max(pxShow[0]-pxWin,0), min(pxShow[0]+pxWin+1,image_.shape[1])
            yy = max(pxShow[1]-pxWin,0), min(pxShow[1]+pxWin+1,image_.shape[2])
            roi = Rectangle(
                        (xx[0]-.5+dd/2,yy[0]-.5+dd/2),
                        width=xx[1]-xx[0]-dd,
                        height=yy[1]-yy[0]-dd,
                        fill=False,
                        edgecolor="C1"
                    )
            ax.add_patch(roi)
            ax.text(xx[0],yy[0],j,va="top",color="C1")
    fig.tight_layout()
    fig,axs = plt.subplots(len(pxShows),1,figsize=(15,2.5*len(pxShows)))
    try: axs[0]
    except:
        axs = [axs]
    for j in range(len(pxShows)):
        ax = axs[j]
        pxShow = pxShows[j]
        xx = max(pxShow[0]-pxWin,0), min(pxShow[0]+pxWin+1,image_.shape[1])
        yy = max(pxShow[1]-pxWin,0), min(pxShow[1]+pxWin+1,image_.shape[2])
        ii = (twSlice, slice(*xx), slice(*yy) )
        profile = image_[ii].mean(axis=(1,2))
        rebProfile = rebin(profile, timeAverageWindow,0)
        
        c = ax.plot(rebTime,rebProfile,lw=1,label="time profile for the roi")[0].get_color()
        if showRaw: ax.plot(times_,profile,"-",lw=.3,ms=2,c=c)
        ax.set_xlim((timeWindow))
        ax.text(0,1,"\n  %i"%j,va="top",transform=ax.transAxes)
        if firstDeriv:
            offset = np.percentile(profile,5)-profile.std()
            c = ax.plot(rebTime[:-1],np.diff(rebProfile)+offset,lw=1,label="first derivative")[0].get_color()
            if showRaw: ax.plot(times_[:-1],np.diff(profile)+offset,lw=.3,c=c)
        ax.legend()
        ax.grid()

def plotTimeProfileForROI(pxShow,pxWin,
                          timeWindow=None,
                          timeAverageWindow=1,
                          times_=None,
                          image_=None,
                          firstDeriv=False,
                          stdDev=False
                          ):
    from matplotlib.patches import Rectangle
    global dimensions
#     pxWin  = 6 # 0 corresponds to single pixel, 1 corresponds to 3 px, 2 -> 5, 3 -> 7 etc
#     pxShow = (180,7)
    if times_ is None: times_ = times
    if image_ is None: image_ = image
    if timeWindow is None:
        timeWindow = (times_.min(), times_.max())
    xx = max(pxShow[0]-pxWin,0), min(pxShow[0]+pxWin+1,image_.shape[1])
    yy = max(pxShow[1]-pxWin,0), min(pxShow[1]+pxWin+1,image_.shape[2])
    twSlice = slice(*(np.where(times_>=tl)[0][0] for tl in timeWindow))

    ii = (twSlice, slice(*xx), slice(*yy) )
    profile = image_[ii].mean(axis=(1,2))
    times_ = times_[twSlice]
    rebTime    = rebin(times_ , timeAverageWindow,0)
    rebProfile = rebin(profile, timeAverageWindow,0)

    # setup figure layout
    figwidth,figheight = 12,3
#     ar = dimensions["Y"]/dimensions["X"]#*figwidth/figheight
#     hspacing = [ar,.3/figwidth]
#     hspacing += [1-sum(hspacing)]
#     fig = plt.figure(figsize=(figwidth,figheight))
    dd = .2
#     axs = []
#     for i in range(0,len(hspacing),2):
#         axs += [fig.add_axes([0,1-sum(hspacing[:i+1]),1,hspacing[i]], 
#     #                          sharex=axs[-1] if len(axs) else None
#                             )]

#     # plot mean of the whole FOV as a reference and inside ROI

#     axs = [
#         fig.add_axes([.01, 1-.1, dimensions["X"]/dimensions["Y"]*.1, .1]),
#         plt.subplot(111)
#     ]
    plt.figure(figsize=(10,4))
    ax = plt.subplot(111)
#     ax = axs[0]
    ax.imshow(image_.mean(axis=(0)).T)
    roi = Rectangle(
                (xx[0]-.5+dd/2,yy[0]-.5+dd/2),
                width=xx[1]-xx[0]-dd,
                height=yy[1]-yy[0]-dd,
                fill=False,
                edgecolor="C1"
            )
    ax.add_patch(roi)
    if stdDev:
        plt.figure(figsize=(10,4))
        ax = plt.subplot(111)
    #     ax = axs[0]
        ax.imshow(image_.std(axis=(0)).T)
        roi = Rectangle(
                    (xx[0]-.5+dd/2,yy[0]-.5+dd/2),
                    width=xx[1]-xx[0]-dd,
                    height=yy[1]-yy[0]-dd,
                    fill=False,
                    edgecolor="C1"
                )
        ax.add_patch(roi)
        

    # plot time profile of the ROI and its first derivative (sometimes informative)
    plt.figure(figsize=(figwidth,figheight))
    ax = plt.subplot(111)
#     ax = axs[1]
    c = ax.plot(times_,profile,"-",lw=.3,ms=2)[0].get_color()
    ax.plot(rebTime,rebProfile,c=c,lw=1,label="time profile for the roi")
    ax.set_xlim((timeWindow))
    if firstDeriv:
        offset = np.percentile(profile,5)-profile.std()
        c = ax.plot(times_[:-1],np.diff(profile)+offset,lw=.3)[0].get_color()
        ax.plot(rebTime[:-1],np.diff(rebProfile)+offset,c=c,lw=1,label="first derivative")
    ax.legend()
    ax.grid()
    return times_,profile

# def printRecur(root, maxLevel=np.inf, stoppingCrit=None):
#     """Recursively prints the tree."""
#     global level, Levels
#     try: Levels[level] += 1
#     except: pass
#     if stoppingCrit is not None:
#         if stoppingCrit(Levels):
#             return None
#     title = root.tag.title().lower().replace(omeTag.lower(),"").replace("{}","")
#     txt = root.text
#     if txt is None: txt = ""
#     else: txt =  ": "+txt
#     count = ".".join(Levels.astype(int).astype(str))#.rstrip(".0")
#     isplit = count.find(".0")
#     if isplit>0:
#         count = count[:isplit]
#     count += (10-len(count))*" "
#     print (count, end="")
#     print (' '*4*level,title,txt)
#     level += 1
#     if level<=maxLevel:
#         for elem in root:
#             printRecur(elem, maxLevel=maxLevel,stoppingCrit=stoppingCrit)
#     level -= 1
#     Levels[level+1:] = 0

def printRecur(root, maxLevel=np.inf, stoppingCrit=None):
    """Recursively prints the tree."""
    global level, Levels
    try: Levels[level] += 1
    except: pass
    if stoppingCrit is not None:
        if stoppingCrit(Levels):
            return None
    title = root.tag.title().lower().replace(omeTag.lower(),"").replace("{}","")
    txt = root.text
    if txt is None: txt = ""
    else: txt =  ": "+txt
    count = ".".join(Levels.astype(int).astype(str))#.rstrip(".0")
    isplit = count.find(".0")
    if isplit>0:
        count = count[:isplit]
    count += (10-len(count))*" "
    print (count, end="")
    print (' '*4*level,title,txt)
    level += 1
    if level<=maxLevel:
        for elem in root:
            printRecur(elem, maxLevel=maxLevel,stoppingCrit=stoppingCrit)
    level -= 1
    Levels[level+1:] = 0

# =============================================================================
# def getSmoothPx(px,timeWindow,avgWindow=30):
#     assert avgWindow%2==0
#     twIndices = slice(*(np.where(times>=tl)[0][0] for tl in timeWindow))
#     return moving_average(LineT[twIndices,px],avgWindow-1)[::avgWindow]
# 
# def getSmoothT(timeWindow,avgWindow=30):
#     assert avgWindow%2==0
#     twIndices = slice(*(np.where(times>=tl)[0][0] for tl in timeWindow))
#     return moving_average(times[twIndices],avgWindow-1)[::avgWindow]
# 
# =============================================================================