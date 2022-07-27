import warnings
from collections import OrderedDict
from itertools import product
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from numba import jit, prange
from scipy.optimize import curve_fit  # ,minimize,basinhopping
from islets.utils import multi_map
from scipy.stats import median_abs_deviation
from scipy.spatial import distance_matrix


mad2std = 1.4826 # https://en.wikipedia.org/wiki/Median_absolute_deviation#Relation_to_standard_deviation

def jitter(x, dist=None, amount=None):
    x = np.array(x)
    if len(x)==0:
        return np.zeros_like(x)
    if dist is None:
        dist = median_abs_deviation(x)
    if amount is None:
        amount = dist
    jtr = np.random.rand(len(x))-.5
    dm = distance_matrix(x.reshape(-1,1),x.reshape(-1,1))
    for j in range(len(x)):
        jtr[j] *= max(0,(dm[j]<dist).sum()-2)
    if np.abs(jtr).max()>0:
        jtr *= amount/np.abs(jtr).max()
    return jtr

def hillCurve(x, h=2, xc=4):
    return (x/xc)**h/(1+(x/xc)**h)*.94+.03

def bspline(cv, n=100, smoothness=3):
    """ Calculate n samples on a period spline
        cv :        Array ov control vertices
        n  :        Number of samples to return
        smoothness: Degree of smoothness
        taken from https://stackoverflow.com/a/33962986/2043500
    """
    from scipy.interpolate import splev, splprep
    cv = np.asarray(cv)
    tck = splprep(cv.T, s=smoothness, per=True, quiet=3)[0]
    # Calculate query range
    u = np.linspace(0,1,n)
    # Calculate result
    return np.array(splev(u, tck)).T


def fast_filter(absdata,
                ironTimeScale,
                gain=1,
                freq=1,
                z_sp=3,
                order=5,
                verbose=False,
                npass=1,
                nreflect="auto",
                dilation=True,
#                 debugAx=None
              ):
    if dilation:
        from cv2 import dilate
    cutFreq = .5/ironTimeScale
    sosFilter_ = sosFilter(cutFreq, freq, order=order)
    if isinstance(nreflect,str):
        if nreflect=="auto":
            nreflect = int(ironTimeScale*freq)
        if nreflect=="half":
            nreflect = absdata.shape[1]//2
    # assuming nreflect is integer
    if verbose: print (f"introducing reflecting boundary {nreflect} points wide.")
    x = np.hstack([absdata[:,:nreflect][:,::-1], absdata, absdata[:,absdata.shape[1]-nreflect:][:,::-1]])
    if verbose: print (f"this changed the shape of the input array to", x.shape)
    absSlow = sosFilter_.run(x)[:,nreflect:x.shape[1]-nreflect]
    absFast = absdata-absSlow
#     if debugAx is not None:
#         plot=True
#         axs,n = debugAx
#         axs[0].plot(absdata[0,:n])
#         axs[0].plot(absSlow[0,:n],"k--")
#         axs[1].plot(absFast[0,:n],lw=.5)
    var = (absSlow*gain)+1
    if z_sp!=0:
        for ipass in range(npass):
            th = z_sp*(var**.5)
#             if plot: axs[0].plot(th[0,:n]+absSlow[0,:n],"k--")
            ff = (absFast > th).astype("uint8")
#             if plot: axs[1].plot(th[0,:n],"k--")
            ii = np.arange(ff.shape[1])
            if ff.any():
                if dilation:
                    dilateKernelSize = int(ironTimeScale*freq*.2)
                    if dilateKernelSize%2==0:
                        dilateKernelSize+=1
                    if dilateKernelSize>=3:
                        if verbose:
                            print ("dilating by", dilateKernelSize)
                        ff = dilate(ff, np.ones((1,dilateKernelSize), dtype = np.uint8)).astype("bool")
#                         plt.plot(ff[0]*10)
#                         plt.imshow(ff, aspect="auto")
#                         plt.show()
#                 dFast = np.zeros_like(absFast)
                for j in range(absdata.shape[0]):
                    if not ff[j].any(): continue
#                     plt.plot(ff[j]*100)
                    y = absdata[j].copy()
                    y[ff[j]] = np.interp( ii[ff[j]], ii[~ff[j]], absdata[j][~ff[j]])
#                     absSlow[j] = np.maximum(0, sosFilter_.run(y))
                    absSlow[j] = sosFilter_.run(y)
#                 absSlow = absSlow + dFast
#                 absSlow[absSlow==0] = absSlow[absSlow>0].min()
                assert (absSlow>=0).all()
                absFast = absdata - absSlow
            var = absSlow*gain+1
    zScore = absFast/(var)**.5
    return absSlow, absFast, zScore

def robust_max(a,Nlast=10,absolute=True):
    if absolute:
        a = np.abs(a)
    a = a[np.isfinite(a)]
    if Nlast<2:
        return a.max()
    else:
        return np.percentile(a,100*(1-Nlast/a.size))

def rebin(a,n,axis=0, dtype="float32", func=np.mean):
    if type(axis)==int and type(n)==int:
        ashape = a.shape
        newShape = ashape[:axis]+(ashape[axis]//n,)+ashape[axis+1:]
        out = np.zeros(newShape,dtype=dtype)
        for i in range(newShape[axis]):
            idx = tuple([slice(None)] * axis + [slice(i*n,(i+1)*n)] + [slice(None)]*(len(ashape)-axis-1))
            x = func(a[idx],axis=axis)
            out[tuple([slice(None)] * axis + [i] + [slice(None)]*(len(ashape)-axis-1))] = x.astype(dtype)
        return out
    else:
        assert len(n)==len(axis)
        out = rebin(a,n[0],axis[0],func=func)
        for i in range(1, len(n)):
            out = rebin(out,n[i],axis[i],func=func,dtype=dtype)
        return out

class sosFilter:
    def __init__(self, cutoff, frequency, order=5, btype="low"):
        from scipy.signal import butter
        nyq = 0.5 * frequency
        try:
            normal_cutoff = cutoff / nyq
        except:
            normal_cutoff = np.array(cutoff) / nyq
        self.sos_pars = butter(order, normal_cutoff, btype=btype, analog=False, output="sos")
        self.frequency = frequency
        
    def run(self,data):
        from scipy.signal import sosfiltfilt
        return sosfiltfilt(self.sos_pars, data)
    
    def get_shape(self, worN=2000):
        from scipy.signal import sosfreqz
        w, h = sosfreqz(self.sos_pars, worN=worN)
        return 0.5*self.frequency*w/np.pi, np.abs(h)

class ManualFilter:
    def __init__(self, cutoff, frequency, percentile=10, average=True):
        self.wIron = int(frequency / cutoff / 2) * 2 + 1
        self.percentile = percentile,
        self.average = average

    def run(self, traces, n_processes=1):
        global iterf
        if self.average:
            def iterf(x_):
                out = lowPass(x_, self.wIron, self.wIron, perc = self.percentile)
                return out
        else:
            def iterf(x_):
                out = lowPass(x_, self.wIron, perc = self.percentile)
                return out
        if len(traces.shape)==1:
            return iterf(traces)
        return multi_map(iterf, traces, processes = n_processes)

def power_spectrum(x, fr, mean=True):
    from scipy.fft import rfft, rfftfreq
    freq = rfftfreq(x.shape[-1], 1/fr)
    FA = rfft(x)
    if mean: 
        ndim = len(x.shape)
        if ndim>1:
            FA = FA.mean(axis=tuple(range(ndim-1)))
    return freq, np.abs(FA)

def get_sep_th(x_,ax=None,plot=False,thMax=None,log=False):
    from scipy.stats import gaussian_kde
    if ax is None:
        fig = plt.figure()
        ax = plt.subplot(111)
    if log:
        x_=np.log10(x_)
    x_ = x_[np.isfinite(x_)]
    gkde = gaussian_kde(x_)
#     if plot:
    ax.hist(x_,100,histtype="step",density=True)
    if thMax is not None:
        x_ = x_[x_<thMax]
    if thMax == "median":
        x_ = x_[x_<np.median(x_)]
    h,edges = np.histogram(x_,100)
    bincenters = (edges[:-1]+edges[1:])/2
#     bincenters = bincenters[:np.argmax(h)]
    gkde_vals = gkde.evaluate(bincenters)
    peaks = np.array(list(crawlDict(gkde_vals.reshape(-1,1), processes=1).keys())).T[0]
    twoHighestPeaks = sorted(np.argsort(-h[peaks])[:2])
    if len(twoHighestPeaks)==1:
        return -np.inf
    lowerPeak, higherPeak = peaks[twoHighestPeaks]
    bincenters = bincenters[lowerPeak:higherPeak]
    if len(bincenters)<3:
        if not plot: plt.close(fig)
        return x_.min()
    gkde_vals = gkde.evaluate(bincenters)
#     if plot:
    ax.plot(bincenters,gkde_vals)
    sinks = np.array(list(crawlDict(-gkde_vals.reshape(-1,1)+gkde_vals.max()+1, processes=1).keys())).T[0]
    sinks = sorted(sinks,key=lambda xi: gkde_vals[xi])
    th = bincenters[sinks[0]]
#     if plot:
        # ax.plot(bincenters,gkde_vals)
    ax.axvline(th,color="r",ls="--")
    if not plot: plt.close()
    if log:
        return 10**(th)
    return th

def decay(time,top,bottom,rate):
#     top,bottom,rate = p
    return bottom+(top-bottom)*np.exp(-time*rate)

def guessDecayPars(y):
    r0 = 1e-6
    t0 = np.percentile(y,60)
    b0 = np.percentile(y,40)
    p0 = (t0,b0,r0)
    return p0

def fit_bleaching(x_, points=100, func=None, lin_corr=True):
    from scipy.interpolate import interp1d
    N = len(x_) 
    t_ = np.arange(len(x_))
    if points!=None:
        nrebin = int(len(x_)//points)
        if nrebin>1:
            if func is None:
                func = np.mean
        else:
            nrebin = 1
    else:
        nrebin = 1
    if nrebin>1:
#         x = rebin(x_, nrebin, func=func)
        x = np.percentile(x_[:len(x_)//nrebin*nrebin].reshape((-1,nrebin)), 10, axis=1)
#         x = rebin(x_, nrebin, func=func)
        t = rebin(t_, nrebin, func=np.mean)
    else:
        x = x_
        t = t_
#     try:
#     out = decayfit(x)
    x = x[len(x)//10:-len(x)//10-1]
    t = t[len(t)//10:-len(t)//10-1]
    out = simple_decayfit(x)
    
    out = interp1d(t, out, fill_value="extrapolate")(t_)
#     except:
#         pass
#     if np.isnan(out).any():
#         out = x_ - x_.mean()
    if lin_corr:
        p = np.polyfit(t_,x_-out,1)
        out += p[1]+p[0]*np.arange(N)
    return out

def decayfit(x,Ntrials=None, outPars=False):
    if Ntrials is None:
        lx = 1
    else:
        lx = Ntrials
    nx = len(x)//(lx+5)
    TT = np.arange(len(x))
    tt = TT.copy()
    for j_ in range(lx):
        try:
            p0 = guessDecayPars(x)
            ff = np.isfinite(x)
            popt = curve_fit(decay,tt[ff],x[ff],p0=p0,bounds=[(None, None),(None, None),(0,None)])[0]
            expDecay = decay(TT,*popt)
            if outPars:
                return expDecay, popt
            else:
                return expDecay
        except:
            if j_%2==0:
                tt=tt[:len(x)-nx]
                x = x[:len(x)-nx]
            else:
                tt=tt[nx:]
                x = x[nx:]
    popt = (np.nan,)*3
    expDecay = decay(TT,*popt)
    if outPars:
        return expDecay, popt
    else:
        return expDecay


def simple_decayfit(x):
    p0 = guessDecayPars(x)
    ff = np.isfinite(x)
    tt = np.arange(len(x))
    try:
        popt = curve_fit(decay,tt[ff],x[ff],p0=p0)[0]
    except:
        warnings.warn("Did not manage to fit the decay.")
        popt = p0
    expDecay = decay(tt,*popt)
    return expDecay

@jit 
def percFilter(x_,perc,filterSize):
    if filterSize%2==0:
        raise ValueError("filter size needs to be odd number")
    delta = filterSize//2
    out = np.zeros_like(x_)
    x_ = np.concatenate((x_[:delta][::-1], x_, x_[-delta:][::-1]))
    for i in prange(len(out)):
        out[i] = np.nanpercentile(x_[i:i+filterSize],perc)
    return out

# @jit 
def runningAverage(x_,filterSize):
    if filterSize%2==0:
        raise ValueError("filter size needs to be odd number")
    delta = filterSize//2
    out = np.zeros_like(x_)
    x_ = np.concatenate((x_[:delta][::-1], x_, x_[-delta:][::-1]))
    for i in range(len(out)):
        out[i] = np.mean(x_[i:i+filterSize],axis=0)
    return out

@jit 
def runningStd(x_,filterSize):
    if filterSize%2==0:
        raise ValueError("filter size needs to be odd number")
    delta = filterSize//2
    out = np.zeros_like(x_)
    x_ = np.concatenate((x_[:delta][::-1], x_, x_[-delta:][::-1]))
    for i in range(len(out)):
        out[i] = np.nanstd(x_[i:i+filterSize],axis=0)
    return out

@jit 
def runningMin(x_,filterSize):
    if filterSize%2==0:
        raise ValueError("filter size needs to be odd number")
    delta = filterSize//2
    out = np.zeros_like(x_)
    x_ = np.concatenate((x_[:delta][::-1], x_, x_[-delta:][::-1]))
    for i in range(len(out)):
        out[i] = np.nanmin(x_[i:i+filterSize],axis=0)
    return out

@jit 
def runningMax(x_,filterSize):
    if filterSize%2==0:
        raise ValueError("filter size needs to be odd number")
    delta = filterSize//2
    out = np.zeros_like(x_)
    x_ = np.concatenate((x_[:delta][::-1], x_, x_[-delta:][::-1]))
    for i in prange(len(out)):
        out[i] = np.nanmax(x_[i:i+filterSize],axis=0)
    return out

def lowPass(x, wIron, wAvg=None, perc=5,npass = 1):
    try: perc[0]
    except: perc = [perc]*npass
    npass = len(perc)
    ironedx = percFilter(x, perc[0], wIron)
    if wAvg is not None:
        ironedx = runningAverage(ironedx,wAvg)
    for i in range(npass-1):
        ironedx += lowPass(x-ironedx,wIron,wAvg,perc[i+1],)
    return ironedx

def getInterestingROIcenters(image_, Nsample=10000, exponent=3, bandwidth=10, n_jobs=4):
    from sklearn.cluster import MeanShift
    p = image_.flatten()**exponent
    p = p/p.sum()
    yd,xd = np.array([(x//image_.shape[1], x%image_.shape[1]) for x in np.random.choice(np.arange(image_.size),p=p,size=Nsample)]).T
    meanshift = MeanShift(bandwidth,n_jobs=n_jobs)
    meanshift.fit(np.array([xd,yd]).T, )
    return meanshift.cluster_centers_.T


def createNNDistanceMatrix(dims):
    from itertools import product
    from scipy.sparse import dok_matrix
    DistMatrix = dok_matrix((np.prod(dims),np.prod(dims)), dtype=np.int0)

    for i,j in product(range(dims[0]),range(dims[1])):
        for di,dj in product([-1,0,1],[-1,0,1]):
            i1 = i+di
            j1 = j+dj
            if i1<0 or j1<0 or i1>=dims[0] or j1>=dims[1]: continue
    #         d = ((np.array(p0)-np.array(p1))**2).sum()**.5
    #         if d<2:
            DistMatrix[i*dims[1]+j,i1*dims[1]+j1] = 1
    return DistMatrix#.tocsr()

def getEventBoundaries(event,spaceDims):
    """in shape: (tmin,tmax),(imin,imax),(jmin,jmax)"""
    nc = np.prod(spaceDims)
    coords = np.array([(coord//nc,(coord%nc)//spaceDims[1],(coord%nc)%spaceDims[1]) for coord in event])
    return np.array([coords.min(0), coords.max(0)+1]).T

def getEventIndices(event,spaceDims):
    """in shape: (tmin,tmax),(imin,imax),(jmin,jmax)"""
    nc = np.prod(spaceDims)
    coords = np.array([(coord//nc,(coord%nc)//spaceDims[1],(coord%nc)%spaceDims[1]) for coord in event])
    return coords

def getEvents(tfMatrix_,DistMatrix_=None, offset=(0,0,0)):
    from scipy.sparse import dok_matrix
    import networkx as nx
    allGraph = dok_matrix((tfMatrix_.size,tfMatrix_.size),dtype=np.int0)
    if DistMatrix_ is None:
        DistMatrix_ = createNNDistanceMatrix(tfMatrix_.shape[1:])
    nc = DistMatrix_.shape[0]
    
    frames = np.where(tfMatrix_.any(axis=(1,2)))[0]
#     frames = frames[:-1][np.diff(frames)==1]

    for iframe in frames:
        curFrame = tfMatrix_[iframe]
        if not curFrame.any(): continue
        hereActiveCells = np.where(curFrame.flatten())[0]
        for ic in hereActiveCells:
            nodeID = iframe*nc+ic
            for diframe in [0,1]:
                iframe1 = iframe+diframe
                if iframe1<0: continue
                if iframe1>=len(tfMatrix_): continue
                for el in DistMatrix_[ic].keys():
                    jc = el[1]
                    if tfMatrix_[iframe1].flatten()[jc]:
                        allGraph[iframe1*nc+jc, nodeID] = 1
                        allGraph[nodeID, iframe1*nc+jc] = 1

    G = nx.from_scipy_sparse_matrix(allGraph)
    G.remove_nodes_from(list(nx.isolates(G)))
    
#     events = [cmp for cmp in nx.connected_components(G) 
#               if #len(cmp)>1 or \
#               allGraph[tuple(cmp)*2]
#              ]
#     events = []
    events = [getEventIndices(ev,tfMatrix_.shape[1:]) for ev in nx.connected_components(G)]
    if any(offset):
        events = [ev+np.repeat([offset],len(ev),axis=0) for ev in events]
              
    return events

def blurInTime(original_, k, s=-1):
    from islets import cmovie
    from cv2 import GaussianBlur
    out = original_.copy()
    out = GaussianBlur(
        out.reshape((len(out),-1)),
#         ksize=(k,1),sigmaX=s, sigmaY=0
        ksize=(1,k),sigmaX=0, sigmaY=s
    ).reshape(out.shape)
    out = cmovie(out)
    return out

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

def agglomerativeLinkage(X,distMatrix):
    Nleaves = len(X)
    idx = list(range(Nleaves))
    Lkg = [[i] for i in range(Nleaves)]
#     dstDict = {}
    Z = []
    while True:
    # for jj in range(11):
        currX = np.array([np.mean([X[el] for el in Lkg[k]],axis=0) for k in idx])
        dsts = distMatrix(currX)
    #     plt.imshow(dsts)
        minDist = np.min(dsts[np.triu_indices(len(dsts),1)])
        i,j = [(i,j) for i,j in np.array(np.where(dsts==minDist)).T if i!=j][0]
        k = len(Lkg)
        Lkg += [Lkg[idx[i]]+Lkg[idx[j]]]
    #     print (idx[i],idx[j], minDist, len(Lkg[k]))
        Z += [[idx[i],idx[j], minDist, len(Lkg[k])]]
#         dstDict[k] = minDist
        if len(Lkg[k])==Nleaves:
            break
        del idx[j], idx[i]
        idx = idx + [k]
    return np.array(Z)

def eventIdx2Rois(ks, evIndices, dims):
    from cv2 import Laplacian, CV_16F
    M = np.zeros(dims)
    for k in ks:
        x,y= evIndices[k]
        # y,x = np.array(evi).T[1:]
        M[x,y] = 1
    return Laplacian(M.T.astype(float),CV_16F,)>0

def get_cluster_color_classes(den, labels = None):
    from collections import defaultdict
    if labels is None:
        labels = np.arange(len(den['leaves']))
    cluster_idxs = defaultdict(list)
    for c, pi in zip(den['color_list'], den['icoord']):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(labels[int(i)])
    return cluster_idxs

# code
def putRois(ks, evIndices, dims):
    M = np.zeros(dims)
    for k in ks:
        x,y= evIndices[k]
        # y,x = np.array(evi).T[1:]
        M[x,y] = 1
    return M # Laplacian(M.T.astype(float),CV_64F,)>0


def clusterCutAndPlot(Xdata,
                      mode = "simple",
                      function=None,
                      functionArgs=None,
                      threshold=4,
                      criterion="maxclust",
                      labels=None,
                      imshow_kw = dict(), #vmin=-vmax,vmax=vmax, cmap="bwr"
                      showClusterBoundaries=False,
                      ):
    if mode == "simple":
        return clusterCutAndPlot(Xdata,
                      mode=None,
                      function=linkage,
                      functionArgs=("complete","correlation"),
                      threshold=threshold,
                      criterion=criterion,
                      labels=labels,
                      imshow_kw = imshow_kw,
                      showClusterBoundaries=showClusterBoundaries,
                      )
    if mode == "agglomerative":
        def distM(xi): return 1-np.corrcoef(xi)
        return clusterCutAndPlot(Xdata,
                      mode=None,
                      function=agglomerativeLinkage,
                      functionArgs=(distM, ),
                      threshold=threshold,
                      criterion=criterion,
                      labels=labels,
                      imshow_kw = imshow_kw,
                      showClusterBoundaries=showClusterBoundaries,
                      )

    from collections import OrderedDict
    if labels is None:
        labels = np.arange(len(Xdata))
    assert len(labels)==len(Xdata)
    Z = function(Xdata,*functionArgs)
    clusterIds = fcluster(Z, threshold, criterion=criterion)
    if criterion=="maxclust":
        threshold = Z[-threshold:-threshold+2,2].mean()
    fig, axs = plt.subplots(1,3,figsize=(18, 8), gridspec_kw={'width_ratios': [1, 1, 2.2]})

    dnd = dendrogram(
        Z,
        orientation="left",
        color_threshold=threshold,
        above_threshold_color="grey",
        leaf_rotation=0,  # rotates the x axis labels
        leaf_font_size=10,  # font size for the x axis labels
        ax = axs[0]
        )
    
    color_classes = get_cluster_color_classes(dnd)
    color_classes = OrderedDict([(c,np.array([labels[j] for j in color_classes[c]])) for c in color_classes])
    axs[0].set_yticks(5+10*np.arange(len(labels)))
    axs[0].set_yticklabels(labels[dnd['leaves']]);
    xt = axs[0].get_xticks()
    axs[0].set_xticklabels(["%3.1f"%n for n in 1-xt])

    ax = axs[1]
    ax.imshow((Xdata)[dnd["leaves"]], origin="lower", **imshow_kw, interpolation="nearest")
    ax.set_aspect("auto")
    axs[0].axvline(threshold,color="red",ls="--")
    if showClusterBoundaries:
        for j in np.where(np.diff(clusterIds[dnd["leaves"]]))[0]:
            axs[1].axhline(j+.5,color="grey")
    axs[1].set_yticks(np.arange(len(labels)))
    axs[1].set_yticklabels(labels[dnd["leaves"]]);

#     try:
#         for c in color_classes.keys():
#             for jEvent in color_classes[c]:
#                 x,y = eventDF["indices"][jEvent]
#                 axs[2].plot(y,x,"s",c=c,ms=.7*500/subdims[0],alpha = .5)
#                 axs[2].text(y.mean(),x.mean(),jEvent,color=c,
#                         va="center",ha="center",
#                         # bbox=dict(facecolor='w', alpha=0.7,edgecolor="none")
#                        )

#         axs[2].set_aspect("equal")
#         axs[2].set_xlim(-.5,subdims[1]-.5)
#         axs[2].set_ylim(subdims[0]-.5,-.5)
#     except:
#         axs[2].remove()
    
    return {"dendrogram":dnd,
            "clusters":color_classes,
            "labels":labels,
            "axes": axs
           }


@jit
def getEvents_(tfMatrix_,DistMatrix_):
    from scipy.sparse import dok_matrix
    allGraph = dok_matrix((tfMatrix_.size,tfMatrix_.size),dtype=np.int0)
    nc = DistMatrix_.shape[0]

    for iframe in prange(len(tfMatrix_)):
        curFrame = tfMatrix_[iframe]
        if not curFrame.any(): continue
        hereActiveCells = np.where(curFrame.flatten())[0]
        for ic in hereActiveCells:
            nodeID = iframe*nc+ic
            for diframe in [0,1]:
                iframe1 = iframe+diframe
                if iframe1<0: continue
                if iframe1>=len(tfMatrix_): continue
                for el in DistMatrix_[ic].keys():
                    jc = el[1]
                    if tfMatrix_[iframe1].flatten()[jc]:
                        allGraph[iframe1*nc+jc, nodeID] = 1
                        allGraph[nodeID, iframe1*nc+jc] = 1
    return allGraph


def crawlDict(image, crawl_th=0, diag=False, processes=10, excludePixels=None, verbose=False):
    global iterf
    if excludePixels is None:
        excludePixels = []
    if verbose:
        print (f"entering crawling dict with {len(excludePixels)} pixels excluded.")
    if np.isfinite(crawl_th):
        excludePixels += list(map(tuple,np.array(np.where(image<crawl_th)).T))
    if verbose:
        print (f"Crawling the image with {len(excludePixels)} pixels excluded.")

    ijs_ = [(i,j) for i,j in product(range(image.shape[0]),range(image.shape[1])) if (i,j) not in excludePixels]
    def iterf(ij):
        return climb((ij[0],ij[1]), image, diag=diag,
                     #excludePixels=excludePixels
                    )
    R_ = multi_map(iterf,ijs_, processes=processes)
    A_ = [ij+r for ij,r in zip(ijs_,R_) if r not in excludePixels]
    A_ = [el for el in A_ if el[-1] is not None]
    B_ = OrderedDict()
    for (i0,j0,i1,j1) in A_:
        if (i1,j1) not in B_:
            B_[(i1,j1)] = []
        B_[(i1,j1)] += [(i0,j0)]
    return B_


def crawl_dict_via_graph(image, validPixelsImage, diag=False, offset=None ):
    gph = nx.DiGraph()
    if offset is None:
        offset = np.zeros(image.ndim, dtype=int)
    for start in zip(*np.where(validPixelsImage)):
        end = climb(start, image, Niter=1, diag=diag)
        gph.add_edge(start, end)
    c = dict()
    for attrs, nodes in zip(nx.algorithms.attracting_components(gph), nx.connected_components(gph.to_undirected())):
        assert len(attrs)==1 # why would it be different? For cyclic graphs perhaps.
        peak = attrs.pop()
        peak = tuple(peak[j]+offset[j] for j in range(image.ndim))
        c[peak] = [ tuple(coords+offset) for coords in nodes]
    return c


def climb(x, blurredWeights,
          diag=True,
          excludePixels=None,
          Niter=1
         ):
    if excludePixels is None:
        excludePixels = []
    dims = np.array(blurredWeights.shape)
    x = tuple(x)
    # x = x+(blurredWeights[tuple(x)],)
    ndim = len(dims)
    v = blurredWeights[x]
    for k in range(Niter):
        vs = []
        xs = []
        for dcoords in product(*([[-1,0,1]]*ndim)):
            if not diag:
                if np.prod(dcoords)!=0: continue
            coords = x + np.asanyarray(dcoords)
            if any(coords<0):
                continue
            if any(coords>=dims):
                continue
            if tuple(coords) in excludePixels:
                continue
            vs += [blurredWeights[tuple(coords)]]
            xs += [coords]
        v1 = max(vs)
        if v1>v:
            x = tuple(xs[np.argmax(vs)])
        else:
            break
    return x