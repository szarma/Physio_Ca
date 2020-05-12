import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import distributions as dst
from scipy.optimize import curve_fit#,minimize,basinhopping
from numba import jit,prange

def get_sep_th(x_,ax=None,plot=False,thMax=None):
    from scipy.stats import gaussian_kde
    from Regions import crawlDict
    if ax is None:
        ax = plt.subplot(111)
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
    peaks = np.array(list(crawlDict(gkde_vals.reshape(-1,1)).keys())).T[0]
    twoHighestPeaks = sorted(np.argsort(-h[peaks])[:2])
    if len(twoHighestPeaks)==1:
        return -np.inf
    lowerPeak, higherPeak = peaks[twoHighestPeaks]
    bincenters = bincenters[lowerPeak:higherPeak]
    if len(bincenters)<3:
        if not plot: plt.close()
        return x_.min()
    gkde_vals = gkde.evaluate(bincenters)
#     if plot:
    ax.plot(bincenters,gkde_vals)
    sinks = np.array(list(crawlDict(-gkde_vals.reshape(-1,1)+gkde_vals.max()+1,).keys())).T[0]
    sinks = sorted(sinks,key=lambda xi: gkde_vals[xi])
    th = bincenters[sinks[0]]
#     if plot:
        # ax.plot(bincenters,gkde_vals)
    ax.axvline(th,color="r",ls="--")
    if not plot: plt.close()
    return th


def decay(time,top,bottom,rate):
#     top,bottom,rate = p
    return bottom+(top-bottom)*np.exp(-time*rate)

def guessDecayPars(y):
    b0 = np.nanpercentile(y,1)
    y  = np.log(y-b0)
#     r0 = np.diff(y)
#     r0 = r0[np.isfinite(r0)]
#     r0 = -np.mean(r0)
    r0 = .1/len(y)#y[1]-y[-2]/len(y)
    t0 = np.nanpercentile(y,99)
    p0 = (np.exp(t0)+b0,b0,r0)
    return p0

def mydebleach(x_):
    try:
        out = decayfit(x_)
    except:
        p = np.polyfit(range(len(x_)),x_,1)
        out = p[1]+p[0]*np.arange(len(x_))
    return out

def decayfit(x,Ntrials=None, outPars=False):
    if Ntrials is None:
        lx = 10
    else:
        lx = Ntrials
    nx = len(x)//(lx+5)
    TT = np.arange(len(x))
    tt = TT.copy()
    for j_ in range(lx):
        try:
            p0 = guessDecayPars(x)
            ff = np.isfinite(x)
            popt = curve_fit(decay,tt[ff],x[ff],p0=p0)[0]
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


# def decayfit(x,Ntrials=None):
#     if Ntrials is None:
#         lx = 10
#     else:
#         lx = Ntrials
#     nx = len(x)//10
#     TT = np.arange(len(x))
#     tt = TT.copy()
#     for j_ in range(lx):
#         try:
#             p0 = guessDecayPars(x)
#             ff = np.isfinite(x)
#             popt = curve_fit(decay,tt[ff],x[ff],p0=p0)[0]
#             expDecay = decay(TT,*popt)
#             return expDecay
#         except:
#             x = x[:-nx]
#             tt = tt[:-nx]
#     return p0

@jit 
def percFilter(x_,perc,filterSize):
    if filterSize%2==0:
        raise ValueError("filter size needs to be odd number")
    delta = filterSize//2
    out = np.zeros_like(x_)
    x_ = np.concatenate((x_[:delta][::-1], x_, x_[-delta:][::-1]))
    for i in prange(len(out)):
        out[i] = np.nanpercentile(x_[i:i+filterSize],perc,axis=0)
    return out
@jit 
def runningMode(x_, filterSize, boundary="reflective"):
    if filterSize%2==0:
        raise ValueError("filter size needs to be odd number")
    delta = filterSize//2
    out = np.zeros_like(x_)
    allowedBoundaries = ["reflective","equal"]
    if boundary not in allowedBoundaries:
        raise ValueError(f"boundary {boundary} not recognized. Allowed values are: "+repr(allowedBoundaries))
    if boundary=="reflective":
        x_ = np.concatenate((x_[:delta][::-1], x_, x_[-delta:][::-1]))
    if boundary=="equal":
        x_ = np.concatenate(([x_[0]]*delta, x_, [x_[-1]]*delta))
    for i in range(len(out)):
        c,v = np.histogram(x_[i:i+filterSize], nbins)
        out[i] = v[np.argmax(c)]
    return out

@jit 
def runningAverage(x_,filterSize):
    if filterSize%2==0:
        raise ValueError("filter size needs to be odd number")
    delta = filterSize//2
    out = np.zeros_like(x_)
    x_ = np.concatenate((x_[:delta][::-1], x_, x_[-delta:][::-1]))
    for i in range(len(out)):
        out[i] = np.nanmean(x_[i:i+filterSize],axis=0)
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

# def getEvents(tfMatrix_,DistMatrix_):
#     import networkx as nx
#     G = nx.from_scipy_sparse_matrix(getEvents_(tfMatrix_,DistMatrix_))
#     events = [cmp for cmp in nx.connected_components(G) 
#               if len(cmp)>1 or allGraph[tuple(cmp)*2]]
#     return events


def blurInTime(original_, k, s=-1):
    from caiman import movie
    from cv2 import GaussianBlur
    out = original_.copy()
    out = GaussianBlur(
        out.reshape((len(out),-1)),
#         ksize=(k,1),sigmaX=s, sigmaY=0
        ksize=(1,k),sigmaX=0, sigmaY=s
    ).reshape(out.shape)
    out = movie(out)
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
    from cv2 import Laplacian,CV_16F
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
    from cv2 import Laplacian,CV_64F
    M = np.zeros(dims)
    for k in ks:
        x,y= evIndices[k]
        # y,x = np.array(evi).T[1:]
        M[x,y] = 1
    return M # Laplacian(M.T.astype(float),CV_64F,)>0

def getEdges(image,rescale=0):
    from scipy.ndimage import binary_fill_holes
    from skimage.feature import canny
    from skimage.morphology import remove_small_objects
    edges = canny(image)
    fill_coins = binary_fill_holes(edges)
    roiEdges = ndi.label(remove_small_objects(fill_coins, 21))[0]
    return mark_boundaries(image*rescale, roiEdges)

def showEventEdges(ks, evIndices, dims, resc=0):
    imm = putRois(ks,evIndices,dims)
    return getEdges(imm)

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
    ax.imshow((Xdata)[dnd["leaves"]], origin="bottom", **imshow_kw)
    ax.set_aspect("auto")
    axs[0].axvline(threshold,color="red",ls="--")
    if showClusterBoundaries:
        for j in np.where(np.diff(clusterIds[dnd["leaves"]]))[0]:
            axs[1].axhline(j+.5,color="grey")
    axs[1].set_yticks(np.arange(len(labels)))
    axs[1].set_yticklabels(labels[dnd["leaves"]]);

    try:
        for c in color_classes.keys():
            for jEvent in color_classes[c]:
                x,y = eventDF["indices"][jEvent]
                axs[2].plot(y,x,"s",c=c,ms=.7*500/subdims[0],alpha = .5)
                axs[2].text(y.mean(),x.mean(),jEvent,color=c,
                        va="center",ha="center",
                        # bbox=dict(facecolor='w', alpha=0.7,edgecolor="none")
                       )

        axs[2].set_aspect("equal")
        axs[2].set_xlim(-.5,subdims[1]-.5)
        axs[2].set_ylim(subdims[0]-.5,-.5)
    except:
        axs[2].remove()
    
    return {"dendrogram":dnd,
            "clusters":color_classes,
            "labels":labels}


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