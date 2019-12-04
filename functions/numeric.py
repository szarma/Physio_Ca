import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import distributions as dst
from scipy.optimize import curve_fit#,minimize,basinhopping
from numba import jit,prange

def decay(time,top,bottom,rate):
#     top,bottom,rate = p
    return bottom+(top-bottom)*np.exp(-time*rate)
from scipy.optimize import curve_fit

def guessDecayPars(y):
    b0 = np.nanmin(y)
    y  = np.log(y-b0)
    r0 = np.diff(y)
    r0 = r0[np.isfinite(r0)]
    r0 = -np.mean(r0)
    t0 = np.nanpercentile(y,99)
    p0 = (np.exp(t0)+b0,b0,r0)
    return p0

def decayfit(x,Ntrials=None):
    if Ntrials is None:
        lx = 10
    else:
        lx = Ntrials
    nx = len(x)//10
    TT = np.arange(len(x))
    tt = TT.copy()
    for j_ in range(lx):
        try:
            p0 = guessDecayPars(x)
            ff = np.isfinite(x)
            popt = curve_fit(decay,tt[ff],x[ff],p0=p0)[0]
            expDecay = decay(TT,*popt)
            return expDecay
        except:
            x = x[:-nx]
            tt = tt[:-nx]
    return p0

@jit 
def percFilter(x_,perc,filterSize):
    if filterSize%2==0:
        raise ValueError("filter size needs to be odd number")
    delta = filterSize//2
    out = np.zeros_like(x_)
    # x_ = np.hstack([[x_[0]]*delta,x_,[x_[-1]]*delta])
    x_ = np.hstack([x_[:delta][::-1], x_, x_[-delta:][::-1]])
#     print (x_.shape, out.shape)
    for i in prange(len(out)):
        out[i] = np.percentile(x_[i:i+filterSize],perc)
    return out
@jit 
def runningAverage(x_,filterSize):
    if filterSize%2==0:
        raise ValueError("filter size needs to be odd number")
    delta = filterSize//2
    out = np.zeros_like(x_)
    # x_ = np.hstack([[x_[0]]*delta,x_,[x_[-1]]*delta])
    x_ = np.hstack([x_[:delta][::-1], x_, x_[-delta:][::-1]])
#     print (x_.shape, out.shape)
    for i in prange(len(out)):
        out[i] = np.nanmean(x_[i:i+filterSize])
    return out
@jit 
def runningStd(x_,filterSize):
    if filterSize%2==0:
        raise ValueError("filter size needs to be odd number")
    delta = filterSize//2
    out = np.zeros_like(x_)
    # x_ = np.hstack([[x_[0]]*delta,x_,[x_[-1]]*delta])
    x_ = np.hstack([x_[:delta][::-1], x_, x_[-delta:][::-1]])
#     print (x_.shape, out.shape)
    for i in prange(len(out)):
        out[i] = np.nanstd(x_[i:i+filterSize])
    return out

def lowPass(x, wIron, wAvg=None, perc=5,npass = 1):
    ironedx = percFilter(x, perc, wIron)
    if wAvg is not None:
        ironedx = runningAverage(ironedx,wAvg)
    for i in range(npass-1):
        ironedx += lowPass(x-ironedx,wIron,wAvg,perc,)
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

def getEvents(tfMatrix_,DistMatrix_):
    from scipy.sparse import dok_matrix
    import networkx as nx
    allGraph = dok_matrix((tfMatrix_.size,tfMatrix_.size),dtype=np.int0)
    nc = DistMatrix_.shape[0]

    for iframe in range(len(tfMatrix_)):
        hereActiveCells = np.where(tfMatrix_[iframe].flatten())[0]
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
    #         break
    #     break

    G = nx.from_scipy_sparse_matrix(allGraph)

    events = [cmp for cmp in nx.connected_components(G) 
              if len(cmp)>1 or allGraph[tuple(cmp)*2]]
    return events


def blurInTime(original_, k, s):
    from caiman import movie
    from opencv import GaussianBlur
    out = original_.copy()
    out = GaussianBlur(
        out.reshape((len(out),-1)),ksize=(k,1),sigmaX=s, sigmaY=0
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
    return Z

def eventIdx2Rois(ks, evIndices, dims):
    from opencv import Laplacian,CV_16F
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

def clusterCutAndPlot(Xdata,
                      function=linkage,
                      functionArgs=("complete","correlation"),
                      threshold=4,
                      criterion="maxclust",
                      labels=None,
                      imshow_kw = dict(), #vmin=-vmax,vmax=vmax, cmap="bwr"
                      showClusterBoundaries=False,
                      ):

    if labels is not None:
        labels = np.arange(len(Xdata))
    assert len(labels)==len(Xdata)
    Z = function(Xdata,*functionArgs)
    clusterIds = fcluster(Z, threshold, criterion=criterion)
    if criterion=="maxclust":
        threshold = Z[-threshold,2]
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
    axs[0].set_yticks(5+10*np.arange(len(labels)))
    axs[0].set_yticklabels(labels[dnd['leaves']]);

    ax = axs[1]
    ax.imshow((Xdata)[dnd["leaves"]], origin="bottom", **imshow_kw)
    ax.set_aspect("auto")
    axs[0].axvline(threshold,color="darkgrey")
    if showClusterBoundaries:
        for j in np.where(np.diff(clusterIds[dnd["leaves"]]))[0]:
            axs[1].axhline(j+.5,color="grey")
    axs[1].set_yticks(np.arange(len(labels)))
    axs[1].set_yticklabels(labels[dnd["leaves"]]);


    for c in color_classes.keys():
    #     if c!="grey": continue
        y,x = np.where(getLines(labels[np.array(color_classes[c])], evIndices))
    #     y,x = np.where(getLines(np.array(color_classes[c]), evIndices))
        plt.plot(x,y,".",ms=.7,ls="none",c=c)
        for k in labels[np.array(color_classes[c])]:
            evi = evIndices[k]
            y,x = np.array(evi).T[1:]
    #         c = axs[2].plot(x,y,".",c=c)[0].get_color()
            axs[2].text(x.mean(),y.mean(),k,color=c,
    #                     va="center",ha="center",
    #                     bbox=dict(facecolor='w', alpha=0.7,edgecolor="none")
                       )
    #     break
    # axs[2].set_xlim(-.5,submovie.shape[1])
    # axs[2].set_ylim(-.5,submovie.shape[2])
    axs[2].set_aspect("equal")