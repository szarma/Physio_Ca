import numpy as np
import pandas as pd
from itertools import product
from collections import OrderedDict
import matplotlib.pyplot as plt
import networkx as nx
from .numeric import rebin
from .utils import multi_map
import plotly.graph_objects as go
from .Regions import Regions as Regions0
from matplotlib._color_data import TABLEAU_COLORS, CSS4_COLORS
from general_functions import getCircularKernel
import pickle
from matplotlib.colors import LogNorm
import os
from sys import exc_info
# MYCOLORS = OrderedDict(TABLEAU_COLORS)
# del MYCOLORS["tab:gray"]
# ks = ["lime",
#     "orangered",
#     "yellowgreen",
#       "lawngreen",
#     # "mediumspringgreen",
#     "aquamarine"]#+list(MYCOLORS.values())
# MYCOLORS.update({k:CSS4_COLORS[k] for k in ks})
# MYCOLORS = list(MYCOLORS.values())
from plotly_express import colors as plc
MYCOLORS = plc.qualitative.Plotly
# MYCOLORS = ["darkred"]

def load_regions(path,
                 mergeDist=1,
                 mergeSizeTh=10,
                 plot=False,
                 verbose=False,
                 calcInterest=True
                ):
    with open(path,"rb") as f:
        regions = pickle.load(f)
    try:
        regions.update()
        pickleDir = os.path.split(path)[0]
        if "examine3" not in dir(regions):
            regions = Regions(regions)
        try:
            protocolFile = os.path.join(pickleDir, [f for f in os.listdir(pickleDir) if "protocol" in f][0])
            regions.import_protocol(protocolFile)
        except:
            pass
        regions.pathToPickle = path
        # regions.sortInOrder()
        regions.detrend_traces()
        regions.infer_gain(plot=plot)
        if plot:
            plt.figure(figsize=(7*20,6))

        ia = 1
        while True:
            size_th = np.percentile(regions.df["size"].values, mergeSizeTh)
            df = getPeak2BoundaryDF(regions.df)
            df = df.query(f"dist<={mergeDist} and size_from<={size_th}")
            if len(df):
                if plot:
                    ax = plt.subplot(1,20,ia)
                    ax.imshow(regions.statImages[regions.mode], cmap="Greys", norm=LogNorm())
                    xl = ax.get_xlim()
                    yl = ax.get_ylim()
                else:
                    ax = None
                suggestGraph = getGraph_of_ROIs_to_Merge(df.iloc[:,:2], regions, plot=plot,ax=ax)
                if plot:
                    ax.set_xlim(xl)
                    ax.set_ylim(yl)
                mergeBasedOnGraph(suggestGraph, regions, verbose=verbose)
            else:
                # print ("No more suggestions.")
                break
            ia += 1
        if plot:
            plt.tight_layout()
        if calcInterest:
            regions.calc_interest()
    except:
        print ("encountered error:", exc_info())
    
    return regions

def crawlDict_restr(image, pixels, diag=False, processes=10, verbose=False):
    global iterf
    def iterf(ij):
        return climb((ij[0],ij[1]), image, diag=diag,)
    R_ = multi_map(iterf,pixels, processes=processes)
    A_ = [ij+r for ij,r in zip(pixels,R_) if r in pixels]
    A_ = [el for el in A_ if el[-1] is not None]
    B_ = OrderedDict()
    for (i0,j0,i1,j1) in A_:
        if (i1,j1) not in B_:
            B_[(i1,j1)] = []
        B_[(i1,j1)] += [(i0,j0)]
    return B_

def climb(x,blurredWeights,
          diag=True,
          excludePixels=[]
         ):
    dims = blurredWeights.shape
    # x = (60,60)
    x = x+(blurredWeights[x[0],x[1]],)
    xs = [x]
    for i in range(1000):
        vs = []
        for di,dj in product([-1,0,1],[-1,0,1]):
            if not diag:
                if di*dj!=0: continue
            i,j = x[0]+di,x[1]+dj
            if i<0 or i>=dims[0] or j<0 or j>=dims[1]:
                continue
            if (i,j) in excludePixels:
                continue
            vs += [(i,j,blurredWeights[i,j])]
        x1 = vs[np.argmax(vs,axis=0)[-1]]
        dx = x1[-1]-x[-1]
        if dx<=0:
            break
        else:
            x = x1
            xs += [x]
    return x[:2]

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


class Regions:
    def __init__(self, movie_,
                 diag=True,
                 debleach=False,
                 gSig_filt=None,
                 mode="highperc+mean",
                 full=True,
                 img_th=0.01,
                 FrameRange=None,
                 processes=7,
                 excludePixels=None,
                 verbose=False,
                 use_restricted=None
                ):
        if isinstance(movie_, Regions0) or isinstance(movie_, Regions):
            for k in movie_.__dict__.keys():
                if verbose:
                    print ("Initiating from another Regions object.")
                setattr(self, k, movie_.__dict__[k])
            return None
            
        self.mode = mode
        if isinstance(movie_, (np.ndarray,)):
            if len(movie_.shape)==2:
                if verbose:
                    print ("Initiating from an image, with a mode", mode)
                self.statImages = {mode:movie_}
            if len(movie_.shape)==3:
                if verbose:
                    print ("Initiating from a movie.")
                self.movie = movie_
                time = np.arange(len(movie_))/movie_.fr
                if FrameRange is None:
                    FrameRange = [0, len(movie_)]
                i0, ie = FrameRange
                self.FrameRange = FrameRange
                self.statImages = getStatImages(movie_[i0:ie], debleach=debleach)
                self.time = time[i0:ie]
                self.Freq = movie_.fr
                self.showTime = {}
        elif isinstance(movie_, dict):
            akey = next(iter(movie_))
            if isinstance(akey,tuple):
                if verbose:
                    print ("Initiating from a crawling dictionary.")
                self.df = pd.DataFrame(OrderedDict([
                    ("peak",  list(movie_.keys())),
                    ("pixels",list(movie_.values()))
                ]))
                del self.mode
            elif isinstance(akey,str):
                if verbose:
                    print ("Initiating from a dictionary assumed to be a dictionary of image stats.")
                self.statImages = movie_
            else:
                raise ValueError("Initializing Regions from a dictionary is only supported for a dictionary of images representing movie statistics, or a pixel crawling dictionary.")
        else:
            raise ValueError("Regions can initialize either from a movie, or an image, or a dictionary. You supplied %s"%str(type(movie_)))
        
        if full and not hasattr(self,"df"):
            self.constructRois(mode=mode, img_th=img_th, diag=diag, gSig_filt=gSig_filt, processes=processes, excludePixels=excludePixels, verbose=verbose, use_restricted=use_restricted)
            
    def constructRois(self, mode, img_th=0, diag=True, gSig_filt=None, processes=5,excludePixels=None, verbose=False,use_restricted=False):
        from .numeric import robust_max
        if mode=="custom":
            image0=self.statImages[mode]
        else:
            k0 = next(iter(self.statImages))
            tmp = np.zeros_like(self.statImages[k0])
            norm = []
            for submode in mode.split("+"):
                im = self.statImages[submode].copy()
                norm += [robust_max(im)]
                im = im/norm[-1]
                tmp += im
            image0 = tmp/len(norm)*np.mean(norm)
            self.statImages[mode] = image0
#             toMin=image0<img_th
            
        from cv2 import GaussianBlur,dilate
        if gSig_filt is None:
            image = image0
#             excludePixels = None
            dks = 3
        else:
            if type(gSig_filt)==int:
                gSig_filt = [gSig_filt]
            dks = max(3,(max(gSig_filt))//2*2+1)
            image = []
            for gSf in gSig_filt:
#                 tmp = high_pass_filter_space(image0,(gSf,)*2)
                tmp  = GaussianBlur(image0,(gSf//2*2+1,)*2,-1)
                tmp -= GaussianBlur(image0,(gSf*2+1,)*2,-1)
                tmp *= gSf
                image += [tmp/robust_max(tmp)]
            image = np.mean(image,axis=0)
        if not use_restricted:
            if excludePixels is None:
                excludePixels = []
        
        if verbose:
            print ("dilating valid pixels by", dks)
        kernel = getCircularKernel(dks)
        ok = dilate((image>img_th).astype(np.uint8), kernel)
        # ok = ok*(self.statImages[self.mode]>0).astype(np.uint8)
        ok = ok.astype(bool)
        self.filterSize = gSig_filt
        self.image = image
        if use_restricted:
            includePixels = list(map(tuple, np.array(np.where(ok)).T))
            if verbose:
                print(f"initiating the cralwing dict on {len(includePixels)} (%.1f%%) pixels only."%(100*len(includePixels)/ok.size))
            B_ = crawlDict_restr(image,
                           includePixels,
                           diag=diag,
                           processes=processes,
                           verbose=verbose
                          )
        else:
            excludePixels += list(map(tuple, np.array(np.where(~ok)).T))
            if verbose:
                print(f"initiating the cralwing dict with {len(excludePixels)} pixels excluded.")
            B_ = crawlDict(image,
                           crawl_th=-np.inf,
                           diag=diag,
                           processes=processes,
                           excludePixels=excludePixels,
                           verbose=verbose
                          )
        if diag:
#             try:
                from .utils import split_unconnected_rois
                B_ = split_unconnected_rois(B_, self.image)
#             except:
#                 print ("Cannot initialize with diagonal crawl. Reverting to diag=False")
#                 diag = False
#                 B_ = crawlDict(image,image.min(),diag=diag,min_gradient=min_gradient, processes=processes)
        
        
        self.df = pd.DataFrame(OrderedDict([
            ("peak",  list(B_.keys())),
            ("pixels",list(B_.values()))
        ]))
#         try:
#             blurKsize = min(self.filterSize)//2*2+1
#             blurKsize = max(3,blurKsize)
#             slightly_blurred_image = GaussianBlur(self.statImages[self.mode],(blurKsize,)*2,-1)
#             self.reassign_peaks(slightly_blurred_image+self.statImages[self.mode])
#         except:
        self.df["peakValue"] = [image[p] for p in B_]
        self.ExcludePixels = excludePixels
        self.update()
    
    def reassign_peaks(self, image, write=True):
        newPeaks = []
        newValues = []
        for i in self.df.index:
            pxs = self.df.pixels[i]
            image_values = [image[px] for px in pxs]
            jchoose = image_values.index(max(image_values))
            newPeaks  += [pxs[jchoose]]
            newValues += [image[pxs[jchoose]]]
        if write:
            self.df["peak"] = newPeaks
            self.df["peakValue"] = newValues
            self.update()
        else:
            return newPeaks, newValues
    
    def get_fov_trace(self, showFreq = 2, pixels=None):
        from .numeric import mydebleach
        from physio_def_1 import rebin
        i0, ie = self.FrameRange
#         n = int(self.movie.fr/showFreq)
        n = int(self.Freq/showFreq)
        if n==0: n=1
        x = rebin(np.arange(i0,ie)/self.Freq,n)
        try:
            y = np.sum([self.df.loc[i,"trace"]*self.df.loc[i,"size"] for i in self.df.index],axis=0)/self.df["size"].sum()
            if n>1:
                y = rebin(y,n)
        except:
            if pixels is None:
                y = self.movie[i0:ie:n].mean(axis=(1,2))
            else:
                y = self.movie[(slice(i0,ie,n),)+pixels].mean(axis=1)
        ydbl = mydebleach(y)
        self.fov_trace = {
                "time": x,
                "raw": y,
                "trend":ydbl
            }

    def update(self, movie_=None):
        self.df["size"] = self.df["pixels"].apply(len)
        self.df["interest"] = [np.sum([self.image[px[0],px[1]] for px in pxs]) for pxs in self.df["pixels"]]
        self.calcEdgeIds()
        self.calcEdges()
        self.df["boundary"] = [edges2nodes(self.df["edges"][j]) for j in self.df.index]
        self.calcNNmap()
        if movie_ is not None:
            self.calcTraces(movie_)
            self.movie = movie_
            self.Freq  = movie_.fr
        else:
            try:
                self.calcTraces()
            except:
                pass
    
    def calcEdgeIds(self):
        dround = np.vstack([(-1,-1),(-1, 1),( 1, 1),( 1,-1),(-1,-1)])
        dedges = []
        for el in zip(dround[:-1],dround[1:]):
            el = np.array(el)+1
            el = np.vstack(sorted(el,key=np.linalg.norm))-1
            dedges += [el*.5]
        dedges = np.stack(dedges)
        edgeID = OrderedDict()
        for k,pixelSet in zip(self.df.peak,self.df.pixels):
            for x,y in pixelSet:
                edges = dedges.copy()
                edges[...,0] += x
                edges[...,1] += y
                for edge in edges.reshape((-1,4)):
                    edge = tuple(edge)
                    if edge not in edgeID:
                        edgeID[edge] = []
                    edgeID[edge] += [k]
                    if len(edgeID[edge])==2 and tuple(edgeID[edge][0])==tuple(edgeID[edge][1]):
                        del edgeID[edge]
        self.edgeIDs = edgeID
    
    def calcEdges(self):
        invEdgeID = OrderedDict()
#         if "edgeIDs" not in locals():
        self.calcEdgeIds()
        for k in self.edgeIDs:
            for p in self.edgeIDs[k]:
                if p not in invEdgeID:
                    invEdgeID[p] = []
                invEdgeID[p] += [k]
        self.df["edges"] = [invEdgeID[p] for p in self.df.peak]
    
    def getEdges(self,ix=None):
        if ix is None:
            out = sum(self.df.edges,[])
        else:
            out = sum(self.df.loc[ix,"edges"],[])
        out = np.unique(out,axis=0)
        return out
    
    def plotEdges(self, ix=None, ax=None, image=True, imkw_args = {}, separate=False, color="k", lw=None,alpha = 1,fill=False):
        if ix is None:
            ix = self.df.index
        if ax is None:
            ax = plt.subplot(111)
        if lw is None:
            lw=.8
        if "cmap" not in imkw_args:
            imkw_args["cmap"] = "Greys"
        if image:
            im = self.statImages[self.mode]
            ax.imshow(im,norm=LogNorm(),**imkw_args)
        if separate:
            for i in ix:
                try:
                    c = self.df.loc[i,"color"]
                except:
                    c = MYCOLORS[i%len(MYCOLORS)]
                y,x = np.array(self.df.loc[i,"boundary"]+self.df.loc[i,"boundary"][:1]).T
                ax.plot(x,y,"-",lw=lw,c=c,alpha=alpha)
                if fill:
                    ax.fill(x,y,c=c,alpha=alpha*.8)
        else:
            tmp = []
            for el in self.df.loc[ix,"boundary"]:
                tmp += el
                tmp += [el[0]]
                tmp += [(np.nan,)*2]

            y,x = np.array(tmp).T
            ax.plot(x,y,color,lw=lw,alpha=alpha)
            
            
        if hasattr(self, "metadata") and "pxSize" in self.metadata:
            lengths = [10,20,50]
            il = np.searchsorted(lengths,self.metadata.pxSize*self.image.shape[1]/10)
            length=lengths[il]
            x0,x1,y0,y1 = np.array([0,length,0,length*3/50])/self.metadata.pxSize + self.image.shape[0]*.02
            ax.fill_between([x0,x1],[y1]*2,[y0]*2, color="k")
            txt = str(length)
            if "pxUnit" in self.metadata:
                txt += self.metadata["pxUnit"]
            ax.text((x0+x1)/2, y1, txt, va="top", ha="center")
            
    def plotPeaks(self, ix=None, ax=None, image=False, ms=1, labels=False,color=None, imkw_args={},absMarker=True):
        if ax is None:
            ax = plt.subplot(111)
        if image:
            im = self.statImages[self.mode]
            ax.imshow(im,norm=LogNorm(),**imkw_args)
        if ix is None:
            ix = self.df.index

        if absMarker:
            sizes = [ms]*len(self.df)
        else:
            sizes = ms * self.df.loc[ix,"size"]**.5
        for i,ms in zip(ix,sizes):
            p = self.df.loc[i,"peak"]
            if color is None:
                try:
                    c = self.df.loc[i,"color"]
                except:
                    c = MYCOLORS[i%len(MYCOLORS)]
            else:
                c = color
            ax.plot(*p[::-1],marker="o",mfc="none",ms=ms,c=c)
            if labels:
                ax.text(*p[::-1],s=str(i),color=c)
    
    def calc_interest(self, zth=4, timescales=[3,10,30,100,300]):
        interesting = np.zeros(len(self.df))
        for ts in [3,10,30,100,300]:
            if 1./self.Freq   > ts/10: continue
            if self.time[-1] < ts*5: continue
            s,f,z = self.fast_filter_traces(ts, write=False)
            interesting += np.mean(z>zth,1)
        self.df["interest"] = interesting
        
    def calcNNmap(self):
        from bidict import bidict
        peak2idx = bidict([(peak,j) for j,peak in zip(self.df.index,self.df.peak)])
        self.peak2idx = peak2idx
        neighborsMap = {k:[] for k in self.df["peak"]}
        for edge in self.edgeIDs:
            if len(self.edgeIDs[edge])>1:
                for e1,e2 in product(self.edgeIDs[edge],self.edgeIDs[edge]):
                    if e1 not in neighborsMap: continue
                    if e2 not in neighborsMap: continue
                    if e1==e2: continue
                    if e2 not in neighborsMap[e1]:
                        neighborsMap[e1] += [e2]
        self.df["neighbors"] = [[peak2idx[pp] for pp in neighborsMap[p]] for p in self.df["peak"]]
        self.df["Nneighbors"] = self.df["neighbors"].apply(len)
        
    def purge_lones(self,min_size=4, verbose=False):
        toDel = []
        for i in self.df.index:
            if self.df.loc[i,"size"]<min_size and self.df.loc[i,"Nneighbors"]==0:
                toDel += [i]
        self.df = self.df.drop(index=toDel)
        if verbose:
            print (f"deleted {len(toDel)} rois. {len(self.df)} remain.")
    
    def calcTraces(self, movie_=None, FrameRange=None):
        if movie_ is None:
            movie_ = self.movie
        if FrameRange is None:
            try: FrameRange = self.FrameRange
            except: FrameRange = [0,len(movie_)]
        i0,ie = FrameRange
        traces = np.ones((len(self.df),(ie-i0)))*np.nan
        for i,ix in enumerate(self.df.index):
            x = self.df.loc[ix,"pixels"]
            x = [ el[0] for el in x ] , [ el[1] for el in x ]
            traces[i] = movie_[i0:ie, x[0], x[1] ].mean(axis=1)
        self.df["trace"] = list(traces)
        time = np.arange(len(movie_))/movie_.fr
        self.time = time[i0:ie]
        
    def detrend_traces(self,fast=True, timescale=200):
#         from .numeric import mydebleach
#         traces = np.vstack(self.df.trace.values)
#         trend = multi_map( mydebleach, traces, processes=processes)
#         self.df["trend"] = trend
#         self.df["detrended"] = list(traces - np.array(trend))
#         traces = np.vstack(self.df.trace.values)
#         print ("Deprecated method. I used to think it makes sense, now I don't think so. The method is still here, not to break your scripts, but it only subtracts mean from the trace.")
        if fast:
            trend = self.df.trace.apply(np.mean)
            self.df["trend"] = trend
            self.df["detrended"] = [self.df.trace[i] - self.df.trend[i] for i in self.df.index]
        else:
            self.slow_filter_traces(timescale)
            self.df["detrended"] = self.df["faster_%g"%timescale]
            self.df["trend"]     = self.df["slower_%g"%timescale]
            del self.df["faster_%g"%timescale], self.df["slower_%g"%timescale]

    def sortFromCenter(self):
        center = np.array(self.image.shape)/2
        self.df["distToCenter"] = [np.sum((np.array(self.df.loc[i,"peak"])-center)**2)**.5 for i in self.df.index]
        self.df.sort_values("distToCenter",inplace=True)
        self.df.index = np.arange(len(self.df))
        self.calcNNmap()

    def sortInOrder(self):
        idxsort = pd.DataFrame(np.vstack(self.df.peak.values), index=self.df.index).sort_values([0,1]).index
        self.df = self.df.loc[idxsort]
        self.df.index = np.arange(len(self.df))
        self.calcNNmap()
    
    def infer_gain(self, plot=False, verbose=False):
        minDt = np.diff(self.time).mean()
        freq = 1/minDt
#         ts = min(50/freq,10)
        ts = 30/freq
        absSlow, absFast, _ = self.fast_filter_traces(ts,write=False, normalize=False,z_sp=0)
        di = 30
        slow_est, fast_vars = [],[]
        for i in range(absFast.shape[0]):
            for j in range(di, absFast.shape[1]-di, absFast.shape[1]//30):
                slow_est  += [absSlow[i,j]]
                fast_vars += [absFast[i,j-di:j+di].var()]
        fast_vars = np.array(fast_vars)
        slow_est = np.array(slow_est)
        logbs = np.log(np.logspace(np.log10(np.percentile(slow_est,2)),np.log10(np.percentile(slow_est,98))))
        d = np.digitize(np.log(slow_est), logbs)
        x = np.array([slow_est[d==i].mean() for i in np.unique(d)])
        y = np.array([np.median(fast_vars[d==i]) for i in np.unique(d)])
        gain = np.mean(y/x)
        slow_est[slow_est<=0] = np.nan
        if plot:
            ax = plt.subplot(111)
            ax.hexbin(slow_est, fast_vars, bins="log",
                      xscale="log",
                      yscale="log",
                      cmap="Greys",
                      mincnt=1
                     )
            c = ax.plot(x,y,"o",mfc="none")[0].get_color()
            ax.plot(x,x*gain,c=c)
            
        if verbose: print ("initial estimate of the gain is", gain)
        
        for _ in range(5):    
            fast_vars[fast_vars>10*gain*slow_est] = np.nan
            if np.isnan(fast_vars).any():
                y = np.array([np.nanmedian(fast_vars[d==i]) for i in np.unique(d)])
                y[y<=0] = y[y>0].min()
                gain = np.nanmean(y/x)
                if verbose: print ("revised estimate of the gain", gain)
                if plot:
                    c = ax.plot(x,y,"o",mfc="none")[0].get_color()
                    ax.plot(x,x*gain,c=c)
        if plot:
            ax.set_title("gain inference")
            ax.set_xlabel("window means")
            ax.set_ylabel("window variances")
            
        self.gain = gain
    
    def slow_filter_traces(self,ironScale, n_processes=10, percentile = [10.], calcStd=False,
                           avg=True, verbose=False, write=True, indices=None, autorebin=True):
        from .numeric import lowPass
        from .utils import multi_map
        global iterf
        if hasattr(self,"Freq"):
            freq = self.Freq
        else:
            freq = 1./np.diff(self.time).mean()
        wIron = int(ironScale*freq)
        if verbose:
            print(f"The movie frequency is {freq:.2f}, so the filter size is {wIron}.")
        nr = 1
        if autorebin:
            if wIron>200:
                nr = int(np.round(wIron/200))
                if nr<1:
                    nr=1
                if nr>1:
                    wIron = wIron//nr
                    if verbose:
                        print (f"Rebinning by {nr}, so the filter size is now {wIron}")
        if wIron%2==0:
            wIron += 1
        if avg:
            def iterf(x_): 
                out = lowPass(x_, wIron, wIron, percentile)
                return out
        else:
            def iterf(x_): 
                out = lowPass(x_, wIron, perc=percentile)
                return out
        if indices is None:
            indices = self.df.index
        assert len(indices)==len(np.unique(indices))
        if not hasattr(self,"gain") or self.gain<=0:
            raise ValueError("Regions need gain to work. Please infer_gain before running this function.")
        traces = np.vstack(self.df.loc[indices,"trace"].values)
        if nr>1:
            traces_ = rebin(traces,nr,1)
        else:
            traces_ = traces
        slow   = np.array(multi_map(iterf,traces_,processes=n_processes,))
        if nr>1:
            from scipy.interpolate import interp1d
            tr = rebin(self.time,nr)
            slow = np.array([interp1d(tr,s, kind="quadratic",fill_value="extrapolate")(self.time) for s in slow])
        fast    = traces-slow
        absFast = np.array([ fast[i]*self.df.loc[ix,"size"] for i,ix in enumerate(indices) ])
        absSlow = np.array([ slow[i]*self.df.loc[ix,"size"] for i,ix in enumerate(indices) ])
        zScore  = absFast/(self.gain*absSlow/nr)**.5
        print(zScore.shape)
        if write:
            for k,v in zip(
                ["slower_%g_"%ironScale, "faster_%g_"%ironScale, "zScore_%g_"%ironScale],
                [slow, fast, zScore]):
                tmp = np.ones((len(self.df),len(self.time)))*np.nan
                tmp[self.df.index.isin(indices)] = v
                self.df[k] = list(tmp)
        else:
            return fast, slow, zScore
            
            
#         def iterf(x_):
#             mad2std = 1.4826
#             out = mad2std*lowPass(np.abs(x_),wIron,wIron*2+1,50.)
#             return out
#         if calcStd:
#             self.df["faster_%g_std"%ironScale] = multi_map(iterf,self.df["faster_%g"%ironScale].values,processes=n_processes)
    
    
    
    def fast_filter_traces(self,
                           ironTimeScale,
                           z_sp=3,
                           order=5,
                           Npoints=None,
                           write=True,
                           verbose=False,
                           usecol="trace",
                           normalize=True
#                            meanSlow2Var=None
                          ):
        from .numeric import sosFilter
        if Npoints is None: Npoints = 15
        minDt = np.diff(self.time).mean()
        freq = 1/minDt
        try:
            self.movie
            if np.abs(freq/self.movie.fr-1)>1e-2:
                print (f"movie frame rate ({self.movie.fr}) and inferred frame ({freq}) rate are different!")
        except:
            pass
        N_dt = ironTimeScale/minDt
        Nrebin = max(1,int(np.round(N_dt/Npoints)))
        if verbose:
            print (f"Nrebin = {Nrebin}")
        C = self.df
        if "trace" in usecol:
            data = np.vstack([C.loc[i,usecol] for i in C.index])
            trend = np.zeros(len(C))
        else:
            data = np.vstack([C.loc[i,"detrended"] for i in C.index])
            trend = C.trend.values
            
        if Nrebin>1:
            freq = freq/Nrebin
            data = rebin(data, Nrebin, axis=1)
            if len(trend.shape)>1:
                trend = rebin(trend, Nrebin, axis=1)
            try: self.showTime
            except: self.showTime = {}
            self.showTime["%g"%ironTimeScale] = rebin(self.time, Nrebin)
        cutFreq = .5/ironTimeScale
        self.sosFilter = sosFilter(cutFreq, freq, order=order)
        dataFilt = self.sosFilter.run(data)

        absSlow = np.vstack([dataFilt[i]*C.loc[ix,"size"] for i,ix in enumerate(C.index)])*Nrebin
        absFast = np.vstack([(data[i]-dataFilt[i])*C.loc[ix,"size"] for i,ix in enumerate(C.index)])*Nrebin

        if z_sp>0:
            from cv2 import dilate
            from .numeric import nan_helper
            var = absSlow
            if hasattr(self,"gain"):
                var = var*self.gain    
            std = var**.5
            ff = (absFast > z_sp*std).astype("uint8")
            if ff.any():
                dilateKernelSize = int(ironTimeScale/minDt/Nrebin*.2)#*.03)
                if dilateKernelSize%2==0:
                    dilateKernelSize+=1
                if dilateKernelSize>=3:
                    if verbose:
                        print ("dilating by", dilateKernelSize)
                    ff = dilate(ff, np.ones(dilateKernelSize, dtype = np.uint8).reshape(1,-1)).astype("bool")
                absFast_tmp = absFast.copy()
                absFast_tmp[ff] = np.nan
                for j in range(C.shape[0]):
                    y = absFast_tmp[j]
                    nans, x= nan_helper(y)
                    if nans.any() and nans.mean()<.5: 
                        y[nans]= np.interp( x(nans), x(~nans), y[~nans] )
#                         y[nans]= np.interp( x(nans), x(~nans), absSlow[j][~nans] )
                dFast = self.sosFilter.run(absFast_tmp)
                absSlow = absSlow + dFast
                absFast = absFast - dFast
        var = absSlow
        if hasattr(self,"gain"):
            var = var*self.gain    
        std = var**.5
        zScore = absFast/std
        if normalize:
            slower = [absSlow[i]/C["size"].iloc[i]/Nrebin for i in range(len(C))]
            faster = [absFast[i]/C["size"].iloc[i]/Nrebin for i in range(len(C))]
        else:
            slower = absSlow
            faster = absFast
        if write:
            C["slower_%g"%ironTimeScale] = slower
            C["faster_%g"%ironTimeScale] = faster
            C["zScore_%g"%ironTimeScale] = list(zScore)
        else:
            return np.array(slower), np.array(faster), zScore

    def calc_raster(self, ts, z_th = 3, Npoints=None, smooth = 0):
        from .numeric import runningAverage
        if "zScore_%g"%ts not in self.df.columns:
            self.fast_filter_traces(ts,Npoints=Npoints)
        zScores = np.vstack(self.df["zScore_%g"%ts])
        if smooth:
            avgSize = 2*smooth+1
            zScores = runningAverage(zScores.T,avgSize).T#*avgSize**.5
        k = "%g"%(ts)
        try:
            self.raster
        except:
            self.raster = {}
        self.raster[k] = zScores>z_th
        
    def peaks2raster(self, ts, npoints = 1000, onlyRaster=True, z_th = 3):
        k = "%g"%ts
        try:
            self.peaks[k]
        except:
            self.calc_peaks(ts, z_th = z_th)
        df = self.peaks[k]
        C = self.df
        rr = np.zeros((len(C),npoints))
        tt = pd.Series(np.linspace(0,self.time.max(),npoints))
        for i,ix in enumerate(C.index):
            ddf = df.query(f"roi=={ix}")
            for _,row in ddf.iterrows():
                rr[i,tt.between(row.t0,row.t0+row.iloc[1])] = 1
        if onlyRaster:
            return rr,
        fig = go.Figure(go.Heatmap(
                x=tt,
                z=rr,
                showscale=False,
                hoverinfo="text",text=[[str(i)] for i in self.df.index]
            ))
        fig.update_yaxes(title_text='roi id')
        fig.update_xaxes(title_text='time [s]')
        fig.update_layout({
            "width":  600,
            "height": 400,
            "margin": dict(l=20, r=10, t=50, b=20),
        })
            # fig.update_xaxes(showticklabels=False)
        return rr, fig
        
    def calc_peaks(self, ts, z_th=3, Npoints=None, smooth=None, verbose=False, save=True, t=None, zScores=None):
        from .numeric import runningAverage
        from scipy.signal import find_peaks, peak_widths
        if zScores is None:
            if "zScore_%g"%ts not in self.df.columns:
                if verbose:
                    print ("filtering...")
                self.fast_filter_traces(ts,Npoints=Npoints)
            zScores = np.vstack(self.df["zScore_%g"%ts])
        if t is None:
            try:
                t = self.showTime["%g"%ts]
            except:
                t = self.time
        dt = np.diff(t).mean()
        if smooth is None:
            smooth = int(ts/dt/5)
            if smooth%2==0: smooth += 1
        if smooth>0:
            if verbose:
                print ("smoothing with kernel", smooth)
            zScores = runningAverage(zScores.T,smooth).T
        peaks = []
        for i,z in zip(self.df.index,zScores):
            pp = find_peaks(z,
                            width=ts/dt/8,
                            height=z_th
                              )
            w,h,x0 = peak_widths(z, pp[0], rel_height=.5)[:3]
            w = w*dt
            x0 = x0*dt + t[0]
            df = pd.DataFrame({"peak height [z-score]":z[pp[0]],"peak half-width [s]":w, "t0":x0})
            df["roi"] = i
            peaks += [df]
        peaks = pd.concat(peaks,ignore_index=True)
        if save:
            k = "%g"%(ts)
            try:
                self.peaks
            except:
                self.peaks = {}
            self.peaks[k] = peaks
        else:
            return peaks
        
    def show_scatter_peaks(self,ts,timeWindows=None):
        from plotly_express import scatter
        peaks = self.peaks["%g"%ts].copy()
        if timeWindows is None:
            timeWindows = [[0,np.inf]]
        peaks["timeWindow"] = [""]*len(peaks)
        for i,tw in enumerate(timeWindows):
#                 print (tw)
            if tw[1]<=tw[0]:continue
            iis = np.where(peaks["t0"].between(*tw) & (peaks["t0"]+peaks[peaks.columns[1]]).between(*tw))[0]
            for j in iis:
                if peaks.loc[j,"timeWindow"]=="":
                    peaks.loc[j,"timeWindow"] = str(i)
                else:
                    j1 = len(peaks)
                    peaks.loc[j1] = peaks.loc[j]
                    peaks.loc[j1,"timeWindow"] = str(i)
#             peaks["timeWindow"][ff] = str(i+1)
        df = peaks.query("timeWindow!=''")    
        from plotly.colors import DEFAULT_PLOTLY_COLORS
        fig = scatter(df,
                      x=df.columns[1],
                      y=df.columns[0],
                      opacity = .2,
                      labels = ["roi"],
                      color=[DEFAULT_PLOTLY_COLORS[int(c)] for c in df["timeWindow"]],
#                       hover_data={k:(k=="roi") for k in df.columns},
#                       hover_data={"roi":True},
                      hover_data=["roi"],
#                       hoverinfo=df["roi"].astype("str"),
                      marginal_x="box",
                      log_y=True,
                      render_mode = "webgl",
                      width=450,
                      height=450,
                     )
#         fig.update_traces(hovertemplate='x: %{x} <br>y: %{y} <br>roi: %{roi}') # 
        fig.update_layout({
        "plot_bgcolor":"white","margin":dict(l=10, r=10, t=20, b=40),"showlegend":False})
        fig.update_layout(
            legend=dict(
                x=.99,
                y=.72,
                traceorder="normal",
                xanchor="right",
            )
        )
        
        fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True, ticks="outside", ticklen=2,
#                          gridcolor="none"
                        )
        fig.update_yaxes(showline=True, linewidth=1, linecolor='black', mirror=True, ticks="outside", ticklen=2,
#                          gridcolor="none"
                        )
        return fig

    def calcIntraCCs(self,movie_,diff=False,indices=None):
        intraCCs = []
        if indices is None:
            indices = self.df.index
        for roi in indices:
            if C.loc[roi,"size"]<1:
                intraCCs += [[np.nan]]
                continue
            trace = C.loc[roi,"trace"]
            pixels = C.loc[roi,"pixels"]
            size = len(pixels)
            indTraces = [movie_[(slice(None),) + px] for px in pixels]
            ccs = []
            for tr in indTraces:
                x,y = trace*size-tr, tr
                if diff:
                    x = np.diff(x)
                    y = np.diff(y)
                ccs += [np.corrcoef(x,y)[0,1]]
            intraCCs += [ccs]

        C.loc[indices,"intraCCs"] = intraCCs

    def import_protocol(self,pathToProtocol):
        protocol = pd.read_csv(pathToProtocol)
        protocol.dropna(how='all', inplace=True)
        if protocol.empty:
            raise ValueError('protocol file contains nothing ')
        try:
            protocol["t_begin"] = pd.to_timedelta(["00:"+el if type(el)==str else "00:00:00" \
                                                           for el in protocol["begin"]]).total_seconds()
            protocol["t_end"] = pd.to_timedelta(["00:"+el if type(el)==str else (self.time[0]+len(self.time)/self.Freq)*1e9 \
                                                         for el in protocol["end"]]).total_seconds()
            self.protocol = protocol
        except:
            pass
        return protocol

    def examine(self, max_rois=10, imagemode=None, debug=False, startShow='',mode="jupyter",name=None,lw=None):
        from .examine import examine
        if imagemode is None:
            imagemode = self.mode
        return examine(self, max_rois=max_rois, imagemode=imagemode, debug=debug, startShow=startShow,mode=mode,name=name)
    
    def examine3(self, max_rois=10, imagemode=None, debug=False, startShow='',mode="jupyter",name=None,lw=None):
        from .examine3 import examine
        if imagemode is None:
            imagemode = self.mode
        return examine(self, max_rois=max_rois, imagemode=imagemode, debug=debug, startShow=startShow,mode=mode,name=name)
    
    def plotTraces(regions, indices, axratios = [1,2], figsize=5, freqShow=2, col="detrended",Offset=5,separate=False):
        xratios = np.array([.1,axratios[0],.1,axratios[1],.1])
        yratios = xratios[:3]
        xr = xratios/sum(xratios)
        yr = yratios/sum(yratios)

        fig = plt.figure(figsize=(xratios.sum()*figsize,yratios.sum()*figsize))
        axs = [
            fig.add_axes([xr[0],yr[0],xr[1],yr[1]]),
            fig.add_axes([xr[:3].sum(),yr[0],xr[3],yr[1]]),
        ]
        regions.plotEdges(ax=axs[0],lw=.5,separate=separate)
        regions.plotEdges(ax=axs[0],ix=indices, separate=True, fill=True, alpha=.5, image=False)
        regions.plotPeaks(ax=axs[0],ix=indices, labels=True)
        nr = int(np.round(regions.Freq/freqShow))
        if nr==0:
            nr=1
        xs = np.vstack(regions.df.loc[indices,col].values)
        if nr>1:
            t = rebin(regions.time, nr)
            xs = rebin(xs,nr,1)
        else:
            t = regions.time

        for i in range(len(xs)):
            xs[i] = xs[i]-np.median(xs[i])
        offset = 0
        if "color" in regions.df.columns:
            colors = regions.df.loc[indices,"color"]
        else:
            colors = [MYCOLORS[ix%len(MYCOLORS)] for ix in indices]
        for x,ix,c in zip(xs, indices, colors):
            axs[1].plot(t,x+offset,lw=.5,color=c)
            axs[1].text(0,offset,str(ix)+" ",color=c,ha="right")
            offset += xs.std()*Offset
        axs[1].set_yticks([])
        axs[1].set_xlabel("time [s]")
        axs[0].set_yticks([])
        axs[0].set_xticks([])
        for sp in ["left","right","top"]: axs[1].spines[sp].set_visible(False)
        if not hasattr(regions,"protocol"):
            return None
        yl = axs[1].get_ylim()[1]
        offset = 2

        for comp, df in regions.protocol.groupby("compound"):
            for ii in df.index:
                t0,t1 = df.loc[ii].iloc[-2:]
                conc = df.loc[ii,"concentration"]
                x,y = [t0,t1,t1,t0,t0],[-1,-1,-2,-2,-1]
                y = np.array(y)-offset
                y = y*yl/20
                plt.fill(x,y,color="grey",alpha =.3)
                plt.text(t0,y[:-1].mean(), " "+conc,va="center", ha="left")
                plt.plot(x,y,color="grey",)

            plt.text(df.t_begin.min(),y[:-1].mean(),comp+" ",va="center", ha="right")
            offset += -1.3    

def getGraph_of_ROIs_to_Merge(df,rreg, plot=False, ax=None,lw=.5,arrow_width=.5):
    C = rreg.df
    Gph = nx.DiGraph()
    ixs = np.unique(df.values.flatten())
    if plot:
        if ax is None:
            plt.figure(figsize=(10,10))
            ax = plt.subplot(111)
        rreg.plotEdges(image=False, ix=ixs, ax=ax, color="k")


    for _,row in df.iterrows():
        i,j = row[["i","j"]]
        l = list(df.query(f"i=={j}")["j"])
        if i in l:
            c="darkgoldenrod"
            iaccept = C.loc[sorted([i,j]),"peakValue"].idxmax()
            if j!=iaccept:
                continue
        else:
            c="r"
        Gph.add_edge(i,j)
        if plot:
            x0,y0 = C.loc[i,"peak"]
            x1,y1 = C.loc[j,"peak"]
            dx = x1-x0
            dy = y1-y0
            ax.arrow(y0,x0,dy,dx,width = arrow_width,
                     linewidth = lw,
                     color="r",
                     zorder=10,
                     length_includes_head=True)
            
            #ax.plot(y1,x1,"o",ms=5,mfc="none",mew=.7,c=c)
            
    if plot:
        plt.gca().set_aspect("equal")
        attractors = np.squeeze([
            list(map(list,nx.attracting_components(Gph.subgraph(el)))) \
                for el in nx.connected_components(Gph.to_undirected())
        ])
        rreg.plotPeaks(ax=ax,ix=attractors,color="c",ms=6)
    return Gph

def plotRoi_to_be_connected(Gph, rreg, nplot=35):
    C = rreg.df
    Gph_ = Gph.to_undirected()
    dd = list(nx.connected_components(Gph_))
    dd = sorted(dd,key=len)[::-1][:nplot]
    nc = 7
    nr = int(np.ceil(len(dd)/nc))

    fig, axs = plt.subplots(nr,nc,figsize=(2*nc,nr*2))
    for i,cl in enumerate(dd):
        try: ax = axs.flat[i]
        except: break
        cl = list(cl)
        gph = Gph.subgraph(nodes=cl)
        pos = np.array([el[::-1] for el in C.loc[cl,"peak"]])
        nx.draw_networkx(gph,
                         ax=ax,
                         node_size=30,
                         node_color="w",
                         pos=dict(zip(cl,pos)),
                         font_size=6
                        )
        rreg.plotEdges(ax=ax,image=False,ix=cl)
        attr = sum(list(map(list,nx.attracting_components(gph))),[])
        rreg.plotEdges(ax=ax,image=False,ix=attr,color="red")
        rreg.plotPeaks(ax=ax,image=False,ix=attr,color="red",ms=1)
        for sp in ax.spines: ax.spines[sp].set_visible(False)
        ax.set_aspect("equal")

    for i in range(i+1,axs.size):
        axs.flat[i].remove()
    plt.subplots_adjust(wspace=0,hspace=0)

def mergeBasedOnGraph(Gph,rreg,verbose=False):
    C = rreg.df
    toDrop = []
    Gph_ = Gph.to_undirected()
    for cl in nx.connected_components(Gph_):
        cl = list(cl)
        gph = Gph.subgraph(nodes=cl)
        attr = sum(list(map(list,nx.attracting_components(gph))),[])
        if len(attr)>2:
            # print ("more than two attractors, not implemented yet, will skip")
            continue
        if len(attr)==2 and len(cl)>2:
            continue
        attr = C.loc[attr,"peakValue"].sort_values().index[-1]
        other = [j for j in cl if j!=attr]
        C.loc[[attr],"pixels"] = [sum(C.loc[cl,"pixels"],[])]
        toDrop += other
    C.drop(index=toDrop,inplace=True)
    if verbose:
        print (f"{len(toDrop)} subsumed into existing ROIs.")
    rreg.update()
#     rreg.sortFromCenter()
    try:
        rreg.calcTraces()
    except:
        pass
    return len(toDrop)

def getPeak2BounAndTraceDF(C):
    peak2bnd = []
    for i in C.index:
        pk = C.loc[i,"peak"]
        if len(C.loc[i,"neighbors"])==0:
            continue
        try: tr_i = np.diff(C.loc[i,"trace"])
        except: pass
        for j in C.loc[i,"neighbors"]:
            bd = C.loc[j,"boundary"]
            x = np.linalg.norm(np.array(bd)-np.repeat([pk],len(bd),axis=0), axis=1)
            out = ( i, j, C.loc[i,"size"], C.loc[j,"size"], x.min(), sum(x==x.min()))
            try:
                tr_j = np.diff(C.loc[j,"trace"])
                cc = np.corrcoef(tr_i,tr_j)[0,1]
                out = out + (cc,)
            except: pass
            peak2bnd += [ out ]
            
    try: peak2bnd = pd.DataFrame(peak2bnd, columns=["i","j","size_i","size_j","dist","nclose","cc"])
    except: peak2bnd = pd.DataFrame(peak2bnd, columns=["i","j","size_i","size_j","dist","nclose"])
    return peak2bnd

def getPeak2BoundaryDF(C):
    peak2bnd = []
    for i in C.index:
        pk = C.loc[i,"peak"]
        if len(C.loc[i,"neighbors"])==0:
            continue
        dists = OrderedDict()
        for j in C.loc[i,"neighbors"]:
            if j not in C.index: continue
            bd = C.loc[j,"boundary"]
            x = np.linalg.norm(np.array(bd)-np.repeat([pk],len(bd),axis=0), axis=1)
            dists[j] = x.min()
        if len(dists):
            jmin = pd.Series(dists).idxmin()
            peak2bnd += [(i,jmin,dists[jmin],C.loc[i,"size"],C.loc[jmin,"size"])]

    peak2bnd = pd.DataFrame(peak2bnd, columns=["i","j","dist","size_from","size_to"])
    return peak2bnd

def getPeak2EdgesDF(C, regions):
    peak2bnd = []
    for i in C.index:
        pk = C.loc[i,"peak"]
        pxi = C.loc[i,"pixels"]
        if len(C.loc[i,"neighbors"])==0:
            continue
        dists = OrderedDict()
        for j in C.loc[i,"neighbors"]:
            pxj = C.loc[j,"pixels"]
            edges = C.loc[j,"edges"]
            emean = np.array(edges).reshape((len(edges),2,2)).mean(axis=1)
            xx = np.linalg.norm(emean-np.repeat([pk],len(emean),axis=0), axis=1)
            barriers = []
            if xx.min()>=1:
                continue
            for k in np.where(xx==xx.min())[0]:
                emin = edges[k]
                if emin[0]==emin[2]:
                    y = int((emin[1]+emin[3])/2)
                    pxs = (int(emin[0]-.5),y),(int(emin[0]+.5),y)
                else:
                    assert emin[1]==emin[3]
                    x = int((emin[0]+emin[2])/2)
                    pxs = (x,int(emin[1]-.5)),(x,int(emin[1]+.5))
                if pxs[0] not in pxi:
                    pxs = pxs[1],pxs[0]
                if pxs[1] not in pxj:
                    print (pxs,i,j)
                    assert pxs[1] in pxj
                barriers += [(k,regions.image[pxs[1]]-regions.image[pxs[0]])]
            barriers = sorted(barriers, key=lambda xi: xi[1])[::-1]
            imin,barrier = barriers[0]
            dists[j] = xx[imin],barrier
        if len(dists)==0:
            continue
        df = pd.DataFrame(dists).T
        jmin = df.sort_values([0,1]).index[0]
        peak2bnd += [(i,jmin)+tuple(df.loc[jmin])]
    peak2bnd = pd.DataFrame(peak2bnd, columns=["i","j","dist","barrier"])
    return peak2bnd


def edges2nodes(x,start=0,direction=1):
    if np.array(x[0]).shape != (2,2):
        x = [(el[:2],el[2:]) for el in x]
    nodes = list(x[start][::direction])

    for i in range(len(x)-1):
        nexts = [edge for edge in x if (edge[0]==nodes[-1] or edge[1]==nodes[-1])]
        for cand in np.unique(sum(nexts,()),axis=0):
            if tuple(cand) not in nodes:
                nodes += [tuple(cand)]
    return nodes

def getStatImages(movie_, debleach=False, downsampleFreq=5):
    if movie_.fr>downsampleFreq:
        n_rebin = int(np.round(movie_.fr/downsampleFreq))
        if n_rebin>=2:
            m_for_image = rebin(movie_,n_rebin)
        else:
            m_for_image = movie_
    else:
        m_for_image = movie_
    statImages = {}
    if debleach:
        m_for_image = m_for_image.astype("float32")
        m_for_image.debleach()

    for f in [np.mean,np.std]:
        statImages[f.__name__] = f(m_for_image,axis=0)
    statImages["highperc"] = np.percentile(m_for_image,100*(1-10/len(m_for_image)), axis=0)
    
    m_for_image = np.diff(m_for_image,axis=0)
    for f in [np.mean,np.std]:
        statImages["diff_"+f.__name__] = f(m_for_image,axis=0)
    statImages["diff_highperc"] = np.percentile(m_for_image,100*(1-10/len(m_for_image)), axis=0)
    
    return statImages
