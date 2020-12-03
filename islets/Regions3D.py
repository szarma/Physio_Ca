import numpy as np
import pandas as pd
from itertools import product
from collections import OrderedDict
import matplotlib.pyplot as plt
import networkx as nx
from .numeric import rebin
from .utils import multi_map
import plotly.graph_objects as go
from matplotlib._color_data import TABLEAU_COLORS, CSS4_COLORS
MYCOLORS = OrderedDict(TABLEAU_COLORS)
del MYCOLORS["tab:gray"]
ks = ["lime",
    "orangered",
    "yellowgreen",
      "lawngreen",
    # "mediumspringgreen",
    "aquamarine"]#+list(MYCOLORS.values())
MYCOLORS.update({k:CSS4_COLORS[k] for k in ks})
MYCOLORS = list(MYCOLORS.values())

MYCOLORS = ["darkred"]

def climb(x, a, diag=True, min_gradient=0, verbose=False):
    dims = a.shape
    xs = [x]
    v = a[tuple(x)]
    for i in range(100):
        ijs = []
        vs = []
        for dij in product(*([[-1,0,1]]*len(dims))):
            if not diag:
                if np.prod(dij)!=0: continue
            ij = x + np.array(dij)
            if any(
                [ij[k]<0 for k in range(len(dims))]+
                [ij[k]>=dims[k] for k in range(len(dims))]
            ):
                continue
            vs += [a[tuple(ij)]]
            ijs += [ij]
        x1 = ijs[np.argmax(vs)]#[:-1]
        dv = max(vs)-v
        if dv<=0 or dv<=min_gradient:
            break
        else:
            x = x1
            v = max(vs)
            xs += [x]
            if verbose:
                print ("=>", xs)
    return tuple(xs[-1])

def crawlDict(image, crawl_th=-np.inf, diag=False, min_gradient=0, processes=1):
    global iterf
    ijs_ = [ij for ij in product(*map(range,image.shape)) if image[ij]>crawl_th]
    def iterf(ij):
        return climb(ij, image, diag=diag, min_gradient=min_gradient)
    R_ = multi_map(iterf,ijs_, processes=processes)
    B_ = OrderedDict()
#     A_ = [ij+(r,) for ij,r in zip(ijs_,R_)]
#     A_ = [el for el in A_ if el[-1] is not None]
#     return A_
    for (ij0,ij1) in zip(ijs_,R_):
        if ij1 not in B_:
            B_[ij1] = []
        B_[ij1] += [ij0]
    return B_


class Regions:
    def __init__(self, movie_, diag=False, min_gradient=0, debleach=False, gSig_filt=None, mode="diff_std", full=True, img_th=-np.inf, FrameRange=None, processes=10):
        self.mode = mode
        if isinstance(movie_, (np.ndarray,)):
            if len(movie_.shape)==2:
                image0 = movie_
                self.statImages = {mode:image0}
            if len(movie_.shape)==3:
                self.movie = movie_
                time = np.arange(len(movie_))/movie_.fr
                if FrameRange is None:
                    FrameRange = [0, len(movie_)]
                i0, ie = FrameRange
                self.FrameRange = FrameRange
                self.statImages = getStatImages(movie_[i0:ie], debleach=debleach, downsampleFreq=5)
                self.time = time[i0:ie]
                self.Freq = movie_.fr
                self.showTime = {}
        elif isinstance(movie_, dict):
            self.statImages = movie_
        else:
            raise ValueError("Regions can initialize either from a movie, or image, or a dictionary of stats images.")
        
        if full:
            self.constructRois(mode=mode, img_th=img_th, diag=diag, min_gradient=min_gradient, gSig_filt=gSig_filt, processes=processes)
            
    def constructRois(self, mode="diff_std", img_th=-np.inf, diag=False, min_gradient=0, gSig_filt=None, processes=5):
        from caiman.motion_correction import high_pass_filter_space
        from pandas import DataFrame
        
        for k0 in self.statImages:
            break
        tmp = np.zeros_like(self.statImages[k0])
        for submode in mode.split("+"):
            im = self.statImages[submode].copy()
            im = im/im.max()
            tmp += im
        self.statImages[mode] = tmp/tmp.max()
        image0 = tmp/tmp.max()
        toMin = image0<=img_th
        if gSig_filt is None:
            image = image0
        else:
            image = []
            if type(gSig_filt)==int:
                gSig_filt = [gSig_filt]
            for gSf in gSig_filt:
                tmp = high_pass_filter_space(image0,(gSf,)*2)
                tmp = np.maximum(tmp,0)
                tmp = tmp/np.percentile(np.abs(tmp),(1-10/tmp.size)*100)
                tmp *= gSf
                image += [tmp]
            image = np.mean(image,axis=0)
        image[toMin] = image.min()
        self.filterSize = gSig_filt
        self.image = image
        B_ = crawlDict(image,image.min(),diag=diag,min_gradient=min_gradient, processes=processes)
        if diag:
            try:
                from scipy.spatial import distance_matrix
                # sort out unconnected pixels
                ks = list(B_.keys())
                for k in ks:
                    pxs = B_[k]
                    if len(pxs)==1:
                        continue
                    dm = distance_matrix(pxs,pxs)
                    gr = nx.Graph(dm<=1)
                    subsets = list(nx.connected_components(gr))
                    if len(subsets)==1:
                        continue
                    for subset in subsets:
                        tmppxs = [pxs[j] for j in subset]
                        vs = [self.image[px] for px in tmppxs]
                        B_[tmppxs[np.argmax(vs)]] = tmppxs
                # sort out "holes" pixels
                allTakenPx = sum(B_.values(),[])
                dims = self.image.shape
                freePx = [(i,j) for i,j in product(range(dims[0]),range(dims[1])) if (i,j) not in allTakenPx]
                dm = distance_matrix(freePx,freePx)
                gg = nx.Graph(dm==1)
                for cc in nx.connected.connected_components(gg):
                    if len(cc)!=1: continue
                    px = freePx[min(cc)]
                    for di in [-1,1]:
                        pxx = px[0]+di
                        if pxx<0: continue
                        if pxx>=dims[0]: continue
                        nnpeak = climb((px[0]+di,px[1]),self.image, diag=diag)
                    B_[nnpeak] += [px]
            except:
                print ("Cannot initialize with diagonal crawl. Reverting to diag=False")
                diag = False
                B_ = crawlDict(image,image.min(),diag=diag,min_gradient=min_gradient, processes=processes)
        
        
        self.df = DataFrame(OrderedDict([
            ("peak",  list(B_.keys())),
            ("pixels",list(B_.values()))
        ]))
        self.df["peakValue"] = [image[p] for p in B_]
        self.update()
#         print(f"Initialized with {len(self.df)} rois.")
    
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
    
    def plotEdges(self, ix=None, ax=None, image =True, imkw_args = {}, separate=False, color="k", lw=None,alpha = 1):
        from matplotlib.colors import LogNorm
        if ix is None:
            ix = self.df.index
        if ax is None:
            ax = plt.subplot(111)
        if lw is None:
            lw=.8
        if image is not False:
            im = self.statImages[self.mode].copy()
            im = im-np.percentile(im,1)
            im = im/np.percentile(im,99.9)
            im = np.minimum(im,1)
            im += .05
            ax.imshow(im,norm=LogNorm(),**imkw_args)
        if separate:
            print ("going in")
            for i in ix:
                c = MYCOLORS[i%len(MYCOLORS)]
                y,x = np.array(self.df.loc[i,"boundary"]).T
                ax.plot(x,y,"-",lw=lw,c=c,alpha=alpha)
        else:
            tmp = []
            for el in self.df.loc[ix,"boundary"]:
                tmp += el
                tmp += [el[0]]
                tmp += [(np.nan,)*2]

            y,x = np.array(tmp).T
            ax.plot(x,y,color,lw=lw,alpha=alpha)
            
    def plotPeaks(self, ix=None, ax=None, image=False, ms=1, labels=False,color=None, imkw_args={},absMarker=True):
        if ax is None:
            ax = plt.subplot(111)
        if image:
            if not len(imkw_args):
                v = np.abs(self.image).max()*3
                ax.imshow(self.image,vmin=-v,vmax=v)
            else:
                ax.imshow(self.image,**imkw_args)
        if ix is None:
            ix = self.df.index

        peaks = self.df.loc[ix,"peak"]
        if absMarker:
            sizes = [ms]*len(self.df)
        else:
            sizes = ms * self.df.loc[ix,"size"]**.5
        for ms,p in zip(sizes,peaks):
            i = self.peak2idx[p]
            if color is None:
                c = MYCOLORS[i%len(MYCOLORS)]
            else:
                c = color
            ax.plot(*p[::-1],marker="o",mfc="none",ms=ms,c=c)
            if labels:
                ax.text(*p[::-1],s=str(i),color=c)
    
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
        
    def purge_lones(self,min_size=4):
        toDel = []
        for i in self.df.index:
            if self.df.loc[i,"size"]<min_size and self.df.loc[i,"Nneighbors"]==0:
                toDel += [i]
        self.df = self.df.drop(index=toDel)
        print (f"deleted {len(toDel)} rois. {len(self.df)} remain.")
    
    def calcTraces(self, movie_=None, FrameRange=None, fr=None):
        if movie_ is None:
            movie_ = self.movie
        if FrameRange is None:
            try: FrameRange = self.FrameRange
            except: FrameRange = [0,len(movie_)]
        i0,ie = FrameRange
        traces = np.ones((len(self.df),(ie-i0)))*np.nan
        for i,ix in enumerate(self.df.index):
            x = np.vstack(self.df.loc[ix,"pixels"])
#             x = [ el[0] for el in x ] , [ el[1] for el in x ]
            traces[i] = movie_[ (slice(i0,ie),) + tuple(x.T) ].mean(axis=1)
        self.df["trace"] = list(traces)
        if fr is None:
            fr = movie_.fr
        time = np.arange(len(movie_))/fr
        self.time = time[i0:ie]
        
    def detrend_traces(self,processes=10):
        from .numeric import mydebleach
        traces = np.vstack(self.df.trace.values)
        trend = multi_map( mydebleach, traces, processes=processes)
        self.df["trend"] = trend
        self.df["detrended"] = list(traces - np.array(trend))

    def sortFromCenter(self):
        center = np.array(self.image.shape)/2
        self.df["distToCenter"] = [np.sum((np.array(self.df.loc[i,"peak"])-center)**2)**.5 for i in self.df.index]
        self.df.sort_values("distToCenter",inplace=True)
        self.df.index = np.arange(len(self.df))
        self.calcNNmap()
    
    def infer_gain(self, plot=False):
        minDt = np.diff(self.time).mean()
        freq = 1/minDt
#         ts = min(50/freq,10)
        ts = 30/freq
        absSlow, absFast = self.fast_filter_traces(ts,write=False)
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
        gain = np.exp(np.mean(np.log(y)-np.log(x)))
        
        if plot:
            ax = plt.subplot(111)
            ax.hexbin(slow_est, fast_vars, bins="log",
                      xscale="log",
                      yscale="log",
                      cmap="hot",
                      mincnt=1
                     )
            ax.plot(x,y,"C0o",mfc="none")
            ax.plot(x,x*gain)
        self.gain = gain
    
    def fast_filter_traces(self, ironTimeScale, z_sp = 2, order=5, Npoints = None, meanSlow2Var=None, write=True):
        from .numeric import sosFilter
        if Npoints is None: Npoints = 15
        minDt = np.diff(self.time).mean()
        freq = 1/minDt
        try:
            self.movie
            if np.abs(freq/self.movie.fr-1)<1e-5:
                print (f"movie frame rate ({self.movie.fr}) and inferred frame ({freq}) rate are different!")
        except:
            pass
        N_dt = ironTimeScale/minDt
        Nrebin = max(1,int(np.round(N_dt/Npoints)))
        print (f"Nrebin = {Nrebin}")
        C = self.df
        useDetrended = "detrended" in C.columns
        if useDetrended:
            data = np.vstack([C.loc[i,"detrended"] for i in C.index])
            trend = np.vstack(C.trend)
        else:
            data = np.vstack([C.loc[i,"trace"] for i in C.index])
        if Nrebin>1:
            freq = freq/Nrebin
            data = rebin(data, Nrebin, axis=1)
            try: self.showTime
            except: self.showTime = {}
            self.showTime["%g"%ironTimeScale] = rebin(self.time, Nrebin)
            if useDetrended:
                trend = rebin(trend,Nrebin, axis=1)
        cutFreq = .5/ironTimeScale
        self.sosFilter = sosFilter(cutFreq, freq, order=order)
        dataFilt = self.sosFilter.run(data)
        if write:
            slowk, fastk = "slower_%g"%ironTimeScale, "faster_%g"%ironTimeScale
        else:
            slowk, fastk = "_slower_%g"%ironTimeScale, "_faster_%g"%ironTimeScale
            
        if useDetrended:
            self.df[slowk] = list(trend+dataFilt)
        else:
            self.df[slowk] = list(dataFilt)
        self.df[fastk] = list(data - dataFilt)

        absFast = np.vstack([C.loc[i,fastk]*C.loc[i,"size"] for i in C.index])*Nrebin
        absSlow = np.vstack([C.loc[i,slowk]*C.loc[i,"size"] for i in C.index])*Nrebin

        if not write:
            del C[slowk], C[fastk]
            return absSlow, absFast
        if z_sp==0:
            return None
        from cv2 import dilate
        from .numeric import nan_helper
        if meanSlow2Var is None:
            var = absSlow
        else:
            var = meanSlow2Var(absSlow)
#         print (hasattr(self,"gain"))
        if hasattr(self,"gain"):
            var = var*self.gain    
        std = var**.5
        ff = (absFast > z_sp*std).astype("uint8")
        if ff.any():
            dilateKernelSize = int(ironTimeScale/minDt/Nrebin*.3)#*.03)
            if dilateKernelSize%2==0:
                dilateKernelSize+=1
            print (dilateKernelSize)
            if dilateKernelSize>=3:
                ff = dilate(ff, np.ones(dilateKernelSize, dtype = np.uint8).reshape(1,-1)).astype("bool")
            absFast_tmp = absFast.copy()
            absFast_tmp[ff] = np.nan
            for j in range(C.shape[0]):
                y = absFast_tmp[j]
                nans, x= nan_helper(y)
                if nans.any(): 
                    y[nans]= np.interp(x(nans), x(~nans), y[~nans])
            dFast = self.sosFilter.run(absFast_tmp)
            absSlow = absSlow + dFast
            absFast = absFast - dFast
        std = absSlow**.5
        zScore = absFast/std
        C["zScore_%g"%ironTimeScale] = list(zScore)
        C["slower_%g"%ironTimeScale] = [absSlow[i]/C["size"].iloc[i]/Nrebin for i in range(len(C))]
        C["faster_%g"%ironTimeScale] = [absFast[i]/C["size"].iloc[i]/Nrebin for i in range(len(C))]

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
        
    def calc_peaks(self, ts, z_th = 3, Npoints=None, smooth=None, tWindow = (None,None)):
        from .numeric import runningAverage
        from scipy.signal import find_peaks, peak_widths
        if "zScore_%g"%ts not in self.df.columns:
            self.fast_filter_traces(ts,Npoints=Npoints)
        zScores = np.vstack(self.df["zScore_%g"%ts])
        try:
            t = self.showTime["%g"%ts]
        except:
            t = self.time
        dt = np.diff(t).mean()
        if smooth is None:
            smooth = int(ts/np.diff(t).mean()/5)
            if smooth%2==0: smooth += 1
        if smooth>0:
            zScores = runningAverage(zScores.T,smooth).T
        peaks = []
        for i,z in zip(self.df.index,zScores):
            pp = find_peaks(z,
                            width=ts/dt/5,
                            height=z_th
                              )
            w,h,x0 = peak_widths(z, pp[0], rel_height=.5)[:3]
            w = w*(t[1]-t[0])
            x0 = x0*(t[1]-t[0])
            df = pd.DataFrame({"peak height [z-score]":z[pp[0]],"peak half-width [s]":w, "t0":x0})
            df["roi"] = i
            peaks += [df]
        peaks = pd.concat(peaks,ignore_index=True)
        k = "%g"%(ts)
        try:
            self.peaks
        except:
            self.peaks = {}
        self.peaks[k] = peaks
        
    def show_scatter_peaks(self,ts,timeWindows = None):
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
        self.protocol = pd.read_csv(pathToProtocol)
        self.protocol["t_begin"] = pd.to_timedelta(["00:"+el if type(el)==str else "00:00:00" \
                                                       for el in self.protocol["begin"]]).total_seconds()
        self.protocol["t_end"] = pd.to_timedelta(["00:"+el if type(el)==str else (self.time[0]+len(self.time)/self.Freq)*1e9 \
                                                     for el in self.protocol["end"]]).total_seconds()

    def examine(self):
        from .examine import examine
        return examine(self)
    def examine1(self, test=False):
        from .examine1 import examine
        return examine(self, test=test)
#         from .examine_old import examine_old
#         return examine_old(self)
    
    

def getGraph_of_ROIs_to_Merge(df,rreg, plot=False, ax=None,lw=.5):
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
            ax.arrow(y0,x0,dy,dx,width = .5,
                     linewidth = lw,
                     color="r",
                     zorder=10,
                     length_includes_head=True)
            
            #ax.plot(y1,x1,"o",ms=5,mfc="none",mew=.7,c=c)
            
    if plot:
        plt.gca().set_aspect("equal")
        
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

def mergeBasedOnGraph(Gph,rreg):
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
    print (f"{len(toDrop)} subsumed into existing ROIs.")
    rreg.update()
    rreg.sortFromCenter()
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
            bd = C.loc[j,"boundary"]
            x = np.linalg.norm(np.array(bd)-np.repeat([pk],len(bd),axis=0), axis=1)
            dists[j] = x.min()
        jmin = pd.Series(dists).idxmin()
        peak2bnd += [(i,jmin,dists[jmin])]

    peak2bnd = pd.DataFrame(peak2bnd, columns=["i","j","dist"])
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

def getStatImages(movie_, debleach=True, downsampleFreq=10):
    if movie_.fr>downsampleFreq:
        n_rebin = int(movie_.fr/downsampleFreq)
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

    m_for_image = np.diff(m_for_image,axis=0)
    for f in [np.mean,np.std]:
        statImages["diff_"+f.__name__] = f(m_for_image,axis=0)

    return statImages
