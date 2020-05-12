from general_functions import multi_map
import numpy as np
import pandas as pd
from itertools import product
from collections import OrderedDict
import matplotlib.pyplot as plt
import networkx as nx


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


def climb(x,blurredWeights,diag=True,min_gradient = 0):
    dims = blurredWeights.shape
    # x = (60,60)
    x = x+(blurredWeights[x[0],x[1]],)
    xs = [x]
    for i in range(100):
        vs = []
        for di,dj in product([-1,0,1],[-1,0,1]):
            if not diag:
                if di*dj!=0: continue
            i,j = x[0]+di,x[1]+dj
            if i<0 or i>=dims[0] or j<0 or j>=dims[1]:
                continue
            vs += [(i,j,blurredWeights[i,j])]
        x1 = vs[np.argmax(vs,axis=0)[-1]]
        dx = x1[-1]-x[-1]
        if dx<=0:
            break
        elif dx<=min_gradient:
            return (None,)
        else:
            x = x1
            xs += [x]
    return x[:2]

# def crawlDict(image, th=-np.inf, diag=False, min_gradient=0):
#     A_ = [(i,j)+climb((i,j),image,diag=diag,min_gradient=min_gradient) for i,j in product(range(image.shape[0]),range(image.shape[1])) if image[i,j]>th]
#     A_ = [el for el in A_ if el[-1] is not None]
#     B_ = OrderedDict()
#     for (i0,j0,i1,j1) in A_:
#         if (i1,j1) not in B_:
#             B_[(i1,j1)] = []
#         B_[(i1,j1)] += [(i0,j0)]
#     return B_


def crawlDict(image, crawl_th=-np.inf, diag=False, min_gradient=0, n_processes=10):
    global iterf
    ijs_ = [(i,j) for i,j in product(range(image.shape[0]),range(image.shape[1])) if image[i,j]>crawl_th]
    def iterf(ij):
        return climb((ij[0],ij[1]), image, diag=diag, min_gradient=min_gradient)
    R_ = multi_map(iterf,ijs_, processes=n_processes)
    A_ = [ij+r for ij,r in zip(ijs_,R_)]
    A_ = [el for el in A_ if el[-1] is not None]
    B_ = OrderedDict()
    for (i0,j0,i1,j1) in A_:
        if (i1,j1) not in B_:
            B_[(i1,j1)] = []
        B_[(i1,j1)] += [(i0,j0)]
    return B_

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
    m_for_image = movie_.astype("float")
    if m_for_image.fr>downsampleFreq:
        from physio_def_1 import rebin
        n_rebin = int(m_for_image.fr/10)
        if n_rebin>=2:
            m_for_image = rebin(m_for_image,n_rebin)
    statImages = {}
    if debleach:
        m_for_image.debleach()

    for f in [np.mean,np.std]:
        statImages[f.__name__] = f(m_for_image,axis=0)

    m_for_image = np.diff(m_for_image,axis=0)
    for f in [np.mean,np.std]:
        statImages["diff_"+f.__name__] = f(m_for_image,axis=0)

    return statImages


class Regions:
    def __init__(self, movie_, crawl_th=0, diag=False, min_gradient=0, debleach=True, gSig_filt=None, mode="diff_std",full=True, img_th = -np.inf):
        self.mode = mode
        if isinstance(movie_, (np.ndarray,)):
            if len(movie_.shape)==2:
                image0 = movie_
                self.statImages = {mode:image0}
            if len(movie_.shape)==3:
                self.movie = movie_
                self.statImages = getStatImages(movie_, debleach=debleach, downsampleFreq=10)
        elif isinstance(movie_, dict):
            self.statImages = movie_
        else:
            raise ValueError("Regions can initialize either from a movie, or image, or a dictionary of stats images.")
            
        for k0 in self.statImages:
            break
        tmp = np.zeros_like(self.statImages[k0])
        for submode in mode.split("+"):
            im = self.statImages[submode].copy()
            im = im/im.max()
            tmp += im
        self.statImages[mode] = tmp/tmp.max()
        image0 = tmp/tmp.max()
        
        if full:
            self.constructRois(image0, crawl_th=crawl_th, img_th=img_th, diag=diag, min_gradient=min_gradient, gSig_filt=gSig_filt)
            
    def constructRois(self, image0, crawl_th=0, img_th=-np.inf, diag=False, min_gradient=0, gSig_filt=None):
        from caiman.motion_correction import high_pass_filter_space
        from pandas import DataFrame
        toMin = image0<img_th
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
        B_ = crawlDict(image,crawl_th,diag=diag,min_gradient=min_gradient)
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
                B_ = crawlDict(image,crawl_th,diag=diag,min_gradient=min_gradient)
            
        self.df = DataFrame(OrderedDict([
            ("peak",  list(B_.keys())),
            ("pixels",list(B_.values()))
        ]))
        self.df["peakValue"] = [image[p] for p in B_]
        self.update()
#         print(f"Initialized with {len(self.df)} rois.")

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
            self.time = np.arange(len(movie_))/movie_.fr
    
#     @property
#     def peak(self):
#         return self.df["peak"]
#     @property
#     def pixels(self):
#         return self.df["pixels"]
#     @property
#     def edges(self):
#         return self.df["edges"]
    
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
    
#     @property
    def getEdges(self,ix=None):
        if ix is None:
            out = sum(self.df.edges,[])
        else:
            out = sum(self.df.loc[ix,"edges"],[])
        out = np.unique(out,axis=0)
        return out
    
####### hidden
    # def plotEdges(self, ix=None, ax=None, image =True, imkw_args = {"cmap":"bwr"}, separate=False, color="k"):
    #     from matplotlib.colors import LogNorm
    #     if ax is None:
    #         ax = plt.subplot(111)
    #     if image:
    #         v = np.abs(self.image).max()*3
    #         ax.imshow(self.image,vmin=-v,vmax=v,**imkw_args)
    #     if separate:
    #         for i in ix:
    #             c = "C%i"%(i%10)
    #             for el in self.df.loc[i,"edges"]:
    #                 y,x = np.array(el).reshape((2,2)).T
    #                 ax.plot(x,y,lw=.8,c=c)
    #     else:
    #         for el in self.getEdges(ix=ix):
    #             y,x = np.array(el).reshape((2,2)).T
    #             ax.plot(x,y,color,lw=.5)
    
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
            
            
    def plotPeaks(self, ix=None, ax=None, image=False, ms=1, labels=False,color=None, imkw_args={},absMarker=False):
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
    
    def calcTraces(self, movie_=None):
        if movie_ is None:
            movie_ = self.movie
        traces = np.ones((len(self.df),len(movie_)))*np.nan
        for i,ix in enumerate(self.df.index):
            x = self.df.loc[ix,"pixels"]
            x = [ el[0] for el in x ] , [ el[1] for el in x ]
            traces[i] = movie_[ :, x[0], x[1] ].mean(axis=1)
        print (traces.shape, )
        self.df["trace"] = list(traces)
        self.time = np.arange(len(movie_))/movie_.fr
        
    def sortFromCenter(self):
        center = np.array(self.image.shape)/2
        self.df["distToCenter"] = [np.sum((np.array(self.df.loc[i,"peak"])-center)**2)**.5 for i in self.df.index]
        self.df.sort_values("distToCenter",inplace=True)
        self.df.index = np.arange(len(self.df))
        self.calcNNmap()

    def filter_traces(self,ironScale, n_processes = 10, percentile = [10.,50.], calcStd=False):
        from general_functions import multi_map
        from numeric import lowPass
        global iterf
        try:
            freq = self.movie.fr
        except:
            freq = 1./np.diff(self.time).mean()
        wIron = int(ironScale*freq)
        if wIron%2==0:
            wIron += 1
        print(f"The movie frequency is {freq:.2f}, so the filter size is {wIron}. This may take some time.")
        def iterf(x_): 
            out = lowPass(x_,wIron,wIron,percentile)
            return out
            
        self.df["slower_%g"%ironScale] = multi_map(iterf,self.df["trace"].values,processes=n_processes)
        self.df["faster_%g"%ironScale] = [self.df.loc[i,"trace"] - \
                                          self.df.loc[i,"slower_%g"%ironScale] for i in self.df.index]
            
        def iterf(x_):
            mad2std = 1.4826
            out = mad2std*lowPass(np.abs(x_),wIron,wIron*2+1,50.)
            return out
        if calcStd:
            self.df["faster_%g_std"%ironScale] = multi_map(iterf,self.df["faster_%g"%ironScale].values,processes=n_processes)

    # def filter_fancy(self,ironScale,n_processes = 10, fancy=False):
    #     from general_functions import multi_map
    #     from numeric import lowPass
    #     global iterf
    #     wIron = int(ironScale*self.movie.fr)
    #     if wIron%2==0:
    #         wIron += 1
    #     print(f"The movie frequency is {self.movie.fr:.2f}, so the filter size is {wIron}. This may take some time.")
    #     filterPars = (wIron,wIron,percentile,2)
    #     if fancy:
    #         print ("Thank you for your trust. I am ignoring the percentile value supplied by you.")
    #         filterPars = (wIron,wIron,10,nRepeats)
    #     def iterf(x_): 
    #         out = lowPass(x_,*filterPars)
    #         return out
    #     self.df["slower_%g"%ironScale] = multi_map(iterf,self.df["trace"].values,processes=n_processes)
    #     self.df["faster_%g"%ironScale] = [self.df.loc[i,"trace"] - \
    #                                       self.df.loc[i,"slower_%g"%ironScale] for i in self.df.index]

    #     if fancy:
    #         from numeric import runningMode
    #         def iterf(x_):
    #             dx = runningMode(x_,wIron*2+1,20)
    #             dx = lowPass(dx,wIron*2+1,wIron*2+1,50,1)
    #             out = x_-dx
    #             return out
    #         self.df["faster_%g"%ironScale] = multi_map(iterf,self.df["faster_%g"%ironScale].values,processes=n_processes)
    #         self.df["slower_%g"%ironScale] = [self.df.loc[i,"trace"] - \
    #                                           self.df.loc[i,"faster_%g"%ironScale] for i in self.df.index]

    #     def iterf(x_):
    #         mad2std = 1.4826
    #         out = mad2std*lowPass(np.abs(x_),wIron,wIron*2+1,50.,1)
    #         return out
    #     self.df["faster_%g_std"%ironScale] = multi_map(iterf,self.df["faster_%g"%ironScale].values,processes=n_processes)
        
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
#         del intraCCs

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
            iaccept = C.loc[[i,j],"peakValue"].idxmax()
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
        nx.draw_networkx(gph,ax=ax,node_size=30,node_color="w",
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
        if len(attr)>1:
            print ("more than two attractors, not implemented yet, will skip")
            continue
        attr = attr[0]
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