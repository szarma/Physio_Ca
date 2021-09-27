import warnings
import json
import os
import pickle
from collections import OrderedDict
from copy import deepcopy
from datetime import date
from itertools import product
from pathlib import Path
from sys import exc_info
from typing import Union
import matplotlib.pyplot as plt
import mgzip
import networkx as nx
import numpy as np
import orjson
import pandas as pd
from scipy.stats import median_abs_deviation
import plotly.graph_objects as go
from matplotlib.colors import LogNorm
from plotly_express import colors as plc

from .general_functions import getCircularKernel
from .numeric import rebin, bspline, crawl_dict_via_graph
from .utils import multi_map

MYCOLORS = plc.qualitative.Plotly
# MYCOLORS = ["darkred"]

def load_regions(path,
                 mergeDist=1,
                 mergeSizeTh=10,
                 plot=False,
                 verbose=False,
                 calcInterest=True,
                 baremin=False
                ):
    with open(path,"rb") as f:
        regions = pickle.load(f)
    try:
        regions.update()
        regions.pathToRois = path
        pickleDir = os.path.split(path)[0]
        regions = Regions(regions)
        try:
            protocolFile = os.path.join(pickleDir, [f for f in os.listdir(pickleDir) if "protocol" in f][0])
            regions.import_protocol(protocolFile)
        except:
            pass
        if not baremin:
            regions.detrend_traces()
            regions.infer_TwoParFit(plot=plot)
            if hasattr(regions,"gain"):
                del regions.gain
        regions.merge_closest(mergeSizeTh=mergeSizeTh, mergeDist=mergeDist, plot=plot, Niter=15, verbose=verbose)
        if calcInterest:
            regions.calc_interest()
    except:
        print ("encountered error:", exc_info())
    
    return regions


class Regions:
    def __init__(self, movie_,
                 diag=False,
                 debleach=False,
                 gSig_filt=None,
                 mode="highperc+mean",
                 full=True,
                 img_th=None,
                 FrameRange=None,
                 #processes=7,
                 #excludePixels=None,
                 verbose=False,
                 #use_restricted=None
                ):
        if hasattr(movie_, "df") and "pixels" in movie_.df.columns:
            for k in movie_.__dict__.keys():
                if verbose:
                    print ("Initiating from another Regions object.")
                setattr(self, k, movie_.__dict__[k])
            self.update()
        elif isinstance(movie_, (np.ndarray,)):
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
                return None
                
            elif isinstance(akey,str):
                if verbose:
                    print ("Initiating from a dictionary assumed to be a dictionary of image stats.")
                self.statImages = movie_
            else:
                raise ValueError("Initializing Regions from a dictionary is only supported for a dictionary of images representing movie statistics, or a pixel crawling dictionary.")
        else:
            raise ValueError("Regions can initialize either from a movie, or an image, or a dictionary. You supplied %s"%str(type(movie_)))
        self.mode = mode
        if gSig_filt is not None:
            if type(gSig_filt) == int:
                gSig_filt = [gSig_filt]
        self.filterSize = gSig_filt
        self.image = self.defineImage(mode=mode, gSig_filt=self.filterSize)

        if full and not hasattr(self,"df"):
            self.constructRois(image=self.image, img_th=img_th, diag=diag, verbose=verbose,
                               #processes=processes,
                               #excludePixels=excludePixels, 
                               #use_restricted=use_restricted
                              )

    def to_json(self,
                output_dir: Union[str, Path],
                file_name: str = '',
                movie=None,
                col: Union[list, None] = None,
                formats: Union[list, None] = None,
                add_date: bool = True,
                use_compression: bool = False) -> None:
        """
        Exports the current regions object to json-format.

        :param output_dir: Directory, where the output will be written to.
        :param file_name: Base filename (will be extended automatically e.g. "5.6" becomes "5.6_rois.json.gz"
        :param movie: Movie, which shall be manually connected to the regions object.
        :param col: Columns, which are required to serialize. Default: ['trace']
        :param formats: Format, which should be used to store json-data. Default: ['vienna']
        :param add_date: If True then the current date will be added to the filename. Default: True
        :param use_compression: If True then the .json file will be compressed to .gz
        :return: None
        """

        if col is None:
            col = ['trace']
        if formats is None:
            formats = ['vienna']
        if type(output_dir) == str:
            output_dir = Path(output_dir)
        if movie is not None:
            self.update(movie)
        file_name = file_name.replace(' ', '_')
        if add_date:
            today = date.today()
            if len(file_name):
                file_name = '_'.join([today.strftime('%Y_%m_%d'), file_name])
            else:
                file_name = today.strftime('%Y_%m_%d')
        if not output_dir.is_dir():
            output_dir.mkdir(parents=True)

        for format in formats:
            if format == 'vienna':
                saving = ['statImages',
                          'mode',
                          'image',
                          'filterSize',
                          'df',
                          'trange',
                          'FrameRange',
                          'analysisFolder',
                          'time',
                          'Freq',
                          'metadata']
                juggle_movie = hasattr(self, 'movie')
                if juggle_movie:
                    movie = self.movie
                    del self.movie
                all_keys = list(self.__dict__.keys())
                subregions = deepcopy(self)
                if juggle_movie:
                    self.movie = movie
                    del movie
                for k in all_keys:
                    if k not in saving:
                        del subregions.__dict__[k]
                for k in self.df.columns:
                    if k not in ['peak', 'pixels', 'peakValue','tag','interest'] + col:
                        del subregions.df[k]

                json_dict = {}
                for k, v in subregions.__dict__.items():
                    value = v
                    if isinstance(v, (pd.DataFrame, pd.Series)):
                        value = json.JSONDecoder().decode(v.to_json(double_precision=15))
                    if isinstance(v, (np.float64, np.int64)):
                        value = str(v)
                    if isinstance(v, np.ndarray):
                        value = v.tolist()
                    if isinstance(v, dict):
                        for k_, v_ in v.items():
                            if isinstance(v_, np.ndarray):
                                # v[k_] = json.dumps(v_.tolist())
                                v[k_] = v_.tolist()
                        value = v

                    json_dict[k] = value

                json_string = orjson.dumps(json_dict).decode()

                roi_file = f'{output_dir.as_posix()}/{file_name}_rois.json'
                if use_compression:
                    if not roi_file.endswith('.gz'):
                        roi_file += '.gz'
                    with mgzip.open(roi_file, 'wt') as file:
                        file.write(json_string)
                else:
                    with open(roi_file, 'wt') as file:
                        file.write(json_string)

                self.pathToPickle = roi_file
                # create_preview_image(self)
                # TODO: Add json to creation process to let it work
            elif format == 'maribor':
                raise NotImplementedError
            else:
                raise NotImplementedError(f'Output format {format} not recognized.')

    @staticmethod
    def from_json(file_path: Union[str, Path], use_compression: Union[bool, None] = None) -> object:
        """
        Static method, which enabled creation of a regions object by reading a (compressed) json file.

        :param file_path: Path to the file, which will be read.
        :param use_compression: Specifies the usage of compression. None means the algorithm decides it by file suffix.
        :return: Regions object.
        :raises FileNotFoundError: Gets raised if file_path does not exist.
        """
        if type(file_path) == str:
            file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError

        if use_compression is None:
            use_compression = True if file_path.suffix == '.gz' else False

        if use_compression:
            with mgzip.open(file_path.as_posix(), 'rt') as file:
                file_content = file.read()
        else:
            with open(file_path.as_posix(), 'rt') as file:
                file_content = file.read()

        # Load JSON-object from file-content:
        json_obj = orjson.loads(file_content.encode())

        # Restoring data types:
        df = pd.DataFrame.from_dict(json_obj['df'])
        df['peak'] = df['peak'].apply(lambda x: tuple(x))
        df['pixels'] = df['pixels'].apply(lambda x: [tuple(y) for y in x])
        del json_obj['df']

        json_obj['metadata'] = pd.Series(json_obj['metadata'])
        json_obj['Freq'] = np.float64(json_obj['Freq'])
        json_obj['image'] = np.array(json_obj['image'])

        for key, value in json_obj['statImages'].items():
            json_obj['statImages'][key] = np.array(value)

        # Create Regions object and set all attributes stored in the JSON file:
        regions = Regions(dict(zip(df['peak'], df['pixels'])))
        for key, value in json_obj.items():
            regions.__setattr__(key, value)
        regions.df = df
        regions.update()
        pickle_dir = file_path.parent
        regions = Regions(regions)
        protocol_files = list(file_path.parent.glob('*protocol*.*'))
        if len(protocol_files) > 0:
            regions.import_protocol(protocol_files[0].as_posix())
        regions.pathToPickle = file_path.as_posix()
        regions.detrend_traces()
        regions.infer_gain(plot=False)
        regions.merge_closest(Niter=15)

        return regions


    def defineImage(self, mode, gSig_filt=None, ):
        from .numeric import robust_max
        if mode == "custom":
            image0 = self.statImages[mode]
        else:
            k0 = next(iter(self.statImages))
            tmp = np.zeros_like(self.statImages[k0])
            norm = []
            for submode in mode.split("+"):
                im = self.statImages[submode].copy()
                norm += [robust_max(im)]
                im = im / norm[-1]
                tmp += im
            image0 = tmp / len(norm) * np.mean(norm)
            if mode not in self.statImages:
                self.statImages[mode] = image0

        from cv2 import GaussianBlur
        if gSig_filt is None:
            image = image0
        else: # gSig_filt is a list of integers
            image = []
            for gSf in gSig_filt:
                tmp = GaussianBlur(image0, (gSf // 2 * 2 + 1,) * 2, -1)
                tmp -= GaussianBlur(image0, (gSf * 2 + 1,) * 2, -1)
                tmp *= gSf
                image += [tmp / robust_max(tmp)]
            image = np.mean(image, axis=0)
        return image


    def constructRois(self, image, img_th=None, dks=3, verbose=False, diag=False):
        from cv2 import dilate,erode
        try:
            dks = max(3, (max(self.filterSize)) // 4 * 2 + 1)
        except:
            pass
        dilation_kernel = getCircularKernel(dks)
        eks = max(3,dks-2)
        erosion_kernel  = getCircularKernel(eks)
        if verbose:
            print ("eroding valid pixels by", eks)
        if img_th is None:
            img_th = median_abs_deviation(image.flat)
        ok = erode((image>img_th).astype(np.uint8), erosion_kernel)
        if verbose:
            print ("dilating valid pixels by", dks)
        ok = dilate(ok, dilation_kernel)
        ok = ok.astype(bool)
        self.validPixels = ok
        B_ = crawl_dict_via_graph(image, ok, diag=diag)
        if diag:
                from .utils import split_unconnected_rois
                B_ = split_unconnected_rois(B_, self.image)
        
        self.df = pd.DataFrame(OrderedDict([
            ("peak",  list(B_.keys())),
            ("pixels",list(B_.values()))
        ]))
        self.df["peakValue"] = [image[p] for p in B_]
        self.update()
        self.merge_closest(mergeDist=.71, mergeSizeTh=100, verbose=verbose)
        
#     def constructRois(self, image, img_th=None, diag=True, dks=3, processes=5, excludePixels=None, verbose=False,use_restricted=True):
#         from cv2 import dilate,erode
#         try:
#             dks = max(3, (max(self.filterSize)) // 4 * 2 + 1)
#         except:
#             pass
#         if not use_restricted:
#             if excludePixels is None:
#                 excludePixels = []
#         dilation_kernel = getCircularKernel(dks)
#         eks = max(3,dks-2)
#         erosion_kernel  = getCircularKernel(eks)
#         if verbose:
#             print ("eroding valid pixels by", eks)
#         if img_th is None:
#           img_th = median_abs_deviation(image.flat)
#         ok = erode((image>img_th).astype(np.uint8), erosion_kernel)
#         if verbose:
#             print ("dilating valid pixels by", dks)
#         ok = dilate(ok, dilation_kernel)
#         # ok = ok*(self.statImages[self.mode]>0).astype(np.uint8)
#         ok = ok.astype(bool)
#         self.validPixels = ok
#         if use_restricted=="break":
#             raise ValueError("Debugging mode")

#         if use_restricted:
#             includePixels = list(map(tuple, np.array(np.where(ok)).T))
#             if verbose:
#                 print(f"initiating the cralwing dict on {len(includePixels)} (%.1f%%) pixels only."%(100*len(includePixels)/ok.size))
#             B_ = crawlDict_restr(image,
#                            includePixels,
#                            diag=diag,
#                            processes=processes,
#                            verbose=verbose
#                           )
#         else:
#             excludePixels += list(map(tuple, np.array(np.where(~ok)).T))
#             if verbose:
#                 print(f"initiating the cralwing dict with {len(excludePixels)} pixels excluded.")
#             B_ = crawlDict(image,
#                            crawl_th=-np.inf,
#                            diag=diag,
#                            processes=processes,
#                            excludePixels=excludePixels,
#                            verbose=verbose
#                           )
#         if diag:
# #             try:
#                 from .utils import split_unconnected_rois
#                 B_ = split_unconnected_rois(B_, self.image)
# #             except:
# #                 print ("Cannot initialize with diagonal crawl. Reverting to diag=False")
# #                 diag = False
# #                 B_ = crawlDict(image,image.min(),diag=diag,min_gradient=min_gradient, processes=processes)
        
        
#         self.df = pd.DataFrame(OrderedDict([
#             ("peak",  list(B_.keys())),
#             ("pixels",list(B_.values()))
#         ]))
# #         try:
# #             blurKsize = min(self.filterSize)//2*2+1
# #             blurKsize = max(3,blurKsize)
# #             slightly_blurred_image = GaussianBlur(self.statImages[self.mode],(blurKsize,)*2,-1)
# #             self.reassign_peaks(slightly_blurred_image+self.statImages[self.mode])
# #         except:
#         self.df["peakValue"] = [image[p] for p in B_]
#         self.ExcludePixels = excludePixels
#         self.update()
    
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
        from .numeric import fit_bleaching, rebin
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
        ydbl = fit_bleaching(y)
        self.fov_trace = {
                "time": x,
                "raw": y,
                "trend":ydbl
            }

    def merge_closest(self, mergeSizeTh=10, mergeDist=1, plot=False, Niter=20, verbose=False, axs=None):
        if plot:
            if axs is None:
                plt.figure(figsize=(7*Niter,6))

        ia = 1
        for ja in range(Niter):
            size_th = np.percentile(self.df["size"].values, mergeSizeTh)
            df = getPeak2BoundaryDF(self.df)
            df = df.query(f"dist<={mergeDist} and size_from<={size_th}")
            if len(df):
                if plot:
                    if axs is None:
                        ax = plt.subplot(1,Niter,ia)
                        ax.imshow(self.statImages[self.mode], cmap="Greys", norm=LogNorm())
                    else:
                        ax = axs[ia%len(axs)]
                    xl = ax.get_xlim()
                    yl = ax.get_ylim()
                else:
                    ax = None
                suggestGraph = getGraph_of_ROIs_to_Merge(df.iloc[:,:2], self, plot=plot,ax=ax)
                if plot:
                    ax.set_xlim(xl)
                    ax.set_ylim(yl)
                mergeBasedOnGraph(suggestGraph, self, verbose=verbose)
            else:
                # print ("No more suggestions.")
                break
            ia += 1
        if plot:
            plt.tight_layout()
    
        
    def update(self, movie_=None):
        self.df["size"] = self.df["pixels"].apply(len)
        self.calcEdges()
        self.df["boundary"] = [edges2nodes(self.df["edges"][j]) for j in self.df.index]
        self.calcNNmap()
        self.df["interest"] = [np.sum([self.image[px[0],px[1]] for px in pxs]) for pxs in self.df["pixels"]]
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
    
    def plotEdges(self,
                  ix=None,
                  ax=None,
                  image=True,
                  imkw_args=None,
                  separate=False,
                  color="darkred",
                  lw=None,
                  alpha=1,
                  fill=False,
                  scaleFontSize=12,
                  norm=LogNorm(vmin=1),
                  spline=True
                  ):
        if imkw_args is None:
            imkw_args = {}
        if ix is None:
            ix = self.df.index
        if ax is None:
            ax = plt.subplot(111)
        if lw is None:
            lw=.5
        if "cmap" not in imkw_args:
            from copy import copy
            mycmap = copy(plt.cm.Greys)
            mycmap.set_bad("lime")
            imkw_args["cmap"] = mycmap
        if image:
            im = self.statImages[self.mode]
            ax.imshow(im,norm=norm,**imkw_args)
        smoothness = min(self.__dict__.get("filterSize",[])+[3] )
        if separate:
            for i in ix:
                try:
                    c = self.df.loc[i,"color"]
                except:
                    c = MYCOLORS[i%len(MYCOLORS)]
                points = self.df.loc[i,"boundary"]+self.df.loc[i,"boundary"][:3]
                if spline:
                    points = bspline(points, smoothness=smoothness)
                points = [(p-points.mean(0))*.9+points.mean(0) for p in points]
                y,x = np.array(points).T
                ax.plot(x,y,"-",lw=lw,c=c,alpha=alpha)
                if fill:
                    ax.fill(x,y,c=c,alpha=alpha*.8)
        else:
            tmp = []
            for el in self.df.loc[ix,"boundary"]:
                if spline:
                    el = list(bspline(el, smoothness=smoothness))
                tmp += el
                tmp += el[:3]
                tmp += [(np.nan,)*2]

            y,x = np.array(tmp).T
            ax.plot(x,y,color,lw=lw,alpha=alpha)
        dim = self.image.shape
        ax.set_xlim(-.5,dim[1]-.5)
        ax.set_ylim(dim[0]-.5, -.5,)
            
        if scaleFontSize<=0: return None
        if hasattr(self, "metadata") and "pxSize" in self.metadata:
            lengths = [10,20,50]
            il = np.searchsorted(lengths,self.metadata.pxSize*self.image.shape[1]/10)
            length=lengths[il]
            x0,x1,y0,y1 = np.array([0,length,0,length*3/50])/self.metadata.pxSize + self.image.shape[0]*.02
            ax.fill_between([x0,x1],[y1]*2,[y0]*2, color="k")
            txt = "\n"*1+str(length)
            if "pxUnit" in self.metadata:
                txt += self.metadata["pxUnit"]
            ax.text((x0+x1)/2, y1+.3*(y1-y0), txt, va="center", ha="center", size=scaleFontSize)
            
    def plotPeaks(self, ix=None, ax=None, image=False, ms=3, labels=False,color=None, imkw_args={},absMarker=True, marker="."):
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
            ax.plot(*p[::-1],marker=marker,ms=ms,c=c)
            if labels:
                ax.text(*p[::-1],s=" "+str(i),color=c)
    
    def calc_interest(self, zth=4, timescales=[3,10,30,100,300]):
        interesting = np.zeros(len(self.df))
        for ts in [3,10,30,100,300]:
            if 1./self.Freq   > ts/10: continue
            if self.time[-1] < ts*5: continue
            s,f,z = self.fast_filter_traces(ts, write=False)
            interesting += np.mean(z>zth,1)
        self.df["interest"] = interesting
        
    def change_frequency(self, fr=2):
        from .movies import movie as cmovie
        traces = np.vstack(self.df.trace)
        fr0 = self.Freq
        trmov = cmovie(traces.T.reshape((-1,len(self.df),1)), fr=fr0)
        trmov = trmov.resize(1,1,fr/fr0)*fr0/fr
        self.df.trace = list(trmov[:,:,0].T)
        self.Freq = fr
        self.time = np.arange(len(trmov))/fr
        
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
        else:
            self.movie = movie_
        if FrameRange is None:
            try: FrameRange = self.FrameRange
            except: FrameRange = [0,len(movie_)]
        else:
            self.FrameRange = FrameRange
        i0,ie = FrameRange
        traces = np.ones((len(self.df),(ie-i0)))*np.nan
        for i,ix in enumerate(self.df.index):
            x = self.df.loc[ix,"pixels"]
            x = [ el[0] for el in x ] , [ el[1] for el in x ]
            traces[i] = movie_[i0:ie, x[0], x[1] ].mean(axis=1)
        self.df["trace"] = list(traces)
        time = np.arange(len(movie_))/movie_.fr
        self.time = time[i0:ie]
        self.Freq = movie_.fr
    
    def reset_filtered(self):
        for col in self.df.columns:
            if "slower" in col or "faster" in col or "zScore" in col:
                del self.df[col]
        self.showTime = {}
        
    def detrend_traces(self,method="debleach", timescale=None):
#         from .numeric import mydebleach
#         traces = np.vstack(self.df.trace.values)
#         trend = multi_map( mydebleach, traces, processes=processes)
#         self.df["trend"] = trend
#         self.df["detrended"] = list(traces - np.array(trend))
#         traces = np.vstack(self.df.trace.values)
#         print ("Deprecated method. I used to think it makes sense, now I don't think so. The method is still here, not to break your scripts, but it only subtracts mean from the trace.")
        if method=="simple":
            trend = self.df.trace.apply(np.mean)
            self.df["trend"] = trend
            self.df["detrended"] = [self.df.trace[i] - self.df.trend[i] for i in self.df.index]
        elif method=="fast":
            self.fast_filter_traces(timescale)
            self.df["detrended"] = self.df["faster_%g"%timescale]
            self.df["trend"]     = self.df["slower_%g"%timescale]
            del self.df["faster_%g"%timescale], self.df["slower_%g"%timescale]
        elif method=="slow":
            self.slow_filter_traces(timescale)
            self.df["detrended"] = self.df["faster_%g"%timescale]
            self.df["trend"]     = self.df["slower_%g"%timescale]
            del self.df["faster_%g"%timescale], self.df["slower_%g"%timescale]
        elif method=="debleach":
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                from .numeric import fit_bleaching
                self.df["trend"] = multi_map(fit_bleaching, self.df.trace, processes=5)
                self.df["detrended"] = [self.df.trace[i] - self.df.trend[i] for i in self.df.index]
        else:
            raise NotImplementedError()

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

    def infer_gain(self, ts=None, plot=False, verbose=False, ret_points=False):
        if ts is None:
            minDt = np.diff(self.time).mean()
            freq = 1 / minDt
            ts = 30 / freq
        absSlow, absFast, _ = self.fast_filter_traces(ts, write=False, normalize=False, filt_cutoff=0)
        di = 30
        slow_est, fast_vars = [], []
        for i in range(absFast.shape[0]):
            for j in range(di, absFast.shape[1] - di, absFast.shape[1] // 30):
                slow_est += [absSlow[i, j]]
                fast_vars += [absFast[i, j - di:j + di].var()]
        fast_vars = np.array(fast_vars)
        slow_est = np.array(slow_est)
        logbs = np.log(np.logspace(np.log10(np.percentile(slow_est, 2)), np.log10(np.percentile(slow_est, 98))))
        d = np.digitize(np.log(slow_est), logbs)
        x = np.array([slow_est[d == i].mean() for i in np.unique(d)])
        y = np.array([np.median(fast_vars[d == i]) if (d == i).sum() > 10 else np.nan for i in np.unique(d)])
        x = x[np.isfinite(y)]
        y = y[np.isfinite(y)]
        if ret_points:
            return x, y
        gain = np.mean(y / x)
        slow_est[slow_est <= 0] = np.nan
        if plot:
            ax = plt.subplot(111)
            ax.hexbin(slow_est, fast_vars, bins="log",
                      xscale="log",
                      yscale="log",
                      cmap="Greys",
                      mincnt=1
                      )
            c = ax.plot(x, y, "o", mfc="none")[0].get_color()
            ax.plot(x, x * gain, c=c)

        if verbose: print("initial estimate of the gain is", gain)

        for _ in range(5):
            fast_vars[fast_vars > 5 * gain * slow_est] = np.nan
            if np.isnan(fast_vars).any():
                x = np.array([slow_est[d == i].mean() for i in np.unique(d)])
                y = np.array([np.nanmedian(fast_vars[d == i]) if (d == i).sum() > 10 else np.nan for i in np.unique(d)])
                x = x[np.isfinite(y)]
                y = y[np.isfinite(y)]
                y[y <= 0] = y[y > 0].min()
                gain = np.nanmean(y / x)
                if verbose: print("revised estimate of the gain", gain)
                if plot:
                    c = ax.plot(x, y, "o", mfc="none")[0].get_color()
                    ax.plot(x, x * gain, c=c)
        if plot:
            ax.set_title("gain inference")
            ax.set_xlabel("window means")
            ax.set_ylabel("window variances")
            ax.grid()
        self.gain = gain

    def infer_TwoParFit(self, ts=None, plot=False, verbose=False, ret_points=False):
        if ts is None:
            minDt = np.diff(self.time).mean()
            freq = 1/minDt
            ts = 30/freq
        if verbose:
            print(f"inferring mean2std parameters for timescale {ts}...")
        absSlow, absFast, _ = self.fast_filter_traces(ts,write=False, normalize=False,filt_cutoff=0)
        di = 30
        slow_est, fast_vars = [],[]
        for i in range(absFast.shape[0]):
            for j in range(di, absFast.shape[1]-di, absFast.shape[1]//30):
                slow_est  += [absSlow[i,j]]
                fast_vars += [absFast[i,j-di:j+di].var()]
        fast_vars = np.array(fast_vars)
        slow_est = np.array(slow_est)
        binedges = np.geomspace(np.percentile(slow_est,1), np.percentile(slow_est,99))
        d = np.digitize( slow_est, binedges)
        irange = range(len(binedges)+1)
        x = np.array([slow_est[d==i].mean() for i in irange])
        y = np.array([np.median(fast_vars[d==i]) if (d==i).sum()>10 else np.nan for i in irange])
        x = x[np.isfinite(y)]
        y = y[np.isfinite(y)]
        if ret_points:
            return x,y
        n = np.mean(y/x)
        k = 1
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
            ax.plot(x,n*x**k,c=c)

        if verbose: print ("initial estimate of the k,n is", (k,n))

        for _ in range(10):
            fast_vars[fast_vars>3*n*slow_est**k] = np.nan
            if np.isnan(fast_vars).any():
                x = np.array([slow_est[d==i].mean() for i in irange])
                y = np.array([np.nanmedian(fast_vars[d==i]) if (d==i).sum()>10 else np.nan for i in irange])
                y[y<=0] = y[y>0].min()
                x = x[np.isfinite(y)]
                y = y[np.isfinite(y)]
                k1, n1 = np.polyfit(np.log(x),np.log(y),1)
                n1 = np.exp(n1)
                if np.abs(k1/k-1)<1e-5 and np.abs(n1/n-1)<1e-5:
                    if verbose: print ("converged")
                    break
                else:
                    k,n = k1,n1
                if verbose:
                    print ("revised estimate of k,n", (k,n))
                if plot:
                    # c = ax.plot(x,y,"o",mfc="none")[0].get_color()
                    c = ax.plot([],"o",mfc="none")[0].get_color()
                    ax.plot(x,n*x**k,c=c)
        if plot:
            ax.set_title("2-parameter mean-to-variance inference")
            ax.set_xlabel("window means")
            ax.set_ylabel("window variances")
            ax.grid()

        self.TwoParFit = k,n
    
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
    
    
    
    def prep_filtering(self,
                       ironTimeScale,
                       Npoints=None,
                       write=True,
                       usecol="trace",
                       order=5,
                       verbose=False,
                       ):
        from .numeric import sosFilter
        if Npoints is None: Npoints = 30
        minDt = np.diff(self.time).mean()
        freq = 1/minDt
        if hasattr(self, "movie"):
            if np.abs(freq/self.movie.fr-1)>1e-2:
                print (f"movie frame rate ({self.movie.fr}) and inferred frame ({freq}) rate are different!")
        N_dt = ironTimeScale/minDt
        Nrebin = max(1,int(np.round(N_dt/Npoints)))
        if verbose:
            print (f"Nrebin = {Nrebin}")
        C = self.df
        if "trace" in usecol:
            data = np.vstack([C.loc[i,usecol] for i in C.index])
            trend = np.zeros(len(C))
        elif usecol=="detrended":
            data = np.vstack([C.loc[i,"detrended"] for i in C.index])
            trend = C.trend.values
        else:
            raise ValueError("usecol can only be a 'trace' (or 'trace_*' or 'detrended'")
            
        if Nrebin>1:
            freq = freq/Nrebin
            data = rebin(data, Nrebin, axis=1)
            if len(trend.shape)>1:
                trend = rebin(trend, Nrebin, axis=1)
        if not hasattr(self,"showTime"):
            self.showTime = {}
        if write:
            if Nrebin>1:
                self.showTime["%g"%ironTimeScale] = rebin(self.time, Nrebin)
            else:
                if "%g"%ironTimeScale in self.showTime:
                    del self.showTime["%g"%ironTimeScale]
        cutFreq = .5/ironTimeScale
        sf = sosFilter(cutFreq, freq, order=order)
        if write:
            self.sosFilter = sf
        return {
            "trend":   trend,
            "detrend": data,
            "Nrebin":  Nrebin,
            "filter":  sf,
        }
    
    def mean2std_funct(self, absmean, reg=3):
        absmean = np.maximum(absmean, reg)
        if hasattr(self, "TwoParFit"):
            k,n = self.TwoParFit
            var = n*absmean**k
        elif hasattr(self, "gain"):
            var = absmean*self.gain
        else:
            var = absmean
        return var**.5
    
    def fast_filter_traces(self,
                            ironTimeScale,
                            filt_cutoff=3,
                            Npoints=None,
                            write=True,
                            verbose=False,
                            usecol="trace",
                            normalize=True,
                            npass=3,
                            test_plotting_kwargs=None,
                            dilate=True,
                            ):
        from cv2 import dilate as cvdilate
        C = self.df
        fdata = self.prep_filtering(ironTimeScale,
                                    Npoints=Npoints,
                                    write=write,
                                    usecol=usecol,
                                    order=5
                                    )
        Nrebin = fdata["Nrebin"]
        if verbose:
            print(f"Nrebin={Nrebin}")
        dataFilt = fdata["filter"].run(fdata["detrend"])
        absSlow = Nrebin * np.vstack([ (fdata[  "trend"][i]+dataFilt[i])*C.loc[ix,"size"] for i,ix in enumerate(C.index) ])
        absFast = Nrebin * np.vstack([ (fdata["detrend"][i]-dataFilt[i])*C.loc[ix,"size"] for i,ix in enumerate(C.index) ])
        z = absFast/self.mean2std_funct(absSlow)
        t = self.time
        if Nrebin>1:
            t = rebin(t,Nrebin)
        if test_plotting_kwargs is not None:
            ax = test_plotting_kwargs["ax"]
            iplot = test_plotting_kwargs["i"]
            ax.plot(t, fdata["detrend"][iplot], lw=.5,c="k")
        if filt_cutoff>0:
            detrend = fdata["detrend"].copy()
            xr = np.arange(len(detrend[0]))
            for i in range(npass):
                if verbose:
                    if i==0:
                        print ("passing through traces: 0", end=" ")
                    else:
                        print (i, end=" ")
                slower = np.array([absSlow[i]/C["size"].iloc[i]/Nrebin - fdata[  "trend"][i] for i in range(len(C))])
                #detrend[z>filt_cutoff] = slower[z>filt_cutoff]
                ff = z > filt_cutoff
                if dilate and ff.any():
                    dt = t[1]-t[0]
                    dilts = ironTimeScale/dt*.1
                    dilateKernelSize = int(dilts)#*.03)
                    if dilateKernelSize%2==0:
                        dilateKernelSize+=1
                    if dilateKernelSize>=3:
                        if verbose:
                            print ("dilating by", dilateKernelSize)
                        ff = cvdilate(ff.astype(np.uint8), np.ones(dilateKernelSize, dtype = np.uint8).reshape(1,-1)).astype("bool")
                    # else:
                        # if verbose:
                        #     print (f"dilation timescale {dilts} is too small")
                ff[:,[0,1,-2,-1]] = False
                for i in range(len(detrend)):
                    ffi = ff[i]
                    detrend[i,ffi] = np.interp(xr[ffi], xr[~ffi], slower[i,~ffi])
                if test_plotting_kwargs is not None:
                    ax.plot(t, slower[iplot] , c="darkred")
                    ax.plot(t, detrend[iplot], c="grey")
                dataFilt = fdata["filter"].run(detrend)
                absSlow = Nrebin * np.vstack([ (fdata[  "trend"][i]+dataFilt[i])*C.loc[ix,"size"] for i,ix in enumerate(C.index) ])
                absFast = Nrebin * np.vstack([ (fdata["detrend"][i]-dataFilt[i])*C.loc[ix,"size"] for i,ix in enumerate(C.index) ])
                z = absFast/self.mean2std_funct(absSlow)
            if verbose:
                print("")
        if normalize:
            slower = [absSlow[i]/C["size"].iloc[i]/Nrebin for i in range(len(C))]
            faster = [absFast[i]/C["size"].iloc[i]/Nrebin for i in range(len(C))]
        else:
            slower = list(absSlow)
            faster = list(absFast)
        if write:
            self.df["slower_%g"%ironTimeScale] = slower
            self.df["faster_%g"%ironTimeScale] = faster
            self.df["zScore_%g"%ironTimeScale] = list(z)
        else:
            return np.array(slower), np.array(faster), z
    
#     def fast_filter_traces(self,
#                            ironTimeScale,
#                            filt_cutoff=3,
#                            order=5,
#                            Npoints=None,
#                            write=True,
#                            verbose=False,
#                            usecol="trace",
#                            normalize=True
# #                            meanSlow2Var=None
#                           ):
#         from .numeric import sosFilter
#         if Npoints is None: Npoints = 15
#         minDt = np.diff(self.time).mean()
#         freq = 1/minDt
#         try:
#             self.movie
#             if np.abs(freq/self.movie.fr-1)>1e-2:
#                 print (f"movie frame rate ({self.movie.fr}) and inferred frame ({freq}) rate are different!")
#         except:
#             pass
#         N_dt = ironTimeScale/minDt
#         Nrebin = max(1,int(np.round(N_dt/Npoints)))
#         if verbose:
#             print (f"Nrebin = {Nrebin}")
#         C = self.df
#         if "trace" in usecol:
#             data = np.vstack([C.loc[i,usecol] for i in C.index])
#             trend = np.zeros(len(C))
#         elif usecol=="detrended":
#             data = np.vstack([C.loc[i,"detrended"] for i in C.index])
#             trend = C.trend.values
#         else:
#             raise ValueError("usecol can only be a 'trace' (or 'trace_*' or 'detrended'")
#
#         if Nrebin>1:
#             freq = freq/Nrebin
#             data = rebin(data, Nrebin, axis=1)
#             if len(trend.shape)>1:
#                 trend = rebin(trend, Nrebin, axis=1)
#         if not hasattr(self,"showTime"):
#             self.showTime = {}
#         if write:
#             if Nrebin>1:
#                 self.showTime["%g"%ironTimeScale] = rebin(self.time, Nrebin)
#             else:
#                 if "%g"%ironTimeScale in self.showTime:
#                     del self.showTime["%g"%ironTimeScale]
#         cutFreq = .5/ironTimeScale
#         self.sosFilter = sosFilter(cutFreq, freq, order=order)
#         dataFilt = self.sosFilter.run(data)
#         for j in range(len(data)):
#             dataFilt[j] = np.maximum(dataFilt[j]+trend[j],.5)
#
#
#         absSlow = np.vstack([dataFilt[i]*C.loc[ix,"size"] for i,ix in enumerate(C.index)])*Nrebin
#         absFast = np.vstack([(data[i]-dataFilt[i])*C.loc[ix,"size"] for i,ix in enumerate(C.index)])*Nrebin
#
#         if filt_cutoff>0:
#             from cv2 import dilate
#             from .numeric import nan_helper
#             var = absSlow
#             if hasattr(self,"gain"):
#                 var = var*self.gain
#             std = var**.5
#             ff = (absFast > filt_cutoff*std).astype("uint8")
#             if ff.any():
#                 dilateKernelSize = int(ironTimeScale/minDt/Nrebin*.2)#*.03)
#                 if dilateKernelSize%2==0:
#                     dilateKernelSize+=1
#                 if dilateKernelSize>=3:
#                     if verbose:
#                         print ("dilating by", dilateKernelSize)
#                     ff = dilate(ff, np.ones(dilateKernelSize, dtype = np.uint8).reshape(1,-1)).astype("bool")
#                 absFast_tmp = absFast.copy()
#                 absFast_tmp[ff] = np.nan
#                 for j in range(C.shape[0]):
#                     y = absFast_tmp[j]
#                     nans, x= nan_helper(y)
#                     if nans.any() and nans.mean()<.5:
#                         y[nans]= np.interp( x(nans), x(~nans), y[~nans] )
# #                         y[nans]= np.interp( x(nans), x(~nans), absSlow[j][~nans] )
#                 dFast = self.sosFilter.run(absFast_tmp)
#                 absSlow = absSlow + dFast
#                 absFast = absFast - dFast
#         var = absSlow
#         if hasattr(self,"gain"):
#             var = var*self.gain
#         std = var**.5
#         zScore = absFast/std
#         if normalize:
#             slower = [absSlow[i]/C["size"].iloc[i]/Nrebin for i in range(len(C))]
#             faster = [absFast[i]/C["size"].iloc[i]/Nrebin for i in range(len(C))]
#         else:
#             slower = absSlow
#             faster = absFast
#         if write:
#             C["slower_%g"%ironTimeScale] = slower
#             C["faster_%g"%ironTimeScale] = faster
#             C["zScore_%g"%ironTimeScale] = list(zScore)
#         else:
#             return np.array(slower), np.array(faster), zScore

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
        
    def events2raster(self, ts, npoints = 1000, onlyRaster=True, z_th = 3):
        k = "%g"%ts
        try:
            self.events[k]
        except:
            self.detect_events(ts, z_th = z_th)
        df = self.events[k]
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
        
    def detect_events(self, ts, z_th=3, Npoints=None, smooth=None, verbose=False, save=True, t=None, zScores=None, min_rel_hw=0.1, debug=False):
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
            smooth = int(ts/dt/10)
            if smooth%2==0: smooth += 1
        if smooth>1:
            if verbose:
                print ("smoothing with kernel", smooth)
            zScores = runningAverage(zScores.T,smooth).T
        if debug:
            self.df["zScore_%g"%ts] = list(zScores)
        events = []
        for i,z in zip(self.df.index,zScores):
            pp = find_peaks(z,
                            width=(ts/dt*min_rel_hw,ts/dt),
                            height=z_th
                              )
            w,h,x0,xe = peak_widths(z, pp[0], rel_height=.5)
            w = w*dt
            x0 = x0*dt + t[0]
            xe = xe*dt + t[0]
            df = pd.DataFrame({"height":z[pp[0]],"halfwidth":w, "t0":x0, "peakpoint":pp[0]*dt+t[0]})
            df["roi"] = i
            events += [df]
        events = pd.concat(events,ignore_index=True)
        events["ts"] = ts
        if save:
            k = "%g"%(ts)
            try:
                self.events
            except:
                self.events = {}
            self.events[k] = events
        else:
            return events
        
    def show_scatter_events(self,ts,timeWindows=None,log_x=False):
        from plotly_express import scatter
        events = self.events[str(ts)].copy()
        if timeWindows is None:
            timeWindows = [[0,np.inf]]
        events["timeWindow"] = [""]*len(events)
        for i,tw in enumerate(timeWindows):
#                 print (tw)
            if tw[1]<=tw[0]:continue
            iis = np.where(events["t0"].between(*tw) & (events["t0"]+events[events.columns[1]]).between(*tw))[0]
            for j in iis:
                if events.loc[j,"timeWindow"]=="":
                    events.loc[j,"timeWindow"] = str(i)
                else:
                    j1 = len(events)
                    events.loc[j1] = events.loc[j]
                    events.loc[j1,"timeWindow"] = str(i)
#             events["timeWindow"][ff] = str(i+1)
        df = events.query("timeWindow!=''")
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
                      log_x=log_x,
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
        C = self.df
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
        protocol = None
        try:
            protocol = pd.read_csv(pathToProtocol, dtype=str)
            protocol.dropna(how='all', inplace=True)
            protocol["t_begin"] = pd.to_timedelta(["00:"+el if type(el)==str else "00:00:00" \
                                                           for el in protocol["begin"]]).total_seconds()
            protocol["t_end"] = pd.to_timedelta(["00:"+el if type(el)==str else (self.time[0]+len(self.time)/self.Freq)*1e9 \
                                                         for el in protocol["end"]]).total_seconds()
            self.protocol = protocol
        except:
            warnings.warn(f"Could not import protocol from {pathToProtocol}.")
        return protocol

    # noinspection PyUnresolvedReferences
    def examine(self, max_rois=10, imagemode=None, debug=False, startShow='',mode="jupyter",name=None,lw=None, test=False):
        if test:
            from .examine_test import examine
        else:
            from .examine import examine
        if imagemode is None:
            imagemode = self.mode
        return examine(self, max_rois=max_rois, imagemode=imagemode, debug=debug, startShow=startShow,mode=mode,name=name)
    
    def examine3(self, max_rois=10, imagemode=None, debug=False, startShow='',mode="jupyter",name=None,lw=None):
        return "examine3 is deprecated. Please, use examine."

    def examine_events(self, df, x, y, debug=False, **otherkwargs):
        from .examine_events import examine_events
        return examine_events(self, df, x, y, debug=debug, **otherkwargs)
    
    def plotTraces(regions, indices, axratios = [1,2], figsize=5, freqShow=2, col="detrended",Offset=5,separate=False):
        if col not in regions.df.columns:
            if col=="detrended":
                regions.detrend_traces()
            else:
                raise ValueError(f"{col} not found in the dataframe.")
            
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
        yl = axs[1].get_ylim()
        dy = yl[1]-yl[0]
        offset = yl[0]/2 - dy/20

        for comp, df in regions.protocol.groupby("compound"):
            for ii in df.index:
                t0,t1 = df.loc[ii].iloc[-2:]
                conc = df.loc[ii,"concentration"]
                x,y = [t0,t1,t1,t0,t0],[-1,-1,-2,-2,-1]
                y = np.array(y)
                y = y*dy/20 + offset
                plt.fill(x,y,color="grey",alpha =.3)
                plt.text(t0,y[:-1].mean(), " "+conc,va="center", ha="left")
                plt.plot(x,y,color="grey",)

            plt.text(df.t_begin.min(),y[:-1].mean(),comp+" ",va="center", ha="right")
            offset -= 1.3*dy/20
            
    def plot_trace(regions, indices, axratios = [1,2.5], figsize=5, freqShow=2, cols=["trace"], Offset=10, separate=False):
        for col in cols:
            if col not in regions.df.columns:
                raise ValueError(f"{col} not found in the dataframe.")
        if hasattr(indices, "__iter__"):
            Offset *= regions.df.loc[indices, cols[0]].apply(median_abs_deviation).mean()
            if "color" in regions.df.columns:
                colors = regions.df.loc[indices, "color"]
            else:
                colors = [MYCOLORS[ix % len(MYCOLORS)] for ix in indices]
            xratios = np.array([.1, axratios[0], .1, axratios[1], .1])
            yratios = xratios[:3]
            xr = xratios / sum(xratios)
            yr = yratios / sum(yratios)
            fig = plt.figure(figsize=(xratios.sum() * figsize, yratios.sum() * figsize))
            axs = [
                fig.add_axes([xr[0], yr[0], xr[1], yr[1]]),
                fig.add_axes([xr[:3].sum(), yr[0], xr[3], yr[1]]),
            ]
            regions.plotEdges(ax=axs[0],lw=.5,separate=separate)
            regions.plotEdges(ax=axs[0],ix=indices, separate=True, fill=True, alpha=.5, image=False)
            regions.plotPeaks(ax=axs[0],ix=indices, labels=True)
            axs[0].set_yticks([])
            axs[0].set_xticks([])
        else:
            Offset = 0
            indices = [indices]
            colors = [None]
            fig, ax = plt.subplots(1,1,figsize=(3*figsize,figsize))
            axs = [None, ax]
        for ic,col in enumerate(cols):
            xs = np.vstack(regions.df.loc[indices,col])
            # if "trace" in col:
            #     t = regions.time
            # else:
            t = regions.time
            if hasattr(regions, "showTime"):
                for k in regions.showTime:
                    if len(t) == xs.shape[1]:
                        break
                    t = regions.showTime[k]
            freq = 1/(t[1]-t[0])
            nr = int(np.round(freq/freqShow))
            if nr>1:
                t = rebin(t, nr)
                xs = rebin(xs,nr,1)
            offset = 0
            for x,ix,c in zip(xs, indices, colors):
                if c is None:
                    axs[1].plot(t,x+offset,lw=1.2 if "slow" in col else .7, color="C%i"%ic, label=col)
                else:
                    axs[1].plot(t,x+offset,lw=.5,color=c)
                    axs[1].text(0,offset+np.median(x[:len(x)//20+1]),str(ix)+" ",color=c,ha="right")
                offset += Offset
        axs[1].set_xlabel("time [s]")
        if c is None:
            axs[1].legend()
        else:
            axs[1].set_yticks([])
            for sp in ["left","right","top"]:
                axs[1].spines[sp].set_visible(False)
        if hasattr(regions,"protocol"):
            yl = axs[1].get_ylim()
            dy = yl[1]-yl[0]
            offset = yl[0]/2 - dy/20
            for comp, df in regions.protocol.groupby("compound"):
                for ii in df.index:
                    t0,t1 = df.loc[ii].iloc[-2:]
                    conc = df.loc[ii,"concentration"]
                    x,y = [t0,t1,t1,t0,t0],[-1,-1,-2,-2,-1]
                    y = np.array(y)
                    y = y*dy/20 + offset
                    plt.fill(x,y,color="grey",alpha =.3)
                    plt.text(t0,y[:-1].mean(), " "+conc,va="center", ha="left")
                    plt.plot(x,y,color="grey",)
                plt.text(df.t_begin.min(),y[:-1].mean(),comp+" ",va="center", ha="right")
                offset -= 1.3*dy/20
        return axs if axs[0] is not None else axs[1]
    
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
        ]).flatten()
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
        unionPixels = sum(C.loc[cl,"pixels"],[])
        unionPixels = list(set(unionPixels))
        C.loc[[attr],"pixels"] = [unionPixels]
        toDrop += other
    C.drop(index=toDrop,inplace=True)
    if verbose:
        print (f"{len(toDrop)} subsumed into existing ROIs.")
    if len(toDrop):
        rreg.update()
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

def getPeak2BoundaryDF(C, verbose=0, distTh=None):
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

def getStatImages(movie_, debleach=False, downsampleFreq=2):
    if movie_.fr>downsampleFreq:
        n_rebin = int(np.round(movie_.fr/downsampleFreq))
        if n_rebin>1:
            m_for_image = movie_.resize(1,1, 1/n_rebin)
        else:
            m_for_image = movie_
    else:
        m_for_image = movie_
    statImages = {}
    # m_for_image = m_for_image.astype("float16")
    if debleach:
        m_for_image.debleach()

    for f in [np.mean,np.std]:
        statImages[f.__name__] = f(m_for_image,axis=0)
    statImages["highperc"] = np.percentile(m_for_image,100*(1-10/len(m_for_image)), axis=0)
    
    # m_for_image = np.diff(m_for_image,axis=0)
    # for f in [np.mean,np.std]:
    #     statImages["diff_"+f.__name__] = f(m_for_image,axis=0)
    # statImages["diff_highperc"] = np.percentile(m_for_image,100*(1-10/len(m_for_image)), axis=0)
    return statImages
