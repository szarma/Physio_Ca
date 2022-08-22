import logging
import warnings
import os
import pickle
from collections import OrderedDict
from copy import deepcopy
from datetime import date
from itertools import product
from pathlib import Path
from sys import exc_info
from typing import Union

import h5py._hl.group
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
from scipy.stats import median_abs_deviation
import plotly.graph_objects as go
from matplotlib.colors import LogNorm
from matplotlib.axes import Axes
from matplotlib.gridspec import GridSpec
import plotly.express.colors as plc
from typing import Dict, List, Optional, Tuple

from .general_functions import getCircularKernel
from .numeric import rebin, bspline, crawl_dict_via_graph
from .utils import multi_map, getGraph_of_ROIs_to_Merge, getPeak2BoundaryDF, getStatImages
from .EventDistillery import define_legs, plot_events
from .serialization import get_serializer

MYCOLORS = plc.qualitative.Plotly
# MYCOLORS = ["darkred"]

def load_regions(path,
                 baremin=False,
                 mergeDist=0,
                 mergeSizeTh=10,
                 plot=False,
                 verbose=False,
                ):
    with open(path,"rb") as f:
        regions = pickle.load(f)
    if not hasattr(regions, "is_3D"):
        regions.is_3D = False
    regions.pathToRois = path
    pickleDir = os.path.split(path)[0]
    # regions = Regions(regions)
    try:
        protocolFile = os.path.join(pickleDir, [f for f in os.listdir(pickleDir) if "protocol" in f][0])
        regions.import_protocol(protocolFile)
    except:
        logging.warning(f"Could not read the protocol file from {pickleDir}")
    try:
        regions.update()
        if not baremin:
            regions.detrend_traces()
            regions.infer_TwoParFit(plot=plot)
            if hasattr(regions,"gain"):
                del regions.gain
        regions.merge_closest(mergeSizeTh=mergeSizeTh, mergeDist=mergeDist, plot=plot, Niter=15, verbose=verbose)
    except:
        print ("encountered error:", exc_info())
    if "interest" in regions.df.columns:
        del regions.df["interest"]
    return regions


class Regions:
    def __repr__(self):
        try:
            return self.df[[c for c in self.df.columns if self.df[c].dtype.kind in "biufc" or c in ["peak"]]].__repr__()
        except:
            warnings.warn("no dataframe created yet.")
            return self
    def __init__(self, movie_,
                 diag=True,
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
            if movie_.ndim>=3:
                if verbose:
                    print ("Initiating from a movie.")
                self.movie = movie_
                self.is_3D = movie_.ndim>3
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

    def to_hdf(self,
               buffer: h5py._hl.group.Group,
               save_cols=None,
               verbose: bool = False,
               overwrite: bool = False,
               ):
        assert hasattr(buffer,'create_group')
        if save_cols is None:
            save_cols = ["trace"]

        from islets.serialization.simple import df2hdf, write_to_hdf, dict2hdf, fig2hdf
        ######### simple attributes
        simples = {k: getattr(self, k) for k in
                   ['mode', "filterSize", "trange", "FrameRange", "is_3D", "Freq", "TwoParFit", ] if
                   hasattr(self, k)}
        dict2hdf(simples, buffer, 'attributes', overwrite = overwrite, verbose = verbose)
        ### dataframe
        columnsToSave = ["peak", "pixels", "peakValue", "tag"] + save_cols
        df_to_save = self.df[[col for col in self.df.columns if col in columnsToSave]]
        df2hdf(df_to_save, buffer, verbose=verbose, overwrite = overwrite)
        ### filtered image
        write_to_hdf(self.image, buffer, "image", overwrite=overwrite, verbose=verbose)
        ### statImages
        if hasattr(self, "statImages"):
            dict2hdf(self.statImages, buffer, "statImages", in_attrs = False, overwrite = overwrite, verbose = verbose)
            for k in self.statImages:
                buffer['statImages'][k].attrs['CLASS'] = "IMAGE"
                buffer['statImages'][k].attrs['IMAGE_SUBCLASS'] = "IMAGE_GRAYSCALE"
                buffer['statImages'][k].attrs['IMAGE_MINMAXRANGE'] = [0, self.statImages[k].max()]
        ### metadata
        if hasattr(self, 'metadata'):
            skip = ["individual Series"]
            md_simples = { ix: self.metadata[ix] for ix in self.metadata.index if ix not in skip}
            for ix in ["Start time", "End time"]:# "Duration"]:
                md_simples[ix]  = md_simples[ix].strftime('%Y-%m-%d %H:%M:%S')
            md_simples["Duration"] = md_simples["Duration"].total_seconds()

            print (md_simples)
            dict2hdf(md_simples, buffer, "metadata", in_attrs = True, overwrite = overwrite, verbose = verbose)
            # if verbose:
            #     print(f"metadata ({list(md_simples.keys())}) saved in 'metadata.")
            for ix in self.metadata.index:
                if ix in md_simples: continue
                if isinstance(self.metadata[ix], pd.DataFrame):
                    df2hdf(self.metadata[ix], buffer['metadata'], ix, verbose = verbose, overwrite = overwrite)
        ### protocol
        if hasattr(self,'protocol'):
            df2hdf(self.protocol, buffer, 'protocol', verbose = verbose, overwrite = overwrite)
        ### preview image
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_axes([0,0,1,1])
        self.plotEdges(ax=ax)
        fig2hdf(fig, buffer, "preview", overwrite = overwrite)
        if verbose:
            print(f"Preview image saved in '{buffer.name}/preview'")

    @staticmethod
    def load(file_path: Union[str, Path],
             file_format: str = 'infer') -> object:
        """Loads a regions object from a supported file type.

        This function is intended to load data from a supported file type into a ready-to-use regions object.
        Supported file types are:
          - HDF5-files
            - HDF5-files need to have two collections. One collections contains the dataframe at "/df". The other
              collection contains the class dictionary at "/__dict__".
          - JSON-files
            - JSON-files need to be encoded in UTF-8.
            - compressed JSON-files are supported to save disk space.

        Parameters
        ----------
        file_path : Union[str, Path]
            Path of the file, which shall be deserialized.
        file_format : {'hdf5', 'json', 'infer'}, default: 'infer'
            Format, which shall be used to load the file. If 'infer' is passed, the function will automatically try to
            use the right format by guessing based on the file suffix.

        Returns
        -------
        Deserialized Regions object.
        """
        serializer = get_serializer(file_path=file_path, file_format=file_format)
        return serializer.load()

#     def save(self,
#              file_path: Union[str, Path],
#              file_format: str = 'infer',
#              movie: np.ndarray = None,
#              columns: Optional[List[str]] = None,
#              add_date: bool = True,
#              **kwargs) -> None:
#         """Stores a regions object into a supported file format.

#         This function is intended to store data into a supported file type into a ready-to-use regions object.
#         Supported file types are:
#           - HDF5-files
#             - HDF5-files will have two collections. One collection contains the dataframe at "/df". The other
#               collection contains the class dictionary at "/__dict__".
#           - JSON-files
#             - JSON-files will be encoded in UTF-8.
#             - compressed JSON-files are supported to save disk space.

#         Parameters
#         ----------
#         file_path : Union[str, Path]
#             Path of the file, which shall be deserialized.
#         file_format : {'hdf5', 'json', 'infer'}, default: 'infer'
#             Format, which shall be used to load the file. If 'infer' is passed, the function will automatically try to
#             use the right format by guessing based on the file suffix.
#         movie : numpy.ndarray
#             Movie, which shall be applied as an update before saving regions.
#         columns : List[str], optional
#             Columns of the dataframe to save except the necessary ones.
#         add_date : bool
#             Determines whether the current date should get added to the name of the saved file.
#         kwargs
#             Arguments to pass, which are special to the serializers. Please consult the save functions of the to check
#             which parameters are supported.

#         Returns
#         -------
#         None
#         """

#         serializer = get_serializer(file_path=file_path, file_format=file_format)
#         serializer.save(self,
#                         movie=movie,
#                         columns_to_save=columns,
#                         add_date=add_date,
#                         **kwargs)

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


    def constructRois(self, image, img_th=None, morph_transform=True, verbose=False, diag=False, merge=True):
        if img_th is None:
            img_th = 0.02  # median_abs_deviation(image.flat)/3
        ok = (image > img_th).astype(np.uint8)
        if morph_transform:
            from cv2 import dilate, erode
            if hasattr(self,"filterSize") and len(self.filterSize):
                erosion_kernel_size  = (max(self.filterSize)) // 8 * 2 + 1
                dilation_kernel_size = (max(self.filterSize)) // 6 * 2 + 1
            else:
                erosion_kernel_size = 3
                dilation_kernel_size = 3
            dilation_kernel = getCircularKernel(dilation_kernel_size)
            erosion_kernel  = getCircularKernel(erosion_kernel_size)
            if erosion_kernel_size>dilation_kernel_size:
                ok = erode(ok, erosion_kernel)
                if verbose:
                    print ("eroding valid pixels by", erosion_kernel_size)
            ok = dilate(ok, dilation_kernel)
            if verbose:
                print ("dilating valid pixels by", dilation_kernel_size)
        ok = ok.astype(bool)
        self.validPixels = ok
        B_ = crawl_dict_via_graph(image, ok, diag=diag)
        if diag:
                from .utils import split_unconnected_rois
                # TODO: this needs review
                B_ = split_unconnected_rois(B_, self.image)

        self.df = pd.DataFrame(OrderedDict([
            ("peak",  list(B_.keys())),
            ("pixels",list(B_.values()))
        ]))
        self.df["peakValue"] = [image[p] for p in B_]
        self.update()
        if merge and not self.is_3D:
            self.merge_closest(mergeDist=.71, mergeSizeTh=100, verbose=verbose, )

    def mergeBasedOnGraph(self, Gph, verbose=0):
        toDrop = []
        Gph_ = Gph.to_undirected()
        connected_components = nx.connected_components(Gph_)
        # if verbose>1:
        #     print(f"There are {len(connected_components)}.")
        for cl in connected_components:
            cl = list(cl)
            if verbose>1:
                print("Connected component:\t", cl)
            gph = Gph.subgraph(nodes=cl)
            attr = sum(list(map(list,nx.attracting_components(gph))),[])
            if verbose>1:
                print("attractor(s):\t\t", attr)
            # if len(attr)>=2:
            #     warnings.warn(f"two or more attractors ({attr}), not implemented yet, will skip")
            #     continue
            # if len(attr)==2 and len(cl)>2:
            #     continue
            attr = self.df.loc[attr,"peakValue"].sort_values().index[-1]
            if verbose>1:
                print("the chosen attractor:\t", attr)
            other = [j for j in cl if j!=attr]
            if verbose>1:
                print("other rois:\t\t", other)
            unionPixels = sum(self.df.loc[cl,"pixels"],[])
            unionPixels = list(set(unionPixels))
            self.df.loc[attr, "Nneighbors"] -= len(other)
            # nns = set([ngh for ngh in self.df.loc[attr,"neighbors"] if ngh not in other])
            # if verbose>1:
            #     print("surviving neighbors:\t", nns)
            for j in other:
                for k in self.df.loc[j, "neighbors"]:
                    if k in other: continue
                    if verbose>2:
                        print (f"\tremoving {j} from neighbor set of {k}")
                    self.df.loc[k, "neighbors"].remove(j)
            survive = []
            for edge in sum(self.df.loc[cl,"edges"],[]):
                if sum([edge in self.df.loc[roi,"edges"] for roi in cl])>1:
                    if edge in self.edgeIDs:
                        del self.edgeIDs[edge]
                else:
                    survive += [edge]
            self.df.loc[[attr],"edges"] = [survive]
            # peak2idx
            self.df.loc[[attr],"boundary"] = [edges2nodes(survive)]
            self.df.loc[[attr],"pixels"] = [unionPixels]
            newSize = len(unionPixels)
            if "trace" in self.df.columns:
                newtrace = np.sum([self.df.loc[j,"trace"]*self.df.loc[j,"size"] for j in cl],axis=0)/newSize
                # newtrace = pd.Series([newtrace],index=[attr])
                # print (newtrace[:10])
                # if not np.isfinite(newtrace).all():
                #     print (attr, other, cl,"!"*10)
                #     return None
                self.df.loc[[attr],"trace"] = pd.Series([newtrace],index=[attr])
            self.df.loc[attr,"size"] = newSize
            toDrop += other
        self.df.drop(index=toDrop,inplace=True)
        if verbose:
            print (f"{len(toDrop)} subsumed into existing ROIs.")
        # if len(toDrop):
        #     rreg.update()
        #     try:
        #         rreg.calcTraces()
        #     except:
        #         pass
        return len(toDrop)

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
            self.calc_peak2idx()
            # self.update()
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

    def merge_closest(self, mergeSizeTh=10, mergeDist=1, plot=False, Niter=20, verbose=0, axs=None):
        if plot:
            if axs is None:
                plt.figure(figsize=(7*Niter,6))

        ia = 1
        for ja in range(Niter):
            size_th = np.percentile(self.df["size"].values, mergeSizeTh)
            df = getPeak2BoundaryDF(self.df, distTh=mergeDist)
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
                suggestGraph = getGraph_of_ROIs_to_Merge(df.iloc[:, :2], self, plot=plot, ax=ax)
                if plot:
                    ax.set_xlim(xl)
                    ax.set_ylim(yl)
                self.mergeBasedOnGraph(suggestGraph, verbose=verbose)
            else:
                # print ("No more suggestions.")
                break
            ia += 1
        if plot:
            plt.tight_layout()


    def update(self, movie_=None):
        self.df["size"] = self.df["pixels"].apply(len)
        if not self.is_3D:
            self.calcEdges()
            self.df["boundary"] = [edges2nodes(self.df["edges"][j]) for j in self.df.index]
            self.calcNNmap()
        # self.df["interest"] = [np.sum([self.image[px[0],px[1]] for px in pxs]) for pxs in self.df["pixels"]]
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
                  color=None,
                  lw=None,
                  alpha=1,
                  fill=False,
                  scaleFontSize=9,
                  spline=True,
                  bound=True,
                  lengths=None,
                  smoothness=None,
                  **kwargs
                  ):
        if ix is None:
            ix = self.df.index
        if ax is None:
            ax = plt.subplot(111)
        if lw is None:
            lw=.5
        if image:
            im = self.statImages[self.mode].copy()
            im[~np.isfinite(im)] = im[np.isfinite(im)].min()
            # im = (1+im)**1.5/(1+np.abs(im))
            # imth = np.percentile(im[im>0],5)
            # im[im<imth] = imth
            if imkw_args is None:
                imkw_args = {}
                from copy import copy
                mycmap = copy(plt.cm.Greys)
                # mycmap.set_bad("lime")
                imkw_args["cmap"] = mycmap
                # imkw_args["norm"] = LogNorm(vmin=imth)
            axim = ax.imshow(im,**imkw_args)
            # plt.colorbar(axim,ax=ax)
        fs = self.__dict__.get("filterSize", [])
        if fs is None:
            fs = []
        if smoothness is None:
            smoothness = int(np.mean(list(fs)+[3]))+1
        if separate:
            for i in ix:
                if color is not None:
                    c = color
                elif "color" in self.df.columns:
                    c = self.df.loc[i,"color"]
                else:
                    c = MYCOLORS[i%len(MYCOLORS)]
                points = self.df.loc[i,"boundary"]+self.df.loc[i,"boundary"][:3]
                if spline:
                    points = bspline(points, smoothness=smoothness)
                else:
                    points = np.array(points)
                # points = [(p-points.mean(0))*.9+points.mean(0) for p in points]
                y,x = np.array(points).T
                ax.plot(x,y,"-",lw=lw,c=c,alpha=alpha,**kwargs)
                if fill:
                    ax.fill(x,y,c=c,alpha=alpha*.8,**kwargs)
        else:
            if color is None:
                color="darkred"
            tmp = []
            for el in self.df.loc[ix,"boundary"]:
                if spline:
                    el = list(bspline(el, smoothness=smoothness))
                tmp += el
                tmp += el[:3]
                tmp += [(np.nan,)*2]

            y,x = np.array(tmp).T
            ax.plot(x,y,color,lw=lw,alpha=alpha, **kwargs)

        if bound:
            dim = self.image.shape
            ax.set_xlim(-.5,dim[1]-.5)
            ax.set_ylim(dim[0]-.5, -.5,)
        else:
            try:
                xlim, ylim = [np.nan] * 2, [np.nan] * 2
                for ln in ax.lines:
                    x, y = ln.get_data()
                    xlim[0] = np.nanmin(list(x) + [xlim[0]])
                    xlim[1] = np.nanmax(list(x) + [xlim[1]])
                    ylim[0] = np.nanmin(list(y) + [ylim[0]])
                    ylim[1] = np.nanmax(list(y) + [ylim[1]])
                ax.set_ylim(ylim)
                ax.set_xlim(xlim)
            except:
                pass
              
        if image and scaleFontSize>=0 and hasattr(self, "metadata") and "pxSize" in self.metadata:
            try:
                float(self.metadata["pxSize"])
            except:
                return None
            if lengths is None:
                lengths = [10,20,50,100,200,500]
            il = np.searchsorted(lengths,self.metadata.pxSize*self.image.shape[1]/10)
            if il>=len(lengths):
                il = len(lengths)-1
            length=lengths[il]

            x0,x1,y0,y1 = np.array([0,length,0,length*3/50])/self.metadata.pxSize + self.image.shape[0]*.02
            ax.fill_between([x0,x1],[y1]*2,[y0]*2, color="k")
            txt = "\n"*1+str(length)
            if "pxUnit" in self.metadata:
                txt += self.metadata["pxUnit"]
            ax.text((x0+x1)/2, y1+.3*(y1-y0), txt, va="center", ha="center", size=scaleFontSize)

    def plotPeaks(self, ix=None, ax=None, image=False, ms=3, labels=False,color=None, imkw_args={},absMarker=True, marker=".", location="peak",labelEdgeDict=None,**kwargs):
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
            if location=="peak":
                p = self.df.loc[i,"peak"]
            elif location=="center":
                p = np.vstack(self.df.loc[i,"pixels"]).mean(0)
            if color is None:
                try:
                    c = self.df.loc[i,"color"]
                except:
                    c = MYCOLORS[i%len(MYCOLORS)]
            else:
                c = color
            ax.plot(*p[::-1],marker=marker,ms=ms,c=c, **kwargs)
            if labels:
                from matplotlib import patheffects
                tx = ax.text(*p[::-1],s=" "+str(i),color=c, ha="center",va="center",**kwargs)
                if labelEdgeDict is not None:
                    tx.set_path_effects([
                        patheffects.Stroke(**labelEdgeDict),
                        patheffects.Normal()
                    ])

    def interpolate_over_breaks(self, rawTraceCol = None):
        if not hasattr(self, "gaps"):
            if rawTraceCol is None:
                rawTraceCol = 'raw_trace' if "raw_trace" in self.df.columns else "trace"
            zeroTraces = np.all([tr == 0 for tr in self.df[rawTraceCol]], 0)
            gaps = np.where(np.abs(np.diff(zeroTraces)))[0].reshape((-1, 2)) + 1
            gaps[:, 1] += 1
            gaps = list(map(tuple, gaps))
            self.gaps = gaps
        else:
            gaps = self.gaps
        if len(gaps)==0:
            warnings.warn("no gaps detected. If you no exactly where they are, you can enter them explicitly",
                          UserWarning)
            return None
        if "raw_trace" not in self.df.columns:
            self.df["raw_trace"] = [tr.copy() for tr in self.df["trace"]]
        self.df["trace"] = [tr.copy() for tr in self.df[rawTraceCol]]
        dix = int(2 * self.Freq)
        for i in self.df.index:
            tr = self.df.loc[i, rawTraceCol]
            for gap in gaps:
                x0 = tr[gap[0] - dix:gap[0]].mean()
                x1 = tr[gap[1]:gap[1] + dix].mean()
                xr = np.arange(gap[0], gap[1] - 1)
                self.df.loc[i, "trace"][xr] = np.interp(xr, gap, [x0, x1])

    def get_activity(self, timescale, timeframe=None, zth=3, saveAs="activity"):
        if hasattr(timescale,"__iter__"):
            activity = np.zeros(len(self.df))
            for ts in timescale:
                activity += self.get_activity(ts,timeframe,zth,None)
        else:
            if 1./self.Freq > timescale/10:
                warnings.warn(f"The requested timescale ({timescale}) is quite small with respect to frequency {self.Freq}. The results might not be the reasonable.")
            duration = self.time[-1]-self.time[0]
            if duration < timescale*5:
                warnings.warn(f"The requested timescale ({timescale}) is large considering that the duration of the recording {duration}. The results might not be the reasonable.")
            s,f,z = self.fast_filter_traces(timescale, write=False)
            nRebin = len(self.time)/z.shape[1]
            if timeframe is None:
                activity = np.mean(z > zth, 1)
            else:
                if nRebin>=2:
                    time = rebin(self.time,int(nRebin))
                else:
                    time = self.time
                fltr = (time>timeframe[0]) & (time<timeframe[1])
                activity = np.mean(z[:,fltr]>zth,1)
        if saveAs is None:
            return activity
        else:
            self.df[saveAs] = activity

    def color_according_to(self,col,cmap="turbo", vmin = None, vmax = None, nan_color=(.7,)*3):
        # check if column type is numeric
        if self.df[col].dtype.kind not in "biufc":
            raise ValueError(f"{col} elements not numeric.")
        x = self.df[col].values.copy().astype("float")
        if vmin is None:
            vmin = np.nanpercentile(x, 1)
        if vmax is None:
            vmax = np.nanpercentile(x, 95)
        x = (x-vmin)/(vmax-vmin)
        rgbs = np.array([nan_color  if np.isnan(xx) else plt.cm.get_cmap(cmap)(xx)[:3] for xx in x])
        rgbs = (256*rgbs).astype(int)
        rgbs = np.minimum(rgbs,255)
        self.df["color"] = ['#%02x%02x%02x' % tuple(rgb) for rgb in rgbs]

    def propose_merging_based_on_shared_edges(self, rois=None, debug=False,breakTies="peakValue"):
        if self.df[breakTies].dtype.kind not in "biufc":
            raise ValueError(f"{breakTies} elements need to be numeric.")
        if rois is None:
            rois = self.df.index
        C = self.df
        # if len(rois)>len(C.index)/2:
        #     warnings.warn("You want to run on more than half of the rois. This will not end well. I'll restrict this to running on the half rois of the smallest size.")
        #     rois = C.loc[rois,"size"].sort_values().index[:len(C)//2]
        # outDF = []
        G = nx.DiGraph()
        for i in rois:
            if debug:
                print("roi:", i)
            edges0 = set(C.loc[i, "edges"])
            if len(C.loc[i, "neighbors"]) == 0:
                continue
            tmp = pd.DataFrame()
            for j in C.loc[i, "neighbors"]:
                if debug:
                    print("neighbor:", j)
                edges = C.loc[j, "edges"]
                tmp.loc[j,"N_common_edges"] = len(edges0.intersection(edges))
                tmp.loc[j,"size"] = C.loc[j,"size"]
            tmp = tmp.sort_values(["N_common_edges", "size"], ascending = [False, False])
            if debug:
                print (tmp)
            jchoose = tmp.index[0]
            # outDF += [(i,jchoose)]
            G.add_edge(i,jchoose,weight=tmp.loc[jchoose,'N_common_edges']/len(edges0))
        a = nx.adj_matrix(G).toarray()
        b = a.T * a
        nodes = list(G.nodes)
        redundantPairs = [[nodes[p] for p in pair] for pair in zip(*np.where(np.triu(b)))]
        for pair in redundantPairs:
            toDel = C.loc[pair,breakTies].sort_values(ascending = False).index
            G.remove_edge(*toDel)
        # outDF = pd.DataFrame(outDF,columns = ["from","to"])
        return G

    def propose_merging_based_on_close_peaks(self, rois=None, debug=False,breakTies="peakValue"):
        if self.df[breakTies].dtype.kind not in "biufc":
            raise ValueError(f"{breakTies} elements need to be numeric.")
        if rois is None:
            rois = self.df.index
        C = self.df
        if len(rois)>len(C.index)/2:
            warnings.warn("You want to run on more than half of the rois. This will not end well. I'll restrict this to running on the half rois of the smallest size.")
            rois = C.loc[rois,"size"].sort_values().index[:len(C)//2]
        # outDF = []
        G = nx.DiGraph()
        for i in rois:
            if debug:
                print("roi:", i)
            if len(C.loc[i, "neighbors"]) == 0:
                continue
            peak0 = np.array(C.loc[i, "peak"])
            peaks = np.vstack(C.loc[C.loc[i, "neighbors"], "peak"])
            dists = np.linalg.norm(peaks - peak, axis = 1)

            tmp = pd.DataFrame()
            for j in C.loc[i, "neighbors"]:
                if debug:
                    print("neighbor:", j)
                peak = C.loc[j, "edges"]


                peak = np.array(C.loc[i, "peak"])


                tmp.loc[j,"N_common_edges"] = len(edges0.intersection(edges))
                tmp.loc[j,"size"] = C.loc[j,"size"]
            tmp = tmp.sort_values(["N_common_edges", "size"], ascending = [False, False])
            if debug:
                print (tmp)
            jchoose = tmp.index[0]
            # outDF += [(i,jchoose)]
            G.add_edge(i,jchoose)
        a = nx.adj_matrix(G).toarray()
        b = a.T * a
        nodes = list(G.nodes)
        redundantPairs = [[nodes[p] for p in pair] for pair in zip(*np.where(np.triu(b)))]
        for pair in redundantPairs:
            toDel = C.loc[pair,breakTies].sort_values(ascending = False).index
            G.remove_edge(*toDel)
        # outDF = pd.DataFrame(outDF,columns = ["from","to"])
        return G

    def change_frequency(self, fr=2):
        from .movies import movie as cmovie
        traces = np.vstack(self.df.trace)
        fr0 = self.Freq
        trmov = cmovie(traces.T.reshape((-1,len(self.df),1)), fr=fr0)
        trmov = trmov.resize(1,1,fr/fr0)*fr0/fr
        self.df.trace = list(trmov[:,:,0].T)
        self.Freq = fr
        self.time = np.arange(len(trmov))/fr

    def calc_peak2idx(self):
        from bidict import bidict
        peak2idx = bidict([(peak,j) for j,peak in zip(self.df.index,self.df.peak)])
        self.peak2idx = peak2idx

    def calcNNmap(self):
        self.calc_peak2idx()
        neighborsMap = {k:[] for k in self.df["peak"]}
        for edge in self.edgeIDs:
            if len(self.edgeIDs[edge])>1:
                for e1,e2 in product(self.edgeIDs[edge],self.edgeIDs[edge]):
                    if e1 not in neighborsMap: continue
                    if e2 not in neighborsMap: continue
                    if e1==e2: continue
                    if e2 not in neighborsMap[e1]:
                        neighborsMap[e1] += [e2]
        self.df["neighbors"] = [[self.peak2idx[pp] for pp in neighborsMap[p]] for p in self.df["peak"]]
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
        traces = np.ones((len(self.df),(ie-i0)), dtype=np.float16)*np.nan
        for i,ix in enumerate(self.df.index):
            x = self.df.loc[ix,"pixels"]
            sl = (slice(i0,ie), ) + tuple([ el[j] for el in x ] for j in range(len(x[0])))
            traces[i] = movie_[sl].mean(axis=1)
        self.df["trace"] = list(traces)
        time = np.arange(len(movie_))/movie_.fr
        self.time = time[i0:ie]
        self.Freq = movie_.fr

    def reset_filtered(self):
        for col in self.df.columns:
            if "slower" in col or "faster" in col or "zScore" in col:
                del self.df[col]
        self.showTime = {}

    def detrend_traces(self,method="debleach", timescale=None, **kwargs):
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
            self.fast_filter_traces(timescale, Npoints=np.inf, **kwargs)
            self.df["detrended"] = self.df["faster_%g"%timescale]
            self.df["trend"]     = self.df["slower_%g"%timescale]
            del self.df["faster_%g"%timescale], self.df["slower_%g"%timescale], self.df["zScore_%g"%timescale]
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
        fast_vars = fast_vars[slow_est>0]
        slow_est = slow_est[slow_est>0]
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
                       mode='fft',
                       ):
        from .numeric import sosFilter
        if Npoints is None: Npoints = 30
        minDt = np.diff(self.time).mean()
        freq = 1/minDt
        if hasattr(self, "movie"):
            if np.abs(freq/self.movie.fr-1)>1e-2:
                print (f"movie frame rate ({self.movie.fr}) and inferred frame ({freq}) rate are different!")
        N_dt = ironTimeScale/minDt
        Nrebin = max(1,int(np.floor(N_dt/Npoints)))
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
        if mode == 'fft':
            sf = sosFilter(cutFreq, freq, order=order)
        elif mode == "manual":
            sf = ManualFilter(ironTimeScale, freq, )
        else:
            raise ValueError("mode can only be 'fft' or 'manual'.")
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
            self.df = self.df.copy()
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
            allFast = np.vstack(self.df["faster_%g"%ts])
            allSlow = np.vstack(self.df["slower_%g"%ts])
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
            allFast = runningAverage(allFast.T,smooth).T
        if debug:
            self.df["zScore_%g"%ts] = list(zScores)
        events = []
        for i,z,f,s in zip(self.df.index,zScores,allFast,allSlow):
            pp = find_peaks(z,
                            width=(ts/dt*min_rel_hw,ts/dt),
                            height=z_th
                              )
            w,h,x0,xe = peak_widths(z, pp[0], rel_height=.5)
            w = w*dt
            t0 = x0*dt + t[0]
            # x0, xe = peak_widths(z, pp[0], rel_height=.9,)[-2:]
            x0 = (x0).astype(int)
            xe = (xe+2).astype(int)
            auc = [np.sum(f[xx0:xxe]/s[xx0:xxe])*dt for xx0,xxe in zip(x0,xe)]

            df = pd.DataFrame({
                "z_max": z[pp[0]],
                "halfwidth": w,
                "t0": t0,
                "peakpoint": pp[0]*dt+t[0],
                "height": f[pp[0]]/s[pp[0]],
                "auc": auc,
                "time":  [ t[xx0:xxe] for xx0, xxe in zip(x0, xe)],
                "trace": [ f[xx0:xxe] + s[xx0:xxe] for xx0,xxe in zip(x0,xe)]
            })
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
            from .protocol import Protocol
            protocol = Protocol.from_file(pathToProtocol, t0 = self.time[0], tend = self.time[-1])
            # protocol = pd.read_csv(pathToProtocol, dtype=str)
            # protocol.dropna(how='all', inplace=True)
            # protocol["t_begin"] = pd.to_timedelta(["00:"+el if type(el)==str else "00:00:00" \
            #                                                for el in protocol["begin"]]).total_seconds()
            # protocol["t_end"] = pd.to_timedelta(["00:"+el if type(el)==str else (self.time[0]+len(self.time)/self.Freq)*1e9 \
            #                                              for el in protocol["end"]]).total_seconds()
            self.protocol = protocol
        except:
            warnings.warn(f"Could not import protocol from {pathToProtocol}.")
        return protocol

    # noinspection PyUnresolvedReferences
    def examine(self, max_rois=10, imagemode=None, debug=False, startShow='',mode="jupyter",name=None,lw=None, test=False,fill=False):
        if test:
            from .examine_test import examine
        else:
            from .examine import examine
        if imagemode is None:
            imagemode = self.mode
        return examine(self, max_rois=max_rois, imagemode=imagemode, debug=debug, startShow=startShow,mode=mode,name=name, fill=fill)

    def examine_events(self, df, x, y, debug=False, **otherkwargs):
        from .examine_events import examine_events
        return examine_events(self, df, x, y, debug=debug, **otherkwargs)

    def generate_plot(self,
                      rois: List[int],
                      events: pd.DataFrame,
                      legs: Dict[str, Tuple[float, float]],
                      offset: Optional[float] = None,
                      halfwidth_scale: Optional[Tuple[float, float]] = None,
                      show_path: bool = False,
                      plot_sum: bool = False,
                      plot_raster: bool = False) -> None:
        """
        Generates a plot, which is ready for publication (or should be), automatically.

        :param rois: List of indices of ROIs, which shall be plotted.
        :param events: Events, which shall be plotted.
        :param legs: Dictionary containing info about the legs of the experiment.
        :param offset: Offset between the traces. If None, it will be calculated automatically.
        :param halfwidth_scale: Sets minimum/maximum of the scale for the halfwidth duration. None: Automatic scale.
        :param show_path: Trigger to show the path of the pickle, where the data is taken from, in the plot.
        :param plot_sum: Trigger to plot the sum of all given rois at the bottom.
        :param plot_raster: Trigger to plot the raster of the given rois at the bottom.
        :return: None
        """

        # original_font_size = plt.rcParams['font.size']
        # plt.rcParams['font.size'] = 25

        fig: Optional[plt.Figure] = None
        spec: Optional[GridSpec] = None
        ax1: Axes
        ax2: Axes
        ax3: Axes
        ax4: Axes
        ax5: Optional[Axes] = None
        ax6: Optional[Axes] = None

        if plot_sum and plot_raster:
            fig = plt.figure(figsize=(16, 14), constrained_layout=True)
            spec = GridSpec(figure=fig,
                            ncols=2, nrows=4,
                            width_ratios=[1, 3],
                            height_ratios=[1, 0.9, 0.5, 0.8])
            ax5 = fig.add_subplot(spec[2, 1])
            ax6 = fig.add_subplot(spec[3, 1])
        elif plot_sum or plot_raster:
            fig = plt.figure(figsize=(16, 12), constrained_layout=True)
            spec = GridSpec(figure=fig,
                            ncols=2, nrows=3,
                            width_ratios=[1, 3],
                            height_ratios=[1, 1, 1])
            if plot_sum:
                ax5 = fig.add_subplot(spec[2, 1])
            elif plot_raster:
                ax6 = fig.add_subplot(spec[2, 1])
        else:
            fig = plt.figure(figsize=(16, 8), constrained_layout=True)
            spec = GridSpec(figure=fig,
                     ncols=2, nrows=2,
                     width_ratios=[1, 3], height_ratios=[1, 1])

        ax1 = fig.add_subplot(spec[0, 0])
        ax2 = fig.add_subplot(spec[0, 1])
        ax3 = fig.add_subplot(spec[1, 0])
        ax4 = fig.add_subplot(spec[1, 1])
        self.plotTraces(rois,axs=[ax1,ax2])

        # self.plotEdges(ax=ax1, lw=.5)
        # self.plotEdges(ax=ax1, ix=rois, separate=True, fill=True, alpha=.5, image=False)
        # self.plotPeaks(ax=ax1, ix=rois, labels=True)
        # ax1.get_xaxis().set_visible(False)
        # ax1.get_yaxis().set_visible(False)
        #
        num_rebin = 10
        time = rebin(self.time, num_rebin)
        # ax2.set_xlim(time[0], time[-1])
        # if offset is None:
        #     offset = determine_offset(num_rebin)
        # off = 0
        # for roi in rois:
        #     y = rebin(self.df.loc[roi, 'detrended'], num_rebin) + off
        #     c = plc.qualitative.Plotly[roi % len(plc.qualitative.Plotly)]
        #     ax2.plot(time, y, c=c)
        #     ax2.text(x=0, y=y.mean(),
        #              s=str(roi),
        #              c=c,
        #              horizontalalignment='right')
        #     off += offset
        #
        # ymin, ymax = ax2.get_ylim()
        # new_ymin = ymin - (ymax - ymin) / 18 * self.protocol['compound'].nunique()
        # vl = self.protocol['t_begin'].append(self.protocol['t_end']).unique()
        # ax2.vlines(vl, new_ymin, ymax, color='black', zorder=1000)
        # ax2.hlines(ymin, 0, self.time[-1], color='black', zorder=1000)
        # ax2.set_ylim(new_ymin, ymax)
        #
        # off = ymin
        # off_diff = (ymin - new_ymin) / self.protocol['compound'].nunique()
        # for comp, df in self.protocol.groupby('compound'):
        #     for ii in df.index:
        #         t0, t1 = df.loc[ii].iloc[-2:]
        #         conc = df.loc[ii, 'concentration']
        #         x, y = [t0, t1, t1, t0, t0], [off - off_diff, off - off_diff, off, off, off - off_diff]
        #         ax2.fill(x, y, color='grey', alpha=.3)
        #         ax2.text(t0, np.mean(y[:-1]), ' ' + conc, va='center', ha='left')
        #         ax2.plot(x, y, color='grey')
        #     ax2.text(df['t_begin'].min(), np.mean([off - off_diff, off - off_diff, off, off]), comp + ' ',
        #              va='center', ha='right')
        #     off -= off_diff
        ax2.get_yaxis().set_visible(False)
        ax2.set_xlabel('time [s]')

        if 'leg' not in events.columns:
            define_legs(events=events, legs=legs)
        bins = np.geomspace(events.halfwidth.min(), events.halfwidth.max())
        binc = (bins[:-1] + bins[1:]) / 2
        for leg, df in events.dropna().sort_values('t0').groupby('leg'):
            h = np.histogram(df['halfwidth'], bins)[0]
            he = h ** .5
            totime = np.diff(legs[leg])[0] / 60
            c = ax3.plot(h / totime, binc, label=leg)[0].get_color()
            ax3.errorbar(h / totime, binc, xerr=(he / totime), c=c, ls='none', marker='.')
        _, ymax = ax3.get_ylim()
        if halfwidth_scale is not None:
            ax3.set_ylim(halfwidth_scale[0], halfwidth_scale[1])
        else:
            ax3.set_ylim(.4, ymax)
        xmin, xmax = ax3.get_xlim()
        if xmin < 0:
            xmin = 0
        ax3.set_xlim(xmax, xmin)
        ax3.set_yscale('log')
        ax3.set_xlabel('average # of events / min')
        ax3.set_ylabel('halfwidth [s]')
        ax3.set_title(f'Number of ROIs = {self.df.shape[0]}')
        ax3.legend()

        ax4.hexbin(
            x=events['t0'],
            y=events['halfwidth'],
            gridsize=350,
            yscale='log',
            mincnt=1,
            cmap='hot'
        )
        ax4.set_xlim(time[0], time[-1])
        if halfwidth_scale is not None:
            ax4.set_ylim(halfwidth_scale[0], halfwidth_scale[1])
        else:
            ax4.set_ylim(.4, ymax)

        if plot_sum and (ax5 is not None):
            sum_trace = np.sum(self.df.loc[rois, 'detrended'].to_numpy(), axis=0)
            y = rebin(sum_trace, num_rebin)
            ax5.set_xlim(time[0], time[-1])
            ax5.plot(time, y)
            ymin, ymax = ax5.get_ylim()
            ax5.vlines(vl, ymin, ymax, color='black')
            ax5.set_ylim(ymin, ymax)
            ax5.get_yaxis().set_visible(False)
            ax5.set_xlabel('time [s]')

        if plot_raster and (ax6 is not None):
            plot_events(events=events, ax=ax6)
            ax6.set_facecolor('black')
            ax6.get_xaxis().set_visible(False)
            ax6.get_yaxis().set_visible(False)
            ax6.set_xlim(time[0], time[-1])

        if show_path:
            if hasattr(self, 'pathToRois'):
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    if plot_raster:
                        ax6.get_xaxis().set_visible(True)
                        ax6.set_xticks([], [])
                        ax6.set_xlabel(self.pathToRois)
                    else:
                        ax4.set_xticks([], [])
                        ax4.set_xlabel(self.pathToRois)
        else:
            ax4.get_xaxis().set_visible(False)

        fig.show()

    def plotTraces(self, indices, axratios=[1, 2],
                   figsize=5,
                   freqShow=2,
                   col="detrended",
                   Offset=5,
                   separate=False,
                   axs=None,
                   protocol=True,
                   roiColor="black",
                   roiAlpha = .6,
                   traceColor = None,
                   labels = True
                   ):
        if col not in self.df.columns:
            if col=="detrended":
                self.detrend_traces()
            else:
                raise ValueError(f"{col} not found in the dataframe.")
        createAxes = axs is None
        if createAxes:
            xratios = np.array([.1,axratios[0],.1,axratios[1],.1])
            yratios = xratios[:3]
            xr = xratios/sum(xratios)
            yr = yratios/sum(yratios)

            fig = plt.figure(figsize=(xratios.sum()*figsize,yratios.sum()*figsize))
            axs = [
                fig.add_axes([xr[0],yr[0],xr[1],yr[1]]),
                fig.add_axes([xr[:3].sum(),yr[0],xr[3],yr[1]]),
            ]
        if hasattr(axs,"__iter__"):
            roiAxes, traceAxes = axs[:2]
            self.plotEdges(ax=roiAxes, lw=.5, separate=separate, color=roiColor)
            self.plotEdges(ax=roiAxes, ix=indices, separate=True, fill=False, alpha=roiAlpha, image=False, lw=1.5)
            self.plotPeaks(ax=roiAxes, ix=indices, labels=labels)
        else:
            traceAxes = axs

        nr = int(np.round(self.Freq / freqShow))
        if nr==0:
            nr=1
        xs = np.vstack(self.df.loc[indices, col].values)
        if len(self.df[col].iloc[0])==len(self.time):
            if nr>1:
                t = rebin(self.time, nr)
                xs = rebin(xs,nr,1)
            else:
                t = self.time
        else:
            t = self.showTime[col.split("_")[-1]]

        for i in range(len(xs)):
            xs[i] = xs[i]-np.median(xs[i])
        offset = 0
        if traceColor is None:
            if "color" in self.df.columns:
                colors = self.df.loc[indices, "color"]
            else:
                colors = [MYCOLORS[ix%len(MYCOLORS)] for ix in indices]
        else:
            colors = [traceColor]*len(indices)
        for x,ix,c in zip(xs, indices, colors):
            traceAxes.plot(t,x+offset,lw=.5,color=c)
            if labels:
                traceAxes.text(t[0],offset+x[:len(x)//20].mean(),str(ix)+" ",color=c,ha="right", va="center")
            offset += xs.std()*Offset
        if createAxes:
            axs[1].set_yticks([])
            axs[1].set_xlabel("time [s]")
            axs[0].set_yticks([])
            axs[0].set_xticks([])
            for sp in ["left","right","top"]: axs[1].spines[sp].set_visible(False)
        if not protocol or not hasattr(self, "protocol"):
            return None
        else:
            from .utils import add_protocol
            add_protocol(traceAxes, self.protocol)


    def plot_trace(regions, indices, axratios = [1,2.5], figsize=5, freqShow=2, cols=["trace"], Offset=10, separate=False,showProtocol=True, spacing=.1):
        for col in cols:
            if col not in regions.df.columns:
                raise ValueError(f"{col} not found in the dataframe.")
        if hasattr(indices, "__iter__"):
            Offset *= regions.df.loc[indices, cols[0]].apply(median_abs_deviation).mean()
            if "color" in regions.df.columns:
                colors = regions.df.loc[indices, "color"]
            else:
                colors = [MYCOLORS[ix % len(MYCOLORS)] for ix in indices]
            xratios = np.array([ spacing, axratios[0], spacing, axratios[1], spacing])
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
            if len(indices):
                xs = np.vstack(regions.df.loc[indices,col])
            else:
                tlength = regions.df[col].iloc[0].shape[0]
                xs = np.zeros(shape=(1,tlength))
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
            c="w"
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
        if hasattr(regions,"protocol") and showProtocol:
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

