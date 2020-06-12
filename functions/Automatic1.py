import numpy as np
import os
from sys import path as syspath
syspath.append(os.path.expanduser("~/srdjan_functs/"))

from physio_def_1 import *
from Regions1 import Regions


from plotly.subplots import make_subplots
import plotly.graph_objects as go
# from plotly.colors import DEFAULT_PLOTLY_COLORS
# DEFAULT_PLOTLY_COLORS[3] = "aqua"
# DEFAULT_PLOTLY_COLORS[-3] = "darkgoldenrod"
from Regions1 import MYCOLORS
from plotly import colors as plxcolors
import numpy as np

def lassoToMask(points, dims):
    from matplotlib.path import Path
    p = Path(points)
    mask = np.array([[p.contains_point((i,j))  for j in range(dims[1])] for i in range(dims[0])])
    return mask

def saveRois(regions,outDir,filename="",movie=None,col="trace",formats=["vienna","maribor"],add_date=True):
        feedback = []
#     try:
        from copy import deepcopy,copy
        from datetime import date
        import pickle
        import pandas as pd
        from os.path import isdir
        from os import makedirs
        regions.sortFromCenter()
        if movie is not None:
            regions.update(movie)
        filename = filename.replace(" ","_")
        if add_date:
            today = date.today()
            if len(filename):
                filename = "_".join([today.strftime("%Y_%m_%d"),filename])
            else:
                filename = today.strftime("%Y_%m_%d")
        if not isdir(outDir):
            makedirs(outDir)
            feedback += f"Output {outDir} directory created."

        traces = pd.DataFrame(np.vstack(regions.df[col]).T)
        try:
            traces["time"] = regions.showTime[col.split("_")[-1]]
        except:
            traces["time"] = regions.time
        traces = traces[["time"]+list(traces.columns[:-1])]
        for format in formats:
            print (format)
            if format=="vienna":
                saving = ['statImages', 'mode', 'image', 'filterSize', 'df', 'trange', "FrameRange", "analysisFolder", "time", "Freq"]
                movie = regions.movie
                del regions.movie
                allAttrs = list(regions.__dict__.keys())
                subRegions = deepcopy(regions)
                regions.movie = movie
                del movie
                for k in allAttrs:
                    if k not in saving:
                        del subRegions.__dict__[k]
                for k in regions.df.columns:
                    if k not in ["peak", "pixels", "trace"]:
                        del subRegions.df[k]
                        
                roifile = f"{outDir}/{filename}_rois.pkl"
                with open(roifile,"wb") as f:
                    pickle.dump(subRegions,f)
                feedback += [f"ROI info saved in {roifile}."]

            elif format=="maribor":
                tracefile = f"{outDir}/{filename}_trace_for_mb.txt"
                np.savetxt(tracefile, traces.values)
                feedback += [f"Traces saved in {tracefile}."]
                coordFile = f"{outDir}/{filename}_coords_for_mb.txt"
                coords = np.array([np.mean(pxs,axis=0) for pxs in regions.df["pixels"]])
                np.savetxt(coordFile, coords)
                feedback += [f"Coordinates saved in {coordFile}."]
            else:
                return "Output format not recognized. Currently, only 'vienna' and 'maribor' are implemented."
#     except:
#         from sys import exc_info
#         feedback += ["ERROR: "+ exc_info().__repr__()]
        return feedback



def plotTraces(regions,indices=None, showCol=["trace"], offset=0, traceRebin=1):
    from numeric import runningAverage
    if indices is None:
        indices = regions.df.sort_values("size",ascending=False).index[:10]
    
    f = go.Figure()
    times = {}
    for col in showCol:
        if "_" in col:
            try:
                time = regions.showTime[col.split("_")[-1]]
                likeTrace = True
            except:
                time = regions.time
                likeTrace = False
        else:
            likeTrace = False
            time = regions.time
        if likeTrace and traceRebin>1:
            time = rebin(time, traceRebin)
        times[col] = time
    ia = 0
    for i in indices:
        if len(indices)==1:
            cl = "black"
        else:
            cl = MYCOLORS[i%len(MYCOLORS)]
        for col in showCol:
            y = regions.df.loc[i,col]
            t = times[col]
            if y.shape!=t.shape:
                y = rebin(y,traceRebin)
            f.add_trace(go.Scatter(
                x=t,
                y=y+ia*offset,
                line=dict(width=1.3 if "slow" in col else .3,color=cl),
                mode="lines",
                name=str(i),
                hovertext=str(i),
                hoverinfo="text"
            ))
        ia += 1
    f.update_layout({
        "height":370,
        "width":600,
        "margin":dict(l=10, r=10, t=50, b=10),
        #"plot_bgcolor":"white",
        "xaxis": {
            "title": 'time [s]',
        },
    })
    return f

    
def fov_trace(region, twoRows=True):
    
    if twoRows:
        f = make_subplots(rows=2, cols=1,
                          row_heights=[0.7, 0.3],
                          shared_xaxes=True, 
                         )
    else:
        f = go.Figure()
    x = region.fov_trace["time"]
    if twoRows:
        tr0 = go.Scatter(
            x = x,
            y = region.fov_trace["raw"],
            mode = "lines",
            line=dict(width=.9,color="navy"),
            showlegend = False,
            )
        tr1 = go.Scatter(
            x = x,
            y = region.fov_trace["trend"],
            mode = "lines",
            line=dict(width=.6,color="navy"),
            showlegend = False,
            )
        f.add_trace(tr0,row=1, col=1)
        f.add_trace(tr1,row=1, col=1)
        
    y = region.fov_trace["raw"] - region.fov_trace["trend"]
    tr = go.Scatter(
        x = x,
        y = y,
        mode = "lines",
        line=dict(width=.9,color="navy"),
        showlegend = False,
        )
    if twoRows:
        f.add_trace(tr, row=2, col=1)
        f.update_layout({
                "xaxis2": {
                    "rangeslider":{"visible":True, "range":(0,x.max()+np.diff(x)[0])},
                    "title": 'time [s]',
                },
            "height":420,
            "width":480,
            "margin":dict(l=10, r=10, t=20, b=20),
            })
    else:
        f.add_trace(tr)
        f.update_layout({
                "xaxis": {
                    "rangeslider":{"visible":True},
                    "title": 'time [s]',
                    "range": (0,x.max()+np.diff(x)[0])
                },
            "height":250,
            "width":500,
            "margin":dict(l=10, r=10, t=20, b=40),
            })
    return f

def showRoisOnly(regions, indices=None, im=None, showall=True):
    if indices is None:
        indices = regions.df.sort_values("size",ascending=False).index[:10]
    f = go.Figure()
    for i in indices:
        cl = MYCOLORS[i%len(MYCOLORS)]

        bds = regions.df.loc[i,"boundary"]
        bds += [bds[0]]
        y,x = np.array(bds).T
        ypts,xpts = np.array(regions.df.pixels[i]).T
        ln = go.Scatter(x=x,y=y,
                        line=dict(width=1,color=cl),
                        #mode="markers+lines",
                        #mode="lines",
                        #marker={"size":2},
                        showlegend = False,
                        name = str(i),
                        hoverinfo='text',
                        hovertext=["%i"%(i)]*len(bds),
                        fill="toself",
                        opacity = .5,
                     )
        f.add_trace(ln)
    if len(indices):    
        y,x = np.vstack([np.mean(regions.df.pixels[i],axis=0) for i in indices]).T
        pts = go.Scatter(x=x,y=y,
                    mode="markers",
                    showlegend = False,
                    # opacity=0.5,
                    # name=list(map(str,indices)),
                    marker=dict(color=[MYCOLORS[i%len(MYCOLORS)] for i in indices],size=3),
                    hovertext=list(map(str,indices)),
                    hoverinfo="text"
                 )
        f.add_trace(pts)
    else:
        f.add_trace(go.Scatter(x=[0],y=[0],
                    mode="markers",
                    marker=dict(color="blue",size=3, opacity=0),
                    hovertext=None,
                 ))
        
    if im!="none":
        # f.add_heatmap(z=im, hoverinfo='skip',showscale=False,colorscale=plxcolors.sequential.Greys)
        imgpointer = createStaticImage(im,regions,showall=showall, separate=True)

        f.add_layout_image(
            dict(
                source=imgpointer,
                xref="x",
                yref="y",
                x=-.5,
                y=(im.shape[0]-.5),
                sizex=im.shape[1],
                sizey=im.shape[0],
                sizing="stretch",
                opacity=1,
                layer="below")
            )
    
    f.update_layout({
        #"title":regions.mode+" (filtered)",
        "height":400,
        "width":360*im.shape[1]/im.shape[0],
        "margin":dict(l=10, r=10, t=50, b=20),
        "xaxis": {
            "linecolor": 'black',
            "linewidth": 1,
            "mirror": True,
            "tickvals": [],
            "range":[-.5,im.shape[1]-.5]
        },
        "yaxis": {
            "linecolor": 'black',
            "linewidth": 1,
            "mirror": True,
            "tickvals": [],
            "range":[-.5,im.shape[0]-.5]
        },
        'clickmode': 'event+select'
    })
        
    return f

def createStaticImage(im,regions,showall=False,color="grey",separate=False):
    if im is None:
        im = regions.image
    from PIL import Image as PilImage
    import matplotlib.pyplot as plt
    bkg_img_file = "/tmp/%i.png"%np.random.randint(int(1e10))
    figsize=np.array(im.shape)[::-1]/30
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0, 0, 1, 1])
    try:
        im = np.clip(im, np.percentile(im,.2), np.percentile(im,(1-20/im.size)*100))
    except:
        pass
    ax.imshow(np.log(im),cmap="Greys",origin="bottom")
    for sp in ax.spines: ax.spines[sp].set_visible(False)
    if showall:
        try:
            regions.plotEdges(ax=ax,color=color,image=False,lw=figsize[0]*.15,separate=separate)
        except:
            pass
    plt.xticks([])
    plt.yticks([])
    plt.savefig(bkg_img_file,dpi=150)
    plt.close(fig)
    return PilImage.open(bkg_img_file)

def getRigidShifts(movie_, gSig_filt):
    # global dview
    from caiman import stop_server, cluster, save_memmap
    from caiman.motion_correction import MotionCorrect
    from caiman.source_extraction.cnmf import params as params
    # subsampling frequency for Motion Correction
    base_name = "/tmp/tmp"
    # try: dview.terminate()
    # except: pass
    
    #if 'dview' in globals():
    #    stop_server(dview=dview)
    c, dview, n_processes = cluster.setup_cluster(
        backend='local', n_processes=None, single_thread=False)

    fnames = save_memmap([movie_], base_name=base_name, order='C', border_to_0=0, dview=dview)
    # motion correction parameters
    opts = params.CNMFParams(params_dict={
        'fnames'              : fnames,
        "max_deviation_rigid" : int(np.ceil(gSig_filt[0]/2)),
        'border_nan'          : True,
        'pw_rigid'            : False,
        'gSig_filt'           : gSig_filt,
        'nonneg_movie'        : True
    }) 
    mc = MotionCorrect(fnames, dview=dview, **opts.get_group('motion'))
    mc.motion_correct(save_movie=False)
    mc.shifts_rig = np.array(mc.shifts_rig)
    try: dview.terminate()
    except: pass
    return mc.shifts_rig

def createRegions(movie_,gSig_filt, drop_weak=True, drop_small=False, mode="diff_std",plot=False, diag=False, full=True, FrameRange=None):
    regions = Regions(movie_, gSig_filt=gSig_filt, mode=mode, diag=diag, full=full, FrameRange=FrameRange)
    C = regions.df
    if drop_weak:
        ## drop very weak Rois
        from numeric import get_sep_th
        img = regions.statImages[regions.mode].copy()
        imgmin = img[img>0].min()
        imth = np.exp(get_sep_th(np.log(imgmin+img.flatten()),plot=plot))
        if np.mean(img<imth)>.3:
            imth = np.percentile(img,30)
        img[img<imth] = 0
        C["medianValue"] = regions.df.pixels.apply(lambda xi: np.median([img[x,y]   for x,y in xi]))
        C["nNonZero"]    = regions.df.pixels.apply(lambda xi: np.sum(   [img[x,y]>0 for x,y in xi]))
#         toDrop = C.query(f"nNonZero<={regions.filterSize**2/3} and medianValue==0").index
        toDrop = C.query(f"medianValue==0").index
#         x = C["peakValue"]
#         x = np.log10(x+np.percentile(x,5))
#         th = 10**get_sep_th(x)
#         pc  = np.mean(C["peakValue"]<th)
#         if pc > .2:
#             th = np.percentile(C["peakValue"],5)
#         toDrop = C.query(f"peakValue<{th}").index
        C.drop(index=toDrop,inplace=True)
        print (f"Dropped {len(toDrop)} as too weak. Left with {len(C)}.")
    
    
    if drop_small:   ## decide on a threshold to remove as too small rois
        tooSmall = cell_halfwidth_in_px**2 # in pixels
        toDrop = C.query(f"size<={tooSmall}").index
        C.drop(index=toDrop,inplace=True)
        print (f"Dropped {len(toDrop)} as too small. Left with {len(C)}.")
    
    regions.sortFromCenter()
    return regions

def resample(movie, newFreq):
    n_rebin = int(np.ceil(movie.fr/newFreq))
    newmovie = rebin(movie,n_rebin)
    newmovie.fr = movie.fr/n_rebin
    return newmovie

def getFigure(w=300,h=300,c="lightgrey"):
    fig = go.Figure(layout={
        "width":w,
        "height":h,
        # "paper_bgcolor":c,
        "plot_bgcolor":c,
        "margin":dict(zip("lrtb",[0]*4)),
        "xaxis":{"range":[0,1],"tickvals":[]},
        "yaxis":{"range":[0,1],"tickvals":[]},

    })
    return fig


