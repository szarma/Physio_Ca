import numpy as np
import os
from sys import path as syspath
syspath.append(os.path.expanduser("~/srdjan_functs/"))

from physio_def_1 import *
from Regions import Regions


from plotly.subplots import make_subplots
import plotly.graph_objects as go
# from plotly.colors import DEFAULT_PLOTLY_COLORS
# DEFAULT_PLOTLY_COLORS[3] = "aqua"
# DEFAULT_PLOTLY_COLORS[-3] = "darkgoldenrod"
from Regions import MYCOLORS
from plotly import colors as plxcolors
import numpy as np

def saveRois(regions,outDir,filename="",movie=None,col="trace",formats=["vienna","maribor"]):
        feedback = []
#     try:
        from copy import deepcopy
        from datetime import date
        import pickle
        import pandas as pd
        from os.path import isdir
        from os import makedirs
        regions.sortFromCenter()
        if movie is not None:
            regions.update(movie)
        filename = filename.replace(" ","_")
        today = date.today()
        if len(filename):
            filename = "_".join([today.strftime("%Y_%m_%d"),filename])
        else:
            filename = today.strftime("%Y_%m_%d")
        if not isdir(outDir):
            makedirs(outDir)
            feedback += f"Output {outdir} directory created."

        traces = pd.DataFrame(np.vstack(regions.df[col]).T)
        try:
            traces["time"] = regions.showTime[col]
        except:
            traces["time"] = regions.time
        traces = traces[["time"]+list(traces.columns[:-1])]
        for format in formats:
            print (format)
            if format=="vienna":
                saving = ['statImages', 'mode', 'image', 'filterSize', 'df']
                allAttrs = list(regions.__dict__.keys())
                subRegions = deepcopy(regions)
                for k in allAttrs:
                    if k not in saving:
                        del subRegions.__dict__[k]
                for k in regions.df.columns:
                    if k not in ["peak", "pixels"]:
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


def plotStatsFig(regions,which=["mean","diff_std"], showRois=False):
    fig1 = make_subplots(rows=1, cols=2, horizontal_spacing=0.02,
#                          shared_yaxes=True,
                    subplot_titles=which)
    fig1.update_layout(
        height = 360,
        margin=dict(l=20, r=10, t=50, b=20),
        width = 600,
        title = {"font":{"size": 7}},
    )
    for k in ["xaxis","yaxis","xaxis2","yaxis2",]:
        fig1.update_layout({
            k : {
            "linecolor": 'black',
            "linewidth": 1,
            "mirror": True,
            "tickvals": [],
        }
        })
    for j,im in enumerate([regions.statImages[which[0]], regions.statImages[which[1]], ]):
        imgpointer = createStaticImage(im,regions,showall=showRois,color="darkgoldenrod")
        fig1.add_layout_image(
            dict(
                source=imgpointer,
                xref="x",
                yref="y",
                x=-.5,
                y=regions.image.shape[0]-.5,
                sizex=im.shape[1],
                sizey=im.shape[0],
                sizing="stretch",
                opacity=1,
                layer="below")
            ,col=j+1,row=1)
    
    for c in [1,2]:
        fig1.update_xaxes({"range":[-.5,im.shape[1]-.5]},col=c,row=1,showticklabels=False)
        fig1.update_yaxes({"range":[-.5,im.shape[0]-.5]},col=c,row=1,showticklabels=False)
    return fig1

def plotRawTraces(regions,indices=None, nAvg=1):
    from physio_def_1 import rebin
    if indices is None:
        indices = regions.df.sort_values("size",ascending=False).index[:10]
    f = go.Figure()
    t = regions.time
    if nAvg > 1:
        t = rebin(t,nAvg)
        print ("rebinned with ",nAvg)
    for i in indices:
        cl = MYCOLORS[i%len(MYCOLORS)]
        y = regions.df.loc[i,"trace"]
        if nAvg>1:
            y = rebin(y,nAvg)
        f.add_trace(go.Scatter(
            x=t,
            y=y,
            line=dict(width=1,color=cl),
            mode="lines",
            name=str(i),
            hovertext=str(i),
            hoverinfo="text"
        ))
    f.update_layout({
        "height":360,
        "width":500,
        "margin":dict(l=10, r=10, t=50, b=10),
        #"plot_bgcolor":"white",
        "xaxis": {
            "title": 'time [s]',
        },
    })
    return f


def plotTraces(regions,indices=None, showCol=["trace"], time=None, offset=0):
    if indices is None:
        indices = regions.df.sort_values("size",ascending=False).index[:10]
    
    f = go.Figure()
    if time is None:
        t = regions.time
    else:
        t = time
    ia = 0
    for i in indices:
        if len(indices)==1:
            cl = "black"
        else:
            cl = MYCOLORS[i%len(MYCOLORS)]
        for col in showCol:
            y = regions.df.loc[i,col]
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

    
def fov_trace(region):
    
    f = make_subplots(rows=2, cols=1,
                      row_heights=[0.7, 0.3],
                      shared_xaxes=True, 
                     )
    
    f.add_trace(
        go.Scatter(
            x = region.fov_trace["time"],
            y = region.fov_trace["raw"],
            mode = "lines",
            line=dict(width=.9,color="navy"),
            showlegend = False,

        ),
            row=1, col=1
    )
    f.add_trace(
        go.Scatter(
            x = region.fov_trace["time"],
            y = region.fov_trace["trend"],
            mode = "lines",
            line=dict(width=.6,color="navy"),
            showlegend = False,

        ),
            row=1, col=1
    )
    y = region.fov_trace["raw"] - region.fov_trace["trend"]
    f.add_trace(
        go.Scatter(
            x = region.fov_trace["time"],
            y = y,
            mode = "lines",
            line=dict(width=.9,color="navy"),
            showlegend = False,
        ), row=2, col=1
    )
    f.update_layout({
#             "height":360,
#             "width":500,
    #         "margin":dict(l=10, r=10, t=50, b=10),
            "xaxis2": {
                "rangeslider":{"visible":True},
                "title": 'time [s]',
            },
        "height":420,
        "width":480,
        "margin":dict(l=10, r=10, t=50, b=20),
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

def createRegions(movie_,gSig_filt, drop_weak=True, drop_small=False, mode="diff_std",plot=False, diag=False, full=True):
    regions = Regions(movie_, gSig_filt=gSig_filt, mode=mode, diag=diag, full=full)
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



# def addRoisToFig(regions,indices, fig1):
#     fig1.data = fig1.data[:2]
#     for i in indices:
#         cl = DEFAULT_PLOTLY_COLORS[i%10]
#         bds = regions.df.loc[i,"boundary"]
#         bds += [bds[0]]
#         y,x = np.array(bds).T
#         ypts,xpts = np.array(regions.df.pixels[i]).T
#         for c in [1,2]:
#             ln = go.Scatter(x=x,y=y,
#                             line=dict(width=1,color=cl),
#                             mode="lines",
#                             showlegend = False,#c==1,
#                             name = str(i),
#                             hoverinfo='text',
#                             hovertext=["%i"%(i)]*len(bds),
#                          )
#             fig1.add_trace(ln,col=c,row=1)
            
#             pts = go.Scatter(x=xpts,y=ypts,
#                     mode="markers",
#                     showlegend = False,
#                     opacity=0.0,
#                     name=str(i),
#                     marker=dict(color=cl),
#                     hovertext=["%i"%(i)]*len(xpts),
#                     hoverinfo="text"
#                  )
#             fig1.add_trace(pts,col=c,row=1)
#     return fig1
# def plotStatsFig(regions,ix=[],which=["mean","diff_std"]):
#     fig1 = make_subplots(rows=1, cols=2, horizontal_spacing=0.02,
#                          shared_yaxes=True,
#                     subplot_titles=which)
#     for j,im in enumerate([regions.statImages[which[0]], regions.statImages[which[1]], ]):
#         fig1.add_heatmap(z=im,col=j+1,row=1, hoverinfo='skip',showscale=False,colorscale=plxcolors.sequential.Greys)
#     fig1.update_layout(
#         height = 300,
#         margin=dict(l=20, r=10, t=50, b=20),
#         width = 500,
#         title = {"font":{"size": 7}},
#     )
#     for k in ["xaxis","yaxis","xaxis2","yaxis2",]:
#         fig1.update_layout({
#             k : {
#             "linecolor": 'black',
#             "linewidth": 1,
#             "mirror": True,
#             "tickvals": [],
#         }
#         })
#     if len(ix): fig1 = addRoisToFig(regions,ix,fig1)
#     for c in [1,2]:
#         fig1.update_xaxes({"range":[-.5,im.shape[1]-.5]},col=c,row=1,showticklabels=False)
#         fig1.update_yaxes({"range":[-.5,im.shape[0]-.5]},col=c,row=1,showticklabels=False)
#     return fig1

# def parseRoiChoices(regions, roi_mode, roi_number, maxNRois=200):
#     try:
#         roi_number = int(roi_number)
#     except:
#         roi_number = 3
#     if roi_number<1:
#         roi_number=1
#     if roi_number>maxNRois:
#         roi_number=maxNRois
#     if roi_mode=="interest":
#         ix = regions.df.sort_values("interest",ascending=False).index[:roi_number]
#     if roi_mode=="central":
#         ix = regions.df.index[:roi_number]
#     if roi_mode=="outer":
#         ix = regions.df.index[::-1][:roi_number]
#     if roi_mode=="largest":
#         ix = regions.df.sort_values("size",ascending=False).index[:roi_number]
#     if roi_mode=="rnd":
#         ix = np.random.choice(regions.df.index, roi_number)
#     return ix



