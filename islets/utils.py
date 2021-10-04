import os
from collections import OrderedDict
import tifffile
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import warnings
import ffmpeg
import pandas as pd
from matplotlib import pyplot as plt

def load_tif_seq(file_names, subindices = None):
    # copied from CaImAn
    with tifffile.TiffSequence(file_names) as tffl:
        input_arr = tffl.asarray()
    return input_arr

def memmap_tif(file_name, subindices = None):
    out = tifffile.memmap(file_name)[subindices]
    return out


def load_tif(file_name, subindices = None):
    # copied from CaImAn
    with tifffile.TiffFile(file_name) as tffl:
        multi_page = True if tffl.series[0].shape[0] > 1 else False
        if len(tffl.pages) == 1:
            warnings.warn('Your tif file is saved a single page' +
                            'file. Performance will be affected')
            multi_page = False
        if subindices is not None:
            # if isinstance(subindices, (list, tuple)): # is list or tuple:
            if isinstance(subindices, list):  # is list or tuple:
                if multi_page:
                    if len(tffl.series[0].shape) < 4:
                        input_arr = tffl.asarray(key=subindices[0])[:, subindices[1], subindices[2]]
                    else:  # 3D
                        shape = tffl.series[0].shape
                        ts = np.arange(shape[0])[subindices[0]]
                        input_arr = tffl.asarray(key=np.ravel(ts[:, None] * shape[1] +
                                                              np.arange(shape[1]))
                                                 ).reshape((len(ts),) + shape[1:])[
                            :, subindices[1], subindices[2], subindices[3]]
                else:
                    input_arr = tffl.asarray()[tuple(subindices)]

            else:
                if multi_page:
                    if len(tffl.series[0].shape) < 4:
                        input_arr = tffl.asarray(key=subindices)
                    else:  # 3D
                        shape = tffl.series[0].shape
                        ts = np.arange(shape[0])[subindices]
                        input_arr = tffl.asarray(key=np.ravel(ts[:, None] * shape[1] +
                                                              np.arange(shape[1]))
                                                 ).reshape((len(ts),) + shape[1:])
                else:
                    input_arr = tffl.asarray()
                    input_arr = input_arr[subindices]

        else:
            input_arr = tffl.asarray()

        input_arr = np.squeeze(input_arr)
    return input_arr

def coltrans(x, vmin=None, vmax=None, tilt=1, offset=0.1):
    from .numeric import robust_max
    iterable = hasattr(x,"__iter__")
    if iterable:
        if vmax is None:
            vmax = robust_max(x)
        if vmin is None:
            vmin = -robust_max(-x)
    else:
        iterable = False
        x = np.array([x])
    y = np.log(np.clip(x,vmin,vmax))**(1./tilt)
    if iterable:
        y -= y.min()
        y /= y.max()
#         y = np.minimum(offset+y,y.max())
        y = np.sqrt(offset**2+y**2)
        y /= y.max()
        return y
    else:
        return y[0]


def save_tiff(movie, movieFilename):
    import PIL
    im = PIL.Image.fromarray(movie[0])
    im.save(movieFilename,
            save_all=True,
            append_images=[PIL.Image.fromarray(movie[i]) for i in range(1,len(movie))],
            compression="tiff_deflate"
           )

def get_video_dimensions(file):
    probe = ffmpeg.probe(file)
    height, width = probe["streams"][0]["height"], probe["streams"][0]["width"]
    return width, height

def show_frame(file, frame=0, ax=None, show=True):
    import ffmpeg
    out, _ = (
        ffmpeg
        .input(file)
        .filter('select', 'eq(n,{})'.format(frame))
        .output("pipe:", format="rawvideo", pix_fmt="rgb24")
    #     .output('pipe:', vframes=1, format='image2', vcodec='mjpeg')
        .run(capture_stdout=True)
    )
    w, h = get_video_dimensions(file)
    values  = np.frombuffer(out,np.uint8).reshape(h, w, 3)
    if show:
        setax = ax is None
        if setax:
            fig = plt.figure(facecolor="lime")
            ax = fig.add_axes([.1,.1,.8,.8])
        ax.imshow(values)
        if setax:
            for sp in ax.spines: ax.spines[sp].set_visible(False)
    return values

def renderLatex(formula, fontsize=12, dpi=200, format='png', file=None):
    from io import BytesIO
    """Renders LaTeX formula into image or prints to file.
    """
    fig = plt.figure(figsize=(0.01, 0.01))
    fig.text(0, 0,
#              s = u'${}$'.format(formula),
             s = formula,
             fontsize=fontsize
            )

    output = BytesIO() if file is None else file
#     with warnings.catch_warnings():
#         warnings.filterwarnings('ignore', category=MathTextWarning)
    fig.savefig(output, dpi=dpi, transparent=True, format=format,
                bbox_inches='tight', pad_inches=0.05, frameon=False)

    plt.close(fig)

    if file is None:
        output.seek(0)
        return output

def hex_to_rgb(value):
    """Return (red, green, blue) for the color given as #rrggbb."""
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb_to_hex(red, green, blue):
    """Return color as #rrggbb for the given color values."""
    return '#%02x%02x%02x' % (red, green, blue)

def get_series_dir(pathToExp, series):
    folder = pathToExp+f"_analysis/"
    if not os.path.isdir(folder):
        return []
    
    relevantsubdirs = [sd for sd in os.listdir(folder) if series == sd or series == "_".join(sd.split("_")[:-1])]
    return relevantsubdirs

def get_filterSizes(px, physSize=5.5):
    base = int(np.ceil(physSize/px))
    wider = int(np.ceil(base*1.1))
    if wider==base: wider += 1
    toComb = int(np.ceil(base*1.3))
    if toComb <= wider: toComb += 1
    return [(base,), (wider,), (base,wider), (base,toComb)]

def split_unconnected_rois(B_, image=None):
    import networkx as nx
    from scipy.spatial import distance_matrix
    ######### sort out unconnected pixels
    ks = list(B_.keys())
    for k in ks:
        pxs = list(B_[k])
        if len(pxs)==1:
            continue
        dm = distance_matrix(pxs,pxs)
        gr = nx.Graph(dm<=1)
        subsets = list(nx.connected_components(gr))
        if len(subsets)==1:
            continue
        for subset in subsets:
            tmppxs = [pxs[j] for j in subset]
            if k in tmppxs:
                k1 = k
            else:
                if image is None:
                    k1 = tmppxs[0]
                else:
                    vs = [image[px] for px in tmppxs]
                    k1 = tmppxs[np.argmax(vs)]
            B_[k1] = tmppxs
    return B_

# def 
#     ######### sort out "holes" pixels
#     allTakenPx = sum(B_.values(),[])
#     dims = self.image.shape
#     freePx = [(i,j) for i,j in product(range(dims[0]),range(dims[1])) if (i,j) not in allTakenPx]
#     dm = distance_matrix(freePx,freePx)
#     gg = nx.Graph(dm==1)
#     for cc in nx.connected.connected_components(gg):
#         if len(cc)!=1: continue
#         px = freePx[min(cc)]
#         for di in [-1,1]:
#             pxx = px[0]+di
#             if pxx<0: continue
#             if pxx>=dims[0]: continue
#             nnpeak = climb((px[0]+di,px[1]),self.image, diag=diag)
#         B_[nnpeak] += [px]
#     return B_

def closeup_movie(regions, indices, movie=None, labels=False):
    if movie is None:
        movie = regions.movie
#     from .utils import show_movie
    allpixels = np.vstack(sum(regions.df.loc[indices,"pixels"],[]))
    i0, j0 = allpixels.min(0)
    ie, je = allpixels.max(0)+1
    i0 = max(i0-10,0)
    j0 = max(j0-10,0)
    ie += 10
    je += 10
    def addplot(ax_):
        regions.plotEdges(ax=ax_, ix=indices, separate=True, image=False, scaleFontSize=0)
        regions.plotPeaks(ax=ax_, ix=indices, labels=labels)
    m = movie[:,i0:ie,j0:je]
    a = show_movie(m, additionalPlot = addplot, offset = (j0,i0), figScale=3, autoadjust=False)
    return a

def mode(l):
    from collections import Counter
    return Counter(l).most_common(1)[0][0]

def autocorr(sett, dtrange, nsplits = 1):
    from numpy import zeros, corrcoef, array, mean, std
    if nsplits == 1:
        ret = zeros(len(dtrange))
        for k,i in enumerate(dtrange):
            if i==0:
                ret[k] = 1.
            else:
                ret[k] = corrcoef(sett[:len(sett)-i],sett[i:])[0,1]
        return ret
    else:
        out = []        
        for j in range(nsplits):
            ret = zeros(len(dtrange))
            ss = sett[j*len(sett)//nsplits : (j+1)*len(sett)//nsplits]
            for k,i in enumerate(dtrange):
                if i==0:
                    ret[k] = 1.
                else:
                    ret[k] = corrcoef(ss[:len(ss)-i],ss[i:])[0,1]
            out += [ret]
        out = array(out)
        return ( mean(out,axis=0), std(out,axis=0) )

def order(testlist):
    import numpy as np
    tmp = sorted([[i,el] for i,el in enumerate(testlist)], key=lambda xi: xi[1])
    return np.array([el[0] for el in tmp])
    
def tally(mylist):
    from collections import Counter
    return sorted(Counter(mylist).most_common(),key=lambda duple: duple[0])

def multi_map(some_function, iterable, processes=1, library="multiprocessing"):
    assert processes==int(processes)
    processes = int(processes)
    if processes==1:
        out = map(some_function, iterable)
    elif processes>1:
        if library=="threading":
            from concurrent.futures import ThreadPoolExecutor
            with ThreadPoolExecutor(max_workers=processes) as executor:
                out = executor.map(some_function, iterable)
        elif library=="multiprocessing":
            from multiprocessing import Pool
            try:
                pool = Pool(processes)
                out = pool.map(some_function, iterable)
            finally:
                pool.close()
                pool.join()
        else:
            print(f"Libraries can only be 'multiprocessing' and 'threading'. '{library}' is not known/implemented.")
            out=None
    else:
        print ("invalid number of processes", processes)
        out = None
    return out

def create_preview_image(regions, filepath=None, show=False):
    from copy import copy
    cmap = copy(plt.cm.Greys)
    cmap.set_bad("lime")
    dims = regions.image.shape
    fig = plt.figure(figsize=(5,5*np.divide(*dims)))
    ax = fig.add_axes([0.01,0.01,.98,.98])
    regions.plotEdges(imkw_args={"cmap":cmap},color="darkred", ax=ax, separate=False)
    text = ax.text(.98,.03,len(regions.df),size=16,transform = ax.transAxes, ha="right",color="goldenrod")
    text.set_bbox(dict(facecolor='black', alpha=0.5))
#     text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'),
#                            path_effects.Normal()])
    if filepath is None:
        if hasattr(regions,"pathToPickle"):
            folder, roiName = os.path.split(regions.pathToPickle)
            filepath = ".image_"+roiName.replace("_rois","").replace(".pkl",".png")
            filepath = os.path.join(folder, filepath)
    fig.savefig(filepath,dpi=75)
    if not show:
        plt.close(fig)

def show_movie(m_show,
               figScale = 1,
               out="jshtml",
               fps = 30,
               saveName=None,
               NTimeFrames=100,
               log=True,
               additionalPlot=None,
               dpi=100,
               tmax=None,
               autoadjust=True,
               cmapArgs=None,
               offset=(0,0),
              ):
    from .numeric import rebin

    from matplotlib import animation
    if tmax is not None:
        pass
    m_show = m_show.copy()
    if NTimeFrames is not None:
        n_rebin = len(m_show)//NTimeFrames
        if n_rebin>1:
            m_show = rebin(m_show, n_rebin)
    if autoadjust:
        m_show += 1
        for p in range(1,5):
            baseline = np.percentile(m_show,p)
            m_show = np.maximum(m_show, baseline)
            if np.all(m_show>0): break
    if log:
        m_show = np.log(m_show)
    figsize = np.array(m_show.shape[1:3][::-1])/100*figScale
    currentBackend = matplotlib.get_backend()
    plt.switch_backend('agg')
    fig = plt.figure(figsize=figsize,dpi=dpi)
    ax = fig.add_axes([0.,0.,1,1])
    extent = (offset[0]-.5, offset[0]-.5+m_show.shape[2], offset[1]-.5+m_show.shape[1], offset[1]-.5, )
    if cmapArgs is None:
        im = ax.imshow(m_show[0], cmap="Greys", vmin=m_show.min(), vmax=m_show.max(), extent=extent)
    else:
        im = ax.imshow(m_show[0], extent=extent, **cmapArgs)
    tx = ax.text(1,0," \n",
                 transform = ax.transAxes,
                 ha="right",
                 family="Monospace",
                 va="center",color="darkgoldenrod")
#     tx.set_path_effects([
#         path_effects.Stroke(linewidth=.7, foreground='black'),
#         path_effects.Normal()
#     ])
    ax.set_xticks([])
    ax.set_yticks([])
    if tmax is not None:
        dt = tmax/len(m_show)
    if additionalPlot is not None:
        additionalPlot(ax)
    plt.close(fig)
    def init():
        im.set_data(m_show[0])
        if tmax is not None:
            tx.set_text("0:00")
#             tx.set_path_effects([path_effects.Stroke(linewidth=.7, foreground='black'),
#                            path_effects.Normal()])
        return (im,)
    def animate(i):
        im.set_data(m_show[i])
        if tmax is not None:
            time = i*dt
            mins = int(time/60)
            sec  = int(time-60*mins)
#             ms   = "%i"%(10*(time-60*mins-sec))
#             tx.set_text(f"{mins}:{sec:02d}.{ms} \n")
            tx.set_text(f"{mins}:{sec:02d} \n")
#         tx.set_path_effects([path_effects.Stroke(linewidth=.7, foreground='black'),
#                        path_effects.Normal()])
        return (im,)
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(m_show),
                                   interval=1000/fps,
                                   blit=True)

    plt.switch_backend(currentBackend)
    if out=="html5":
        from IPython.display import HTML
        return HTML(anim.to_html5_video())
    elif out=="jshtml":
        from IPython.display import HTML
        return HTML(anim.to_jshtml(default_mode="reflect"))
    elif out=="save" or saveName is not None:
#         try:
        anim.save(saveName, extra_args=['-vcodec', 'libx264'])
#         except:
#             saveName = input("please enter a valid filename. Otherwise, I'll save it as 'video.mp4'.")
#             try: anim.save(saveName, extra_args=['-vcodec', 'libx264'])
#             except:
#                 saveName = "video.mp4"
#                 anim.save(saveName, extra_args=['-vcodec', 'libx264'])
#         return None
    else:
        raise ValueError("out can only be one of the following: 'html5, jshtml, save'")


def getFigure(w=300,h=300,c="lightgrey"):
    import plotly.graph_objects as go
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


def showRoisOnly(regions, indices=None, im=None, showall=True, lw=None):
    import plotly.graph_objects as go
    from .Regions import MYCOLORS
    if indices is None:
#         indices = regions.df.sort_values("size",ascending=False).index
        indices = regions.df.index
    f = go.Figure()
    if "color" in regions.df.columns:
        colors = regions.df.loc[indices,"color"]
    else:
        colors = [MYCOLORS[i%len(MYCOLORS)] for i in indices]
    
#     for i in indices:
#         try:
#             cl = regions.df.loc[i, "color"]
#         except:
#             cl = MYCOLORS[i%len(MYCOLORS)]
#         bds = regions.df.loc[i,"boundary"]
#         bds += [bds[0]]
#         y,x = np.array(bds).T
#         ypts,xpts = np.array(regions.df.pixels[i]).T
#         ln = go.Scatter(x=x,y=y,
#                         line=dict(width=.7,color=cl),
#                         #mode="markers+lines",
#                         mode="lines",
#                         #marker={"size":2},
#                         hoveron = 'points+fills',
#                         showlegend = False,
#                         name = str(i),
#                         hoverinfo='text',
#                         hovertext=["%i"%(i)]*len(bds),
#                         fill="toself",
#                         #opacity = .5,
#                         fillcolor='rgba(255, 0, 0, 0.05)',
#                      )
#         f.add_trace(ln)
    if len(indices):    
#         y,x = np.vstack([np.mean(regions.df.pixels[i],axis=0) for i in indices]).T
        y,x = np.vstack([regions.df.loc[i,"peak"] for i in indices]).T
        pts = go.Scatter(x=x,y=y,
                    mode="markers",
                    showlegend = False,
                    # opacity=0.5,
                    # name=list(map(str,indices)),
                    marker=dict(color=colors,size=4),
                    hovertext=list(map(str,indices)),
                    hoverinfo="text"
                 )
        f.add_trace(pts)
#     else:
#         f.add_trace(go.Scatter(x=[0],y=[0],
#                     mode="markers",
#                     marker=dict(color="blue",size=3, opacity=0),
#                     hovertext=None,
#                  ))
        
    if im!="none":
        # f.add_heatmap(z=im, hoverinfo='skip',showscale=False,colorscale=plxcolors.sequential.Greys)
        imgpointer = createStaticImage(None,
                                       regions,
                                       showall=showall,
#                                        separate=bool(len(MYCOLORS)-1),
#                                        origin="lower",
                                       lw=lw
                                      )

        f.add_layout_image(
            dict(
                source=imgpointer,
                xref="x",
                yref="y",
                x=-.5,
#                 y=(im.shape[0]-.5),
                y=-.5,
                sizex=im.shape[1],
                sizey=-im.shape[0],
                sizing="stretch",
                opacity=1,
                layer="below")
            )
    h,w = 360,360*im.shape[1]/im.shape[0]
    if w>500:
        h = 500/w*h
        w = 500
    h += 70
    w += 20
    f.update_layout({
        #"title":regions.mode+" (filtered)",
        "height":h,
        "width":w,
        "margin":dict(l=10, r=10, t=50, b=20),
        "xaxis": {
            "zeroline" : False,
            "showgrid" : False,
            "linecolor": 'black',
            "linewidth": 1,
            "mirror": True,
            "tickvals": [],
            "range":[-.5,im.shape[1]-.5]
        },
        "yaxis": {
            "zeroline" : False,
            "showgrid" : False,
            "linecolor": 'black',
            "linewidth": 1,
            "mirror": True,
            "tickvals": [],
#             "range":[-.5,im.shape[0]-.5],
            "range":[im.shape[0]-.5,-.5],
        },
        'clickmode': 'event+select',
        "dragmode":'lasso'
    })
    f.update_yaxes(
        scaleanchor = "x",
        scaleratio = 1,
      )
    try:
        lengths = [10,20,50]
        il = np.searchsorted(lengths,regions.metadata.pxSize*regions.image.shape[1]/10)
        length=lengths[il]
        x0,x1,y0,y1 = np.array([0,length,0,length*3/50])/regions.metadata.pxSize + regions.image.shape[0]*.02
        f.add_shape(
                    type="rect",
                    x0=x0,y0=y0,x1=x1,y1=y1,
                    line=dict(width=0),
                    fillcolor="black",
                    xref='x', yref='y'
                )
        f.add_trace(go.Scatter(
            x=[(x0+x1)/2],
            y=[y1*1.2],
            text=[f"<b>{length}µm</b>"],
            mode="text",
            showlegend=False,
            textposition='bottom center',
            textfont={"color":"black"},
            hoverinfo="skip",
        ))
    except:
        pass
    
    return f
    
    

def createStaticImage(im,regions,showall=True,color="grey",separate=True, returnPath=False, cmap=None,origin="lower",lw=None):
    if im is None:
        im = regions.statImages[regions.mode]
    if lw is None:
        lw = .1
    from PIL import Image as PilImage
    currentBackend = matplotlib.get_backend()
    plt.switch_backend('agg')
    if cmap is None:
        from copy import copy
        cmap = copy(plt.cm.Greys)
        cmap.set_bad("lime")
    bkg_img_file = "/tmp/%i.png"%np.random.randint(int(1e10))
    figsize=np.array(im.shape)[::-1]/30
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0, 0, 1, 1])
    im[im==0] = np.nan
    try:
        im = np.clip(im, np.percentile(im,1), np.percentile(im,(1-20/im.size)*100))
    except:
        pass
    ax.imshow(np.log(im+1),cmap=cmap,origin=origin)
    for sp in ax.spines: ax.spines[sp].set_visible(False)
    if showall:
        try:
            regions.plotEdges(ax=ax,color=color,image=False,lw=figsize[0]*lw,separate=separate,scaleFontSize=0)
        except:
            pass
    plt.xticks([])
    plt.yticks([])
    # plt.ylim(plt.ylim()[::-1])
    plt.savefig(bkg_img_file,dpi=150)
    plt.close(fig)
    plt.switch_backend(currentBackend)
    if returnPath:
        return bkg_img_file
    
    return PilImage.open(bkg_img_file)


def saveRois(regions,outDir,filename="",movie=None,col=["trace"],formats=["vienna"],add_date=True):
        feedback = []
#     try:
        from copy import deepcopy
        from datetime import date
        import pickle
        import pandas as pd
        from os.path import isdir
        from os import makedirs
#         regions.sortFromCenter()
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
            feedback += [f"Output {outDir} directory created."]

        for format in formats:
            if format=="vienna":
                saving = ['statImages', 'mode', 'image', 'filterSize', 'df', 'trange', "FrameRange", "analysisFolder", "time", "Freq","metadata"]
                juggleMovie = hasattr(regions, "movie")
                if juggleMovie:
                    movie = regions.movie
                    del regions.movie
                allAttrs = list(regions.__dict__.keys())
                subRegions = deepcopy(regions)
                if juggleMovie:
                    regions.movie = movie
                    del movie
                for k in allAttrs:
                    if k not in saving:
                        del subRegions.__dict__[k]
                for k in regions.df.columns:
                    if k not in ["peak", "pixels", "peakValue","tag","interest"]+col:
                        del subRegions.df[k]
                        
                roifile = f"{outDir}/{filename}_rois.pkl"
                with open(roifile,"wb") as f:
                    pickle.dump(subRegions,f)
                regions.pathToPickle = roifile
                create_preview_image(regions)
                feedback += [f"ROI info saved in {roifile}."]

            elif format=="maribor":
                warnings.warn(f"you requested to save multiple traces, but maribor format supports only one, so only first ({col[0]}) will be used.")
                col = col[0]
                traces = pd.DataFrame(np.vstack(regions.df[col]).T)
                try:
                    traces["time"] = regions.showTime[col.split("_")[-1]]
                except:
                    traces["time"] = regions.time
                traces = traces[["time"]+list(traces.columns[:-1])]
                tracefile = f"{outDir}/{filename}_{col}.txt"
                np.savetxt(tracefile, traces.values)
                feedback += [f"Traces saved in {tracefile}."]
                coordFile = f"{outDir}/{filename}_coords.txt"
                coords = np.array([np.mean(pxs,axis=0) for pxs in regions.df["pixels"]])
                np.savetxt(coordFile, coords)
                feedback += [f"Coordinates saved in {coordFile}."]
            else:
                return "Output format not recognized. Currently, only 'vienna' and 'maribor' are implemented."
#     except:
#         from sys import exc_info
#         feedback += ["ERROR: "+ exc_info().__repr__()]
        return feedback


def getGraph_of_ROIs_to_Merge(df,rreg, plot=False, ax=None,lw=.5,arrow_width=.5):
    if plot:
        ixs = np.unique(df.values.flatten())
        if ax is None:
            plt.figure(figsize=(10,10))
            ax = plt.subplot(111)
        rreg.plotEdges(image=False, ix=ixs, ax=ax, color="k")

    C = rreg.df
    Gph = nx.DiGraph()

    for _,row in df.iterrows():
        i,j = row[["i","j"]]
        # from_,to_ = C.loc[[i,j],"peakValue"].sort_values().index
        l = list(df.query(f"i=={j}")["j"])
        c = 'red'
        if i in l:
            if C.loc[i,"peakValue"]>C.loc[j,"peakValue"]:
                c = "darkgoldenrod"
            else:
                Gph.add_edge(i, j)
        else:
            Gph.add_edge(i, j)

        if plot:
            x0,y0 = np.mean(C.loc[i,"pixels"],0)
            x1,y1 = C.loc[j,"peak"]
            dx = x1-x0
            dy = y1-y0
            ax.arrow(y0,x0,dy,dx,width = arrow_width,
                     linewidth = lw,
                     color=c,
                     zorder=10,
                     length_includes_head=True)

            #ax.plot(y1,x1,"o",ms=5,mfc="none",mew=.7,c=c)

    if plot:
        plt.gca().set_aspect("equal")
        attractors = sum([list(attr) for attr in nx.attracting_components(Gph)],[])
        # print (attractors)
        rreg.plotPeaks(ax=ax,ix=attractors,color="c",ms=6, zorder=10)

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
        rreg.plotEdges(ax=ax,image=False,ix=cl, spline=False)
        attr = sum(list(map(list,nx.attracting_components(gph))),[])
        rreg.plotEdges(ax=ax,image=False,ix=attr,color="red",spline=False)
        rreg.plotPeaks(ax=ax,image=False,ix=attr,color="red",ms=1)
        for sp in ax.spines: ax.spines[sp].set_visible(False)
        ax.set_aspect("equal")

    for i in range(i+1,axs.size):
        axs.flat[i].remove()
    plt.subplots_adjust(wspace=0,hspace=0)


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
        if distTh is not None:
            bd = C.loc[i, "boundary"]
            x = np.linalg.norm(np.array(bd) - np.repeat([pk], len(bd), axis=0), axis=1)
            if x.min()>distTh:
                continue
        if verbose:
            print(i,"distances:",)
        # dists = OrderedDict()
        for j in C.loc[i,"neighbors"]:
            if j not in C.index:
                warnings.warn(f"{j} is listed as a neighbor of {i}, but it does not exist")
                continue
            bd = C.loc[j,"boundary"]
            x = np.linalg.norm(np.array(bd)-np.repeat([pk],len(bd),axis=0), axis=1)
            d = x.min()
            if distTh is not None and d<=distTh:
                peak2bnd += [(i,j,d,C.loc[i,"size"],C.loc[j,"size"])]
            # dists[j] = x.min()
            if verbose:
                print(f"\tto {j}:", d)
        # if len(dists):
        #     jmin = pd.Series(dists).idxmin()
        #     peak2bnd += [(i,jmin,dists[jmin],C.loc[i,"size"],C.loc[jmin,"size"])]

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