import logging
import os
from collections import OrderedDict
import tifffile
import matplotlib
import networkx as nx
import numpy as np
import ffmpeg
import pandas as pd
from tqdm import tqdm
from matplotlib import pyplot as plt

def load_tif_seq(file_names):
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
            logging.warning('Your tif file is saved a single page' +
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
                        # noinspection PyUnresolvedReferences
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
                        # noinspection PyUnresolvedReferences
                        input_arr = tffl.asarray(key=np.ravel(ts[:, None] * shape[1] +
                                                              np.arange(shape[1]))
                                                 ).reshape((len(ts),) + shape[1:])
                else:
                    input_arr = tffl.asarray()
                    input_arr = input_arr[subindices]

        else:
            input_arr = tffl.asarray()

        # noinspection PyUnresolvedReferences
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
    # noinspection PyUnresolvedReferences
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
    '''This is apparently a horrible way to write a tiff, Use Tifffile instead.
    Need to remove.'''
    #TODO: remove
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
        .global_args('-loglevel', 'error')
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

def rgb_to_hex(x):
    """Return color as #rrggbb for the given color values."""
    red, green, blue = (np.array(x[:3])*255).astype(int)
    return '#%02x%02x%02x' % (red, green, blue)

def get_series_dir(pathToExp, series):
    folder = pathToExp+f"_analysis/"
    if not os.path.isdir(folder):
        return []
    #relevantsubdirs = [sd for sd in os.listdir(folder) if series == sd or series == "_".join(sd.split("_")[:-1])]
    relevantsubdirs = [sd for sd in os.listdir(folder) if series in sd and os.path.isdir(os.path.join(folder,sd))]
    return relevantsubdirs

def get_filterSizes(px, physSize=5.5):
    base = int(np.ceil(physSize/px))
    wider = int(np.ceil(base*1.1))
    if wider==base: wider += 1
    return [(wider,), (base, wider)]
    # toComb = int(np.ceil(base*1.3))
    # if toComb <= wider: toComb += 1
    # return [(base,), (wider,), (base,wider), (base,toComb)]

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
                    # noinspection PyUnresolvedReferences
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

def closeup_movie(regions, indices, movie=None, labels=False,**kwargs):
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
        regions.plotEdges(ax=ax_, ix=indices, separate=True, image=False, scaleFontSize=0, bound=False)
        regions.plotPeaks(ax=ax_, ix=indices, labels=labels)
    m = movie[:,i0:ie,j0:je]
    a = show_movie(m, additionalPlot = addplot, offset = (j0,i0), figScale=3, autoadjust=False, **kwargs)
    return a

def mode(l):
    from collections import Counter
    return Counter(l).most_common(1)[0][0]

def autocorr(sett, dtrange, nsplits = 1):
    # noinspection PyUnresolvedReferences,PyUnresolvedReferences
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
#     cmap = copy(plt.cm.Greys)
#     cmap.set_bad("lime")
    dims = regions.image.shape
    fig = plt.figure(figsize=(5,5*np.divide(*dims)))
    ax = fig.add_axes([0.01,0.01,.98,.98])
    regions.plotEdges(color="darkred", ax=ax, separate=False)
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
               fps = 60,
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
    if m_show.ndim==4:
        logging.warning("Assuming this a 3D recording.")
        t,z,h,w = m_show.shape
        m_show = np.concatenate(np.transpose(m_show, (1, 0, 2, 3)), axis = -1, )
        for i in range(2):
            m_show = np.insert(m_show, i + np.arange(w, m_show.shape[-1], w + i),
                          np.ones((h, 1)) * m_show.max() / 2, axis = 2)
    else:
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


def boil_down(regions, timescales=None, verbose=False, minNevents = 10, time_ranges=None, skip_sequential_filtering=False):
    if not skip_sequential_filtering:
        from .EventDistillery import sequential_filtering
        if timescales is None:
            timescales = 2. ** np.arange(-1, 20)
        timescales = timescales[timescales < regions.time[-1] / 3]
        timescales = timescales[timescales > 5 / regions.Freq]
        sequential_filtering(regions, timescales=timescales, verbose=verbose)
    timescales = regions.timescales
    if time_ranges is None:
        time_ranges = {"all": (regions.time[0], regions.time[-1])}
    output = {}
    for tr in time_ranges:
        tbegin, tend = time_ranges[tr]
        median_ccs = {}
        for ts in timescales:
            k = "%g" % ts
            time = regions.showTime.get(k, regions.time)
            time_mask = (time>tbegin) & (time<tend)
            zs = np.vstack(regions.df["zScore_%g" % ts].values)[:,time_mask]
            meanZ = zs.sum(0)
            median_ccs[ts] = np.median([np.corrcoef(meanZ, z)[0, 1] for z in zs])
        t_sync=pd.Series(median_ccs).sort_values().idxmax()
        median_cc = median_ccs[t_sync]
        events = regions.events["%g" % t_sync]
        events["tend"] = events["t0"] + events["halfwidth"]
        periods = []
        halfwidths = []
        for roi, evroi in events.query(f"t0>{tbegin} and tend<{tend}").groupby("roi"):
            if len(evroi) < minNevents:
                continue
            periods += [evroi.t0.sort_values().diff().median()]
            halfwidths += [evroi.halfwidth.median()]
        output[tr] = dict(t_sync=t_sync, halfwidth=np.median(halfwidths), period=np.median(periods), cc2sum=median_cc)
    return pd.DataFrame(output).T

def get_ccs(regions, timescales=None, verbose=False, time_ranges=None, skip_sequential_filtering=False, mode="toMean", col="zScore"):
    validModes = ["toMean", "cross"]
    if mode not in validModes:
        raise ValueError("mode can only accept one of the following: "+str(validModes))
    if not skip_sequential_filtering:
        from .EventDistillery import sequential_filtering
        if timescales is None:
            timescales = 2. ** np.arange(-1, 20)
        timescales = timescales[timescales < regions.time[-1] / 10]
        timescales = timescales[timescales > 5 / regions.Freq]
        sequential_filtering(regions, timescales=timescales, verbose=verbose)
    timescales = regions.timescales
    if time_ranges is None:
        time_ranges = {"all": (regions.time[0], regions.time[-1])}
    output = {}
    for tr in time_ranges:
        if mode=="toMean":
            ccs = np.ones((len(timescales), len(regions.df)))
        if mode=="cross":
            ccs = np.ones((len(timescales), len(regions.df), len(regions.df)))
        tbegin, tend = time_ranges[tr]
        for it,ts in enumerate(timescales):
            k = "%g" % ts
            time = regions.showTime.get(k, regions.time)
            time_mask = (time>tbegin) & (time<tend)
            zs = np.vstack(regions.df[col+"_%g" % ts].values)[:,time_mask]
            if mode=="toMean":
                meanZ = zs.sum(0)
                ccs[it] = [np.corrcoef(meanZ, z)[0, 1] for z in zs]
            if mode=="cross":
                ccs[it] = np.corrcoef(zs)
                
        output[tr] = ccs
    return output

def showRoisOnly(regions, indices=None, im=None, showall=True, lw=None,fill=False):
    import plotly.graph_objects as go
    from .Regions import MYCOLORS
    if indices is None:
#         indices = regions.df.sort_values("size",ascending=False).index
        indices = regions.df.index
    if im is None:
        im = regions.statImages[regions.mode]
    f = go.Figure()
    if "color" in regions.df.columns:
        colors = regions.df.loc[indices,"color"]
    else:
        colors = [MYCOLORS[i%len(MYCOLORS)] for i in indices]
    if len(indices):
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
        
    if im!="none":
        im = np.maximum(im, np.percentile(im,5))
        imgpointer = createStaticImage(regions,
                                       im = im,
                                       showall=showall,
                                       lw=lw,
                                       fill=fill
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
    
    return f
    
    

def createStaticImage(regions,im=None,showall=True,returnPath=False,origin="lower",lw=None,fill=False,**imkwargs):
    if im is None:
        im = regions.statImages[regions.mode]
    if lw is None:
        lw = .1
    from PIL import Image as PilImage
    currentBackend = matplotlib.get_backend()
    plt.switch_backend('agg')
    # if cmap is None:
        # from copy import copy
        # cmap = copy(plt.cm.Greys)
        # cmap.set_bad("lime")
    bkg_img_file = "/tmp/%i.png"%np.random.randint(int(1e10))
    figaspect = im.shape[0]/im.shape[1]
    figsize=(8, 8*figaspect)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0, 0, 1, 1])
    im[im==0] = np.nan
    try:
        im = np.clip(im, np.percentile(im,1), np.percentile(im,(1-20/im.size)*100))
    except:
        pass
    if "cmap" not in imkwargs:
        imkwargs["cmap"] = plt.cm.Greys
    ax.imshow( np.log(im+1), origin=origin, **imkwargs )
    for sp in ax.spines: ax.spines[sp].set_visible(False)
    if showall:
        try:
            regions.plotEdges(ax=ax, image=False, lw=figsize[0]*lw, separate=True, scaleFontSize=25, fill=fill )
        except:
            pass
    plt.xticks([])
    plt.yticks([])
    # plt.ylim(plt.ylim()[::-1])
    plt.savefig(bkg_img_file,dpi=100)
    plt.close(fig)
    plt.switch_backend(currentBackend)
    if returnPath:
        return bkg_img_file
    
    return PilImage.open(bkg_img_file)

def motion_correct(movie, m_rshifted, freqMC=0.5, max_dev=(5, 5), plot_name="shifts.png", pinpoint_template = 0, Niter = 4, mode = "full", verbose = True, template = None, shifts = None):
    freq = movie.fr
    if shifts is None:
        n_rebin = int(np.round(freq / freqMC))
        if n_rebin < 1:
            n_rebin = 1
        if verbose: (f"The original frequency is {freq}.")
        if n_rebin > 1:
            if verbose: print(f"The movie will be rebinned in time by {n_rebin} for shifts extraction.")
        reb_movie = movie.resize(1, 1, 1. / n_rebin).astype("float32")
        reb_movie = reb_movie.gaussian_blur_2D(max_dev[0] * 2 + 1, max_dev[1] * 2 + 1, -1, -1)

        if verbose: print(f'Extracting shifts started. The output will be saved in {plot_name} file.')
        #     fig, axs = plt.subplots(2, 1, figsize=(8, 8), sharex=True, sharey=True, dpi=150)
        fig, axs = plt.subplots(1, 1, figsize = (8, 3.5), sharex = True, sharey = True, dpi = 150)
        axs = [axs, axs]
        #     axs[0].set_title("shifts (in pixels)")
        axs[0].set_title("vertical: sold, horizontal: dashed")
        axs[0].set_ylabel("shifts (in pixels)")
        #     axs[0].set_ylabel("vertical")
        #     axs[1].set_ylabel("horizontal")
        axs[1].set_xlabel("time [s]")
        max_shift_vert, max_shift_hor = max_dev

        shifts = np.zeros(shape = (len(reb_movie), 2), dtype = "float32")
        if mode in ["pairwise","full"]:
            dshifts = [[0, 0]]
            for iframe in range(reb_movie.shape[0]-1):
                mtmp = reb_movie[iframe:iframe + 2]
                dshifts += [mtmp.extract_shifts(max_shift_vert//2+1, max_shift_hor//2+1, template = mtmp[0])[0][1]]
            shifts += np.cumsum(np.array(dshifts), axis = 0)

            shifts[:, 0] -= np.median(shifts[:, 0])
            shifts[:, 1] -= np.median(shifts[:, 1])

            reb_movie = reb_movie.apply_shifts(shifts)

            c = axs[0].plot([])[0].get_color()
            axs[0].plot(np.arange(len(shifts)) / reb_movie.fr, shifts[:, 0], label = "pairwise inference", c = c)
            axs[1].plot(np.arange(len(shifts)) / reb_movie.fr, shifts[:, 1], "--", c = c)
            plt.tight_layout();fig.savefig(plot_name)

        if mode in ["template-based", "full"]:
            if template is None:
                ichoose = int(pinpoint_template * len(reb_movie))
                template = reb_movie[ichoose]


            for i_extract in range(Niter):
                dshifts = reb_movie.extract_shifts(
                    max_shift_w = max_shift_hor,
                    max_shift_h = max_shift_vert,
                    template = template
                )[0]
                shifts += np.array(dshifts)
                reb_movie.apply_shifts(dshifts)
                maxshifts = np.abs(dshifts).max(axis = 0)
                if verbose:
                    print("maximal shifts are ", maxshifts)
                # print ("all dshifts ", dshifts)
                c = axs[0].plot([])[0].get_color()
                axs[0].plot(np.arange(len(shifts)) / reb_movie.fr, shifts[:, 0], label = i_extract, c = c)
                axs[1].plot(np.arange(len(shifts)) / reb_movie.fr, shifts[:, 1], "--", c = c)
                plt.tight_layout();fig.savefig(plot_name)
                i_extract += 1
                if (maxshifts[0] < max_dev[0]) and (maxshifts[1] < max_dev[1]):
                    break

        axs[0].legend()
        maxyl = np.abs([ax.get_ylim() for ax in axs]).max()
        for ax in axs:
            ax.set_ylim(-maxyl, maxyl)
        plt.tight_layout();fig.savefig(plot_name)
        if verbose: print(f'Extracting shifts finished.')

        #### if rebinned, then need to expand shift to fit the whole movie
        if n_rebin > 1:
            from . import cmovie
            shifts = cmovie(shifts.reshape((1,) + shifts.shape)).resize(n_rebin, 1, 1)[0]
    # pad the remaining frames if exist:
    npad = len(movie) - len(shifts)
    if npad > 0:
        shifts = np.vstack([shifts, [shifts[-1]] * npad])

    # prepare for saving
    chunkSize = 3e9 # gigabytes
    singleFrameSize32 = movie[0].astype("float32").nbytes
    di = int(chunkSize/singleFrameSize32)+1
    range_ = range(0, len(shifts), di)
    if verbose:
        print(f'Applying shifts...', flush=True)
        range_ = tqdm(range_)

    for i in range_:
        sl = slice(i, min(len(shifts), i + di))
        tmpm = movie[sl].astype("float32").apply_shifts(shifts[sl])
        tmpm[tmpm < 0] = 0
        # tmpm = np.round(tmpm)
        # tmpm = dst.poisson.rvs(mu=tmpm)
        # tmpm[tmpm > maxv] = maxv
        m_rshifted[sl] = tmpm.astype(m_rshifted.dtype)
    if verbose: print("Done.")
    return shifts


def correct_phase(movie, m_phasecor, freqMC=.5, max_dev=5, plot_name="phase_shift.png", verbose = False, shifts = None):
    if shifts is None:
        from . import cmovie
        freq = movie.fr
        n_rebin = int(np.round(freq / freqMC))
        if n_rebin < 1:
            n_rebin = 1
        if verbose: print(f"The original frequency is {freq}.")
        if n_rebin > 1:
            if verbose: print(f"The movie will be rebinned in time by {n_rebin} for shifts extraction.")
        reb_movie = movie.resize(1, 1, 1. / n_rebin).astype("float32")

        fig, ax = plt.subplots(1, 1, figsize = (8, 3.5))
        ax.set_title("shifts for odd lines")
        ax.set_ylabel("shifts (in pixels)")
        ax.set_xlabel("time [s]")

        # first pass on the pairs of neighboring frames:
        dshifts = []
        for iframe in range(reb_movie.shape[0]):
            evenFrame = reb_movie[iframe,::2]
            evenFrame = (evenFrame[1:]+evenFrame[:-1])/2
            oddFrame = reb_movie[iframe,1::2][:-1]
            tmpshift = cmovie.extract_shifts(oddFrame.reshape((1,)+oddFrame.shape), max_dev, 0, template=evenFrame)[0][0]
            dshifts += [tmpshift]
        shifts = np.array(dshifts).astype("float32")

        ax.plot(np.arange(len(shifts)) / reb_movie.fr, shifts[:,1], label = "pairwise inference", c = "darkgrey")
        # ax.plot(np.arange(len(shifts)) / reb_movie.fr, np.array(dshifts), label = "pairwise inference", c = "darkgrey")
        shifts[:,0] = 0

        if verbose: print(f'Extracting shifts finished.')
        plt.tight_layout()
        fig.savefig(plot_name)

        if n_rebin > 1:
            from . import cmovie
            shifts = cmovie(shifts.reshape((1,) + shifts.shape)).resize(n_rebin, 1, 1)[0]
        # pad the remaining frames if exist:
        npad = len(movie) - len(shifts)
        if npad > 0:
            shifts = np.vstack([shifts, [shifts[-1]] * npad])

    # prepare for saving
    chunkSize = 3e9 # gigabytes
    singleFrameSize32 = movie[0].astype("float32").nbytes
    di = int(chunkSize/singleFrameSize32)+1
    range_ = range(0, len(shifts), di)
    if verbose:
        print(f'Applying shifts...', flush=True)
        range_ = tqdm(range_)

    for i in range_:
        sl = slice(i, min(len(shifts), i + di))
        tmpm = movie[sl,1::2].astype("float32").apply_shifts(shifts[sl])
        tmpm[tmpm < 0] = 0
        # tmpm[tmpm > maxv] = maxv
        m_phasecor[sl,1::2] = tmpm.astype(m_phasecor.dtype)
    if verbose: print("Done.")
    return shifts

def load_json(path):
    from json import loads
    from .Regions import Regions
    with open( path ) as f:
        txt = f.read()
    data = loads(txt)
    if "Freq" in data:
        data["Freq"] = float(data["Freq"])
    data["df"] = pd.DataFrame(data["df"])
    data["df"]["trace"] = data["df"]["trace"].apply(np.array)
    for k in data["statImages"]:
        data["statImages"][k] = np.array(data["statImages"][k])
    regions = Regions(dict(zip(data["df"]["peak"].apply(tuple), data["df"]["pixels"])))
    regions.df.index = [int(i) for i in regions.df.index]
    regions.df["pixels"] = [[tuple(px) for px in pxs] for pxs in regions.df["pixels"]]
    for k in list(data.keys()):
        setattr(regions, k, data[k])
        del data[k]
    regions.df["peak"] = regions.df["peak"].apply(tuple)
    if "interest" in regions.df.columns:
        regions.df["activity"] = regions.df["interest"]
        del regions.df["interest"]
    regions.df["peak"] = regions.df["peak"].apply(tuple)
    regions.df["pixels"] = [[tuple(px) for px in pxs] for pxs in regions.df["pixels"]]
    regions.df.index=list([int(j) for j in regions.df.index])
    regions.image = np.array(regions.image)
    regions.time = np.array(regions.time)
    if hasattr(regions,"metadata"):
        regions.metadata = pd.Series(regions.metadata)
    regions.update()
    regions.pathToRois = path
    try:
        folder = os.path.split(path)[0]
        protocolFile = os.path.join(folder, [f for f in os.listdir(folder) if "protocol" in f][0])
        regions.import_protocol(protocolFile)
    except:
        pass
    return regions

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
                saving = ['statImages', 'mode', 'image', 'filterSize', 'df', 'trange', "FrameRange", "analysisFolder", "time", "Freq","metadata","TwoParFit"]
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
                    if k not in ["peak", "pixels", "peakValue","tag","activity"]+col:
                        del subRegions.df[k]
                        
                roifile = f"{outDir}/{filename}_rois.pkl"
                with open(roifile,"wb") as f:
                    pickle.dump(subRegions,f)
                regions.pathToPickle = roifile
                if regions.is_3D:
                    feedback += ["Could not create a preview image for 3D regions. Not implemented yet."]
                else:
                    create_preview_image(regions)
                feedback += [f"ROI info saved in {roifile}."]

            elif format=="maribor":
                logging.warning(f"you requested to save multiple traces, but maribor format supports only one, so only first ({col[0]}) will be used.")
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

def add_protocol(ax, protocol, color="lightgrey", location="bottom"):
    yl = ax.get_ylim()
    dy = yl[1]-yl[0]
    if location not in ["top","bottom","all"]:
        raise ValueError("location can be either 'top' or 'bottom'")
    if location=="all":
        offset = 0
    else:
        dy = dy/20
    if location=="top":
        offset = yl[1] + dy
    if location=="bottom":
        offset = yl[0] - dy

    for comp, df in protocol.groupby("compound"):
        for ii in df.index:
            t0,t1 = df.loc[ii,["t_begin","t_end"]]
            conc = df.loc[ii,"concentration"]
            x,y = [t0,t1,t1,t0,t0],[1,1,0,0,1]
            y = np.array(y)
            y = y*dy + offset
            c = df.loc[ii].get("color",color)
            ax.fill(x,y,color=c,)
            ax.text(t0,y.mean()-dy*0.15, " "+conc,va="center", ha="left")
            ax.plot(x,y,color="black",lw=1)
        ax.text(df.t_begin.min(),y[:-1].mean(),comp+" ",va="center", ha="right")
        offset += 1.3*dy*(-1)**int(location=="bottom")


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
                         font_size=10
                        )
        rreg.plotEdges(ax=ax,image=True,ix=cl, spline=False,scaleFontSize=0,bound=False)
        attr = sum(list(map(list,nx.attracting_components(gph))),[])
        rreg.plotEdges(ax=ax,image=False,ix=attr,color="red",spline=False,scaleFontSize=0,bound=False)
        rreg.plotPeaks(ax=ax,image=False,ix=attr,color="red",ms=1)
        for sp in ax.spines: ax.spines[sp].set_visible(False)
        ax.set_aspect("equal")

    for i in range(i+1,axs.size):
        axs.flat[i].remove()
    plt.subplots_adjust(wspace=3e-2,hspace=3e-2)


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
                logging.warning(f"{j} is listed as a neighbor of {i}, but it does not exist")
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
    if len(peak2bnd):
        peak2bnd = pd.concat([dfi.sort_values("size_to").iloc[[-1]] for i, dfi in peak2bnd.groupby("i")])
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


def _init_logger():
    """This is so that Javabridge doesn't spill out a lot of DEBUG messages
    during runtime.
    Originally from CellProfiler/python-bioformats.
    We copied it from https://github.com/pskeshu/microscoper
    """
    from bioformats import javabridge as jb

    rootLoggerName = jb.get_static_field("org/slf4j/Logger",
                                         "ROOT_LOGGER_NAME",
                                         "Ljava/lang/String;")

    rootLogger = jb.static_call("org/slf4j/LoggerFactory",
                                "getLogger",
                                "(Ljava/lang/String;)Lorg/slf4j/Logger;",
                                rootLoggerName)

    logLevel = jb.get_static_field("ch/qos/logback/classic/Level",
                                   "ERROR",
                                   "Lch/qos/logback/classic/Level;")

    jb.call(rootLogger,
            "setLevel",
            "(Lch/qos/logback/classic/Level;)V",
            logLevel)

def getStatImages(movie_, debleach=False, downsampleFreq=1):
    if movie_.fr>downsampleFreq:
        n_rebin = int(np.round(movie_.fr/downsampleFreq))
        if n_rebin>3:
            # m_for_image = movie_.resize(1,1, 1/n_rebin)
            from .movies import movie as cmovie
            from .numeric import rebin
            m_for_image = cmovie(rebin(movie_,n_rebin), fr=movie_/n_rebin)
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


def saveMovie(movie, filename, maxFreq=1, maxdim=256, showtime=True, **kwargs):
    from .utils import show_movie
    from .numeric import rebin
    if maxFreq<movie.fr:
        nrebinT = int(np.ceil((movie.fr / maxFreq)))
    else:
        nrebinT = 1
    nrebinV = int(np.round(movie.shape[1]/maxdim))
    nrebinV = max(nrebinV,1)
    nrebinH = int(np.round(movie.shape[2]/maxdim))
    nrebinH = max(nrebinH,1)
    if (nrebinT>1) or (nrebinV>1) or (nrebinH>1):
        showMovie = rebin(movie, (nrebinT,nrebinV,nrebinH), axis=(0,1,2))
        from . import cmovie
        showMovie = cmovie(showMovie, fr=movie.fr / nrebinT)
    else:
        showMovie = movie+1

    if filename=="embed":
        return show_movie(showMovie,
                          NTimeFrames=len(showMovie),
                          tmax=len(showMovie)/showMovie.fr if showtime else None,
                          **kwargs
                          )
    else:
        show_movie(showMovie,
                   out="save",
                   saveName=filename,
                   NTimeFrames=len(showMovie),
                   tmax=len(showMovie)/showMovie.fr if showtime else None,
                   **kwargs
                   )
        return 0


def autocorr2d(sett, dxrange, dyrange):
    from numpy import zeros, corrcoef
    Nx, Ny = sett.shape
    ret = zeros((len(dxrange), len(dyrange)))
    for kx,dx in enumerate(dxrange):
        for ky,dy in enumerate(dyrange):
            ret[kx,ky] = corrcoef(sett[  :Nx-dx,   :Ny-dy].flatten(),
                                  sett[dx:     , dy:     ].flatten())[0,1]
    return ret

def get_trace(row):
    return tuple([np.array(row[k].replace("\n"," ").strip(" []").split(), dtype="float") for k in ["time","trace"]])

def get_evshape(row, regions,):
    ts = row.ts
    sts = "%g"%ts
    roi = row.roi
    t = regions.showTime.get(sts, regions.time)
    it0 = np.searchsorted(t,row.t0)
    dt = t[1]-t[0]
    dhw = int(np.ceil(row.halfwidth/dt))
    dix = max(5,(dhw//5*5))
    slc = slice(it0-dix, it0+dix*4 )
    x = regions.df.loc[roi,"faster_%g"%ts][slc]# + regions.df.loc[roi,"slower_%g"%ts][slc]
    return x


def import_data(mainFolder, constrain="", forceMetadataParse=False, verbose=0):
    from .general_functions import td2str
    from tqdm import tqdm
    from .Recording import Recording, parse_leica
    # tqdm().pandas()
    recordings = []
    for cur, ds, fs in os.walk(mainFolder):
        #### if you wish to restrict to only certain folders: ####
        for f in fs:
            if not (f.endswith(".lif") or f.endswith(".nd2") or f.endswith(".czi")):
                continue
            if any([constr.strip() in cur + f for constr in constrain.split(",")]):
                path = os.path.join(cur, f)
                recordings += [path]
    recordings = sorted(recordings)

    status = []
    ilifs = 0
    for pathToRecording in tqdm(recordings):
        if verbose >= 1:
            print("#" * 20, pathToRecording)
        try:
            rec = Recording(pathToRecording)
        except:
            logging.warning("Could not import %s" % pathToRecording)
            continue
        recType = "Leica" if pathToRecording.endswith(".lif") else "nonLeica"
        if forceMetadataParse:
            rec.parse_metadata()
            rec.save_metadata()
        if recType == "Leica":
            try:
                sers = parse_leica(rec)
            except:
                sers = list(rec.metadata['Name'])
                
        else:
            sers = [rec.Experiment.split(".")[0]]

        analysisFolder = os.path.join(rec.folder, rec.Experiment + "_analysis")
        if not os.path.isdir(analysisFolder):
            continue
        existingSeries = []
        for fs in os.listdir(analysisFolder):
            if not os.path.isdir(os.path.join(analysisFolder, fs)): continue
            if len(os.listdir(os.path.join(analysisFolder, fs))) == 0: continue
            if fs[0] == ".": continue
            ############
            fssplit = fs.split('_')
            if len(fssplit) == 1:
                existingSeries += [fs]
                continue
            trange = fssplit[-1].split("-")
            if len(trange) > 1:
                fs = "_".join(fssplit[:-1])
            existingSeries += [fs]
        if verbose > 1:
            print(existingSeries)
        sers = np.unique(sers + existingSeries)
        for series in sers:
            subdirs = get_series_dir(pathToRecording, series)
            if verbose >= 2:
                print("series=", series, ", with subdirs:", subdirs)
            for ser in subdirs:
                if verbose >= 2:
                    print("ser=", ser)
                if recType != "Leica":
                    series = "all"
                try:
                    rec.import_series(series, onlyMeta = True)
                except:
                    continue
                if not hasattr(rec, "Series") or series not in rec.Series:
                    continue
                md = pd.Series()
                md["path to exp"] = pathToRecording
                md["experiment"] = os.path.split(pathToRecording)[-1]
                md["series"] = series

                saveDir = os.path.join(analysisFolder, ser)
                if len(rec.Series) == 0 or "metadata" not in rec.Series[series]: continue
                for k, v in rec.Series[series]["metadata"].items():
                    md[k] = v
                if "_" in ser:
                    fssplit = ser.split("_")
                    trange = fssplit[-1].split("-")
                    if verbose > 2:
                        print(trange)
                    if len(trange) >= 2:
                        try:
                            t0, t1 = [float(t.strip("s")) for t in fssplit[-1].split("-", maxsplit = 1)]
                            md["Time Range"] = "%i-%i" % (t0, t1)
                            md["Duration [s]"] = t1 - t0
                        except:
                            print("Oops, having problems parsing ", ser)
                            continue
                    else:
                        md["Time Range"] = "all"
                        md["Duration [s]"] = md["SizeT"] / md["Frequency"]
                else:
                    md["Time Range"] = "all"
                    md["Duration [s]"] = md["SizeT"] / md["Frequency"]
                fs = get_filterSizes(md.pxSize)
                movies = [os.path.join(saveDir, fn) for fn in os.listdir(saveDir) if fn.endswith(".mp4")]
                if len(movies):
                    movies = sorted(movies, key = lambda xi: os.stat(xi).st_mtime)
                    md["movie done"] = True
                    movieFilename = movies[-1]
                    md["path to movie"] = movieFilename
                else:
                    md["path to movie"] = "None"
                    md["movie done"] = False

                if md["movie done"]:
                    md["movie size [MB]"] = np.round(os.path.getsize(movieFilename) / 10 ** 6, 1)
                md["date"] = md["Start time"].date().__str__()
                for k in ["bit depth", "End time", "Name", "frame_range"]:  # , "individual Series"
                    try:
                        del md[k]
                    except:
                        pass
                times = ["00:00"] + [td2str(el) for el in md["individual Series"]["Duration"].cumsum()]
                md["Duration"] = times[-1]
                md["Series Durations"] = " \n".join(["%s [%s-%s]" % (name.lstrip("Series0"), t0, t1) for name, t0, t1 in
                                                     zip(md["individual Series"]["Name"], times[:-1], times[1:])])
                del md["individual Series"]
                pklsDone = {}
                for fsize in fs:
                    pickleFile = os.path.join(saveDir, ".".join(map(str, fsize)) + "_rois.pkl")
                    pickleThere = os.path.isfile(pickleFile)
                    pklsDone[fsize] = pickleThere
                md["pickles done"] = pklsDone
                if md["movie done"]:
                    pathToProtocol = movieFilename.replace(".mp4", "_protocol.txt").replace("_corrected", "").replace(
                        "_original", "")
                else:
                    pathToProtocol = os.path.join(saveDir, "%s_%s_protocol.txt" % (md["experiment"], md["series"]))
                md["path to protocol"] = pathToProtocol
                md["protocol done"] = False
                try:
                    protocol = pd.read_csv(pathToProtocol)
                    if len(protocol):
                        md["protocol done"] = True
                        protocol = " ".join(np.unique([
                            "%s:%s" % (
                            row["compound"].capitalize() if "glu" in row["compound"].lower() else row["compound"],
                            row["concentration"].replace(" ", "")) for _, row in protocol.iterrows()]))
                        protocol = protocol.replace("Glucose", "Glu")
                        md["protocol"] = protocol
                except:
                    pass
                pathToAddInfo = os.path.split(pathToProtocol)[0]
                pathToAddInfo = os.path.join(pathToAddInfo, "additional_info.txt")
                md["path to add_info"] = pathToAddInfo
                md["add_info done"] = os.path.isfile(pathToAddInfo)
                if md["add_info done"] and os.path.getsize(pathToAddInfo) > 10:
                    # print (ser, )
                    try:
                        addInfo = pd.read_csv(pathToAddInfo, sep = ":", header = None, index_col = 0).T
                    except:
                        md["add_info done"] = False
                        continue
                    if len(addInfo) == 0:
                        md["add_info done"] = False
                        continue
                    for kk in addInfo.columns:
                        # print ("%s:%s"%(kk, addInfo[kk].iloc[0]), end=" ")
                        md[str(kk).strip()] = str(addInfo[kk].iloc[0]).strip()
                status += [dict(md.items())]

        ilifs += 1
    #     if ilifs>3:
    #         break
    status = pd.DataFrame(status)
    if "protocol" not in status.columns:
        status["protocol"] = [""] * len(status)
    return status


def import_data_new(mainFolder, constrain="", forceMetadataParse=False, verbose=0, extensions=None):
    from islets.general_functions import td2str
    from tqdm import tqdm
    from islets.Recording import Recording

    if extensions is None:
        extensions = ["tif", "tiff", "nd2", "lif"]
    recordings = []
    for cur, ds, fs in os.walk(mainFolder):
        for f in fs:
            if not any([f.endswith(ext) for ext in extensions]):
                continue
            if any([constr.strip() in cur + f for constr in constrain.split(",")]):
                path = os.path.join(cur, f)
                recordings += [path]
    recordings = sorted(recordings)

    status = []
    for pathToRecording in tqdm(recordings):
        if verbose >= 1:
            print("#" * 20, pathToRecording)
        try:
            rec = Recording(pathToRecording)
        except:
            print ("Could not import %s" % pathToRecording)
            continue
        # if pathToRecording.endswith(".nd2"):
        #     recType = "Nikon"
        # if pathToRecording.endswith(".lif"):
        #     recType = "Leica"
        # if pathToRecording.endswith("tif") or pathToRecording.endswith("tiff"):
        #     recType = "tif"
        if forceMetadataParse:
            rec.parse_metadata()
            rec.save_metadata()

        analysisFolder = os.path.join(rec.folder, rec.Experiment + "_analysis")
        if not os.path.isdir(analysisFolder):
            continue
        subDirs = []
        for fs in os.listdir(analysisFolder):
            if fs[0] == ".": continue
            fullpath = os.path.join(analysisFolder, fs)
            if os.path.isdir(fullpath) and len(os.listdir(fullpath)) > 0:
                subDirs += [fullpath]

        for subdir in subDirs:
            md = pd.Series()
            md["folder"] = subdir
            md["path to exp"] = pathToRecording
            md["experiment"] = os.path.split(pathToRecording)[-1]

            fs = os.path.split(subdir)[1]
            if "_c" in fs:
                fs, channel = fs.split("_c")
                md["channel"] = int(channel)
            else:
                md["channel"] = 0
            if "_" in fs:
                fs, trange = fs.split('_', maxsplit=1)
                md["time_range"] = trange
            else:
                trange = None
            md["series"] = fs

            rec.import_series(fs, onlyMeta=True, restrict=trange, channel=md["channel"])

            for k, v in rec.Series[md["series"]]["metadata"].items():
                md[k] = v

            md["Duration [s]"] = md["SizeT"] / md["Frequency"]

            movies = [os.path.join(subdir, fn) for fn in os.listdir(subdir) if fn.endswith(".mp4")]
            md["movies"] = sorted(movies, key=lambda xi: os.stat(xi).st_mtime)[::-1]
            images = [os.path.join(subdir, fn) for fn in os.listdir(subdir) if fn.endswith(".png")]
            md["images"] = sorted(images, key=lambda xi: os.stat(xi).st_mtime)[::-1]

            md["date"] = md["Start time"].date().__str__()
            for k in ["End time", "Name", "frame_range"]:  # , "individual Series"
                try:
                    del md[k]
                except:
                    pass
            times = ["00:00"] + [td2str(el) for el in md["individual Series"]["Duration"].cumsum()]
            md["Duration"] = times[-1]
            md["Series Durations"] = " \n".join(["%s [%s-%s]" % (name.lstrip("Series0"), t0, t1) for name, t0, t1 in
                                                 zip(md["individual Series"]["Name"], times[:-1], times[1:])])
            del md["individual Series"]

            protocols = [os.path.join(subdir, fn) for fn in os.listdir(subdir) if fn.endswith("protocol.txt")]
            if len(protocols) == 0:
                pathToProtocol = os.path.join(subdir, "%s_%s_protocol.txt" % (md["experiment"], md["series"]))
            elif len(protocols) == 1:
                pathToProtocol = os.path.join(subdir, protocols[0])
            else:
                pathToProtocol = os.path.join(subdir, sorted(protocols, key=lambda xi: os.path.getsize(xi))[-1])

            md["path to protocol"] = pathToProtocol
            try:
                protocol = pd.read_csv(pathToProtocol)
                if len(protocol):
                    protocol = " ".join(set([
                        "%s:%s" % (
                        row["compound"].capitalize() if "glu" in row["compound"].lower() else row["compound"],
                        row["concentration"].replace(" ", "")) for _, row in protocol.iterrows()]))
                    protocol = protocol.replace("Glucose", "Glu")
                    md["protocol"] = protocol
            except:
                pass

            pathToAddInfo = os.path.join(subdir, "additional_info.txt")
            md["path to add_info"] = pathToAddInfo

            try:
                addInfo = pd.read_csv(pathToAddInfo, sep=":", header=None, index_col=0).T
                addInfo = addInfo.iloc[0]
                for kk, vv in addInfo.items:
                    md[str(kk).strip()] = str(vv).strip()
            except:
                pass

            status += [dict(md.items())]

    status = pd.DataFrame(status)
    if "protocol" not in status.columns:
        status["protocol"] = [""] * len(status)
    return status
