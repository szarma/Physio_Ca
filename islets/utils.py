import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib
import networkx as nx
import pandas as pd

def coltrans(x, vmin=None, vmax=None, tilt=1, offset=0.1):
    from .numeric import robust_max
    try:
        x[0]
        iterable = True
        if vmax is None:
            vmax = robust_max(x)
        if vmin is None:
            vmin = -robust_max(-x)
    except:
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

def filter_spikes_per_roi(roiSpikes, regions, plot=False, halfwidth_toll=.33, freqShow=5):
    roi = roiSpikes["roi"].unique()
    assert len(roi)==1
    roi = roi[0]
    origtrace = regions.df.loc[roi,"trace"]
    tsslows = {ts: regions.df.loc[roi,"slower_%g"%ts] for ts in roiSpikes.ts.unique()}
    roiSpikes = pd.DataFrame(roiSpikes)
    t = regions.time
    roiSpikes["status"] = "ok"
    
    if plot:
#         currentBackend = matplotlib.get_backend()
#         if plot not in ["interactive","static"]:
#             raise ValueError("plot can take only three values: None, 'interactive' and 'static'.")
#         if plot=="interactive" and currentBackend!="nbAgg":
#             plt.switch_backend("nbAgg")
#         if plot=="static" and currentBackend!='module://ipykernel.pylab.backend_inline':
#             plt.switch_backend('module://ipykernel.pylab.backend_inline')
#         plot = plot is not None
        fig, axs = plt.subplots(3,1,figsize=(13,14), sharex=True, sharey=True)
        
        def draw_spike(row, ax):
            ax.fill(
                row.t0+row.halfwidth*np.array([0,0,1,1,0]),
                row.its+np.array([0,1,1,0,0])*.8-.5,
                color=row.color,
                #linewidth=.7,
            )
        for ix, row in roiSpikes.iterrows():
            draw_spike(row, axs[1])
    colorDict = {}
    for ix, row in roiSpikes.iterrows():
        it0 = np.searchsorted(t,row.t0)-1
        ite = np.searchsorted(t,row.t0+row.halfwidth)+1
        spikeLength = ite-it0
        if spikeLength<3:
            status = "too short"
            roiSpikes.loc[ix,"status"] = status
            if status not in colorDict:
                colorDict[status] = axs[0].plot([],label=status)[0].get_color()
            row = pd.Series(row)
            row.color = colorDict[status]
            if plot: draw_spike(row,axs[0])
            continue
        slowtrace = tsslows[row.ts][it0:ite].copy()
        if row.ts*4 in tsslows:
            slowtrace-=tsslows[row.ts*4][it0:ite]
        else:
            if row.ts*2 in tsslows:
                slowtrace-=tsslows[row.ts*2][it0:ite]
        ddslow = np.diff(np.diff(slowtrace))
        ddorig = np.diff(np.diff(origtrace[it0:ite]))
        # if ddorig.mean()>0:
        if spikeLength>20 and ddslow.mean()>0.1 and ddslow.mean()>-ddorig.mean()*.5:
            print (ddslow.mean(), ddorig.mean())
            status = "spurious"
            roiSpikes.loc[ix,"status"] = status
            if status not in colorDict:
                colorDict[status] = axs[0].plot([],label=status)[0].get_color()
            row = pd.Series(row)
            row.color = colorDict[status]
            if plot: draw_spike(row, axs[0])
            if plot: axs[0].text(row.t0+row.halfwidth/2, row.its+.3, "x", va="bottom", ha="center")
            continue
            
        tcur = t[it0:ite]
        bkg = tsslows[row.ts][it0]+np.arange(spikeLength)*(tsslows[row.ts][ite-1]-tsslows[row.ts][it0])/(spikeLength-1)
        purePeak = origtrace[it0:ite]-bkg
        amplitude = np.percentile(purePeak, 95)
        if amplitude/bkg.mean()<.1:
#         if amplitude<1.*np.abs(bkg[-1]-bkg[0]):
            status = "too steep"
            roiSpikes.loc[ix,"status"] = status
            if plot:
                if status not in colorDict:
                    colorDict[status] = axs[0].plot([],label=status)[0].get_color()
                row = pd.Series(row)
                row.color = colorDict[status]
                draw_spike(row,axs[0])
                axs[0].text(row.t0+row.halfwidth/2, row.its+.3, ("/" if bkg[-1]>bkg[0] else "\\"), va="bottom", ha="center")
            
    if len(colorDict):
        axs[0].legend()
    g = nx.DiGraph()
    for isp, spike in roiSpikes.query("status=='ok'").sort_values("halfwidth").iterrows():
        conflicts = roiSpikes.query(  f"t0>{spike.t0   - spike.halfwidth*halfwidth_toll}")
        conflicts = conflicts.query(  f"t0<{spike.t0   + spike.halfwidth*halfwidth_toll}")
        conflicts = conflicts.query(f"tend>{spike.tend - spike.halfwidth*halfwidth_toll}")
        conflicts = conflicts.query(f"tend<{spike.tend + spike.halfwidth*halfwidth_toll}")
        conflicts = conflicts.query(f"index!={spike.name}")
        g.add_node(isp)
        for other in list(conflicts.index):
            g.add_edge(isp, other)
    if plot:
        ax = axs[1]
        nx.draw_networkx(g,
                         ax=ax,
                         with_labels=False,
                         pos={node: (roiSpikes.loc[node, "t0"]+.5*roiSpikes.loc[node, "halfwidth"], roiSpikes.loc[node, "its"]) for node in g.nodes},
                         node_size=0,
                         edge_color="grey"
                        )
        #ax.set_yticks(np.arange(len(timescales)));
        #ax.set_yticklabels(timescales);
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    ######################################################
    df_filt = []
    for ixs in nx.connected_components(g.to_undirected()):
        row = pd.Series(roiSpikes.loc[list(ixs)].sort_values("height").iloc[-1])
        for col in ["t0","tend","height"]:
            row.col = roiSpikes.loc[list(ixs),col].mean()
        row.halfwidth = row.tend-row.t0
        df_filt += [row]
    df_filt = pd.DataFrame(df_filt)
    for col in ["loghalfwidth","score","status"]:
        if col in df_filt.columns:
            del df_filt[col]
    if plot:
#         if len(df_filt):
#             df_filt["color"] = list(plt.cm.hot(coltrans(df_filt.height, vmax=10)))
        for ix, row in df_filt.iterrows():
            draw_spike(row, axs[2])
        nr = np.round(regions.Freq/freqShow)
        if nr>1:
            nr = int(nr)
            from .numeric import rebin
            x = rebin(t, nr)
            y = rebin(origtrace, nr)
        else:
            x,y = t, origtrace
#         y = y - np.percentile(y,1)
#         y = y / np.percentile(y,99)
#         y = y * roiSpikes.its.max() + roiSpikes.its.max()/3
        for ax in axs:
            ax.set_facecolor("k")
            axx = ax.twinx()
            axx.set_ylabel("light intensity")
            axx.plot(x, y, c="grey", lw=.7)
            ax.set_ylabel("filtering timescale [s]")
            ax.set_xlabel("time [s]]")
            ax.set_ylim(ax.get_ylim()[::-1])
        axs[0].set_title("discarded")
        axs[1].set_title("all")
        axs[2].set_title("filtered")
        ax.set_yticks(roiSpikes.its.unique());
        ax.set_yticklabels(["%g"%ts for ts in roiSpikes.ts.unique()]);
        ax.set_xlim(roiSpikes.t0.min()-1, roiSpikes.tend.max()+1)
        plt.tight_layout()
        df_filt, axs
    return df_filt, None

def save_tiff(movie, movieFilename):
    import PIL
    im = PIL.Image.fromarray(movie[0])
    im.save(movieFilename,
            save_all=True,
            append_images=[PIL.Image.fromarray(movie[i]) for i in range(1,len(movie))],
            compression="tiff_deflate"
           )


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
    import numpy as np
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
        import matplotlib.patheffects as path_effects
        from pandas import Timedelta
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
    ax = fig.add_axes([0.01,0.01,.98,.98])
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
        return HTML(anim.to_jshtml())
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
    plt.ylim(plt.ylim()[::-1])
    plt.savefig(bkg_img_file,dpi=150)
    plt.close(fig)
    plt.switch_backend(currentBackend)
    if returnPath:
        return bkg_img_file
    
    return PilImage.open(bkg_img_file)


def saveRois(regions,outDir,filename="",movie=None,col=["trace"],formats=["vienna"],add_date=True):
        feedback = []
#     try:
        from copy import deepcopy,copy
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
                
                traces = pd.DataFrame(np.vstack(regions.df[col]).T)
                try:
                    traces["time"] = regions.showTime[col.split("_")[-1]]
                except:
                    traces["time"] = regions.time
                traces = traces[["time"]+list(traces.columns[:-1])]
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
