import numpy as np

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

def multi_map(some_function, iterable, processes=1):
    assert type(processes) == int
    if processes==1:
        out = map(some_function, iterable)
    elif processes>1:
        from multiprocessing import Pool
        pool = Pool(processes)
        out  = pool.map(some_function, iterable)
        pool.close()
        pool.join()
    else:
        print ("invalid number of processes", processes)
        quit()
    return out


def showMovie(m_show, figsize = (6,6), out="jshtml",fps = 30, saveName=None, NTimeFrames=100,log=True,additionalPlot=None):
    import matplotlib.pyplot as plt
    from matplotlib import animation
    if NTimeFrames is not None:
        n_rebin = len(m_show)//NTimeFrames
        if n_rebin>1:
            m_show = rebin(m_show, n_rebin)
    if log:
        for p in range(1,5):
            baseline = np.percentile(m_show,p)
            m_show = np.maximum(m_show, baseline)
            if np.all(m_show>0): break
        m_show = np.log(m_show)
    fig, ax = plt.subplots(figsize=figsize,dpi=150)
    im = ax.imshow(m_show[0].T, cmap="Greys", vmin=0, vmax=m_show.max())
    if additionalPlot is not None:
        additionalPlot(ax)
    plt.close(fig)
    def init():
        im.set_data(m_show[0].T)
        return (im,)
    def animate(i):
        im.set_data(m_show[i].T)
        return (im,)
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(m_show),
                                   interval=1000/fps,
                                   blit=True)
    if out=="html5":
        from IPython.display import HTML
        return HTML(anim.to_html5_video())
    if out=="jshtml":
        from IPython.display import HTML
        return HTML(anim.to_jshtml())
    if out=="save" or saveName is not None:
        try:
            anim.save(saveName)
        except:
            saveName = input("please enter a valid filename. Otherwise, I'll save it as 'video.mp4'.")
            try: anim.save(saveName)
            except:
                saveName = "video.mp4"
                anim.save(saveName)
        return None
    
def show_movie(m_show, figScale = 1, out="jshtml",fps = 30, saveName=None, NTimeFrames=100,log=True,additionalPlot=None, dpi=100):
    import matplotlib.pyplot as plt
    from matplotlib import animation
    if NTimeFrames is not None:
        n_rebin = len(m_show)//NTimeFrames
        if n_rebin>1:
            m_show = rebin(m_show, n_rebin)
    if log:
        for p in range(1,5):
            baseline = np.percentile(m_show,p)
            m_show = np.maximum(m_show, baseline)
            if np.all(m_show>0): break
        m_show = np.log(m_show)
    figsize = np.array(m_show.shape[1:])/dpi*figScale
    fig = plt.figure(figsize=figsize,dpi=dpi)
    ax = fig.add_axes([0.01,0.01,.98,.98])
    im = ax.imshow(m_show[0], cmap="Greys", vmin=m_show.min(), vmax=m_show.max())
    if additionalPlot is not None:
        additionalPlot(ax)
    plt.close(fig)
    def init():
        im.set_data(m_show[0])
        return (im,)
    def animate(i):
        im.set_data(m_show[i])
        return (im,)
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(m_show),
                                   interval=1000/fps,
                                   blit=True)
    if out=="html5":
        from IPython.display import HTML
        return HTML(anim.to_html5_video())
    if out=="jshtml":
        from IPython.display import HTML
        return HTML(anim.to_jshtml())
    if out=="save" or saveName is not None:
        try:
            anim.save(saveName, extra_args=['-vcodec', 'libx264'])
        except:
            saveName = input("please enter a valid filename. Otherwise, I'll save it as 'video.mp4'.")
            try: anim.save(saveName, extra_args=['-vcodec', 'libx264'])
            except:
                saveName = "video.mp4"
                anim.save(saveName, extra_args=['-vcodec', 'libx264'])
        return None
    
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


def showRoisOnly(regions, indices=None, im=None, showall=True):
    import plotly.graph_objects as go
    from .Regions import MYCOLORS
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

def createStaticImage(im,regions,showall=False,color="grey",separate=False, returnPath=False):
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
    if returnPath:
        return bkg_img_file
    
    return PilImage.open(bkg_img_file)
