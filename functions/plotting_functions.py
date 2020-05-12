from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.colors import DEFAULT_PLOTLY_COLORS
from plotly import colors as plxcolors
import numpy as np


def plotStatsFig(regions,ix=[]):
    fig1 = make_subplots(rows=1, cols=2, horizontal_spacing=0.02,
                         shared_yaxes=True,
                    subplot_titles=("Standard Deviation","StDev: Filtered Contrast-Enhanced"))
    for j,im in enumerate([regions.std_image, regions.image]):
        fig1.add_heatmap(z=im,col=j+1,row=1, hoverinfo='skip',showscale=False,colorscale=plxcolors.sequential.Greys)
    fig1.update_layout(
        height = 250,
        margin=dict(l=10, r=10, t=50, b=10),
        width = 500,
        title = {"font":{"size": 7}},
        xaxis = {
            "linecolor": 'black',
            "linewidth": 1,
            "mirror": True,
            "tickvals": [],
#             "range":[-.5,im.shape[1]-.5]
        },
        yaxis =  {
            "linecolor": 'black',
            "linewidth": 1,
            "mirror": True,
            "tickvals": [],
#             "range":[-.5,im.shape[1]-.5]
        }
    )
    if len(ix): fig1 = addRoisToFig(regions,ix,fig1)
    for c in [1,2]:
        fig1.update_xaxes({"range":[-.5,im.shape[1]-.5]},col=c,row=1,showticklabels=False)
        fig1.update_yaxes({"range":[-.5,im.shape[0]-.5]},col=c,row=1,showticklabels=False)
    return fig1

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
#                             showlegend = c==1,
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

def plotRawTraces(regions,indices=None, nAvg = 1):
    from physio_def_1 import rebin
    if indices is None:
        indices = regions.df.sort_values("size",ascending=False).index[:10]
    f = go.Figure()
    t = regions.time
    if nAvg > 1:
        t = rebin(t,nAvg)
        print ("rebinned with ",nAvg)
    for i in indices:
        cl = DEFAULT_PLOTLY_COLORS[i%10]
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
        "height":400,
        "width":800,
        "margin":dict(l=10, r=10, t=50, b=10),
    #     "paper_bgcolor":"grey",
        "xaxis": {
            "title": 'time [s]',
        },
    })
    return f