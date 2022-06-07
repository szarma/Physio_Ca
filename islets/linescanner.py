import json
from sys import exc_info
import numpy as np
import plotly.graph_objects as go
from dash.dependencies import Input, Output
try:
    from jupyter_plotly_dash import JupyterDash
except ModuleNotFoundError:
    from jupyter_dash import JupyterDash
try:
    from dash import dcc, html
except ImportError:
    import dash_core_components as dcc
    import dash_html_components as html

from .numeric import rebin
from .utils import getFigure


# import pandas as pd
# import matplotlib.pyplot as plt


def take_part(ls, key, extent):    
    if extent is None:
        extent = (None,)*4
    x0,x1,t0,t1 = extent
    if x0 is None: x0 = 0
    if x1 is None: x1 = ls.data.shape[0]
    if t0 is None: 
        it0 = 0
    else:
        it0 = np.searchsorted(ls.time, t0)
    if t1 is None: 
        it1 = ls.data.shape[1]
    else:
        it1 = np.searchsorted(ls.time, t1)
    x0 = int(x0)
    x1 = min(int(x1)+1, ls.data.shape[0])
    if x1-x0<=1: x1=x0+1
    try:
        getattr(ls,key)
    except:
        timescale = float(key.split("_")[-1])
        ls.fast_filter_traces(timescale,z_sp=0)
    return x0,x1,it0,it1

def plot_heatmap(ls, key="zScore_2", extent=None, NTimepoints=300):
    from .numeric import robust_max
    x0,x1,it0,it1 = take_part(ls,key,extent)
    nt = it1-it0
    nr = int(nt/NTimepoints)
    Z = getattr(ls,key)[x0:x1, it0:it1]
    tr = ls.time[it0:it1]
    if nr>1:
        z = rebin(Z,nr,1)
        tr = rebin(tr, nr)
    else:
        z=Z
    fig = go.Figure(
        data=go.Heatmap(z=z,
                        showscale=False, 
                        zmin=-robust_max(z,min(100,int(z.size*.1))) if ("zScore" in key or key=="detrended") else .5,
                        zmax=robust_max(z,min(100,int(z.size*.1))),
                        colorscale="RdBu_r" if ("zScore" in key or key=="detrended") else "hot",
                        x=tr,
                        y=np.arange(x0,x1), 
                        hoverinfo="skip",
                        
                       ),
        layout={
            "dragmode":'select',
            "width":350,
            "height":380,
            "margin":{"t":20,**dict(zip("lbr",[3]*3))},
            "xaxis":{"title":"time [s]"},
            "yaxis":{"title":"line index"}
        }
    )
    fig.add_traces(go.Scatter(x=[ls.time[[it0,it1-1]].mean()],y=[(x0+x1)/2],opacity=0, hoverinfo='skip'))
    
#     fig.add_annotation(
#                 x=0.05,
#                 y=0.95,
#                 text=f"rebinned by {nr}",
#                 xref="paper",
#                 yref="paper",
#                 showarrow=False,
#                 font_size=20
#     )
    return fig


def plot_trace(ls, key="data", extent=None, nr=1):
    x0,x1,it0,it1 = take_part(ls,key,extent)
    Z = getattr(ls,key)[x0:x1, it0:it1]
    z = Z.mean(0)
    tr = ls.time[it0:it1]
    if nr>1:
        z = rebin(z,nr)
        tr = rebin(tr, nr)
    
    fig = go.Figure(
        data=go.Scatter(x=tr, y=z,line_width=1,line_color="darkgrey"),
        layout={
            "width":500,
            "height":390,
            "margin":{"t":30,**dict(zip("lbr",[3]*3))},
            "xaxis":{"title":"time [s]"}
        }
                    )
    return fig


# overviewGraph = 
# zoomGraph     = dcc.Graph(id='zoomed-fig',  figure=getFigure(),)
styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll',
        'overflowY': 'scroll',
        'height':"300px",
        'width':"400px",
    },
    "flex":{
        "padding":"10px","max-width":"550px","max-height":"800px","border":'thin lightgrey solid',
    }
}

def examine(linescan, debug=False):
    marks = np.array([i*10**j for i in [1,2,5] for j in range(5)])
    marks.sort()
    marks = marks[marks<linescan.Freq]
    app = JupyterDash(__name__, width=1500, height=1200)

    app.layout = html.Div([
        html.H1("LineScanner"),
        html.Div([
            html.Div([
                html.H3(graphID.capitalize()+" figure"),
                dcc.RadioItems(id=f"{graphID}-radio",
                               options=[{"label":"raw data","value":"data"}, {"label":"filtered","value":"zScore_2"}, {"label":"detrended","value":"detrended"}], 
                               value= "detrended", 
    #                            value= {"zoom":None,"overview":"zScore_2"}[graphID]
                              ),
                dcc.Graph(id=f'{graphID}-fig',figure=getFigure(),)
            ]+[
                html.Div(
                    ["rebin by",
                     html.Div(dcc.Slider(id="rebin-slider", 
                                         min=0, max=len(marks)-1,
                                         updatemode="drag",
                                         value=np.searchsorted(marks,linescan.Freq/10)-1
    #                                      marks = {i:{"label":str(marks[i])} for i in range(len(marks))}
                                        ), style={"width":"300px"}),
                     html.Div(id="slider-out",children=""),
    #                  dcc.Input(id="rebin-value",type=int)
                    ],style={"display":"flex", "flex-wrap":"wrap","align-items":"bottom","flex-shrink":3}),
              ]*int(graphID=="trace"),style=styles['flex']) for graphID in ["overview", "zoom","trace"]
        ],style={"display":"flex", "flex-wrap":"wrap","align-items":"top","flex-shrink":3}),
        html.Div([
            html.Pre(id='relayout-data', style={"display":"block" if debug else "none"}),
            html.Pre(id="output",        style={"display":"block" if debug else "none"}),
            html.Pre(id="output2",       style={"display":"block" if debug else "none"}),
        ],style={"display":"flex", "flex-wrap":"wrap","align-items":"top","flex-shrink":3}),
        ])

    #### callbacks
    @app.callback(
        Output('slider-out', 'children'),
        [Input("rebin-slider","value")]
    )
    def show_rebin(value):
        return str(marks[value])

    @app.callback(
        Output('overview-fig', 'figure'),
        [Input("overview-radio","value")]
    )
    def show_overview(value):
        return plot_heatmap(linescan, value )


    @app.callback(
        Output('relayout-data', 'children'),
        [Input('overview-fig', 'selectedData')])
    def display_data(data_):
        return json.dumps(data_["range"], indent=2)


    @app.callback(
        [Output('zoom-fig', 'figure'),
        Output('output', 'children')],
        [Input('overview-fig', 'selectedData'), 
         Input("zoom-radio","value")])
    def replot(data_,radioValue):
        try:
            t0,t1 = data_["range"]["x"]
            x0,x1 = data_["range"]["y"]
            return plot_heatmap(linescan, radioValue, (x0,x1,t0,t1), ), "ok"
        except:
            fig = getFigure(350,380,"white")
            text="Choose a box <br> select tool <br> in the overview plot <br> to select area <br> to zoom in"
            fig.add_annotation(
                        x=0.5,
                        y=0.5,
                        text=text,
                        xref="paper",
                        yref="paper",
                        showarrow=False,
                        font_size=20
            )
            return fig, str(exc_info())

    @app.callback(
        [Output('trace-fig', 'figure'),
        Output('output2', 'children')],
        [Input('overview-fig', 'selectedData'), 
         Input('zoom-fig', 'selectedData'), 
         Input("trace-radio","value"),
         Input("rebin-slider","value")
        ])
    def plot_trace_callback(data_ov, data_zoom, radioValue, i_rebin):
        try:
            if data_zoom is None:
                if data_ov is None:
                    data_ = None
                else:
                    data_ = data_ov
            else:
                data_ = data_zoom
            if data_ is None:
                extent = None
            else:
                t0,t1 = data_["range"]["x"]
                x0,x1 = data_["range"]["y"]    
                extent = (x0,x1,t0,t1)
            nr = int(marks[i_rebin])
            try:
                fig = plot_trace(linescan, radioValue, extent, nr=nr)
            except:
                fig = getFigure()
    #         text = f"rebinned by {nr}"
    #         if extent is not None:
    #             text += "<br>extent=(%.1f,%.1f,%.1f,%.1f)"%extent
    #         fig.add_annotation(
    #                     x=0.05,
    #                     y=0.95,
    #                     text=text,
    #                     xref="paper",
    #                     yref="paper",
    #                     showarrow=False,
    #                     font_size=20
    #         )
            return fig,"ok"
        except:
            return getFigure(), str(exc_info())


    return app