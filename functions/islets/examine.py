import numpy as np
import pandas as pd
# from itertools import product
# from collections import OrderedDict
import matplotlib.pyplot as plt
# import networkx as nx
from .numeric import rebin
# from .utils import multi_map
import plotly.graph_objects as go


def examine(self):
    from .utils import showRoisOnly
    outputStyle = {
        "color":"navy",
        "font-family":"Courier New",
        "font-size":"80%",
        }
    infoStyle = {
        "font-size":"80%",
        "color":"grey",
        }
    from jupyter_plotly_dash import JupyterDash
    from dash.dependencies import Input, Output, State
    import dash_core_components as dcc
    import dash_html_components as html
    import plotly.graph_objects as go
    from dash import no_update
    roisImage = showRoisOnly(self,indices=self.df.index, im=self.statImages[self.mode])
    roisImage.update_layout({"dragmode":'lasso'},)
    SelectedRois = html.Div([
        html.Div([
               "Selected ROIs:",
                dcc.Input(id="selected-rois",
                    type="text",
                    debounce=False,
                    size=6,
                    value="",
                 ),
            ],style={"display": "inline-box"}),

            html.Button('Discard unselected', id='discard-button', style={"display":"inline-box"},
                        n_clicks=1),
            html.Div(id="discard-feedback",children="",
                 style={"display":"inline-box",**outputStyle,}
                )])

    FilterBox = html.Div([
        html.Div([
           "Filter Timescale",
            dcc.Input(id="filter-input",
                type="str",
                debounce=True,
                size=6,
                value="30",
             ),
            html.Button("Filter",id="filter-button",n_clicks=1),
            html.Div("[it will accept also 'trace' or 'detrended']",style={"display": "inline-box",**infoStyle})
            ],style={"display": "inline-box"}),

        html.Div(id="filter-feedback",children="",
             style={"display":"inline-box",**outputStyle,}),
        html.Div([
           "rebin",
            dcc.Input(id="rebin-input",
#                     type="int",
                debounce=True,
                size=3,
                value=max(1,int(self.Freq/2)),
                style={"display":"inline-box"}
             ),
            html.Div([""],style={"margin-left":"20px"}),
            dcc.Checklist(id="sum-checklist",
                options = [{"label":"sum","value":"sum"}],
                value=[],)
        ]),

    ])

    APP_LAYOUT = [html.Div([
            html.Div([
                dcc.Graph(id="roi-selector",figure = roisImage),
                SelectedRois
            ],style={"max-width":"550px","max-height":"550px","border":"thin grey solid"}),
            html.Div([
                dcc.Graph(id="trace-show",),
                FilterBox
            ],style={"max-width":"550px","max-height":"800px","border":"thin grey solid"})
        ],
            style={"display":"flex", "flex-wrap":"wrap","align-items":"top","flex-shrink":3}
        ) for ks in [["roi-selector","range-pickers",
    #                   "roi-hover"
                     ]]]
    app = JupyterDash(__name__,
                      width=1000,
    #                   height=3000,
                     )
    @app.callback(
        Output("selected-rois", "value"),
        [Input("roi-selector", "selectedData")],
        )
    def showSelected(selData):
        if selData is None:
            return "all"
        ix = np.array( [p["hovertext"] for p in selData["points"]]).astype(int)
        ix = np.unique(ix)
        return ",".join(ix.astype(str))

    @app.callback(
        [Output("discard-feedback", "children"),
         Output("roi-selector","figure")],
        [Input("discard-button",   "n_clicks")],
        [State("selected-rois", "value")]
                 )
    def discard_callback(n_clicks,selected):
        if n_clicks <= 0:
            return no_update
        if selected == "":
            out = ""
        else:
            try:
                selectedIndices = np.unique(list(eval(selected)))
                nremoved = len(self.df)-len(selectedIndices)
                self.df = self.df.loc[selectedIndices]
                out = "%i rois removed."%(nremoved)
                fig = showRoisOnly(self,indices=self.df.index, im=self.statImages[self.mode])
                fig.update_layout({"dragmode":'lasso'},)
            except:
                out = "something's off."
                fig = go.Figure()
        return out, fig

    @app.callback(
        [Output("filter-feedback", "children"),
         Output("trace-show","figure")],
        [Input("selected-rois", "value"),
         Input("filter-button","n_clicks"),
         Input("rebin-input","value"),
         Input("sum-checklist","value")],
        [State("filter-input","value"),
         State("trace-show","relayoutData"),
        ]
                 )
    def filter_and_plot_callback(selected,n_clicks,nRebin,checklist,filtTs,rlod):
        if n_clicks <= 0:
            return no_update
        if selected in ["all",""]:
            selectedIndices = list(self.df.index)
        else:
            selectedIndices = list(eval(selected+","))
        nRebin=int(nRebin)
#             import json
#             out = json.dumps(rlod)
#             {"xaxis.range[0]": 322.98541153852864, "xaxis.range[1]": 528.4209918884629, "yaxis.range[0]": -9.125982292630475, "yaxis.range[1]": 34.41639176490383}
        out = ""#+" | ".join([str(nRebin),str(type(nRebin)),str(nRebin>1)])
#             out = "start with "+str(selected)
#             out += "| "+str(selectedIndices)
        try:
            filtTs = float(filtTs)
            k = "faster_%g"%filtTs
            if k not in self.df:
                self.fast_filter_traces(filtTs,Npoints=np.inf)
                out += "filtering for features shorter than ~%gs done."%filtTs
            y = np.sum([self.df.loc[ix,k]*self.df.loc[ix,"size"] for ix in selectedIndices],0)/self.df.loc[selectedIndices,"size"].sum()
        except:
            y = np.sum(self.df.loc[selectedIndices,filtTs],0)
        try:
            t = self.showTime["%g"%filtTs]
        except:
            t = self.time
        if nRebin>1:
            try:
                y = rebin(y,nRebin)
                t = rebin(t,nRebin)
            except:
                out+= "rebinning didnt work"
        from plotly.subplots import make_subplots
        try:
            treatments = self.protocol["compound"].unique()
        except:
            treatments = []
        h,w = 370+20*len(treatments),600
        margin = dict(zip("tblr",[30,10,20,20]))
        fg = make_subplots(rows = 2, shared_xaxes=True, row_heights = [1-.05*len(treatments)-.01,.05*len(treatments)],
                           vertical_spacing=0,
                           start_cell="bottom-left"

                          )
        fg.add_trace(go.Scatter(
            x=t,
            y=y,
            line_width=.7,
            line_color="darkred"),row=1,col=1
                     )
        it = 0
        treatments = [tr for tr in treatments if tr[:3]!="glu"] + [tr for tr in treatments if tr[:3]=="glu"]
        for treat in treatments:
            fg.add_trace(go.Scatter(
                x = [self.protocol.t_begin.min()],
                y = [0.4+it],
                mode="text",
                text=[treat[:3]+" "],
                textposition="middle left",
                showlegend=False,
                        ),row=2,col=1)
            for _,row in self.protocol.query(f"compound=='{treat}'").iterrows():
                t1,t2 = row.t_begin,row.t_end
                fg.add_trace(go.Scatter(
                    x = [t1,t2,t2,t1,t1],
                    y = np.array([0,0,1,1,0])*.8+it,
                    mode="lines",
                    line_color="grey",
                    showlegend=False,
                    fill="toself",
                    opacity = .4),row=2,col=1)
                fg.add_trace(go.Scatter(
                    x = [t1],
                    y = [0.4+it],
                    mode="text",
                    text=[" "+row.concentration],
                    textposition="middle right",
                    showlegend=False,
                            ),row=2,col=1)
                fg.update_yaxes({"tickvals":[]},row=2,col=1)

            it += 1
        fg.update_layout({
            "height":h+margin["t"]+margin["b"],
            "margin":margin,
            "width":w,
            "xaxis":{"title":"time [s]"},
            "plot_bgcolor":"white",
            "showlegend":False
        })

        fg.update_xaxes(showline=True, linewidth=1, linecolor='black',mirror=True,
                        ticks="outside", ticklen=2, row=1, col=1,
                        showgrid=False)
        fg.update_yaxes(showline=True, linewidth=1, linecolor='black',mirror=True,
                        ticks="outside", ticklen=2, row=1, col=1,
                        showgrid=False)
        try:
            fg.update_xaxes(range=[rlod["xaxis.range[0]"], rlod["xaxis.range[1]"]])
            fg.update_yaxes(range=[rlod["yaxis.range[0]"], rlod["yaxis.range[1]"]])
        except:
            pass
        return out, fg


    app.layout = html.Div(children=APP_LAYOUT,
                          style={"family":"Arial"}
                         )
    return app

