import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .numeric import rebin
import plotly.graph_objects as go
from .utils import getFigure
import dash
# from dash import callback_context as ctx


def examine(self, test=False, max_rois=10, imagemode="diff_std", debug=False, startShow="all"):
    if type(startShow)!=str:
        startShow = str(list(startShow))
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
    image = self.statImages[imagemode]
    roisImage = getFigure(300,300)
#     if test:
#         roisImage = getFigure(300,300)
#     else:
#         roisImage = showRoisOnly(self,indices=self.df.index, im=image)
#         roisImage.update_layout({"dragmode":'lasso'},)
    SelectedRois = html.Div([
            html.Div([
                   "Selected ROIs:",
                    dcc.Input(id="selected-rois",
                        type="text",
                        debounce=False,
                        size=6,
                        value=",".join(list(map(str,self.df.index[:max_rois]))),
                        ),
                ],style={"visibility":"visible" if debug else "hidden"}),
            html.Button('Discard unselected', id='discard-button', style={"display":"inline-block"},
                        n_clicks=0),
            html.Div(id="discard-feedback",children="_",
                 style={"display":"inline-block",**outputStyle,}
                ),
    ])
    if not "detrended" in self.df.columns:
        try:
            self.detrend_traces()
        except:
            pass
    
    FilterBox = html.Div([
        html.Div([
           "Choose columns",
            dcc.Dropdown(id="cols-input",
#                 debounce=False,
#                 size=6,
                value=["detrended" if "detrended" in self.df.columns else "trace"],
                 options=[{"value":c,"label":c} for c in self.df.columns if \
                          hasattr(self.df[c].iloc[0],"shape") \
                          and len(self.df[c].iloc[0].shape)   \
                         ],
                 multi=True
             ),
            ],
#             style={"display": "inline-block"}
        ),
        html.Div([
           "Filter Timescale",
            dcc.Input(id="filter-input",
                type="str",
                debounce=True,
                size=6,
                placeholder="30",
#                 value="30",
              style={"width":"20px","margin-right":"5px","margin-left":"3px"},
             ),
            html.Button("Filter",id="filter-button",n_clicks=0),
            ],style={"display": "inline-block"}),

        html.Div(id="filter-feedback",children="",
             style={"display":"inline-block",**outputStyle,}),
        html.Div([
           "rebin",
            dcc.Input(id="rebin-input",
                debounce=True,
                size=3,
                value=max(1,int(self.Freq/2)),
                style={"display":"inline-block","margin-right":"20px","margin-left":"3px"}
             ),
           "offset",
            dcc.Input(id="offset-input",
                debounce=True,
                size=3,
                value=np.round(self.df["detrended"].apply(np.std).mean()*5) if "detrended" in self.df.columns else 0 ,
                style={"display":"inline-block","margin-right":"20px","margin-left":"3px"}
             ),
            html.Div(dcc.Checklist(id="sum-checklist",
                options = [{"label":"sum","value":"sum"}],
                value=[]),style={"display":"inline-block"}),
        ]),
    ])

    APP_LAYOUT = [html.Div([
            html.Div([
                html.Div([
                       "Plot Rois:",
                        dcc.Input(id="selectspec-rois",
                            type="text",
                            debounce=True,
                            size=25,
                            value=startShow,
                            ),
                    ],),
                html.Div(id="hidden",children="0",
                     style={"display":"block",**outputStyle,"visibility":"visible" if debug else "hidden"}
                    ),
                dcc.Graph(id="roi-selector",figure = roisImage),
                SelectedRois
            ],style={"max-width":"550px","max-height":"550px",
#                      "border":"thin grey solid"
                    }),
            html.Div([
                dcc.Graph(id="trace-show",figure=getFigure(400,300)),
                FilterBox
            ],style={"max-width":"550px","max-height":"800px",
#                      "border":"thin grey solid"
                    })
        ],
            style={"display":"flex", "flex-wrap":"wrap","align-items":"top","flex-shrink":3}
        ) for ks in [["roi-selector","range-pickers",
    #                   "roi-hover"
                     ]]]
    app = JupyterDash(__name__,
                      width=1200,
                      height=700,
                     )
#     @app.callback(
#         Output("roi-selector", "figure"),
#         [Input("selectspec-rois", "value")],
#         )
#     def roi_overview_callback(selspec):
#         roisImage = showRoisOnly(self,
#                                  indices=self.df.index,
#                                  im=image
#                                 )
#         roisImage.update_layout({"dragmode":'lasso'},)
#         return roisImage
#         return no_update
        
    
#     @app.callback(
#         Output("selected-rois", "value"),
#         [Input("roi-selector", "selectedData")],
#         [State("selectspec-rois", "value")]
#         )
#     def showSelected(selData, shownData):
#         if selData is None:
#             ix = self.df.index
#         else:
#             ix = np.array( [p["hovertext"] for p in selData["points"]]).astype(int)
#             ix = np.unique(ix)
#         if np.abs(ix-self.df.index).max()==0:
#             out = "all"
#         else:
#             out = ",".join(ix.astype(str))
#         return out
# #         import json
# #         return json.dumps(selData["points"])

    @app.callback(
        Output("selected-rois", "value"),
        [Input("roi-selector", "selectedData"),],
        )
    def showSelected(selData):
        if selData is None:
            ix = self.df.index
        else:
            ix = np.array( [p["hovertext"] for p in selData["points"]]).astype(int)
        ix = np.unique(ix)
        return ",".join(ix.astype(str))

#     @app.callback(
#         Output("selected-rois", "value"),
#         [Input("roi-selector", "selectedData")],
#         )
#     def showSelected(selData):
#         if selData is None:
#             return "all"
#         ix = np.array( [p["hovertext"] for p in selData["points"]]).astype(int)
#         ix = np.unique(ix)
#         try:
#             more = figure["data"][-1].selectedpoints
#             ix = np.unique(list(ix)+list(more))
#         except:
#             pass
#         return ",".join(ix.astype(str))

    @app.callback(
        [Output("discard-feedback", "children"),
         Output("hidden", "children"),
         Output("roi-selector","figure")],
        [Input("discard-button","n_clicks"),
         Input("selectspec-rois","value"),],
        [State("selected-rois", "value"),
         State("hidden","children")]
                 )
#     def discard_callback(n_clicks,selected):
    def discard_callback(n_clicks,selspec,selected,curnc):
        curnc = int(curnc)
        if n_clicks<=curnc:
            mode="write_selected"
        else:
            mode="discard"
        if selspec=="all" or selspec=="":
            seeIndices = self.df.index
        else:
            seeIndices = list(eval(selspec))
            #seeIndices = np.where(self.df.index.isin(seeIndices))[0]
        out = ""
        if mode=="discard" and selected not in ["all",""]:
            try:
                selectedIndices = np.unique(list(eval(selected)))
                nremoved = len(self.df)-len(selectedIndices)
                self.df = self.df.loc[selectedIndices]
                out = "%i rois removed."%(nremoved)
                seeIndices = np.intersect1d(seeIndices, selectedIndices)
            except:
                out = "something's off."
                
        
        fig = showRoisOnly(self, indices=seeIndices, im=image, showall=False)
        fig.update_layout({"dragmode":'lasso'},)
#         showTraces = seeIndices[:3]
#         fig.update_traces(selectedpoints=np.where(self.df.index.isin(showTraces))[0])
        
            
        return out, str(n_clicks), fig#, ",".join(showTraces.astype(str))

    @app.callback(
        Output("trace-show","figure"),
        [Input("selected-rois", "value"),
         Input("selectspec-rois", "value"),
         Input("cols-input","value"),
         Input("rebin-input","value"),
         Input("offset-input","value"),
         Input("sum-checklist","value")],
        [State("trace-show","relayoutData")]
                 )
    def plot_callback(selected,shownRois,cols,nRebin,offset,checklist,rlod):
        from sys import exc_info
        try:
            if cols is None or cols==[]:
                return no_update
            if selected in ["all",""]:
                selectedIndices = list(self.df.index)
            else:
                selectedIndices = list(eval(selected+","))
            if shownRois in ["all",""]:
                shownIndices = list(self.df.index)
            else:
                shownIndices = list(eval(shownRois+","))
            selectedIndices = np.intersect1d(selectedIndices, shownIndices)
            nRebin=int(nRebin)
            offset=float(offset)
            toSum = bool(len(checklist))
            if not toSum:
                selectedIndices = selectedIndices[:max_rois]
            from plotly.subplots import make_subplots
            from .Regions import MYCOLORS
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
            dt0 = self.time[1]-self.time[0]
            for col in cols:
                if toSum:
                    ys = np.array([np.sum([self.df.loc[ix,col]*self.df.loc[ix,"size"] for ix in \
                                           selectedIndices],0)/self.df.loc[selectedIndices,"size"].sum()])
                    ixlbl = ["avg"]
                else:
                    ys = np.vstack(self.df.loc[selectedIndices,col])
                    ixlbl = selectedIndices
                try:    t = self.showTime[col.split("_")[-1]].copy()
                except: t = self.time.copy()
#                 t = self.time.copy()
                dt = t[1]-t[0]
                nrb = int(nRebin*dt0/dt)
                if nrb>1:
                    ys = rebin(ys,nrb,axis=1)
                    if "zScore" in col:
                        ys *= nrb**.5
                    t = rebin(t,nrb)
                for iy,y,ix in zip(np.arange(len(ys)),ys, ixlbl):
                    if toSum:
                        cl="navy"
                    else:
                        cl = MYCOLORS[ix%len(MYCOLORS)]
                    fg.add_trace(go.Scatter(x=t,y=y+iy*offset,line_width=.7,line_color=cl,hovertext=f"{col} [ix={ix}]",hoverinfo="text", ),row=1,col=1)
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
#                 fg.update_xaxes(range=[rlod["xaxis.range[0]"], rlod["xaxis.range[1]"]], row=1, col=1,)
#                 fg.update_yaxes(range=[rlod["yaxis.range[0]"], rlod["yaxis.range[1]"]], row=1, col=1,)
            except:
                pass

            return fg
        except:
            out += exc_info().__repr__().replace(",","<br>")
            fg = getFigure(300,300)
            fg.add_trace(go.Scatter(
                    x = [0],
                    y = [0],
                    mode="text",
                    text=[out],
                    textposition="top right",
                    showlegend=False,))
            return fg

#         return out, fg

    @app.callback(
        [Output("filter-feedback", "children"),
         Output("cols-input", "options")
        ],
        [Input("filter-button","n_clicks")],
        [State("filter-input","value")]
                 )
    def filter_callback(n_clicks,filtTs):
        from sys import exc_info
        out = ""
        try:
            if n_clicks <= 0:
                return no_update
            try:
                filtTs = float(filtTs)
            except:
                out += "Please input timescale in seconds."
                return out, no_update
            k = "faster_%g"%filtTs

            if k not in self.df:
                self.fast_filter_traces(filtTs, saveIntermed=test)
                out += f"filtering for features shorter than ~{filtTs}s done."
            else:
                out += f"Traces filtered at {filtTs}s already exist"

            options = [{"value":c,"label":c} for c in self.df.columns if \
                              hasattr(self.df[c].iloc[0],"shape") \
                              and len(self.df[c].iloc[0].shape)   \
                             ]
        except:
            out += exc_info().__repr()
            options = no_update
        return out, options


    app.layout = html.Div(children=APP_LAYOUT,
                          style={"family":"Arial"}
                         )
    return app

