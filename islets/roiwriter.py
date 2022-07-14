import json
from sys import exc_info
import numpy as np
from .numeric import rebin
from .utils import getFigure
from .utils import saveRois
from pandas import Series


def specify_rois(self,
                 movie=None,
                 debug=False,
                 mode="jupyter",
                 name=None,
                 lw=None,
                 fill=False,
                 **imkwargs
                 ):
    if name is None:
        name = __name__
    from .utils import showRoisOnly
    outputStyle = {
        "color":"navy",
        "font-family":"Courier New",
        "font-size":"80%",
        }
    preStyle = {
        "font-family":"Courier New",
        "font-size":"80%",
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll',
        'overflowY': 'scroll',
        }

    if mode=="jupyter":
        from jupyter_plotly_dash import JupyterDash as Dash
    if mode=="dash":
        from dash import Dash
    if mode=="jupyter-dash":
        from jupyter_dash import JupyterDash as Dash
    from dash.dependencies import Input, Output, State
    from dash import dcc
    from dash import html
    import plotly.graph_objects as go
    from dash import no_update
    
    def get_pixels(theshape):
        if theshape["type"] == "circle":
            x0 = min(theshape["x0"], theshape["x1"])
            x1 = max(theshape["x0"], theshape["x1"])
            y0 = self.image.shape[0] - min(theshape["y0"], theshape["y1"])
            y1 = self.image.shape[0] - max(theshape["y0"], theshape["y1"])
            x_center = (x0 + x1) / 2
            y_center = (y0 + y1) / 2
            radius_x = (x1 - x0) / 2
            radius_y = (y1 - y0) / 2
            pixels = []
            for x in range(int(x0), int(x1)+1):
                for y in range(int(y1)-1, int(y0)):
                    if (x-x_center)**2/radius_x**2 + (y-y_center)**2/radius_y**2 < 1 :
                        pixels += [(y,x)]
        if theshape["type"] == "rect":
            x0 = min(theshape["x0"], theshape["x1"])
            x1 = max(theshape["x0"], theshape["x1"])
            y0 = min(theshape["y0"], theshape["y1"])
            y1 = max(theshape["y0"], theshape["y1"])
            x0 = int(np.ceil(x0))
            x1 = int(np.ceil(x1))
            y0 = self.image.shape[0] - int(np.ceil(y0))
            y1 = self.image.shape[0] - int(np.ceil(y1))
            pixels = [(y,x) for x in range(x0,x1) for y in range(y1,y0)]
        return pixels
        
    roisImage = go.Figure()#showRoisOnly(self,indices=self.df.index, im=self.statImages[imagemode], lw=lw)
    roisImage.add_trace(
        go.Scatter(
            x = [0],
            y = [0],
            line_width = 0,
            name = ""
        ),
    )
    from .utils import createStaticImage
    im = self.image
    imgpointer = createStaticImage(None, im[::-1], **imkwargs)
    roisImage.add_layout_image(
        dict(
            source = imgpointer,
            xref = "x",
            yref = "y",
            x = -.5,
            y = im.shape[0]-.5,
            sizex = im.shape[1],
            sizey = im.shape[0],
            sizing = "stretch",
            opacity = 1,
            layer = "below")
    )
    h, w = 700, 700 * im.shape[1] / im.shape[0]
    if w > 500:
        h = 500 / w * h
        w = 500
    h += 70
    w += 20
    roisImage.update_layout({
        # "title":regions.mode+" (filtered)",
        "height": h,
        "width": w,
        "margin": dict(l = 10, r = 10, t = 50, b = 20),
        "xaxis": {
            "zeroline": False,
            "showgrid": False,
            "linecolor": 'black',
            "linewidth": 1,
            "mirror": True,
            "tickvals": [],
            "range": [-.5, im.shape[1] - .5]
        },
        "yaxis": {
            "zeroline": False,
            "showgrid": False,
            "linecolor": 'black',
            "linewidth": 1,
            "mirror": True,
            "tickvals": [],
            "range": [-.5, im.shape[0] - .5,],
        },
        # 'clickmode': 'event+select',
        "dragmode": 'drawcircle',
        "newshape": dict(line_width = 2)
    })
    # roisImage.layout.update(modebar_add = "select2d")

    # for mode in ["drawline", "drawopenpath", "drawclosedpath", "drawcircle", "drawrect", "eraseshape"]:
    #     roisImage.update_layout(modebar_add = mode)
    roisImage.update_yaxes(
        scaleanchor = "x",
        scaleratio = 1,
    )
    roisImage.update_layout(
        modebar_add = [
            # 'drawline',
            # 'drawopenpath',
            # 'drawclosedpath',
            'drawcircle',
            'drawellipse',
            'drawrect',
            'eraseshape',
            # "lasso",
            # "select2d"
        ]
    )


    FilterBox = html.Div([
        html.Div([
           "Choose columns",
            dcc.Dropdown(id="cols-input",
                # value=["detrended" if "detrended" in self.df.columns else "trace"],
                value=["trace"],
                 options=[{"value":c,"label":c} for c in ["trace"]],
                 multi=True
             ),
            ],
        ),
        html.Div([
           "Filter Timescale",
            dcc.Input(id="filter-input",
                type="str",
                debounce=True,
                size=10,
                placeholder="30",
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
                value=max(1,int(np.round(self.__dict__.get("Freq",1)/2))),
                style={"display":"inline-block","margin-right":"20px","margin-left":"3px"}
             ),
           "offset",
            dcc.Input(id="offset-input",
                debounce=True,
                size=3,
                value=np.round(self.df["detrended"].apply(np.std).mean()*5) if "detrended" in self.df.columns else 0 ,
                style={"display":"inline-block","margin-right":"20px","margin-left":"3px"}
             ),
            html.Div([
                dcc.Checklist(id="sum-checklist",
                    options = [{"label":"sum","value":"sum"},{"label":"remove median","value":"remove"}],
                    value=["remove"]),
            ],style={"display":"inline-block"}),
        ]),
    ])

    APP_LAYOUT = [
        html.Div([
            html.Div([
                dcc.Graph(id="roi-selector",figure = roisImage),
                html.Pre(id="selected-catcher",children="-"*10,
                     style={"width":"200px","height":"400px","display":"inline-block" if debug else "none",**preStyle}
                    ),
                html.Pre(id="relayout-catcher",children="-"*10,
                     style={"width":"200px","height":"400px","display":"inline-block" if debug else "none",**preStyle}
                    ),
                html.Pre(id="pixels-choice",children="-"*10,
                     style={"width":"200px","height":"400px","display":"inline-block" if debug else "none",**preStyle},
                    ),
            ],style={"max-width":"800px","max-height":"800px","align-items":"top"
#                      "border":"thin grey solid"
                    }),
            html.Div([
                dcc.Graph(id="trace-show",figure=getFigure(400,300)),
                FilterBox,
                # html.Details(title="Movie Closeup",children=movieCloseup, open=debug)
                html.Br(),
            ],style={"max-width":"550px","max-height":"800px",
#                      "border":"thin grey solid"
                    })
        ],
            style={"display":"flex", "flex-wrap":"wrap","align-items":"top","flex-shrink":3}
        ) for ks in [["roi-selector","range-pickers",
    #                   "roi-hover"
                     ]
                     
                    ]
    ]
    if mode=="jupyter":
        app = Dash(name,
                  width=1500,
                  height=1000,
                 )
    else:
        app = Dash(name)


    app.layout = html.Div(
        children=APP_LAYOUT,
        style={"family":"Arial"}
        )


    @app.callback(
        Output("selected-catcher","children"),
        # [Input("roi-selector", "relayoutData"),]
        [Input("roi-selector", "selectedData"),]
                 )
    def catch_sel(rlout):
        return "selectedData:\n"+json.dumps(rlout, indent = 2)

#     @app.callback(
#         Output("pixels-choice","children"),
#         [Input("roi-selector", "relayoutData"),
#          Input("roi-selector", "selectedData"),]
#                  )
#     def catch_pixels(rlout, selData):
#         output = ""
#         try:
#             if "shapes" in rlout:
#                 shapes = rlout["shapes"]
#                 theshape = shapes[-1]
#                 pixels = get_pixels(theshape)
#                 whereIsPeak = np.argmax([im[px] for px in pixels])
#                 peak = pixels[whereIsPeak]
#                 trace = np.mean([ movie[:,v,h] for v,h in pixels ], axis = 0 )
#                 if "tmp" in self.df.index:
#                     self.df.drop(index=['tmp'], inplace=True)
#                 self.df.loc["tmp"] = Series({
#                     "peak":peak,
#                     "pixels":pixels,
#                     "trace":trace,
#                     "size":len(pixels)
#                 })
#                 if not hasattr(self,"time"):
#                     self.time = np.arange(len(movie))/movie.fr
#                     self.movie = movie
#                     self.Freq = movie.fr
#                 output += json.dumps({
#                     "pixels":pixels,
#                     # "center": (x_center, y_center),
#                     # "radii": (radius_x, radius_y)
#                 })
#         except:
#             output += repr(exc_info()).replace(",","<br>")
#         return output

    @app.callback(
        Output("relayout-catcher","children"),
        [Input("roi-selector", "relayoutData"),]
                 )
    def catch_rlo(rlout):
        return "relayoutData:\n"+json.dumps(rlout, indent = 2)



    @app.callback(
        Output("trace-show","figure"),
        [Input("cols-input","value"),
         Input("rebin-input","value"),
         Input("roi-selector", "relayoutData")],
        [State("trace-show","relayoutData")]
                 )
    def plot_callback(cols,nRebin,roi_rlod,rlod):
        out = ""
        try:
            if "shapes" in roi_rlod:
                shapes = roi_rlod["shapes"]
                theshape = shapes[-1]
                pixels = get_pixels(theshape)
                whereIsPeak = np.argmax([im[px] for px in pixels])
                peak = pixels[whereIsPeak]
                if not hasattr(self,"time"):
                    self.time = np.arange(len(movie))/movie.fr
                    self.movie = movie
                    self.Freq = movie.fr
                trace = np.mean([ self.movie[:,v,h] for v,h in pixels ], axis = 0 )
                if "tmp" in self.df.index:
                    self.df.drop(index=['tmp'], inplace=True)
                self.df.loc["tmp"] = Series({
                    "peak":peak,
                    "pixels":pixels,
                    "trace":trace,
                    "size":len(pixels)
                })

            if cols is None or cols==[]:
                return no_update
#             if selected in ["all",""]:
#                 selectedIndices = list(self.df.index)
#             else:
#                 selectedIndices = list(eval(selected+","))
#             if shownRois in ["all",""]:
#                 shownIndices = list(self.df.index)
#             else:
#                 shownIndices = list(eval(shownRois+","))
            selectedIndices = ["tmp"]#np.intersect1d(selectedIndices, shownIndices)
            nRebin=int(nRebin)

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
                ys = np.vstack(self.df.loc[selectedIndices,col])
                ixlbl = selectedIndices
                try:    t = self.showTime[col.split("_")[-1]].copy()
                except: t = self.time.copy()
                dt = t[1]-t[0]
                nrb = int(nRebin*dt0/dt)
                if nrb>1:
                    ys = rebin(ys,nrb,axis=1)
                    if "zScore" in col:
                        ys *= nrb**.5
                    t = rebin(t,nrb)
                for iy,y,ix in zip(np.arange(len(ys)),ys, ixlbl):
                    #cl="black"
                    fg.add_trace(
                        go.Scatter(
                            x=t,
                            y=y,
                            line_width=.7,
                            #line_color=cl,
                            hovertemplate = "t = %{x}",
                            name=f"{col}[{ix}]"
                        ),
                        row=1,
                        col=1
                    )
            it = 0
            treatments = [tr for tr in treatments if tr[:3].lower()!="glu"] + [tr for tr in treatments if tr.lower()[:3]=="glu"]
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
            # fg.update_layout(hovermode = 'x')
            fg.update_xaxes(showspikes = True,
                            spikethickness = 1,
                            spikedash = "solid",
                            # spikecolor = "green",
                            spikesnap = "cursor",
                            spikemode = "across"
                            )
            try:
                fg.update_xaxes(range=[rlod["xaxis.range[0]"], rlod["xaxis.range[1]"]])
                fg.update_yaxes(range=[rlod["yaxis.range[0]"], rlod["yaxis.range[1]"]])
#                 fg.update_xaxes(range=[rlod["xaxis.range[0]"], rlod["xaxis.range[1]"]], row=1, col=1,)
#                 fg.update_yaxes(range=[rlod["yaxis.range[0]"], rlod["yaxis.range[1]"]], row=1, col=1,)
            except:
                pass

            return fg
        except:
            out += str(exc_info()).replace(",","<br>")
            fg = getFigure(300,300)
            fg.add_trace(go.Scatter(
                    x = [0],
                    y = [0],
                    mode="text",
                    text=[out],
                    textposition="top right",
                    showlegend=False,))
            return fg


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
                self.fast_filter_traces(filtTs,)
                # if filtTs<=10:
                #     self.fast_filter_traces(filtTs,Npoints=np.inf,)
                # else:
                #     self.fast_filter_traces(filtTs,Npoints=np.inf,filt_cutoff=0)
                out += f"filtering for features shorter than ~{filtTs}s done."
            else:
                out += f"Traces filtered at {filtTs}s already exist"

            options = [{"value":c,"label":c} for c in self.df.columns if \
                              hasattr(self.df[c].iloc[0],"shape") \
                              and len(self.df[c].iloc[0].shape)   \
                             ]
        except:
            out += repr(exc_info())
            options = no_update
        return out, options


    return app

