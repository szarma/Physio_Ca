import json
import traceback
import numpy as np
from .numeric import rebin
from .utils import getFigure
from .utils import saveRois

def trcbk():
    return repr(traceback.format_tb())

def examine(self, 
            max_rois=10, 
            imagemode=None, 
            debug=False, 
            startShow="all",
            mode="jupyter",
            name=None,
            lw=None,
            fill=False
           ):
    if name is None:
        name = __name__
    if type(startShow)!=str:
        startShow = str(list(startShow))
    if imagemode is None:
        imagemode = self.mode
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
    if mode=="jupyter":
        from jupyter_plotly_dash import JupyterDash as Dash
    if mode=="dash":
        from dash import Dash
    if mode=="jupyter-dash":
        from jupyter_dash import JupyterDash as Dash
    from dash.dependencies import Input, Output, State
    try:
        from dash import dcc, html
    except ImportError:
        import dash_core_components as dcc
        import dash_html_components as html
    import plotly.graph_objects as go
    from dash import no_update
    from networkx import adjacency_matrix
    # from importlib import import_module
    # module = import_module(self.__module__)
    from .utils import getPeak2BoundaryDF, getGraph_of_ROIs_to_Merge
    # getPeak2BoundaryDF = getattr(module, "getPeak2BoundaryDF")
    # getGraph_of_ROIs_to_Merge = getattr(module, "getGraph_of_ROIs_to_Merge")
    # mergeBasedOnGraph = getattr(module, "mergeBasedOnGraph")
    
    roisImage = getFigure()#showRoisOnly(self,indices=self.df.index, im=self.statImages[imagemode], lw=lw)
    roisImage.update_layout({"dragmode":'lasso'},)
    if hasattr(self,"time"):
        if not hasattr(self,"TwoParFit"):
            self.infer_TwoParFit()
    SelectedRois = html.Div([
        html.Div([
               "Selected ROIs:",
                dcc.Input(id="selected-rois",
                    type="text",
                    debounce=False,
                    size=16,
                    value="",#.join(list(map(str,self.df.index[:max_rois]))),
                    ),
            ],style={"display":"block" if debug else "none"}),
        html.Div([
            html.Button('Discard unselected',
                        id='discard-button',
                        n_clicks=0),
            html.Button('Discard selected',
                        id='discard-button-selected',
                        n_clicks=0),],
            style={"display":"inline-block","width":"200px"},
        ),
        html.Div([
            html.Button('Mark for merging', id='mark-button', n_clicks=0,style={"display":"block"}),
            html.Button('Merge', id='merge-button', n_clicks=0),], 
            style={"display":"inline-block","width":"200px"},
        ),
        html.Pre(id="discard-feedback",children="_",
             style={"display":"block",**outputStyle,}
            ),
        html.Div([
            html.Button(id="tag-button", children="Tag cells",
                    n_clicks=0,style={"display":"inline-block"}),
            dcc.Dropdown(id="tag-drop", multi=True,value=[],options = [
                {"label":opt,"value":opt} for opt in [
                    "islet","non-islet","ductal","beta","non-beta","interesting"
                ]],style={"width":"200px",}),
        ],style={'align-items': 'center', 'display': 'flex'}),
        html.Pre(id="tag-feedback",children="_",
             style={**outputStyle}
            ),
        html.Button(id="save-button", children="Save", n_clicks=0),
        html.Pre(id="save-feedback",children="_",
             style={"display":"inline-block",**outputStyle,"padding":"5px"}
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
                # value=["detrended" if "detrended" in self.df.columns else "trace"],
                value=["trace"],
                 options=[{"value":c,"label":c} for c in self.df.columns if \
                          hasattr(self.df[c].iloc[0],"shape") \
                          and len(self.df[c].iloc[0].shape)   \
                         ],
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
    initNcs = {"discard_unsel":0,"discard_sel":0,"mark":0,"merge":0}
    
    movieCloseup = [
        # html.H3("Movie closeup", style={"display":"inline-block"}),
        html.Button("generate closeup", id="create-closeup", n_clicks=0,title="[!experimental feature!] Works only if the regions have the associated movie."),
        html.Div(id="closeup-output")
                   ]
    
    APP_LAYOUT = [
        html.Div([
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
                dcc.Graph(id="roi-selector",figure = roisImage),
                SelectedRois,
                html.Pre(id="hidden",children=json.dumps(initNcs, indent=2),
                     style={"display":"block" if debug else "none",**outputStyle}
                    )
            ],style={"max-width":"550px","max-height":"550px",
#                      "border":"thin grey solid"
                    }),
            html.Div([
                dcc.Graph(id="trace-show",figure=getFigure(400,300)),
                FilterBox,
                # html.Details(title="Movie Closeup",children=movieCloseup, open=debug)
                html.Br(),
                html.Div(movieCloseup)
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
                  width=1200,
                  height=700,
                 )
    else:
        app = Dash(name)


    app.layout = html.Div(children=APP_LAYOUT,
                          style={"family":"Arial"}
                         )
    @app.callback(
        Output("closeup-output","children"),
        [Input("create-closeup", "n_clicks"),],
        [State("selected-rois", "value")],
        )
    def create_closeup(n_clicks,selectedData):
        if debug:
            output = [f"Start n_clicks: {n_clicks}; {selectedData}"]
        else:
            output = []
        try:
            from .utils import closeup_movie
            if n_clicks>0:
                if hasattr(self,"movie"):
                    selected = eval(selectedData)
                    if hasattr(selected, '__getitem__'):
                        selected = list(selected)
                    else:
                        selected = [selected]
                    output += [f"selected={selected}"]
#                     output += [str(type(m))]
                    m = closeup_movie(self, indices = selected)
                    output += [html.Iframe(srcDoc=m._repr_html_(),
                                           width = 700, height = 800
                                           )]
                else:
                    output += ["Oops, you first need to associate a movie to the regions."]
        except:
            output+=[trcbk()]
        return output
        
    @app.callback(
        Output("tag-feedback","children"),
        [Input("tag-button", "n_clicks"),],
        [State("tag-drop", "value"),
         State("selected-rois", "value")]
        )
    def tag(n_clicks,tags,selectedData):
        try:
            if n_clicks>0:
                if "tag" not in self.df.columns:
                    self.df["tag"] = [""]*len(self.df)
                selected = list(eval(selectedData))
                for i in selected:
                    for tag in tags:
                        if tag in self.df.loc[i,"tag"]:
                            continue
                        else:
                            tagExists = bool(len(self.df.loc[i,"tag"]))
                            self.df.loc[i,"tag"] += "|"*int(tagExists)+tag
                return f"Successfully applied tags for {len(selected)} rois."
        except:
            return trcbk()
        
        
        
    @app.callback(
        Output("save-feedback","children"),
        [Input("save-button", "n_clicks"),],
        )
    def save(n_clicks):
        from os.path import split
        try:
            if n_clicks>0:
                folder, fname = split(self.pathToRois)
                return saveRois(self,folder,fname.replace("_rois.pkl",""),formats=["vienna"])
        except:
            return trcbk()
            
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


    @app.callback(
        [Output("discard-feedback", "children"),
         Output("hidden", "children"),
         Output("roi-selector","figure")],
        [Input("discard-button","n_clicks"),
         Input("discard-button-selected","n_clicks"),
         Input("mark-button","n_clicks"),
         Input("merge-button","n_clicks"),
         Input("selectspec-rois","value"),],
        [State("selected-rois", "value"),
         State("hidden","children"),
         State("roi-selector", "relayoutData")]
                 )
    def mark_discard_callback(discard_unsel_clicks,discard_sel_clicks,mark_clicks,merge_clicks,selspec,selected,curnc,rlout):
        out = ["-"*40]
        fig = getFigure()
        outcurns = "nothing"
        try:
#             out += [html.Br(),"what now?"]
            d = {
                 "discard_unsel": discard_unsel_clicks,
                 "discard_sel":  discard_sel_clicks,
                 "mark":         mark_clicks,
                 "merge":        merge_clicks
                }
            outcurns = json.dumps(d, indent=2)
            curnc = json.loads(curnc)
            ##########
            if discard_unsel_clicks>curnc["discard_unsel"]:
                mode="discard_unsel"
            elif discard_sel_clicks>curnc["discard_sel"]:
                mode="discard_sel"
            elif mark_clicks>curnc["mark"]:
                mode="mark"
            elif merge_clicks>curnc["merge"]:
                mode="merge"
            else:
                mode="plot"
            #out += [html.Br(),f"mode = {mode}"]
            ##########
            if selspec=="all" or selspec=="":
                seeIndices = self.df.index
            else:
                seeIndices = list(eval(selspec))
            if selected in ["all",""]:
                selectedIndices = self.df.index
            else:
                try:
                    selectedIndices = np.unique(list(eval(selected)))
                except:
                    selectedIndices = [int(selected)]
            ##########
            if mode=="discard_unsel":
                nremoved = len(self.df)-len(selectedIndices)
                self.df = self.df.loc[selectedIndices]
                out += [html.Br(), "%i rois removed."%(nremoved)]
                self.update()
            if mode=="discard_sel":
                nremoved = len(selectedIndices)
                if nremoved>len(self.df)/2:
                    out += [html.Br(), "Can not remove max 50% of existing rois in one go. Sorry, this is for your own safety :-)"]
                else:
                    # out += [html.Br(), "removing: "+",".join(selectedIndices.astype(str))]
                    self.df.drop(index=selectedIndices, inplace=True)
                    out += [html.Br(), "%i rois removed."%(nremoved)]
                    self.update()
            if mode in ["discard_sel","discard_unsel","plot"]:
                seeIndices = [isee for isee in seeIndices if isee in self.df.index]
                # out += [html.Br(),seeIndices.__repr__()]
                fig = showRoisOnly(self, im=self.statImages[imagemode], showall=True, lw=lw, indices=seeIndices,fill=fill)
        
            ##########
#             if mode=="plot":
#                 fig = showRoisOnly(self, indices=self.df.index, im=self.statImages[imagemode], showall=False)
            ##########
            if mode=="mark" or (mode=="merge" and (not hasattr(self,"mergeGraph"))):
                p2b = getPeak2BoundaryDF(self.df.loc[selectedIndices], distTh=np.inf)
                p2b = p2b[p2b[["i","j"]].isin(selectedIndices).all(axis=1)]
                gph = getGraph_of_ROIs_to_Merge(p2b[["i","j"]],self)
                self.mergeGraph = gph
            if mode=="mark":
                out += [ html.Br(), "If you now merge, this is how rois will be merged"]
                fig = showRoisOnly(self, indices=selectedIndices, im=self.statImages[imagemode], showall=True, lw=lw)
                for j,i in self.mergeGraph.edges:
                    fig.add_annotation(
                      x=self.df.loc[i,"peak"][1],  # arrows' head
                      y=self.df.loc[i,"peak"][0],  # arrows' head
                      ax=self.df.loc[j,"peak"][1],  # arrows' tail
                      ay=self.df.loc[j,"peak"][0],  # arrows' tail
                      xref='x',
                      yref='y',
                      axref='x',
                      ayref='y',
                      text='',  # if you want only the arrow
                      showarrow=True,
                      arrowhead=2,
                      arrowsize=1,
                      arrowwidth=1.5,
                      arrowcolor='red'
                    )
            ##########
            if mode=="merge":
                dn = self.mergeBasedOnGraph(self.mergeGraph)
                del self.mergeGraph
                out += [ html.Br(), f"{dn} rois merged into existing roi(s)"]
                fig = showRoisOnly(self, im=self.statImages[imagemode], showall=True, lw=lw)


            try:
                fig.update_xaxes(range = [rlout["xaxis.range[0]"], rlout["xaxis.range[1]"]])
                fig.update_yaxes(range = [rlout["yaxis.range[0]"], rlout["yaxis.range[1]"]])
            except:
                pass

        except:
            out += [html.Br(),"  "+trcbk()]
            
        
        return out, outcurns, fig
    if hasattr(self,"time"):
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
            out = ""
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
                toSum = "sum" in checklist
                removeMedian = "remove" in checklist
                if not toSum:
                    if "activity" in self.df.columns:
                        selectedIndices = list(
                            self.df.loc[selectedIndices].sort_values(
                            "activity", ascending=False
                            ).index)[:max_rois]
                    else:
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
                            cl="black"
                        else:
                            try:
                                cl = self.df.loc[ix,"color"]
                            except:
                                cl = MYCOLORS[ix%len(MYCOLORS)  ]
                        if removeMedian:
                            y = y-np.median(y)
                        fg.add_trace(
                            go.Scatter(
                                x=t,
                                y=y+iy*offset,
                                line_width=.7,
                                line_color=cl,
                                hovertemplate = "t = %{x}",
                                # hovertext=f"{col} [ix={ix}]",
                                # hoverinfo="text",
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
                out += trcbk().replace(",","<br>")
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

