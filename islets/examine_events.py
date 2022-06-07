import json
import traceback
from sys import exc_info

import numpy as np
import pandas as pd
import plotly.express.colors as plc
from scipy.signal import peak_widths, find_peaks
from .numeric import rebin
from .fitting import generic_function

MYCOLORS = plc.qualitative.Plotly

def examine_events(self, spikeDF, x, y,
    debug=False,
    mode="jupyter",
    name=None,
    fitDF=None,
    show_diagonal=False,
    maxpoints = 1000
    ):
    if name is None:
        name = __name__
        
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

    if hasattr(self,"hb"):
        del self.hb
    from .plotly_hexbin import HBclass
    self.hb = HBclass(spikeDF,x=x,y=y)
    scatterImage = self.hb._hexbin(gridsize=50,return_figure=True)
    
    if fitDF is None:
        if "fshape" in spikeDF.columns:
            fitDF = spikeDF
    if show_diagonal:
        maxv = np.nanmax(self.hb.df.iloc[:,:2].values)
        minv = np.nanmin(self.hb.df.iloc[:,:2].values)
        scatterImage.add_trace(go.Scatter(x=[minv,maxv], y=[minv,maxv], mode="lines",line_color="grey", line_width=.7,showlegend=False,name="y=x"))

    if mode=="jupyter":
        app = Dash(name,
                  width=2000,
                  height=600,
                 )
    else:
        app = Dash(name)


    app.layout = html.Div([
        html.Div([
            dcc.Graph(id="scatterplot", figure = scatterImage),
            html.Div(id="selected-indices-feedback",style = {"width": "%ipx"%scatterImage.layout.width}),
            dcc.Input(id="selected-indices",
                     style={"border": "thin navy solid",
                            "width":"99%",
                            "height":"40px",
                            "display":"block" if debug else "none"
            })
        ], style={"display":"inline-block"}),
        html.Div([
            html.Div([
                dcc.Slider(id="tracer-slider",
                           min=0,
                           max=1,
                           value=None,
                           marks={i:{"label":""} for i in [0,1]},
                           updatemode="drag"
                          ),
            ], style={"width":'600px'}),
            html.Div([
                html.Button('<',id="previous-button", n_clicks=None, style={"display":"inline-block"}),
                html.Pre("", style={"display":"inline-block","margin-left":"5px","margin-right":"5px"},id="roi-indicator"),
                html.Button('>',id="next-button", n_clicks=0, style={"display":"inline-block"}),
            ],style={"margin-left":"150px"}),
            html.Div(id="trace-output",style = {"width": "1000px"})
            ])
    ],
        style={"family":"Arial","display":"flex"}
    )
    
            
    @app.callback(
#         Output("selected", "children"),
        [Output("selected-indices", "value"),
         Output("selected-indices-feedback", "children")],
        [Input("scatterplot", "selectedData"),],
        )
    def showSelected(selData):
        output = None
        feedback = ""
        try:
            if selData is None:
                output = None
                feedback = "Use the graph tools (lasso, rectangle) from the upper right corner to select points to show on the right plot"
            else:
                bin_ids = [pt["pointIndex"] for pt in selData["points"]]
                # return json.dumps(bin_ids, indent=2)   
                orig_indices = list(self.hb.df[
                    self.hb.df.bin_id.isin(bin_ids)
                ].index)
                if len(orig_indices)>maxpoints:
                    show_indices = orig_indices[::(len(orig_indices)//maxpoints+1)]
                    feedback = f"Uh, there is {len(orig_indices)} here. Showing them all will probably crash your browser. I'll show you a subset of {len(show_indices)}."
                else:
                    show_indices = orig_indices
                    feedback = f"{len(show_indices)} data point(s) selected."
                output = json.dumps(show_indices)

        except:
            exc_type, exc_value, exc_traceback = exc_info()
            feedback += "\nError: "+str(exc_type)
            feedback += "\nvalue: " + str(exc_value)
            for el in traceback.format_tb(exc_traceback):
                for el_ in str(el).splitlines():
                    feedback += "\n"+el_
        return output, feedback
        
    
    @app.callback(
        [Output("tracer-slider", "max"),
         Output("tracer-slider", "marks"),
         Output("previous-button", "n_clicks"),
#          Output("roi-indicator","children")
        ],
        [Input("selected-indices", "value"),],
        )
    def setTraces(selData):
        if selData is None:
            return no_update
        if "error" in str(selData).lower():
            return no_update
        
        spIndices = eval(selData)
        df = spikeDF.loc[spIndices]
        df["abundance"] = [(df.roi==roi).sum() for roi in df.roi]
        df = df.sort_values("abundance", ascending=False)
        color_id = (df.roi.diff()!=0).cumsum().values
        marks = {i:{"label":"▮", "style":{"color":MYCOLORS[color_id[i]%len(MYCOLORS)]}} for i in range(len(df))}

        return len(df)-1, marks, len(df)-1#, "loading..."
    
        
    @app.callback(
        [Output("trace-output", "children"),
         Output("roi-indicator","children"),],
        [Input("tracer-slider", "value"),   
         State("selected-indices", "value"),],
        )
    def plotSelected(v, selData):
        out = []
        roiInd = ""
#         out += [str((selData, v))]
        try:
            if selData is not None and "error" not in str(selData).lower():
                spIndices = eval(selData)
            else:
                return out, roiInd
            df = spikeDF.loc[spIndices]
            df["abundance"] = [(df.roi==roi).sum() for roi in df.roi]
            df = df.sort_values(["abundance", "t0"], ascending=[False, True])
            row = df.iloc[v]
            roi = int(row.roi)
            roiInd = "roi: %i @ %.1fs"%(roi, row.t0)
            fig = go.Figure()
            cols2show = ["trace"]
            if "slower_%g"%row.ts in self.df.columns:
                cols2show += ["slower_%g"%row.ts]
            for col in cols2show:
                dim = self.df[col][roi].shape[0]
                for kt in [""] + list(self.showTime.keys()):
                    time = self.showTime.get(kt, self.time)
                    if dim == len(time):
                        break 
                i0 = np.searchsorted(time, row.t0 - row.halfwidth*2. - .5)
                ie = np.searchsorted(time, row.t0 + row.halfwidth*4. + .5)
                t = time[i0:ie]
                x = self.df[col][roi][i0:ie]
                mode = "lines+markers" if (ie-i0<50 and col=="trace") else "lines"
                name=col
                if "slow" in col:
                    mode="lines"
                    nr=0
                else:
                    nr = int(row.halfwidth*self.Freq/50)
                    if nr>1:
                        t = rebin(t, nr)
                        x = rebin(x, nr)
                        name += f" (@ %.3gHz)"%(self.Freq/nr)
                fig.add_trace(go.Scatter(
                    x=t,
                    y=x,
                    mode=mode,
                    marker_size=4,
                    line=dict(width=.7 if col=="trace" else 1,
                              # color=color
                             ),
                    name=name
                ))
                if col=="trace":
                    dh = (x.max()-x.min())*1.2
                    hmin = x.min()-dh/1.2*.1
                    fig.add_trace(go.Scatter(
                        # x=[row.t0]*2,
                        # y=[x.min(), x.max()],
                        x=np.array([0,0,1,1,0])*row.halfwidth+row.t0,
                        y=np.array([0,1,1,0,0])*dh+hmin,
                        mode="lines",
                        line_width=0,
                        opacity=.3,
                        fill="toself",
                        showlegend=False,
                        name=None,
                    ))
                    fig.add_trace(go.Scatter(
                        x=[row.t0],
                        y=[dh+hmin],
                        mode="text",
                        text=["hw0 = %.3gs"%row.halfwidth],
                        textposition="top right",
                        showlegend=False,
                        name=None,
                    ))

            if "time" in row and "trace" in row:
                try:
                    fig.add_trace(go.Scatter(
                        x = row.time,
                        y = row.trace,
                        mode = "lines",
                        line = dict(width = 1),
                        name = "trace smoothed"
                    ))
                except:
                    pass
            index = int(row.name)
            for dummy in [None]:
                if fitDF is None: break
                if index not in fitDF.index: break
                fitrow = fitDF.loc[index]
                if not isinstance(fitrow.fshape,str): break
                if "time" in fitDF.columns and isinstance(fitrow.time, np.ndarray):
                    t = fitrow.time
                    xx = fitrow.yf
                    ll, rl = fitrow.edges
                else:
                    x_fit, bkg = generic_function(fitrow.fshape, fitrow.pars, t, fillnan=False,bkgsum=False)
                    xx = np.ones_like(x_fit)*np.nan
                    # largish = np.abs(np.diff(x_fit))>x_fit.max()/100
                    largish = np.abs(x_fit)>x_fit.max()/100
                    if not largish.any():
                        out += [str(np.round(x_fit,5))]
                        assert False
                    ifit0, ifit_end = np.where(largish)[0][[0,-1]]
                    ifit_end += 2
                    ifit0 = max(0, ifit0-1)
                    xx[ifit0:ifit_end] = x_fit[ifit0:ifit_end]
                    _, w, ll, rl = np.squeeze(peak_widths(xx,find_peaks(xx)[0]))*(t[1]-t[0])
                    ll += t[0]
                    rl += t[0]
                    xx += bkg
                    x_fit += bkg
                name = fitrow.fshape
                if "chi2" in fitrow.index:
                    name +=  " (chi2=%.2g)"%fitrow.chi2
                fig.add_trace(go.Scatter(
                    x=t,
                    y=xx,
                    mode="lines",
                    line=dict(dash="dot",width=1,color="black"),
                    name=name,
                    showlegend=True))
                    # fig.add_trace(go.Scatter(
                    #     x=t[[ifit0,ifit_end-1]],
                    #     y=x_fit[[ifit0,ifit_end-1]],
                    #     mode="lines",
                    #     line=dict(width=3,color="red",),
                    #     opacity = .3,
                    #     showlegend=False
                    # ))
                    # fig.add_trace(go.Scatter(
                    #     t[[ifit0,ifit_end-1]],
                    #     x_fit[[ifit0,ifit_end-1]],
                    #     mode="lines",
                    #     line=dict(width=.5,color="red"),
                    #     fill="toself"
                    # ))
                for l in [ll,rl]:
                    fig.add_trace(go.Scatter(
                        x=[l]*2,
                        y=np.array([0,1])*dh+hmin,
                        mode="lines",
                        line_width=1,
                        line_color="red",
                        showlegend=False,
                    ))
                try:
                    txt = "hwf = %.3gs"%(rl-ll)
                except:
                    txt = "error"
                fig.add_trace(go.Scatter(
                    x=[ll],
                    y=[hmin],
                    mode="text",
                    text=[txt],
                    textposition="bottom right",
                    showlegend=False,
                    name=None,
                ))
                #     fig.add_trace(go.Scatter(
                #         x=t[0]+np.array([ll,ll,rl,rl,ll]),
                #         y=np.array([0,1,1,0,0])*dh+hmin,
                #         mode="lines",
                #         line_width=0,
                #         opacity=.3,
                #         fill="toself",
                #         showlegend=False,
                #     ))
                #
                # except:
                #     pass
    
            fig.update_layout(
                xaxis_title='time [s]',
                margin=dict(l=10, r=10, t=30),
                # width=500
            )
            out += [html.Div(dcc.Graph(figure=fig), style= {'width':"700px","display":"inline-block"})]
            if debug:
                row = row.copy()
                for k,v in list(row.items()):
                    if  hasattr(v, '__iter__') and len(v) > 3:
                        del row[k]
                out += [html.Pre(repr(row), style={ "width":"200px", "overflowX":"scroll", "overflowY":"scroll","display":"inline-block", "border":"thin grey solid"})]

        except:
            exc_type, exc_value, exc_traceback = exc_info()
            out += [html.Br(), str(exc_type)]
            out += [html.Br(), str(exc_value)]
            for el in traceback.format_tb(exc_traceback):
                for el_ in str(el).splitlines():
                    out += [html.Br(), el_]
#             out += [str((exc_type, exc_value)), repr()]
        return out,roiInd
    
    
        
    @app.callback(
         Output("tracer-slider", "value"),
        [Input("previous-button", "n_clicks"),
         Input("next-button", "n_clicks"),],
        [State("tracer-slider", "value"),
         State("tracer-slider", "max"),]
        )
    def clickPrevious(click_previous, click_next, sliderVal, sliderMax):
        from dash import callback_context as ctx
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
        if sliderVal is None:
            sliderVal = 0
        elif "previous" in button_id:
            sliderVal = max(0,sliderVal-1)
        else:
            sliderVal = min(sliderMax,sliderVal+1)
        return sliderVal
        
    return app

