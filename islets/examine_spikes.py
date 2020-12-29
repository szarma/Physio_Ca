import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .numeric import rebin
import plotly.graph_objects as go
from .utils import getFigure
import dash
import json
from .utils import saveRois
from sys import exc_info


def examine_spikes(self, timescale,
                debug=False, 
                startShow="all",
                mode="jupyter",
                name=None,
           ):
    if name is None:
        name = __name__
    if type(startShow)!=str:
        startShow = str(list(startShow))
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
    import dash_core_components as dcc
    import dash_html_components as html
    import plotly.graph_objects as go
    from dash import no_update
    
    time = self.showTime.get("%g"%timescale, self.time)    
    scatterImage = self.show_scatter_spikes(timescale)
    spikeDF = self.spikes["%g"%timescale]
    scatterImage.update_layout({"dragmode":'lasso'},)
    
    
    
    if mode=="jupyter":
        app = Dash(name,
                  width=1200,
                  height=700,
                 )
    else:
        app = Dash(name)


    app.layout = html.Div([
        html.Div([
            dcc.Graph(id="scatterplot", figure = scatterImage),
            html.Pre(id="selected")
        ], style={"display":"inline-block"}),
        html.Div(id="tracer-div")
    ],
        style={"family":"Arial","display":"flex"}
    )
    
            
#     @app.callback(
#         Output("selected", "children"),
#         [Input("scatterplot", "selectedData"),],
#         )
#     def showSelected(selData):
#         if selData is None:
#             return "nothing yet"
#         else:
#             return json.dumps([pt["pointIndex"] for pt in selData["points"]], indent=2)   
        
    @app.callback(
        Output("tracer-div", "children"),
        [Input("scatterplot", "selectedData"),],
        )
    def setTraces(selData):
        if selData is None:
            return "nothing yet"
        spIndices = [pt["pointIndex"] for pt in selData["points"]]
        df = spikeDF.loc[spIndices]
        return html.Div([
            dcc.Slider(id="tracer-slider", min=0, max=len(df)-1,value=0,
               marks={i:{'label':""} for i in range(len(df))},
                       updatemode="drag"
                      ),
            html.Div(id="tracer-output")
        ], style={"width":'600px'})
    
        
    @app.callback(
        Output("tracer-output", "children"),
        [Input("tracer-slider", "value"),],
        [State("scatterplot", "selectedData")]
        )
    def plotSelected(v, selData):
        out = [""]
        try:
            spIndices = [pt["pointIndex"] for pt in selData["points"]]
            df = spikeDF.loc[spIndices]
            df["abundance"] = [(df.roi==roi).sum() for roi in df.roi]
            df = df.sort_values("abundance", ascending=False)
            row = df.iloc[v]
            out += [html.Pre(f"roi: {row.roi}")]
            i0 = np.searchsorted(time, row.t0-1)
            ie = np.searchsorted(time, row.t0 + timescale*1.5)
            t = time[i0:ie]
            x = self.df.trace[row.roi][i0:ie]
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=t,y=x))
            fig.add_trace(go.Scatter(x=[row.t0]*2,y=[x.min(), x.max()]))
            out += [dcc.Graph(figure=fig)]
        except:
            out += [str(exc_info())]
        return out
    
    
    
        
    return app

