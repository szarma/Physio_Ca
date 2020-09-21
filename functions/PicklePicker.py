from jupyter_plotly_dash import JupyterDash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go
import os
import numpy as np
from IPython.core.display import HTML

from sys import path as syspath
syspath.append(os.path.expanduser("~/srdjan_functs/"))
from sys import exc_info
from caiman import movie as cmovie
import json
import pickle

import pandas as pd
from dash import no_update

from islets.Recording import Recording
from islets.Regions1 import Regions
from islets.utils import getFigure, showRoisOnly
from islets.numeric import rebin

from collections import OrderedDict
outputStyle = {
    "color":"navy",
    "font-family":"Courier New",
    "font-size":"80%",
    }
infoStyle = {
    "font-size":"80%",
    "color":"grey",
    }
bodyStyle = {
    }

regions = None

def PicklePicker(pathToRec="",appWidth=800,debug=False):
    allRecs = sorted([os.path.join(cur,f) for cur,ds,fs in os.walk("/data") for f in fs if f.endswith(".lif") or f.endswith(".nd2")])
    allRecs = [f for f in allRecs if os.path.isfile(f.replace(os.path.split(f)[1],"."+os.path.split(f)[1]+".meta"))]
#     if pathToRec is not None and pathToRec in allRec and ser is not None:
#         initOpts = [
#             {"value":el, "label":el} \
#             for el in os.listdir(pathToRec+"_analysis") if ser in el
#                    ]
#     else:
    initOpts = [{"value":None,"label":None}]
    app = JupyterDash(__name__,
                  width=appWidth,
#                   height=1000,
                 )
    APP_LAYOUT = [
        html.Div(
            'Select the experiment',
            style={**bodyStyle, "display":"inline-box"}
        ),
        dcc.Dropdown(
            options = [{"value":path,"label":path} for path in allRecs], 
            id="filepath-dropdown", 
            value = pathToRec if pathToRec in allRecs else '',
            style={**bodyStyle,"width":"70%"}
        ),
        html.Div(
            'Select the series',
            style={**bodyStyle, "display":"inline-box"}
        ),
        dcc.Dropdown(
            options = initOpts,
            id="series-dropdown",
            value=initOpts[0]["value"],
#             options=[{"value":"","label":""}],
            style={**bodyStyle,"width":"70%"}),
        html.Div([
            html.Div([
                'Select filer size',
                dcc.RadioItems(
                    options = [{"value":None,"label":None}],
                    id="pickle-radio",
                    labelStyle={"display":"block"}
                )
                ],style={**bodyStyle,
                         "display":"inline-block",
                         "vertical-align":"text-top"
                        }
            ),
            html.Div(
                id="pickle-previews",
                children = ["None"],
                style={
                    "display":"inline-block",
                    "padding-left":"20px",
                    "vertical-align":"text-top",
                }),
        ],style={"padding-top":"20px"}),
        html.Div(id="pickle-feedback",style={**outputStyle}, children=""),
    ]
    app.layout = html.Div(children=APP_LAYOUT,
                      style={"family":"Arial"}
                     )
    # callbacks
    @app.callback(
        [Output("series-dropdown","value"),
         Output("series-dropdown","options"),
#          Output("pickle-feedback","children"),
        ],
        [Input("filepath-dropdown","value"),]
    )
    def serve_series(pathToRec_):
        try:
            if pathToRec is None:
                return no_update
            if len(pathToRec)==0:
                return no_update
            analysisFolder = pathToRec+"_analysis"
            opts = [{"value":os.path.join(analysisFolder,el), "label":el} for el in sorted(os.listdir(analysisFolder))[::-1] if os.path.isdir(os.path.join(analysisFolder,el)) and el[0]!="."]
            pathToSer = (opts[0]["value"] if len(opts) else None)
            return pathToSer, opts#, ""
        except:
            return no_update#, no_update, str(exc_info())
        
        
    @app.callback(
        [Output("pickle-radio","value"),
         Output("pickle-radio","options"),
         Output("pickle-previews","children")
        ],
        [Input("series-dropdown","value")]
    )
    def serve_pickles(pathToSer,width=300, height=300):
        import base64
        from collections import OrderedDict 
        if pathToSer is None:
            return no_update
        if len(pathToSer)==0:
            return no_update
        preview = OrderedDict()
        options = []
        for f in sorted(os.listdir(pathToSer))[::-1]:
            if f.endswith("pkl"):
                k = f.split("_")[0].replace(".","+")
                pathToRoi = os.path.join(pathToSer,f)
                options += [{ "label": k, "value": pathToRoi }]
                previewFile = pathToRoi.replace(f,".image_"+f.replace("_rois.pkl",".png"))
                if os.path.isfile(previewFile):
                    encoded_image = base64.b64encode(open(previewFile, 'rb').read())
                    preview[k] = html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode()), width="%ipx"%(width*.9),height="%ipx"%(height*.8+20))

        if len(options):
            val = options[0]["value"]
        else:
            val = None
        previews = [html.Div(children=[html.Div(k,style={"padding-left":"20px"}), preview[k]], style=dict(
            height="%ipx"%height,
            width="%ipx"%width,
            display="inline-block",
            border='thin lightgrey solid' 
        )) for k in preview]
        return val, options, html.Div(previews,style={"width":"%ipx"%(width*2.1)})



    @app.callback(
        Output("pickle-feedback","children"),
        [Input("pickle-radio","value")]
    )
    def import_pickle(path):
        try:
            if debug:
                globpath = path
            if path is None:
                return no_update
            if len(path)==0:
                return no_update
            global regions
#             from islets.numeric import mydebleach
            from islets.utils import multi_map
            with open(path,"rb") as f:
                regions = pickle.load(f)
            regions.update()
            regions.detrend_traces()
            regions.infer_gain()
            pickleDir = os.path.split(path)[0]
            try:
                protocolFile = os.path.join(pickleDir, [f for f in os.listdir(pickleDir) if "protocol" in f][0])
                regions.import_protocol(protocolFile)
            except:
                pass
            feedback = [
                "Regions imported successfully.",
                "Original movie:",
                html.Div(children=[
                    "dimension (T,X,Y): (%i, %i, %i)"%((len(regions.time), )+regions.image.shape),
                    html.Br(),
                    "duration (h:m:s): "+str(pd.Timedelta(regions.time.max(), unit="s")).split()[-1].split(".")[0],
                    html.Br(),
                    "frequency (Hz): %.3g"%regions.Freq,
                ],style={'padding-left': '30px'}),
                "Number of ROIs: %i"%len(regions.df),
            ]
        #     feedback = sum([[el, html.Br()] for el in feedback],[])
        #     roisImage = showRoisOnly(regions,indices=regions.df.index, im=regions.statImages[regions.mode])
        #     return feedback,roisImage,1
        except:
            feedback = str(exc_info())
        return feedback
    
    return app