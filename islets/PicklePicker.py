try:
    from jupyter_plotly_dash import JupyterDash
except ModuleNotFoundError:
    from jupyter_dash import JupyterDash
try:
    from dash import dcc, html
except ImportError:
    import dash_core_components as dcc
    import dash_html_components as html

from dash.dependencies import Input, Output, State
import os

from sys import exc_info
import pandas as pd
from dash import no_update

from .Recording import Recording
from .Regions import load_regions
from .LineScan import LineScan

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

class PicklePicker:
    def __init__(self, pathToRec=None, series=None, appWidth=1500,
                 debug=False, appHeight=1200, showExamine=True, max_rois=10):
        self.allRecs = sorted([os.path.join(cur, f) for cur, ds, fs in os.walk("/data") for f in fs
                               if f.endswith(".lif") or f.endswith(".nd2") or f.endswith(".tif")])
        self.allRecs = [f for f in self.allRecs if os.path.isfile(
            f.replace(os.path.split(f)[1], "." + os.path.split(f)[1] + ".meta"))]
        self.pathToRec = pathToRec
        self.series = series
        self.appWidth = appWidth
        self.debug = debug
        self.appHeight = appHeight
        self.showExamine = showExamine
        self.max_rois = max_rois

        self.regions = None
        self.linescan = None
        self.show()

    def show(self):
        initOpts = [{"value": "", "label": ""}]
        app = JupyterDash(__name__,
                          width=self.appWidth,
                          height=self.appHeight,
                          #                   height=1000,
                          )
        APP_LAYOUT = [
            html.Div(
                'Select the experiment',
                style={**bodyStyle, "display": "inline-box"}
            ),
            dcc.Dropdown(
                options=[{"value": path, "label": path} for path in self.allRecs],
                id="filepath-dropdown",
                value=self.pathToRec if self.pathToRec in self.allRecs else '',
                style={**bodyStyle, "width": "70%"}
            ),
            html.Div(
                'Select the series',
                style={**bodyStyle, "display": "inline-box"}
            ),
            dcc.Dropdown(
                options=initOpts,
                id="series-dropdown",
                value=initOpts[0]["value"],
                style={**bodyStyle, "width": "70%"}),
            html.Div(id="series-feedback", children="-----", style=outputStyle),
            html.Div([
                html.Div([
                    html.Div(id="frase-container"),
                    dcc.RadioItems(
                        options=[{"value": None, "label": None}],
                        id="pickle-radio",
                        labelStyle={"display": "inline-block"}
                    )
                ], style={**bodyStyle,
                          "display": "inline-block",
                          "vertical-align": "text-top"
                          }
                ),
                html.Div(
                    id="pickle-previews",
                    children=["None"],
                    style={
                        "display": "inline-block",
                        "padding-left": "20px",
                        "vertical-align": "text-top",
                    }),
            ], style={"padding-top": "20px"}),
            html.Div(id="pickle-feedback", style={**outputStyle}, children=""),
            html.Iframe(id="examiner", width="99%", height="700px",
                        style={"display": "inline-block" if self.showExamine else "none"})
        ]
        app.layout = html.Div(children=APP_LAYOUT,
                              style={"family": "Arial"}
                              )

        # callbacks
        @app.callback(
            [Output("series-dropdown", "value"),
             Output("series-dropdown", "options"),
             Output("series-feedback", "children"),
             ],
            [Input("filepath-dropdown", "value"), ]
        )
        def serve_series(pathToRec_):
            feedback = []
            opts = initOpts
            val = ""
            try:
                if pathToRec_ is not None and len(pathToRec_) != 0:
                    analysisFolder = pathToRec_ + "_analysis"
                    if os.path.isdir(analysisFolder):
                        subdirs = sorted(os.listdir(analysisFolder))[::-1]
                        series_ = [el for el in subdirs if
                                   os.path.isdir(os.path.join(analysisFolder, el)) and el[0] != "."]
                        opts = [{"value": os.path.join(analysisFolder, el), "label": el} for el in series_]
                        try:
                            rec = Recording(pathToRec_)
                        #                         feedback += [str(rec.metadata), str("rec" in locals())]
                        except:
                            feedback += ["Could not load recording metadata."]
                        if "rec" in locals():
                            try:
                                rec.metadata = rec.metadata.set_index("Name")
                                for opt in opts:
                                    lbl = opt["label"]
                                    lbl = lbl.split("-")[0]
                                    if lbl not in rec.metadata.index: continue
                                    if rec.metadata.loc[lbl, "line scan"] != "none":
                                        opt["label"] += "*"
                                opts = opts[::-1]
                            except:
                                feedback += [str(rec.metadata.index), str(exc_info())]
                        #                     val = opts[0]["value"]
                        val = [opt for opt in opts if self.series in opt["label"]]
                        if len(val) >= 1:
                            val = val[0]["value"]
                        else:
                            val = opts[0]["value"]
                        # feedback += ["options: "+str(opts)]
                    else:
                        feedback += [f"directory {analysisFolder} apparently does not exist."]

            except:
                feedback += [str(exc_info())]
            return val, opts, feedback

        @app.callback(
            [
             Output("pickle-radio", "options"),
             Output("pickle-previews", "children"),
             Output("frase-container", "children")
             ],
            [Input("series-dropdown", "value")]
        )
        def serve_pickles(pathToSer, width=300, height=300):
            import base64
            from collections import OrderedDict
            if pathToSer is None:
                return no_update
            if len(pathToSer) == 0:
                return no_update
            preview = OrderedDict()
            options = []
            frase = ""
            for f in sorted(os.listdir(pathToSer))[::-1]:
                if f.endswith("pkl"):
                    k = f.split("_rois.pkl")[0].replace(".", "+")
                    pathToRoi = os.path.join(pathToSer, f)
                    options += [{"label": k, "value": pathToRoi}]
                    previewFile = pathToRoi.replace(f, ".image_" + f.replace("_rois.pkl", ".png"))
                    if os.path.isfile(previewFile):
                        encoded_image = base64.b64encode(open(previewFile, 'rb').read())
                        preview[k] = html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode()),
                                              width="%ipx" % (width * .9), height="%ipx" % (height * .8 + 20))
            if len(options) > 0:
                frase = 'Select filer size'
            else:
                for f in sorted(os.listdir(pathToSer))[::-1]:
                    if not f.endswith("png"): continue
                    k = f.split("_")[-1].split(".")[0]
                    options += [{"label": k, "value": os.path.join(pathToSer, k)}]
                    previewFile = os.path.join(pathToSer, f)
                    encoded_image = base64.b64encode(open(previewFile, 'rb').read())
                    preview[k] = html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode()),
                                          width="%ipx" % (width * .7), height="%ipx" % (height * .8 + 20))
                if len(options):
                    frase = "Select series"

            # if len(options):
            #     val = options[0]["value"]
            # else:
            #     val = None
            previews = [html.Div(children=[html.Div(k, style={"padding-left": "20px"}), preview[k]], style=dict(
                height="%ipx" % height,
                width="%ipx" % width,
                display="inline-block",
                border='thin lightgrey solid'
            )) for k in preview]
            return options, html.Div(previews, style={"width": "%ipx" % (width * 4.1)}), frase

        @app.callback(
            [Output("pickle-feedback", "children"),
             Output("examiner", "src")
             ],
            [Input("pickle-radio", "value")]
        )
        def import_pickle(path):
            link2app = ""
            feedback = []
            try:
                if self.debug:
                    globpath = path
                if path is None:
                    return no_update
                if len(path) == 0:
                    return no_update
                if path.endswith("pkl"):  ####### pickles ################
                    self.regions = load_regions(path)
                    self.regions.get_activity(max(10, 1/self.regions.Freq*5))
                    feedback += [
                        "Regions imported successfully.",
                        "Original movie:",
                        html.Div(children=[
                            "dimension (T,X,Y): (%i, %i, %i)" % ((len(self.regions.time),) + self.regions.image.shape),
                            html.Br(),
                            "duration (h:m:s): " +
                            str(pd.Timedelta(self.regions.time.max(), unit="s")).split()[-1].split(".")[0],
                            html.Br(),
                            "frequency (Hz): %.3g" % self.regions.Freq,
                        ], style={'padding-left': '30px'}),
                        "Number of ROIs: %i" % len(self.regions.df),
                    ]
                    if self.showExamine:
                        japp = self.regions.examine(max_rois=self.max_rois)
                        japp._repr_html_()
                        link2app = "https://ctn.physiologie.meduniwien.ac.at" + japp.get_app_root_url()
                else:  ####### line scan ################
                    self.pathToRec, ser = os.path.split(path)
                    self.pathToRec = os.path.split(self.pathToRec)[0].split("_analysis")[0]
                    feedback += [f"Loading the line scan {ser} from {self.pathToRec}"]
                    rec = Recording(self.pathToRec)
                    md = rec.metadata.set_index("Name")
                    #                 feedback += [str(md.index)]
                    line_scan = md.loc[ser, "line scan"]
                    #                 feedback += [line_scan]
                    rec.import_series(ser, isLineScan=(line_scan == "single"))
                    #                 feedback += [html.Br(),"data imported"]
                    data = rec.Series[ser]["data"].sum(1)
                    #                 feedback += [html.Br(),str(data.shape)]
                    self.linescan = LineScan(
                        data=data.T,
                        metadata=rec.Series[ser]["metadata"],
                        name="%s: %s" % (rec.Experiment[:-4], ser),
                    )
                    japp = self.linescan.examine()
                    japp._repr_html_()
                    link2app = "https://ctn.physiologie.meduniwien.ac.at" + japp.get_app_root_url()
            #         except ImportError as err:
            #             feedback += [str(exc_info())+f" in {err.path}"]
            #             feedback += [f" in {dir(err)}"]
            #             tb = err.__traceback__
            #             feedback += [f" in {dir(tb)}"]
            #             tb = tb.__dict__
            #             for kk in ['tb_frame', 'tb_lasti', 'tb_lineno', 'tb_next']:
            #                 feedback += [f"{kk}: {getattr(tb,kk)}"]
            #             feedback += [f"{err.msg}"]

            except:
                feedback += [str(exc_info())]
            return feedback, link2app

        return app
