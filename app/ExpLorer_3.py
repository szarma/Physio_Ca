# -*- coding: utf-8 -*-
import dash
import dash_table
from dash_table.Format import Format, Scheme#, Sign, Symbol
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State

import json
import os
from sys import exc_info
from sys import path as syspath
syspath.append("/home/jupyter-srdjan/srdjan_functs/")

from general_functions import td2str#, td2secs


def import_data(mainFolder, constrain, forceMetadataParse=False):
    recordings = []
    for cur,ds,fs in os.walk(mainFolder):
        #### if you wish to restrict to only certain folders: ####
        if constrain not in cur: continue
        for f in fs:
            if not (f.endswith(".lif") or f.endswith(".nd2")):
                continue
            path = os.path.join(cur,f)
            recordings += [path]
    recordings = sorted(recordings)
    
    from islets.Recording import Recording, parse_leica
    from islets.utils import get_series_dir, get_filterSizes
    import numpy as np
    import pandas as pd
    
    status = []
    ilifs = 0
    for pathToRecording in recordings:
        print ("#"*20, pathToRecording)
        try:
            rec = Recording(pathToRecording)
        except:
            print ("Could not import ", pathToRecording)
            continue
        recType = "Nikon" if pathToRecording.endswith(".nd2") else "Leica"
        if forceMetadataParse:
            rec.parse_metadata()
            rec.save_metadata()
        if recType=="Leica":
            sers = parse_leica(rec)
        else:
            print ("Nikon not yet supported. Bug me to enable it.")
            continue
            

        analysisFolder = os.path.join(rec.folder, rec.Experiment+"_analysis")
        if not os.path.isdir(analysisFolder):
            os.makedirs(analysisFolder)

        for series in sers:
#             if recType=="Leica":
            subdirs = get_series_dir(pathToRecording, series)
#             else:
#                 subdirs = os.listdir(analysisFolder)
            for ser in subdirs:
                md = pd.Series()
                md["path to exp"] = pathToRecording
                md["experiment"] = os.path.split(pathToRecording)[-1]
                md["series"] = series
                try:
                    rec.import_series(series, onlyMeta=True)
                except:
                    print (f"could not import {series}")
                    status += [md]
                    continue
                
                saveDir = os.path.join(analysisFolder, ser)
                for k,v in rec.Series[series]["metadata"].items(): md[k] = v
                if "_" in ser:
                    t0,t1 = [float(t) for t in ser.rstrip("s").split("_")[-1].split("-")]
                    md["Time Range"] = "%i-%i"%(t0,t1)
                    md["Duration [s]"] = t1-t0
                else:
                    md["Time Range"] = "all"
                    md["Duration [s]"] = md["SizeT"]/md["Frequency"]
                fs = get_filterSizes(md.pxSize)
                movieFilename = os.path.join(saveDir, rec.Experiment+"_"+series+".mp4")
                md["path to movie"] = movieFilename
                md["movie done"] = os.path.isfile(movieFilename)
                if md["movie done"]:
                    md["movie size [MB]"] = np.round(os.path.getsize(movieFilename)/10**6,1)
                md["date"] = md["Start time"].date().__str__()
                for k in ["bit depth", "Start time", "End time","Name","frame_range"]: # , "individual Series"
                    try:    del md[k]
                    except: pass
                times = ["00:00"]+[td2str(el) for el in md["individual Series"]["Duration"].cumsum()]
                md["Series Durations"] = " \n".join(["%s [%s-%s]"%(name.lstrip("Series0"), t0, t1) for name, t0, t1 in zip(md["individual Series"]["Name"], times[:-1], times[1:])])
#                 md["Series Durations"] = ["%s [%s]"%(r["Name"].replace("Series","S"), td2str(td_nanfloor(r["Duration"]))) for _,r in md["individual Series"].iterrows()]
#                 md["Series Durations"] = "\n".join(md["Series Durations"])
                del md["individual Series"]
                pklsDone = {}
                for fsize in fs:
                    pickleFile = os.path.join(saveDir, ".".join(map(str,fsize))+"_rois.pkl")
                    pickleThere = os.path.isfile(pickleFile)
                    pklsDone[fsize] = pickleThere
                md["pickles done"] = pklsDone
                pathToProtocol = movieFilename.replace(".mp4","_protocol.txt").replace("_"+md["Time Range"]+"s","")
                md["path to protocol"] = pathToProtocol
                md["protocol done"] = False
                try:
                    protocol = pd.read_csv(pathToProtocol)
                    if len(protocol):
                        md["protocol done"] = True
                        protocol = " ".join(np.unique([
                            "%s:%s"%(row["compound"].capitalize() if "glu" in row["compound"].lower() else row["compound"], row["concentration"].replace(" ","")) for _,row in protocol.iterrows()]))
                        protocol = protocol.replace("Glucose","Glu")
                        md["protocol"] = protocol
                except:
                    pass
            
                status += [dict(md.items())]
        ilifs +=1
    #     if ilifs>3:
    #         break

    return pd.DataFrame(status)


cols_editable = []
orderedCols = ['date', 'experiment', 'series',
               "Time Range",
               'Frequency',
               'Duration',
               'protocol',
               'pxSize',
               'SizeT','SizeX', 'SizeY', 'SizeZ',
               'comments',
               'Series Durations',
               "movie size [MB]",
               'Duration [s]',
               "path to movie",
               "path to exp",
               "path to protocol",
               "protocol done"
#                'Slice number',
#                'mouse/human',
#                'sex', 'microscope',
#                'dye', 'part of pancreas','filename', 'ext', 'bit depth','Start time', 'End time',  'Parsing Error', 'path', "parse protocols"
     ]

showCols = ['date', 'experiment', 'series', 'Frequency', 'pxSize', 'Duration','SizeX', 'SizeY','protocol',"comments"]         

def prepareDF(mainFolder, constrain=""):
    status = import_data(mainFolder, constrain)
    if "comments" not in status:
        status["comments"] = [""]*len(status)

    df = status[
        [c for c in orderedCols if c in status.columns]#+[c for c in df.columns if c not in orderedCols]
               ].copy()
    del df["Duration"]
    df["Duration"] = status["Duration"].apply(td2str).values

    return df

df = prepareDF("/data/Sandra/2020/2020_07_08/")
df = df.iloc[:0]
exceptCols = ["path to movie", "path to protocol"]

Nrows = 3
from math import ceil
allCols = [{"label": col, "value":col} for col in df.columns if col not in exceptCols]+[{"label": "show all", "value":"all"}]
Ncols = int(ceil(len(allCols)/Nrows))
colLists = [allCols[i*Nrows:(i+1)*Nrows] for i in range(Ncols)]

checkboxes = []
for i,cols in enumerate(colLists):
    a = dcc.Checklist(
            id="cols_chooser_%i"%i,
            options = cols,
            value = [ col["value"] for col in cols if col["value"] in showCols],
            style = {
              "padding": "5px",
              "display":"inline-block",
            }
        )
    checkboxes += [a]
checkboxes = html.Div(checkboxes)

presentationCols = []
for i in df.columns:
    singleCol = {
        "name": i,
        "id": i,
        "editable": i in cols_editable,
     } 
    if i=="pxSize":
        singleCol["type"] = "numeric"
        singleCol["format"] = Format(precision=3, scheme=Scheme.fixed,)
    if i in ["Duration [s]", "Frequency"]:
        singleCol["type"] = "numeric"
        singleCol["format"] = Format(precision=1, scheme=Scheme.fixed,)
    presentationCols += [singleCol]


conditionalFormats = [
    {'if': {'column_id': col}, 'width': '100px'} \
    for col in ["experiment","date","series","Series Durations","comments","protocol"]
]
                
app = dash.Dash(__name__,suppress_callback_exceptions=True)

app.layout = html.Div(children=[
    html.H2(children='Exp!Lorer'),
#     html.H4("the experiments explorer"),
    dcc.Markdown(" ~ _the experiments explorer_ ~"),
    html.Div([
        "Folder to parse:",
        dcc.Input(id="main_folder", size="33%", value="/data/Sandra/2020/"),
    ],style={"width":"35%","display":"inline-block"}),
    html.Div([
        "Constrain:",
        dcc.Input(id="constrain",   size="33%", value="2020_07_08"),
    ],style={
        "width":"35%",
#         "display":"inline-block"
            }),   
    html.Button(id="update-table",children="Update table",n_clicks=1),
    checkboxes,
    html.Div(id="parser",children="",style={
            "font-family":"Courier New",
            "font-size":"80%",
            'border': 'thin lightgrey solid',
            "width":"80%",
            "height":"100px",
            'overflowX': 'scroll',
            'overflowY': 'scroll',
            'display':'none',
    }),
    
    html.Div(id="video", children="nothing yet", style={"font-family":"Courier New",}),
    html.Div(
        dash_table.DataTable(
            id='table',
            columns=[el for el in presentationCols if el["name"] in showCols],
            row_selectable='single',
            data=df.to_dict('records'),
            filter_action="native",
            sort_action="native",
            css=[{'selector': 'table', 'rule': 'table-layout: fixed'}],
            style_table={
                'overflowX': 'auto'
            },
            style_cell={
                'width': '50px',
                'overflow': 'hidden',
                'textOverflow': 'ellipsis',
            },
            tooltip_duration=None,
            style_cell_conditional=conditionalFormats,
            ),
    ),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Br(),


])

@app.callback(
     Output("table","columns"),
    [Input("cols_chooser_%i"%i,"value") for i in range(Ncols)]
)
def cols_sel(*collist):
    collist = sum(collist,[])
    if "all" in collist:
        collist = df.columns
    return [el for el in presentationCols if el["name"] in collist]

@app.callback(
    [Output("table","data"),
     Output("table","selected_rows"),
     Output("table","tooltip_data")
    ],
    [Input("update-table","n_clicks"),],
    [State("main_folder","value"),
     State("constrain","value")]
)
def update_table(n_clicks, main, restrict):
    if n_clicks>0:
        df = prepareDF(main,restrict)
        tooltip_data=[{column: {'value': str(value).replace("\n","\n\n"), 'type': 'markdown'} for column, value in row.items() if column in ["protocol","Series Durations"]} for row in df.to_dict('rows')]
        return df.to_dict('records'), [0], tooltip_data
    

@app.callback(
    Output("protocol-save-out","children"),
    [Input("protocol-save-button","n_clicks")],
    [State("protocol-text","value"),
     State("protocol-path","children"),]
)
def save_protocol(n_clicks, txt, path):
    if n_clicks>0:
        path = path.split(")")[0].split("(")[-1]
        with open(path,"w") as f:
            f.write(txt)
        txt1 = open(path).read()
        if txt==txt1:
            out = "File saved ✔"
        else:
            out = "File not saved ✘"
        return out



@app.callback(
     Output("parser","children"),
    [
        Input("table","active_cell"),
        Input("table","selected_rows"),
     ],
    [State("table","data"),]
)
def show(active_cell, selected_rows, data):
    out = {
        "active_cell":   active_cell,
        "selected_rows": selected_rows,
#         "data": data
           }
    i = active_cell["row"]
    col = active_cell["column"]
    if selected_rows is not None:
        out["selected_row"] = selected_rows[-1]
    return json.dumps(out)

@app.callback(
     Output("video","children"),
    [Input("table","selected_rows")],
    [State("table","data"),]
)
def serve_stuff(selected_rows, data):
    try:
        if selected_rows is None:
            return "Select a row to show details." #######################
        ix = selected_rows[-1]
        moviePath = data[ix]["path to movie"]
        out = []
        if moviePath is not None:
            out += [
                html.Div([
#                     moviePath,
                    "Video",
                    html.Br(),
                    html.Video(
                        src=app.get_asset_url(moviePath.lstrip("/")),
                        width = 350,
                        controls=True,
                        loop=True
                        )
                ],style={"display":"inline-block","vertical-align":"text-top","width":"400px"})
            ]
        protocolPath = data[ix]["path to protocol"]
        if not os.path.isfile(protocolPath):
            os.system(f"touch {protocolPath}")
        txt = open(protocolPath).read()
        out += [
            html.Div([
                "Protocol file",
                html.Br(),
                html.Div(
                    protocolPath,
                    id="protocol-path",
                    style={"display":"none"}
                ),
                dcc.Textarea(
                    id="protocol-text",
                    value = txt,
                    style={
                        "font-family":"Courier New",
#                         "font-size":"80%",
                        'border': 'thin lightgrey solid',
                        "width":"95%",
                        "height":"250px",
                        'overflowX': 'scroll',
                        'overflowY': 'scroll',}
                ),
                html.Button(
                    'Save protocol',
                    id='protocol-save-button',
                    n_clicks=0,
                ),
                html.Div(id="protocol-save-out")
            ],
            style={"width":"400px","height":"400px","display":"inline-block","vertical-align":"text-top"})
        ]
        
        return out
    except:
        return str(exc_info())



if __name__ == '__main__':
    app.run_server(debug=True)