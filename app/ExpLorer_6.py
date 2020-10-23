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
# from islets.Recording import Recording
from general_functions import td2str#, td2secs

addinfoFeatures = """comments
sex
strain
species
dye
slice number
part of pancreas
microscope
""".splitlines()

def serve_examiner(username="srdjan", pickleFilename = "/data/Sandra/2019/2019_08_09/Experiment27.lif_analysis/Series014/13.15_rois.pkl"):
    with open(f"/data/.tmp/empty_roi_examiner.ipynb","r") as f:
        nb_text = f.read()
    with open(f"/data/.tmp/{username}_roi_examiner.ipynb","w") as f:
#         f.write(nb_text%pickleFilename)
        f.write(nb_text.replace("#####", pickleFilename))
    return f"https://ctn.physiologie.meduniwien.ac.at/user/{username}/notebooks/local_data/.tmp/{username}_roi_examiner.ipynb"



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
                pathToAddInfo = os.path.split(pathToProtocol)[0]
                pathToAddInfo = os.path.join(pathToAddInfo, "additional_info.txt")
                md["path to add_info"] = pathToAddInfo
                md["add_info done"] = os.path.isfile(pathToAddInfo)
                if md["add_info done"] and os.path.getsize(pathToAddInfo)>10:
                    # print (ser, )
                    addInfo = pd.read_csv(pathToAddInfo, sep=":", header=None, index_col=0).T
                    if len(addInfo)==0:
                        md["add_info done"] = False
                        continue
                    for kk in addInfo.columns:
                        # print ("%s:%s"%(kk, addInfo[kk].iloc[0]), end=" ")
                        md[str(kk).strip()] = str(addInfo[kk].iloc[0]).strip()
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
               "protocol done",
               "comments",
               "sex",
               "strain",
               "species",
               "dye",
               "slice number",
               "part of pancreas",
               "microscope",
               "path to add_info",
               "add_info done"
     ]

showCols = ['date', 'experiment', 'series', 'Frequency', 'pxSize', 'Duration','SizeX', 'SizeY','protocol',"comments","Time Range", "movie size [MB]",]         

def prepareDF(mainFolder, constrain="",ishow=None):
    status = import_data(mainFolder, constrain)
    for addinfo in addinfoFeatures:
        if addinfo not in status:
            status[addinfo] = [""]*len(status)
    df = status[
        [c for c in orderedCols if c in status.columns]#+[c for c in df.columns if c not in orderedCols]
               ].copy()
    del df["Duration"]
    df["Duration"] = status["Duration [s]"].apply(td2str).values
    if ishow is not None:
        df = df.iloc[:ishow]
    return df

startFolder = "/data/Sandra/2020/2020_07_08/"
exceptCols = ["path to movie", "path to protocol", "path to add_info"]
Nrows = 4
from math import ceil

conditionalFormats = [
    {'if': {'column_id': col}, 'width': '100px'} \
    for col in ["date","Series Durations","comments","protocol"]
]
conditionalFormats += [
    {'if': {'column_id': col}, 'width': '80px'} \
    for col in ["series"]
]
conditionalFormats += [
    {'if': {'column_id': col}, 'width': '120px'} \
    for col in ["experiment"]
]

users = ["srdjan","johannes","sandra","marjan","arianna"]
app = dash.Dash(__name__,suppress_callback_exceptions=True)

# if __name__=="__main__":
df = prepareDF(startFolder,"",0)
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

sortedCols=df.columns
dfrec = df.to_dict('records')
del df

presentationCols = []
for i in sortedCols:
    if i in exceptCols: continue
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



app.layout = html.Div(children=[
    html.H2(children='Exp!Lorer'),
#     html.H4("the experiments explorer"),
    dcc.Markdown(" ~ _the experiments explorer_ ~"),
    html.Div([
        "Who are you?",
        dcc.Dropdown(id="username",
                     options = [{"value":un, "label":un} for un in users],
                     value="srdjan",style={"width":"30%","display":"inline-block"}),
    ],style={'align-items': 'center', 'display': 'flex'}),
    html.Div([
        html.Div("Folder to parse:",style={"width":"130px","display":"inline-block"}),
        dcc.Input(id="main_folder", value="/data/Sandra/2020/2020_07_08/"),
        html.Br(),
        html.Div("Constrain:",style={"width":"130px","display":"inline-block"}),
        dcc.Input(id="constrain",   size="33%", value="",placeholder="e.g. 2019_08_"),
        html.Button(id="update-table",children="Update table",n_clicks=1),
    ],),   
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

        html.Div(id="explorer-link", children="nothing yet", style={"font-family":"Courier New"}),
    html.Div(id="video", children="nothing yet", style={"font-family":"Courier New",}),
    html.Div(
        dash_table.DataTable(
            id='table',
            columns=[el for el in presentationCols if el["name"] in showCols],
            row_selectable='single',
            data=dfrec,
            style_header={'whiteSpace':'normal'},
            filter_action="native",
            sort_action="native",
            css=[{'selector': 'table', 'rule': 'table-layout: fixed'}],
            style_table={
                'overflowX': 'auto'
            },
            style_cell={
                'width': '70px',
                'overflow': 'hidden',
                'textOverflow': 'ellipsis',
            },
            tooltip_duration=None,
            style_cell_conditional=conditionalFormats,
            ),
    ),
] + [html.Br()]*10)


@app.callback(
     Output("table","columns"),
    [Input("cols_chooser_%i"%i,"value") for i in range(Ncols)]
)
def cols_sel(*collist):
    collist = sum(collist,[])
    if "all" in collist:
        collist = sortedCols
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
    Output("addinfo-save-out","children"),
    [Input("addinfo-save-button","n_clicks")],
    [State("addinfo-text","value"),
     State("addinfo-path","children"),]
)
def save_addinfo(n_clicks, txt, path):
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
     Output("explorer-link","children"),
    [Input("table","selected_rows"),
     Input("username","value")],
    [State("table","data")]
)
def serve_explorer_link(selected_rows, username, data):
    try:
        if selected_rows is None:
            return "Select a row to show details." #######################
        else:
            return serve_examiner(username=username)
    except:
        return str(exc_info())

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
        ###### additional information
        
        addinfoPath = data[ix]["path to add_info"]
        if not os.path.isfile(addinfoPath):
            with open(addinfoPath,"w") as ff:
                ff.write(":\n".join(addinfoFeatures))
        with open(addinfoPath,"r") as ff:
            txt = ff.read()
        out += [
            html.Div([
                "Addinfo file",
                html.Br(),
                html.Div(
                    addinfoPath,
                    id="addinfo-path",
                    style={"display":"none"}
                ),
                dcc.Textarea(
                    id="addinfo-text",
                    value = txt,
                    wrap="off",
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
                    'Save additional information',
                    id='addinfo-save-button',
                    n_clicks=0,
                ),
                html.Div(id="addinfo-save-out")
            ],
            style={"width":"400px","height":"400px","display":"inline-block","vertical-align":"text-top"})
        ]

        return out
    except:
        return str(exc_info())



if __name__ == '__main__':
    app.run_server(debug=True)