# -*- coding: utf-8 -*-
import sys
debug = sys.argv[0]=="explorer1.py"
if debug:
    print ("debug mode")
import subprocess
import traceback
import dash
from dash import html, dcc, dash_table
from dash.dash_table.Format import Format, Scheme  # , Sign, Symbol
from dash.dependencies import Input, Output, State
from dash import no_update
import json
import os
from islets.general_functions import td2str
from islets.utils import import_data

addinfoFeatures = """comments
sex
strain
species
dye
slice number
part of pancreas
microscope
""".splitlines()

def trcbk():
    return repr(traceback.format_tb())

users = sorted(["srdjan","johannes","sandra","marjan","nastja","ya-chi","dean","lidija","anita","natalia","jan_k","tito"])

def save_meta(n_clicks, txt, path, filetype):
    if n_clicks > 0:
        path = path.split(")")[0].split("(")[-1]
        folder = os.path.split(path)[0]
        nConflicts = len([f for f in os.listdir(folder) if filetype in f]) - 1
        if nConflicts:
            out = html.Div(
                f"There is at least one other {filetype} file in the same folder. Please remove it before saving.",
                style = {"font-size": "small", "font-color": "darkred"})
        else:
            # with open(path, "w") as f:
            #     f.write(txt)
            proc = subprocess.Popen(
                f'''sudo sh -c "echo -n '{txt}' > '{path}'"''',
                shell = True,
                stdout = subprocess.PIPE,
                text=True
            )
            outs, errs = proc.communicate()
            txt1 = open(path).read()
            if txt == txt1:
                out = "File saved ✔"
            else:
                out = "File not saved ✘"
                # out += "\nOutput: "+str(outs)
                # out += "\nError: "+str(errs)
                # output_list = [li for li in difflib.ndiff(txt, txt1) if li[0] != ' ']
                # out += "\ndifference: > "+ repr(output_list) + " <"
        return out

def serve_examiner(username="srdjan",
                   rec="",
                   ser="",
                   ):
    with open(f"../notebooks/new_empty_roi_examiner.ipynb", "r") as f:
        nb_text = f.read()
    nb_text = nb_text.replace("##pathToRec##", rec)
    nb_text = nb_text.replace("##series##", ser)
    with open(f"/data/.tmp/{username}_roi_examiner.ipynb","w") as f:
        f.write(nb_text)
    return f"https://ctn.physiologie.meduniwien.ac.at/user/{username}/notebooks/local_data/.tmp/{username}_roi_examiner.ipynb"


cols_editable = []
orderedCols = ['date', 'experiment', 'series',
               'line scan',
               "Time Range",
               'Frequency',
               'Duration',
               'protocol',
               'pxSize',
               'SizeT', 'SizeX', 'SizeY', 'SizeZ',
               'comments',
               'Series Durations',
               "movie size [MB]",
               'Duration [s]',
               "path to movie",
               "path to exp",
               "path to protocol",
               "protocol done",
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

showCols = ['date', 'experiment', 'series', 'Frequency', 'pxSize', 'Duration', 'SizeX', 'SizeY', 'protocol', "comments",
            "Time Range", "movie size [MB]", "line scan"]


def prepareDF(mainFolder, constrain="", ishow=None):
    status = import_data(mainFolder, constrain)
    for addinfo in addinfoFeatures:
        if addinfo not in status:
            status[addinfo] = [""] * len(status)
    df = status[
        [c for c in orderedCols if c in status.columns]  # +[c for c in df.columns if c not in orderedCols]
    ].copy()
    del df["Duration"]
    df["Duration"] = status["Duration [s]"].apply(td2str).values
    if ishow is not None:
        df = df.iloc[:ishow]
    return df


startFolder, constr = "/data/Sandra/", "125"
exceptCols = ["path to movie", "path to protocol", "path to add_info"]
Nrows = 4
from math import ceil

conditionalFormats = [
    {'if': {'column_id': col}, 'width': '100px'} \
    for col in ["date", "Series Durations", "comments", "protocol"]
]
conditionalFormats += [
    {'if': {'column_id': col}, 'width': '100px'} \
    for col in ["series"]
]
conditionalFormats += [
    {'if': {'column_id': col}, 'width': '120px'} \
    for col in ["experiment"]
]


app = dash.Dash(__name__, suppress_callback_exceptions = True)

df = prepareDF(startFolder, constr, 0)
allCols = [{"label": col, "value": col} for col in df.columns if col not in exceptCols] + [
    {"label": "show all", "value": "all"}]
Ncols = int(ceil(len(allCols) / Nrows))
colLists = [allCols[i * Nrows:(i + 1) * Nrows] for i in range(Ncols)]

checkboxes = []
for i, cols in enumerate(colLists):
    a = dcc.Checklist(
        id = "cols_chooser_%i" % i,
        options = cols,
        value = [col["value"] for col in cols if col["value"] in showCols],
        style = {
            "padding": "5px",
            "display": "inline-block",
        }
    )
    checkboxes += [a]
checkboxes = html.Div(checkboxes)

sortedCols = df.columns
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
    if i == "pxSize":
        singleCol["type"] = "numeric"
        singleCol["format"] = Format(precision = 3, scheme = Scheme.fixed, )
    if i in ["Duration [s]", "Frequency"]:
        singleCol["type"] = "numeric"
        singleCol["format"] = Format(precision = 1, scheme = Scheme.fixed, )
    presentationCols += [singleCol]

buttonStyle = {"font-family": "Courier New", "background-color": "navy", "display": "inline-block", "color": "white",
               "padding": "5px"}
app.layout = html.Div(
    children = [
                   html.H2(children = 'Exp!Lorer'),
                   #     html.H4("the experiments explorer"),
                   dcc.Markdown(" ~ _the experiments explorer_ ~"),
                   html.Div([
                       html.Div("Who are you?",
                                style = {"width": "130px", "display": "inline-block"}),
                       dcc.Dropdown(id = "username",
                                    options = [{"value": un, "label": un} for un in users],
                                    #                      value="srdjan",
                                    placeholder = "Please choose your username",
                                    style = {"width": "250px", "display": "inline-block"}
                                    ),
                   ], style = {'align-items': 'center', 'display': 'flex'}),
                   html.Div([
                       html.Div("Folder to parse:",
                                style = {"width": "130px", "display": "inline-block"}),
                       dcc.Input(
                           id = "main_folder",
                           value = startFolder,
                           size = "47"
                       ),
                       html.Br(),
                       html.Div("Constrain:", style = {"width": "130px", "display": "inline-block"}),
                       dcc.Input(
                           id = "constrain",
                           value = constr if "constr" in globals() else "",
                           placeholder = "e.g. _08_, lif",
                           size = "30"
                       ),
                       html.Button(id = "update-table", children = "Update table", n_clicks = 1),
                       html.Pre(id = "table-feedback", style = {"font-color":"lightgrey"} ),
                   ], ),
                   checkboxes,
                   html.Div(id = "parser", children = "", style = {
                       "font-family": "Courier New",
                       "font-size": "80%",
                       'border': 'thin lightgrey solid',
                       "width": "80%",
                       "height": "100px",
                       'overflowX': 'scroll',
                       'overflowY': 'scroll',
                       'display': 'none',
                   }),

                   html.Div(id = "explorer-link", children = "nothing yet", style = buttonStyle),
                   html.Div(id = "video", children = "nothing yet",
                            style = {"font-family": "Courier New", }),
                   html.Div(
                       dash_table.DataTable(
                           id = 'table',
                           columns = [el for el in presentationCols if el["name"] in showCols],
                           row_selectable = 'single',
                           data = dfrec,
                           #             selected_rows=[min(len(dfrec),3)],
                           style_header = {'whiteSpace': 'normal'},
                           filter_action = "native",
                           sort_action = "native",
                           css = [{'selector': 'table', 'rule': 'table-layout: fixed'}],
                           style_table = {
                               'overflowX': 'auto'
                           },
                           style_cell = {
                               'width': '70px',
                               'overflow': 'hidden',
                               'textOverflow': 'ellipsis',
                           },
                           tooltip_duration = None,
                           style_cell_conditional = conditionalFormats,
                       ),
                   ),
               ] + [html.Br()] * 10)


@app.callback(
    Output("table", "columns"),
    [Input("cols_chooser_%i" % i, "value") for i in range(Ncols)]
)
def cols_sel(*collist):
    collist = sum(collist, [])
    if "all" in collist:
        collist = sortedCols
    return [el for el in presentationCols if el["name"] in collist]


@app.callback(
    [Output("table", "data"),
     Output("table", "selected_rows"),
     Output("table", "tooltip_data"),
     Output("table-feedback", "children")
     ],
    [Input("update-table", "n_clicks"), ],
    [State("main_folder", "value"),
     State("constrain", "value")]
)
def update_table(n_clicks, main, restrict):
    if n_clicks > 0:
        try:
            df = prepareDF(main, restrict)
            tooltip_data = [
                {column: {'value': str(value).replace("\n", "\n\n"), 'type': 'markdown'} for column, value in
                 row.items() if column in ["protocol", "Series Durations"]} for row in df.to_dict('rows')]
            return df.to_dict('records'), [0], tooltip_data, "ok"
        except:
            return no_update, no_update, no_update, trcbk()#str(exc_info())


@app.callback(
    Output("protocol-save-out", "children"),
    [Input("protocol-save-button", "n_clicks")],
    [State("protocol-text", "value"),
     State("protocol-path", "children"), ]
)
def save_protocol(n_clicks, txt, path):
    return save_meta(n_clicks, txt, path, filetype="protocol")

@app.callback(
    Output("addinfo-save-out", "children"),
    [Input("addinfo-save-button", "n_clicks")],
    [State("addinfo-text", "value"),
     State("addinfo-path", "children"), ]
)
def save_addinfo(n_clicks, txt, path):
    return save_meta(n_clicks, txt, path, filetype = "additional_info")


@app.callback(
    Output("parser", "children"),
    [
        Input("table", "active_cell"),
        Input("table", "selected_rows"),
    ],
    [State("table", "data"), ]
)
def show(active_cell, selected_rows, data):
    out = {
        "active_cell": active_cell,
        "selected_rows": selected_rows,
        #         "data": data
    }
    i = active_cell["row"]
    col = active_cell["column"]
    if selected_rows is not None:
        out["selected_row"] = selected_rows[-1]
    return json.dumps(out)


@app.callback(
    Output("explorer-link", "children"),
    [Input("table", "selected_rows"),
     Input("username", "value")],
    [State("table", "data")]
)
def serve_explorer_link(selected_rows, username, data):
    try:
        if selected_rows is None:
            return "Select a row to show details."  #######################
        if username is None:
            return "Please choose your username from the dropdown menu above."
        else:
            ix = selected_rows[-1]
            rec = data[ix]["path to exp"]
            ser = data[ix]["series"]
            link = serve_examiner(username = username, rec = rec, ser = ser)
            out = html.A('Examine', href = link, target = '_blank')
            #             out = dcc.Link("Examine", href=link, )
            return out
    except:
        return trcbk()


@app.callback(
    Output("video", "children"),
    [Input("table", "selected_rows"), Input("username", "value")],
    [State("table", "data"), ]
)
def serve_stuff(selected_rows, username, data):
    try:
        if selected_rows is None:
            return "Select a row to show details."  #######################
        ix = selected_rows[-1]
        protocolPath = data[ix]["path to protocol"]
        out = []
        try:
            folder = os.path.split(protocolPath)[0]
            hrefFolder = f"https://ctn.physiologie.meduniwien.ac.at/user/{username}/tree/local_{folder[1:]}"
            out += [html.Div(html.A("open folder", href = hrefFolder, target = '_blank'), style = buttonStyle),
                    html.Br()]
        except:
            pass
        line_scan = data[ix].get("line scan","none")
        if line_scan == "none":
            moviePath = data[ix]["path to movie"]
            if moviePath is not None:
                hrefMovie = f"https://ctn.physiologie.meduniwien.ac.at/user/{username}/files/local_{moviePath[1:]}?download=1"
                out += [
                    html.Div([
                        #                         "Path to video:",
                        "Video",
                        html.Br(),
                        html.Video(
                            src = app.get_asset_url(moviePath.lstrip("/")),
                            width = 350,
                            controls = True,
                            loop = True
                        ),
                        # dcc.Link(
                        #     children = str(moviePath),
                        #     style={"font-size":"50%"},
                        #     href=hrefMovie
                        #           ),
                        html.A(str(moviePath), href = hrefMovie, target = '_blank', style = {"font-size": "50%"}),
                    ], style = {"display": "inline-block",
                                "vertical-align": "text-top",
                                "width": "400px"}),
                ]
        else:
            saveDir = os.path.split(protocolPath)[0]
            allPngs = [os.path.join(saveDir, el) for el in os.listdir(saveDir) \
                       if el.endswith("png") and el[0] != "."
                       ]
            imageContainer = html.Div([
                html.Div(
                    html.Img(src = app.get_asset_url(imgPath.lstrip("/")), width = 300),
                    style = {"display": "inline-block"}
                ) for imgPath in allPngs
            ], style = {"width": 303 * len(allPngs)})
            out += [
                html.Div("Use horizontal scroll to see all series."),
                html.Div(
                    imageContainer,
                    style = {
                        "display": "inline-block",
                        "vertical-align": "text-top",
                        'overflowX': 'scroll',
                        "width": "320px",
                        "padding": "5px",
                        "height": "500px"
                    })
            ]
        ###### protocol #################
        try:
            if not os.path.isfile(protocolPath):
                with open(protocolPath, "w") as f:
                    f.write("")
            txt = open(protocolPath).read()
            protocolFolder = os.path.split(protocolPath)[0]
            noProts = len([ff for ff in os.listdir(protocolFolder) if "protocol" in ff])
            protocolContent = [
                "Protocol file",
                html.Br(),
                html.Div(
                    protocolPath,
                    id = "protocol-path",
                    style = {"display": "none"}
                ),
                dcc.Textarea(
                    id = "protocol-text",
                    value = txt,
                    style = {
                        "font-family": "Courier New",
                        #                         "font-size":"80%",
                        'border': 'thin lightgrey solid',
                        "width": "95%",
                        "height": "250px",
                        'overflowX': 'scroll',
                        'overflowY': 'scroll', }
                ),
                html.Button(
                    'Save protocol',
                    id = 'protocol-save-button',
                    n_clicks = 0,
                ),
                html.Div(id = "protocol-save-out")
            ]
        #             if not os.path.isfile(protocolPath):
        #                 with open(protocolPath,"w") as f: f.write("")
        #             txt = open(protocolPath).read()
        #             link = "https://ctn.physiologie.meduniwien.ac.at/user/srdjan/edit/local_"+protocolPath.lstrip("/")
        #             protocolContent = [
        #                 html.A(  'Edit protocol file', href=link, target='_blank')
        #             ]
        except:
            protocolContent = [trcbk()]
        protocolContent += [html.Div(str(protocolPath), style = {"font-size": "50%"}), ]
        out += [html.Div(protocolContent,
                         style = {"width": "400px",
                                  "height": "400px",
                                  "display": "inline-block",
                                  "vertical-align": "text-top"})
                ]

        ###### additional information #################
        try:
            addinfoPath = data[ix]["path to add_info"]
            if not os.path.isfile(addinfoPath):
                with open(addinfoPath,"w") as ff:
                    ff.write(":\n".join(addinfoFeatures+[""]))
            with open(addinfoPath, "r") as ff:
                txt = ff.read()
            addInfoContent = [
                "Addinfo file",
                html.Br(),
                html.Div(
                    addinfoPath,
                    id = "addinfo-path",
                    style = {"display": "none"}
                ),
                dcc.Textarea(
                    id = "addinfo-text",
                    value = txt,
                    wrap = "off",
                    style = {
                        "font-family": "Courier New",
                        #                         "font-size":"80%",
                        'border': 'thin lightgrey solid',
                        "width": "95%",
                        "height": "250px",
                        'overflowX': 'scroll',
                        'overflowY': 'scroll', }
                ),
                html.Button(
                    'Save additional information',
                    id = 'addinfo-save-button',
                    n_clicks = 0,
                ),
                html.Div(id = "addinfo-save-out")
            ]
        except:
            addInfoContent = [trcbk()]
        addInfoContent += [html.Div(str(addinfoPath), style = {"font-size": "50%"}), ]

        out += [html.Div(addInfoContent,
                         style = {"width": "400px",
                                  "height": "400px",
                                  "display": "inline-block",
                                  "vertical-align": "text-top"})
                ]

    except:
        out += [trcbk()]
    return out


if __name__ == '__main__':
    app.run_server(debug = True)
