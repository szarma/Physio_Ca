# -*- coding: utf-8 -*-
import difflib
import subprocess

import dash
from dash import html, dcc, dash_table
from dash.dash_table.Format import Format, Scheme  # , Sign, Symbol
from dash.dependencies import Input, Output, State
from dash import no_update
import json
import os
from sys import exc_info
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


# startFolder = "/data/Sandra/2020/2020_07_08/"
constr = None
startFolder, constr = "/data/Ariana", ""
# startFolder = "/data/Sandra/2019/2019_09_03/"
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

users = sorted(["srdjan", "johannes", "sandra", "marjan", "nastja", "ya-chi", "dean", "lidija", "anita", "natalia"])
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

buttonStyle = {"font-family": "Courier New",
               "background-color": "navy",
               "display": "inline-block",
               "color": "white",
               "padding": "5px"}

preStyle = {"font-family": "Courier New",
            "font-size": "80%",
            'border': 'thin lightgrey solid',
            "width": "80%",
            "height": "100px",
            'overflowX': 'scroll',
            'overflowY': 'scroll',
            }
app.layout = html.Div(
    children = [
        html.H2(children = '(Re)Processer'),
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
               html.Div(id = "table-feedback"),
           ], ),
           checkboxes,
           html.Div(id = "parser", children = "", style = {'display': 'none', **preStyle}),
           html.Br(),
           html.Div("Please choose what you wish to reprocess."),
           html.Br(),
           html.Div("At the bottom of the page you can see how the script will look like."),
           html.Div(
               dash_table.DataTable(
                   id = 'table',
                   columns = [el for el in presentationCols if el["name"] in showCols],
                   row_selectable = 'multi',
                   data = dfrec,
                   selected_rows = [],
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

           # html.Pre(id = "script-sketch",
           #          children = "nothing yet",
           #          style = preStyle
           #          ),
           html.Div(id = "script-sketch",
                    children = "nothing yet",
                    ),
        ] +
        [html.Br()] * 10
)

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
            return df.to_dict('records'), [], tooltip_data, "ok"
        except:
            return no_update, no_update, no_update, str(exc_info())

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
    Output("script-sketch", "children"),
    [Input("table", "selected_rows"),
     Input("username", "value")],
    [State("table", "data"), ]
)
def serve_stuff(selected_rows, username, data):
    out = []
    try:
        out +=  [dcc.Markdown(
            """Select rows you wish to (re)process.  
            When you are done clicking, you can also manually add optional arguments like `--notify` or `--restrict`.  
            Finally, go to your home folder, create a new text file, call it e.g. `process_script.sh`, and copy the content of this text box into it.  
            Save. Then open a terminal and type `bash process_script.sh` and hit enter."""
            , style={"font-size":"70%"})]

        commands = []
        for ix in selected_rows:
            pathToExp = data[ix]["path to exp"]
            if pathToExp.endswith("lif"):
                series = data[ix]["series"]
            else:
                series = "all"
            command = f"""/data/useful/scripts/full_process.py '{pathToExp}' --series='{series}'"""
            restrict = data[ix]["Time Range"]
            # saveDir  = f"{pathToExp}_analysis/{series}"
            if restrict != "all":
                restrict = restrict.replace("-","_", count=1)
                # saveDir += restrict
                command += f" --restrict='{restrict}'"
            commands += [command+"\n"]
        out += [dcc.Textarea(
            id = "script",
            value = "".join(commands),
            style= preStyle
        ),
            # html.Button("Save",id="save-script",style = {"display": "block",buttonStyle})
        ]


    except:
        out += [html.Dic(str(exc_info()))]
    return out


if __name__ == '__main__':
    app.run_server(debug = True)