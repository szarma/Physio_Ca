# -*- coding: utf-8 -*-
import dash
import dash_table
from dash_table.Format import Format, Scheme#, Sign, Symbol
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

import json
import os
from sys import exc_info
from sys import path as syspath
syspath.append("/home/jupyter-srdjan/srdjan_functs/")
import pandas as pd
from islets.Recording import Recording, parse_leica
from general_functions import td2str, td2secs
import parsing_3 as parser
# parser.update_protocols()
df = parser.import_protocols(parseMetadata=True)
# df = df.iloc[:4]


def get_series_file(row, ext):
    folder, out = get_series_dir(row)
    fldr_content = os.listdir(folder)
    files = [el for el in fldr_content if el.endswith(ext)]
    if len(files)==0:
        out += [f"Sorry, no {ext} file found in {folder}."]
        return None, out #######################
    
    fileSizes = pd.Series({afile: os.path.getsize(folder+afile)/1e6 for afile in files})
    fileSizes = fileSizes.sort_values()
    if len(files)>1:
        out += [html.Br(), f"WARNING: More than one {ext} file found in {folder}: "]
        for afile in fileSizes.index:
            out += [html.Br(), "%s [%.3gMB]"%(afile, fileSizes[afile])]
        out += [html.Br(), "Serving you the larger."]
    else:
        afile = fileSizes.index[-1]
        out += [html.Br(), "Serving %s [%.3gMB]"%(folder+afile, fileSizes[afile])]
    filepath = folder+afile
    return filepath, out

def get_series_dir(row):
    folder = row.path+f"_analysis/"
    if not os.path.isdir(folder):
        return None, ["%s analysis directory not found."%folder ]

    relevantsubdirs = [sd for sd in  os.listdir(folder) if row.series in sd]
    if len(relevantsubdirs)==0:
        return None, [f"{row.series} not found in {folder}"]
    out = []
    subdir = sorted(relevantsubdirs, key=len)[0]
    if len(relevantsubdirs)>1:
        out += [html.Br(), f"Multiple subdirecories found in {folder}. Taking {subdir}"]
    folder = folder+subdir+"/"
    return folder, out


protocols = []
parse_prot = []
for _,row in df.iterrows():
    protocolFile, output = get_series_file(row,"protocol.txt")
    try:
        prot = pd.read_csv(protocolFile)
        x = ", ".join(pd.np.unique([f"{row['compound']}:{row['concentration']}".replace(" ","") for _,row in prot.iterrows()]))
    except:
        x = "N/A"
    protocols  += [x]
    parse_prot += [output]

df["protocol"] = protocols
df["parse protocols"] = parse_prot


orderedCols = ['date', 'experiment', 'series',
               'Frequency',
               'Duration',
               'Duration [s]',
               'protocol',
               'pxSize',
               'SizeX', 'SizeY',
               'SizeT',
               'comments',
               'Slice number',
               'mouse/human',
               'sex', 'microscope',
               'dye', 'part of pancreas','filename', 'ext', 'bit depth','Start time', 'End time',  'Parsing Error', 'path', "parse protocols"]
cols_editable = ['comments','Slice number','mouse/human','sex','dye', 'part of pancreas',]
df = df[[c for c in orderedCols if c in df.columns]+[c for c in df.columns if c not in orderedCols]]

df["Duration [s]"] = df["Duration"].apply(td2secs)
df["Duration"] = df["Duration"].apply(td2str)

showCols = ['date', 'experiment', 'series', 'Frequency', 'pxSize', 'Duration','SizeX', 'SizeY',
               'protocol',]

# df["pxSize"] = df["pxSize"].round(3)

Nrows = 5
from math import ceil
allCols = [{"label": col, "value":col} for col in df.columns]+[{"label": "show all", "value":"all"}]
Ncols = int(ceil(len(allCols)/Nrows))
colLists = [allCols[i*Nrows:(i+1)*Nrows] for i in range(Ncols)]
# print (Ncols, Nrows, colLists)




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
    print ("-"*20,i)
    print (a)
    checkboxes += [
        a
    ]
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

app = dash.Dash(__name__)

app.layout = html.Div(children=[
    html.H2(children='ExpLorer'),
#     html.H4("the experiments explorer"),
    dcc.Markdown(" ~ _the experiments explorer_ ~"),
    checkboxes,
    html.Div(id="output",children="",style={
            "font-family":"Courier New",
            "font-size":"80%",
            'border': 'thin lightgrey solid',
            "width":"80%",
            "height":"100px",
            'overflowX': 'scroll',
            'overflowY': 'scroll',
#             'display':'none',
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
#             column_editable=cols_editable,
#             editable=cols_editable,
#             fixed_rows={'headers': True},
#             style_data={
#                 'whiteSpace': 'normal',
#                 'height': 'auto',
#                 'lineHeight': '15px',
#             },
            style_cell={
#                 'height': 'auto',
#                 'minWidth': '50px', 'width': '50px', 'maxWidth': '50px',
#                 'whiteSpace': 'normal',
#             },
#             style_cell={
#                 'overflow': 'hidden',
#                 'textOverflow': 'ellipsis',
#                 'maxWidth': 0,
            },
            style_cell_conditional=[
                {'if': {'column_id': 'date'},
                 'width': '30px'},
                {'if': {'column_id': 'Duration [s]'},
                 'format': '%.1f'},
                {'if': {'column_id': 'experiment'},
                 'width': '60px'},
                {'if': {'column_id': 'series'},
                 'width': '100px'},
                {'if': {'column_id': 'comments'},
                 'width': '200px'},
#                 {'if': {'column_id': 'protocol'},
#                  'width': '200px', 'minWidth': '200px', 'maxWidth': '200px', },
            ],
#             tooltip_data=[
#                 {
#                     column: {'value': str(value), 'type': 'markdown'}
#                     for column, value in row.items()
#                 } for row in df.to_dict('rows')
#             ],
#             tooltip_duration=None
            style_table={
                'overflowX': 'auto'
            },

            ),
    ),


])

@app.callback(
     Output("table","columns"),
    [Input("cols_chooser_%i"%i,"value") for i in range(Ncols)]
)
def cols_sel(*collist):
    collist = sum(collist,[])
    if "all" in collist:
        collist = df.columns
#     return [{"name": i, "id": i, "editable": i in cols_editable} for i in df if i in collist]
    return [el for el in presentationCols if el["name"] in collist]



@app.callback(
     Output("output","children"),
    [
        Input("table","active_cell"),
        Input("table","selected_rows"),
     ],
)
def show(active_cell, selected_rows):
    out = {
        "active_cell":   active_cell,
        "selected_rows": selected_rows,
           }
    if selected_rows is not None:
        out["selected_row"] = selected_rows[0]
    return json.dumps(out)

@app.callback(
     Output("video","children"),
    [Input("table","selected_rows")]
)
def show_video(selected_rows):
    try:
        if selected_rows is None:
            return "Select a row to show video." #######################
        ix = selected_rows[0]
        row = df.loc[ix]
        moviePath, out = get_series_file(row, ".mp4")
        if moviePath is not None:
            out += [
                html.Br(),
                html.Video(
                    src=app.get_asset_url(moviePath.lstrip("/")),
                    width = 400,
                    controls=True,
                    loop=True,
                    )
            ]
        return out
    except:
        return str(exc_info())



if __name__ == '__main__':
    app.run_server(debug=True)