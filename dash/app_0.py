# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
# import os

# external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__,
                # external_stylesheets=external_stylesheets
                )

app.layout = html.Div(children=[
    html.H1(children='Hello Dash'),
    # html.H2(children='Hello Dash'),
    dcc.Input(debounce=True, id="series-folder-input", value=""),
    dcc.Dropdown(
        id='rois-list-dropdown',
        options=[],
        # value='NYC'
    ),
])

@app.callback(
    [Output("rois-list-dropdown","options")],
    [Input("series-folder-input","value")],
)
def get_rois_from_folder(serdir):
    if serdir!="":
        if os.path.isdir(serdir):
            labels = [{"label":f, "value":f} for f in os.listdir(serdir)]
            return labels
        else:
            return []
    else:
        return []

if __name__ == '__main__':
    app.run_server(
        debug=True,
        # use_reloader=True
    )


# from dash import Dash
# import dash_html_components as html
# from dash.dependencies import Input, Output
# import dash_core_components as dcc
#
# app = Dash()
#
# app.layout = html.Div(children=[
#     html.H1("Title"),
#     # html.H1("Title"),
#     "Input:",
#     dcc.Input(id="input_box", placeholder="some input", value=""),
#     html.Div(id="output_box", children="nothing yet")
#     ])
#
# # @app.callback(
# #     Output("output_box","children"),
# #     [Input("input_box","value")]
# # )
# # def io_callback(v):
# #     # from dash import no_update
# #     if v is not None:
# #         if v=="":
# #             return ""#no_update
# #         return "Ouput: "+v
#
#
# if __name__ == '__main__':
#     app.run_server(debug=True)
