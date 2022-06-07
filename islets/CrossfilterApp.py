try:
    from jupyter_plotly_dash import JupyterDash
except ModuleNotFoundError:
    from jupyter_dash import JupyterDash
try:
    from dash import dcc, html
except ImportError:
    import dash_core_components as dcc
    import dash_html_components as html

from dash.dependencies import Input, Output
import plotly.graph_objects as go

from collections import OrderedDict


def crossfilterApp(df,appWidth=1000, panelWidth=300, ms=10):
    baseFig = go.Figure(layout={
        "width":panelWidth/2,
        "height":panelWidth/2,
        "paper_bgcolor":"blue",
        "xaxis":{"tickvals":[]},
        "yaxis":{"tickvals":[]},
        "margin":dict(zip("lrbt",[3]*4))
    })
    def get_figure(df, ix, iy, selectedpoints, selectedpoints_local):
        x_col = df.columns[ix]
        y_col = df.columns[iy]
        if selectedpoints_local and selectedpoints_local['range']:
            ranges = selectedpoints_local['range']
            selection_bounds = {'x0': ranges['x'][0], 'x1': ranges['x'][1],
                                'y0': ranges['y'][0], 'y1': ranges['y'][1]}
        else:
            selection_bounds = {'x0': df[x_col].min(), 'x1': df[x_col].max(),
                                'y0': df[y_col].min(), 'y1': df[y_col].max()}
        return {
            'data': [{
                'x': df[x_col],
                'y': df[y_col],
                'text': df.index,
                'textposition': 'top',
                'selectedpoints': selectedpoints,
                'customdata': df.index,
                'type': 'scatter',
                'mode': 'markers',
                'marker': { 'color': 'rgba(0, 116, 217, 0.7)', 'size': ms },
                'unselected': {
                    'marker': { 'opacity': 0.3 },
                    # make text transparent when not selected
                    'textfont': { 'color': 'rgba(0, 0, 0, 0)' }
                }
            }],
            'layout': {
                "width":panelWidth,
                "height":panelWidth,
    #             'margin': {'l': 20, 'r': 0, 'b': 15, 't': 5},
                'dragmode': 'select',
                'hovermode': False,
                "xaxis":{"title":x_col},
                "yaxis":{"title":y_col},
                "margin":dict(zip("tblr",[50]*3+[10])),
                # Display a rectangle to highlight the previously selected region
                'shapes': [dict({
                    'type': 'rect',
                    'line': { 'width': 1, 'dash': 'dot', 'color': 'darkgrey' }
                }, **selection_bounds
                )]
            }
        }
    N = df.shape[1]
    graphIDs = OrderedDict([('g_%i_%i'%(i,j),(i,j)) for i in range(N) for j in range(i+1,N)])
    app = JupyterDash(__name__,width=appWidth)
    app.layout = html.Div([
        html.Div(dcc.Graph(id=grID,figure=baseFig),style={"display":"inline-block"})
            for grID in graphIDs]+[
        html.Div(id="output",children="",style={
            "font-family":"Courier New",
            "font-size":"80%",
            'border': 'thin lightgrey solid',
            "width":"80%",
            "height":"100px",
            'overflowX': 'scroll',
            'overflowY': 'scroll',})
    ])

    OUTPUTS = [Output(grID, 'figure') for grID in graphIDs]+[Output("output","children")]
    INPUTS  = [Input (grID, 'selectedData') for grID in graphIDs]
    @app.callback(OUTPUTS,INPUTS)
    def callback(*args):
        from sys import exc_info
        try:
            global SELECTED_POINTS
            from numpy import intersect1d
            selectedpoints = df.index
            for selected_data in args:
                if selected_data and selected_data['points']:
                    selectedpoints = intersect1d(selectedpoints,
                        [p['customdata'] for p in selected_data['points']])
            SELECTED_POINTS = selectedpoints
            outputs = [get_figure(df, graphIDs[grID][0], graphIDs[grID][1], selectedpoints, selection) for grID,selection in zip(graphIDs,args)]
            feedback = "selected: "+", ".join([str(p) for p in selectedpoints])
            return outputs+[feedback]
        except:
            return [baseFig]*len(graphIDs)+["ERROR:"+str(exc_info())]


    return app