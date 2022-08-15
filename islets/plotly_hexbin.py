import json

try:
    from dash import dcc, html
except ImportError:
    import dash_core_components as dcc
    import dash_html_components as html
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from dash.dependencies import Input, Output
from jupyter_plotly_dash import JupyterDash as  Dash
from scipy.spatial import distance_matrix


# def mpl_to_plotly(cmap, N):
#     h = 1.0/(N-1)
#     pl_colorscale = []
#     for k in range(N):
#         C = list(map(np.uint8, np.array(cmap(k*h)[:3])*255))
#         pl_colorscale.append([round(k*h,2), f'rgb({C[0]}, {C[1]}, {C[2]})'])
#     return pl_colorscale

def make_hexagon(prototypical_hex, offset, fillcolor, linecolor="grey"):   
    new_hex_vertices = [vertex + offset for vertex in prototypical_hex]
    vertices = np.asarray(new_hex_vertices[:-1])
    # hexagon center
    center=np.mean(vertices, axis=0)
    if linecolor is None:
        linecolor = fillcolor
    #define the SVG-type path:    
    path = 'M '
    for vert in new_hex_vertices:
        path +=  f'{vert[0]}, {vert[1]} L' 
    return  dict(type='path',
                 line=dict(color=linecolor, 
                           width=0.5),
                 path=  path[:-2],
                 fillcolor=fillcolor, 
                ), center 

def get_hexbin_attributes(hb):
    out = hb.get_offsets(), hb.get_facecolors(), hb.get_array()
    paths = hb.get_paths()
    points_codes = list(paths[0].iter_segments())
    prototypical_hexagon = [item[0] for item in points_codes]
    return (np.array(prototypical_hexagon), )+ out

# def pl_cell_color(mpl_facecolors):
#     return [ f'rgb({int(R*255)}, {int(G*255)}, {int(B*255)})' for (R, G, B, A) in mpl_facecolors]

class HBclass:
    def __init__(self, df, x, y):
        self.df = df#[[x,y]].copy()
    
#     @classmethod
    def _hexbin(self, df=None, x=None, y=None,
               debug=False,
               log_x=False,
               log_y=False,
               sizeScaling = 1,
               return_figure=False,
               **hexbin_kwargs
        ):
        if df is not None:
            self.df = df[[x,y]].copy()
        if x is None:
            x = self.df.columns[0]
        if y is None:
            y = self.df.columns[1]
        hexbin_kwargs.update(mincnt=1)
        hexbin_kwargs["gridsize"] = hexbin_kwargs.get("gridsize",(50,50))
        hexbin_kwargs["bins"] = hexbin_kwargs.get("bins","log")

        if log_x or log_y: raise NotImplementedError()
        if log_x:
            xshow = x+"_logged"
            df[xshow] = np.log10(df[x])
        else:
            xshow = x

        if log_y:
            yshow = y+"_logged"
            df[yshow] = np.log10(df[y])
        else:
            yshow = y
        currentBackend = matplotlib.get_backend()
        plt.switch_backend('Agg')
        fig = plt.figure(figsize=(0.1,0.1))
        HB = plt.hexbin(self.df[xshow], self.df[yshow], **hexbin_kwargs)
        plt.show()
        plt.switch_backend(currentBackend)

        hexagon_vertices, offsets, mpl_facecolors, counts = get_hexbin_attributes(HB)
        shapes = []
        centers = []
        for k in range(len(offsets)):
            if counts[k]:
                shape, center = make_hexagon(
                    hexagon_vertices,
                    offsets[k],
                    "blue"#cell_color[k]
                )
                shapes.append(shape)
                centers.append(center)

        # plcmap = mpl_to_plotly(plt.cm.hot, 101)

        X, Y = zip(*centers)
        xrescale = max(X)-min(X)
        yrescale = max(Y)-min(Y)
        centers = np.array(centers)
        rescaled_centers = np.transpose([centers.T[0]/xrescale, centers.T[1]/yrescale])
        rescaled_points = np.transpose([self.df[xshow].values/xrescale, self.df[yshow].values/yrescale])
        self.df["bin_id"] = np.argmin(distance_matrix(rescaled_points,rescaled_centers),1)
        text = [f'x: {round(X[k],2)}<br>y: {round(Y[k],2)}<br>counts: {int(counts[k])}<br>bin id: {k}' for k in range(len(X))]
        maxgsize = hexbin_kwargs.get("gridsize",100)
        if not isinstance(maxgsize, int):
            maxgsize = max(maxgsize)
        
        trace = go.Scatter(
                     x=list(X), 
                     y=list(Y), 
                     mode='markers',
                     showlegend=False,
                     marker=dict(size=350/maxgsize*sizeScaling, 
                                 # symbol="hexagon",
                                 color=np.log(1+counts) if hexbin_kwargs.get("bins")=="log" else counts, 
                                 colorscale="hot",
                                #  showscale=True,
                                #  colorbar=dict(
                                #              thickness=20,  
                                #              ticklen=4
                                #              )
                                ),
                   text=text, 
                   hoverinfo='text'
                  )

        layout = go.Layout(
                           width=550, height=500,
                           xaxis=dict(title=self.df.columns[0],showgrid=True),
                           yaxis=dict(title=self.df.columns[1]),
                           hovermode='closest',
                           dragmode='lasso',
                          )
        fig = go.Figure(data=[trace], layout=layout)
        fig.update_layout({"margin":dict(l=10, r=10, t=20, b=40)})
        if return_figure:
            return fig

        app = Dash(__name__, width=1000, height=1000)
        app.layout = html.Div([
            dcc.Graph(figure=fig, id="hexbin-graph"),
            html.Div(id="hexbin-selected",
                     children="Nothing yet",
                     style={"display":"block" if debug else "none"}
                    )
        ])
        @app.callback(
            Output("hexbin-selected", "children"),
            [Input("hexbin-graph","selectedData")]
        )
        def hexbin_selector(selData):
            out = []
            try:
                out = [html.Pre(
                        json.dumps(selData["points"], indent=4),
                         style={
                             "overflowX":"scroll",
                             "overflowY":"scroll",
                             "border": "thin navy solid",
                             "width":"30%",
                             "height":"200px",
                             "display":"inline-block"
                               }
                       )]
                out += [html.Pre(
                        json.dumps([pt["pointIndex"] for pt in selData["points"]], indent=4),
                         style={
                             "overflowX":"scroll",
                             "overflowY":"scroll",
                             "border": "thin green solid",
                             "width":"30%",
                             "height":"200px",
                             "display":"inline-block"
                               }
                       )]
                subdf = self.df[self.df.bin_id.isin([pt["pointIndex"] for pt in selData["points"]])]
                subdf = list(subdf.index)
                out += [html.Pre(
                        json.dumps(subdf, indent=4),
                        id="hexbin-selected-indices",
                         style={
                             "overflowX":"scroll",
                             "overflowY":"scroll",
                             "border": "thin red solid",
                             "width":"30%",
                             "height":"200px",
                             "display":"inline-block"
                               }
                       )]
                return out
            except:
                out +=["Oops"]
            return out
        return app
    
def hexbin(df=None, x=None, y=None,
           debug=False,
           log_x=False,
           log_y=False,
           sizeScaling = 1,
           **hexbin_kwargs
    ):
    hb_ = HBclass(df,x,y,)
    return hb_._hexbin(
               debug=debug,
               log_x=log_x,
               log_y=log_y,
               sizeScaling=sizeScaling,
               **hexbin_kwargs,
    )