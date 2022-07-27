import numpy as np
import pandas as pd
from collections import OrderedDict
import matplotlib.pyplot as plt
import plotly.graph_objects as go


# noinspection PyUnresolvedReferences
def climb(x, a, diag=True, min_gradient=0, verbose=False):
    dims = a.shape
    xs = [x]
    v = a[tuple(x)]
    for i in range(100):
        ijs = []
        vs = []
        for dij in product(*([[-1,0,1]]*len(dims))):
            if not diag:
                if np.prod(dij)!=0: continue
            ij = x + np.array(dij)
            if any(
                [ij[k]<0 for k in range(len(dims))]+
                [ij[k]>=dims[k] for k in range(len(dims))]
            ):
                continue
            vs += [a[tuple(ij)]]
            ijs += [ij]
        x1 = ijs[np.argmax(vs)]#[:-1]
        dv = max(vs)-v
        if dv<=0 or dv<=min_gradient:
            break
        else:
            x = x1
            v = max(vs)
            xs += [x]
            if verbose:
                print ("=>", xs)
    return tuple(xs[-1])


from .Regions import Regions

class Regions3D(Regions):

    def other_function(self):
        return "exists"
