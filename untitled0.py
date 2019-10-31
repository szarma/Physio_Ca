### Preface
import numpy as np
import pandas as pd
from os.path import isfile, isdir
import matplotlib.pyplot as plt

from sys import path as syspath
syspath.append("./functions/")
from physio_def_1 import getApparentFreq, importFrames, getTimes

from plotFirst_1 import plotImage

import pandas as pd

from collections import OrderedDict

npzFile = "/Volumes/physio/team/slakrupnik/project/experiments/NIKON/20190806/512X512_image_slice3_16la_2ca.npz"

isfile(npzFile)

datafile = np.load(npzFile)

datafile.files

metadata = pd.read_csv(npzFile.replace("npz","txt"), index_col=0).iloc[0]
metadata