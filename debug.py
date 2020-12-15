#%% Imports
from copy import copy
import islets
import json
import mgzip
import gzip
import numpy as np

#%% Load test regions
regions = islets.load_regions('/home/johannes/Downloads/5.6_rois.pkl')

#%% Function implementation:
from copy import deepcopy
import pandas as pd

col = ['trace']

saving = ['statImages',
          'mode',
          'image',
          'filterSize',
          'df',
          'trange',
          "FrameRange",
          "analysisFolder",
          "time",
          "Freq",
          "metadata"]

juggleMovie = hasattr(regions, "movie")
if juggleMovie:
    movie = regions.movie
    del regions.movie
allAttrs = list(regions.__dict__.keys())
subRegions = deepcopy(regions)
if juggleMovie:
    regions.movie = movie
    del movie
for k in allAttrs:
    if k not in saving:
        del subRegions.__dict__[k]
for k in regions.df.columns:
    if k not in ["peak", "pixels", "peakValue", "tag", "interest"] + col:
        del subRegions.df[k]

json_dict = {}
for k, v in subRegions.__dict__.items():
    value = v
    if isinstance(v, (pd.DataFrame, pd.Series)):
        value = json.JSONDecoder().decode(v.to_json(double_precision=15))
    if isinstance(v, (np.float64, np.int64)):
        value = str(v)
    if isinstance(v, np.ndarray):
        value = v.tolist()
    if isinstance(v, dict):
        for k_, v_ in v.items():
            if isinstance(v_, np.ndarray):
                # v[k_] = json.dumps(v_.tolist())
                v[k_] = v_.tolist()
        value = v

    json_dict[k] = value

json_string = json.dumps(json_dict)

#%% Save file
with mgzip.open('/home/johannes/Downloads/test.json.gz', 'wt') as file:
    file.write(json_string)


#%% Open file
with mgzip.open('/home/johannes/Downloads/test.json.gz', 'rt') as file:
    s = file.read()

#%% Load JSON-dict
json_obj = json.loads(s)

#%% Recreate Regions:
stat_images = copy(json_obj['statImages'])
df = pd.DataFrame.from_dict(json_obj['df'])
df['peak'] = df['peak'].apply(lambda x: tuple(x))
df['pixels'] = df['pixels'].apply(lambda x: [tuple(y) for y in x])
metadata = pd.Series(json_obj['metadata'])

for k, v in stat_images.items():
    stat_images[k] = np.array(v)




# %% Create regions
regions = islets.Regions(dict(zip(df['peak'],
                                  df['pixels'])))
