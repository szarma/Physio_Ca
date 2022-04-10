### functions


```python
import pandas as pd
import numpy as np
np.corrcoef([0,1], [0,2])
# import bioformats as bf
import os
import islets
import matplotlib.pyplot as plt
%config InlineBackend.figure_format = 'retina'
from scipy.stats import distributions as dst
plt.rcParams["animation.embed_limit"] = 200
from copy import deepcopy
from scipy.stats import median_abs_deviation
```

### Import pickle
__NB__: in this notebook you __do not__ merge, discard, or change rois in any way.


```python
pathToRoi = "/data/Sandra/2020/2020_11_05/Experiment74b.lif_analysis/Series004-16/2021_11_08_9_srdjan_rois.pkl"
regions = islets.load_regions(pathToRoi)
```


```python
regions.get_activity([10,100])
```


```python
plt.figure()
regions.plotEdges()
```


```python
regions.examine()
```

### Extract events


```python
rawEvents = islets.EventDistillery.sequential_filtering(regions,)
```


```python
ExampleRoi = 187
```


```python
singleRaw = rawEvents.query(f"roi=={ExampleRoi}").copy()
```


```python
islets.EventDistillery.plot_events(singleRaw,
                                   modify=True # to add color to events
                                  )
```


```python
%matplotlib notebook
```


```python
distilled, axs, origs = islets.EventDistillery.distill_events_per_roi(singleRaw, regions, plot=True)
```

This is where you inspect if things make sense. In most cases they should and you need to continue to distilling with default parameters. 

In this case, the default way may not be the best way.

In particular, it looks as like we cannot rely on the `EventDistillery` to catch the big jump(s) around 3300s.
Luckily, by looking at the examiner, it seems that there are not many traces with these jumps.

I will choose to perform default distillation up to a timescale of ~128s, and deal with large events in a different way.


```python
Events = islets.EventDistillery.distill_events(rawEvents.query("ts<=128"), regions, )
```


```python
Events
```

### Make sure it makes sense


```python
Events["log10 halfwidth"] = np.log10(Events.halfwidth)
Events["log10 z_max"    ] = np.log10(Events.z_max    )
Events["log10 height"   ] = np.log10(Events.height   )
```


```python
Events.columns
```


```python
regions.examine_events(Events,y="log10 z_max", x="log10 halfwidth")
```


```python
Events = Events.query("z_max>4")
```

### Save


```python
for col in Events.columns:
    if "log10" in col:
        del Events[col]
```


```python
eventSave = pathToRoi.split("_rois")[0]+"_events.csv"
eventSave
```


```python
Events.to_csv(eventSave)
```
