```python
# Importing 
%load_ext autoreload
%autoreload 2
import pandas as pd
import numpy as np
import os
import islets
import matplotlib.pyplot as plt
%config InlineBackend.figure_format = 'retina'
import islets._manuscript_functions as mf
import statsmodels.api as sm
```

## Import


```python
pathToRegions = '/data/Sandra/2021/2021_04_27/Experiment112b.lif_analysis/Series001-4/2022_01_06_5_rois.pkl'
regions = islets.load_regions(pathToRegions)
```


```python
Events = pd.read_csv('/data/Sandra/2021/2021_04_27/Experiment112b.lif_analysis/Series001-4/2022_01_06_5_auto_events.csv', index_col=0)
```

There is a submodule called `_manuscript_functions`, I imported it at the beginning as `mf`

Let's plot a hexbin of time vs halfwidth


```python
fig, axs = mf.plot_events(Events, regions)
```

    /home/srdjan/github/Physio_Ca/islets/protocol.py:121: MatplotlibDeprecationWarning: Adding an axes using the same arguments as a previous axes currently reuses the earlier instance.  In a future version, a new instance will always be created and returned.  Meanwhile, this warning can be suppressed, and the future behavior ensured, by passing a unique label to each axes instance.
      axp = fig.add_axes([pos.x0, pos.y0+pos.height*1.02, pos.width, axHeight], label="protocol")



![png](temp_files/temp_6_1.png)


If you use this function on a shorter experiment, you will see that the figure size changes depending on the duration.

There is also a function which produces a figure with a bit more information content.

## Basic manipulations of hexbin

Adding color to your protocol plot:


```python
# make a copy of the original, as we do not want mess it up
protocol = regions.protocol.copy()
protocol
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>compound</th>
      <th>concentration</th>
      <th>begin</th>
      <th>end</th>
      <th>t_begin</th>
      <th>t_end</th>
      <th>conc</th>
      <th>unit</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>glucose</td>
      <td>8mM</td>
      <td>00:00</td>
      <td>20:00</td>
      <td>0.0</td>
      <td>1200.000000</td>
      <td>8</td>
      <td>mM</td>
    </tr>
    <tr>
      <th>1</th>
      <td>glucose</td>
      <td>6mM</td>
      <td>20:00</td>
      <td>28:58</td>
      <td>1200.0</td>
      <td>1738.000000</td>
      <td>6</td>
      <td>mM</td>
    </tr>
    <tr>
      <th>2</th>
      <td>glucose</td>
      <td>8mM</td>
      <td>28:58</td>
      <td>48:58</td>
      <td>1738.0</td>
      <td>2938.000000</td>
      <td>8</td>
      <td>mM</td>
    </tr>
    <tr>
      <th>3</th>
      <td>glucose</td>
      <td>6mM</td>
      <td>48:58</td>
      <td>57:57</td>
      <td>2938.0</td>
      <td>3477.000000</td>
      <td>6</td>
      <td>mM</td>
    </tr>
    <tr>
      <th>4</th>
      <td>glucose</td>
      <td>8mM</td>
      <td>57:57</td>
      <td>77:57</td>
      <td>3477.0</td>
      <td>4677.000000</td>
      <td>8</td>
      <td>mM</td>
    </tr>
    <tr>
      <th>5</th>
      <td>glucose</td>
      <td>6mM</td>
      <td>77:57</td>
      <td>NaN</td>
      <td>4677.0</td>
      <td>5133.939394</td>
      <td>6</td>
      <td>mM</td>
    </tr>
    <tr>
      <th>6</th>
      <td>pH</td>
      <td>7.4</td>
      <td>00:00</td>
      <td>28:58</td>
      <td>0.0</td>
      <td>1738.000000</td>
      <td>7.4</td>
      <td></td>
    </tr>
    <tr>
      <th>7</th>
      <td>pH</td>
      <td>7.1</td>
      <td>28:58</td>
      <td>48:58</td>
      <td>1738.0</td>
      <td>2938.000000</td>
      <td>7.1</td>
      <td></td>
    </tr>
    <tr>
      <th>8</th>
      <td>pH</td>
      <td>7.4</td>
      <td>48:58</td>
      <td>NaN</td>
      <td>2938.0</td>
      <td>5133.939394</td>
      <td>7.4</td>
      <td></td>
    </tr>
  </tbody>
</table>
</div>




```python
# add color
protocol["color"] = [
    "lightgrey",
    "whitesmoke",
    "lightgrey",
    "whitesmoke",
    "lightgrey",
    "whitesmoke",
    "whitesmoke",
    "lightsalmon",
    "whitesmoke",
]
protocol
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>compound</th>
      <th>concentration</th>
      <th>begin</th>
      <th>end</th>
      <th>t_begin</th>
      <th>t_end</th>
      <th>conc</th>
      <th>unit</th>
      <th>color</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>glucose</td>
      <td>8mM</td>
      <td>00:00</td>
      <td>20:00</td>
      <td>0.0</td>
      <td>1200.000000</td>
      <td>8</td>
      <td>mM</td>
      <td>lightgrey</td>
    </tr>
    <tr>
      <th>1</th>
      <td>glucose</td>
      <td>6mM</td>
      <td>20:00</td>
      <td>28:58</td>
      <td>1200.0</td>
      <td>1738.000000</td>
      <td>6</td>
      <td>mM</td>
      <td>whitesmoke</td>
    </tr>
    <tr>
      <th>2</th>
      <td>glucose</td>
      <td>8mM</td>
      <td>28:58</td>
      <td>48:58</td>
      <td>1738.0</td>
      <td>2938.000000</td>
      <td>8</td>
      <td>mM</td>
      <td>lightgrey</td>
    </tr>
    <tr>
      <th>3</th>
      <td>glucose</td>
      <td>6mM</td>
      <td>48:58</td>
      <td>57:57</td>
      <td>2938.0</td>
      <td>3477.000000</td>
      <td>6</td>
      <td>mM</td>
      <td>whitesmoke</td>
    </tr>
    <tr>
      <th>4</th>
      <td>glucose</td>
      <td>8mM</td>
      <td>57:57</td>
      <td>77:57</td>
      <td>3477.0</td>
      <td>4677.000000</td>
      <td>8</td>
      <td>mM</td>
      <td>lightgrey</td>
    </tr>
    <tr>
      <th>5</th>
      <td>glucose</td>
      <td>6mM</td>
      <td>77:57</td>
      <td>NaN</td>
      <td>4677.0</td>
      <td>5133.939394</td>
      <td>6</td>
      <td>mM</td>
      <td>whitesmoke</td>
    </tr>
    <tr>
      <th>6</th>
      <td>pH</td>
      <td>7.4</td>
      <td>00:00</td>
      <td>28:58</td>
      <td>0.0</td>
      <td>1738.000000</td>
      <td>7.4</td>
      <td></td>
      <td>whitesmoke</td>
    </tr>
    <tr>
      <th>7</th>
      <td>pH</td>
      <td>7.1</td>
      <td>28:58</td>
      <td>48:58</td>
      <td>1738.0</td>
      <td>2938.000000</td>
      <td>7.1</td>
      <td></td>
      <td>lightsalmon</td>
    </tr>
    <tr>
      <th>8</th>
      <td>pH</td>
      <td>7.4</td>
      <td>48:58</td>
      <td>NaN</td>
      <td>2938.0</td>
      <td>5133.939394</td>
      <td>7.4</td>
      <td></td>
      <td>whitesmoke</td>
    </tr>
  </tbody>
</table>
</div>



Many named colors you can find here
https://matplotlib.org/2.0.2/examples/color/named_colors.html


```python
# And now use this instance in plotting
fig, axs = mf.plot_events(Events, regions, protocol=protocol)
```

    /home/srdjan/github/Physio_Ca/islets/protocol.py:121: MatplotlibDeprecationWarning: Adding an axes using the same arguments as a previous axes currently reuses the earlier instance.  In a future version, a new instance will always be created and returned.  Meanwhile, this warning can be suppressed, and the future behavior ensured, by passing a unique label to each axes instance.
      axp = fig.add_axes([pos.x0, pos.y0+pos.height*1.02, pos.width, axHeight], label="protocol")



![png](temp_files/temp_14_1.png)


## Define legs

Say now we also define legs in this experiment

To guide you, you can have a look a the following numbers:


```python
regions.protocol.legs
```




    {(0.0, 1200.0): {'glucose': '8mM', 'ph': '7.4'},
     (1200.0, 1738.0): {'glucose': '6mM', 'ph': '7.4'},
     (1738.0, 2938.0): {'glucose': '8mM', 'ph': '7.1'},
     (2938.0, 3477.0): {'glucose': '6mM', 'ph': '7.4'},
     (3477.0, 4677.0): {'glucose': '8mM', 'ph': '7.4'},
     (4677.0, 5133.939393939394): {'glucose': '6mM', 'ph': '7.4'}}




```python
legDict = {
    "leg0": ( 500,1300),
    "leg1": (2500,3000),
    "leg2": (4000,4900),
}
# I intentially made them unequal durations
```

lets also define the colors of these legs


```python
legColorDict = {
    "leg0": "blue",
    "leg1": "red",
    "leg2": "darkgoldenrod",
}
```


```python
fig, axs = mf.plot_events(Events, regions, protocol=protocol)
axhex = axs[0]
for leg in legDict:
    mf.emphasize_region(axhex, legDict[leg],[.4,200],color=legColorDict[leg])
```

    /home/srdjan/github/Physio_Ca/islets/protocol.py:121: MatplotlibDeprecationWarning: Adding an axes using the same arguments as a previous axes currently reuses the earlier instance.  In a future version, a new instance will always be created and returned.  Meanwhile, this warning can be suppressed, and the future behavior ensured, by passing a unique label to each axes instance.
      axp = fig.add_axes([pos.x0, pos.y0+pos.height*1.02, pos.width, axHeight], label="protocol")



![png](temp_files/temp_22_1.png)


Let's annotate Events with these legs.


```python
islets.EventDistillery.define_legs(Events, legDict)
```

We can now plot only defined events:


```python
fig, axs = mf.plot_events(Events, regions, protocol=protocol, plottype="scatter", legColorDict=legColorDict)
```

    /home/srdjan/github/Physio_Ca/islets/protocol.py:121: MatplotlibDeprecationWarning: Adding an axes using the same arguments as a previous axes currently reuses the earlier instance.  In a future version, a new instance will always be created and returned.  Meanwhile, this warning can be suppressed, and the future behavior ensured, by passing a unique label to each axes instance.
      axp = fig.add_axes([pos.x0, pos.y0+pos.height*1.02, pos.width, axHeight], label="protocol")



![png](temp_files/temp_26_1.png)


When the events have legs, you can also use a function that outputs a figure with more information content:


```python
mf.big_plot(regions, Events, )
```

    100%|██████████| 1/1 [00:03<00:00,  3.77s/it]


    There are 244 rois with more than 10 events (72%).
    Out of them, 122 are sampled 20 times to estimate mean and std of the firing rate.



![png](temp_files/temp_28_2.png)


You might also wish to see the cumulative density function (see explanation below), or the average trace of all rois (raw, unfiltered).


```python
mf.big_plot(regions, Events, cdf=True, plot_sum=True)
```

    There are 244 rois with more than 10 events (72%).
    Out of them, 122 are sampled 20 times to estimate mean and std of the firing rate.



![png](temp_files/temp_30_1.png)


### basic plots of legs


```python
# We can also have their histograms:
hwBinEdges = np.geomspace(Events.halfwidth.min(), Events.halfwidth.max())
hwBinCenters = (hwBinEdges[:-1]*hwBinEdges[1:])**.5

ax = plt.subplot(111)
for leg in legDict:
    ev = Events.query(f"leg=='{leg}'")
    ax.hist(ev["halfwidth"],hwBinEdges, color=legColorDict[leg], label=leg, histtype="step")
ax.set_xscale("log")
ax.legend(loc = (1.01,0.01))
```




    <matplotlib.legend.Legend at 0x7fe7a2fd5da0>




![png](temp_files/temp_32_1.png)



```python
# But, I like to smooth them
ax = plt.subplot(111)
for leg in legDict:
    ev = Events.query(f"leg=='{leg}'")
    h = np.histogram(ev['halfwidth'],hwBinEdges)[0]
    ax.plot(hwBinCenters, h, color=legColorDict[leg], label=leg)
ax.set_xscale("log")
ax.legend(loc = (1.01,0.01))
```




    <matplotlib.legend.Legend at 0x7fe7a22fda58>




![png](temp_files/temp_33_1.png)



```python
# or their CDF
ax = plt.subplot(111)
for leg in legDict:
    ev = Events.query(f"leg=='{leg}'")
    ax.hist(ev["halfwidth"],hwBinEdges, color=legColorDict[leg], label=leg, histtype="step", cumulative=True, density=True)
ax.set_xscale("log")
ax.legend(loc = (1.01,0.01))
```




    <matplotlib.legend.Legend at 0x7fe7a2fcea20>




![png](temp_files/temp_34_1.png)



```python
# and full CDF 
ax = plt.subplot(111)
for leg in legDict:
    ev = Events.query(f"leg=='{leg}'")
    x = ev["halfwidth"].values.copy()
    x.sort()
    ax.plot(x, np.linspace(0,1,len(x)), color=legColorDict[leg], label=leg)
ax.set_xscale("log")
ax.legend(loc = (1.01,0.01))
```




    <matplotlib.legend.Legend at 0x7fe7a01c4a20>




![png](temp_files/temp_35_1.png)


Looks like the median halfwidth increases in `leg1` and then decreases in `leg2`. Let's investigate this.

## Investigate halfwidths
say we are interested only in "bursting", i.e. events shorter than say 10 seconds.

### Intro


```python
events_short = Events.query("halfwidth<10")
```


```python
# and full CDF 
ax = plt.subplot(111)
for leg in legDict:
    ev = events_short.query(f"leg=='{leg}'")
    x = ev["halfwidth"].values.copy()
    x.sort()
    ax.plot(x, np.linspace(0,1,len(x)), color=legColorDict[leg], label=leg)
ax.set_xscale("log")
ax.legend(loc = (1.01,0.01))
```




    <matplotlib.legend.Legend at 0x7fe7a3a7fe48>




![png](temp_files/temp_40_1.png)


We could also plot boxplots


```python
events_short.boxplot("halfwidth", by="leg", showfliers=False, whis = (2.5,97.5))
plt.grid(False)
```


![png](temp_files/temp_42_0.png)


Looks like the events prolong in the leg1 and shrink back a bit in leg2
Let's try to quantify this effect.


```python
result = mf.get_hw_effect( events_short, ref_leg="leg0", control_for_roi=False )
# ref_leg is reference leg, the leg which we want to compare against.
```


<style type="text/css">
            /* Fix details summary arrow
               not shown in Firefox
               due to bootstrap
               display: block;
             */
            summary {
                display: list-item;
            }
            </style><div>reference leg: leg0</div><table><tr><td style="text-align:left">Obtained effect for leg '<span style='font-family:monospace'>leg1</span>':</td><td>+10%,</td><td>with the p-value of 2e-31</td></tr><tr><td style="text-align:left">Obtained effect for leg '<span style='font-family:monospace'>leg2</span>':</td><td>+7%,</td><td>with the p-value of 1.9e-19</td></tr></table><details><summary style='color:navy'>full output</summary><table class="simpletable">
<caption>OLS Regression Results</caption>
<tr>
  <th>Dep. Variable:</th>        <td>log10_hw</td>     <th>  R-squared:         </th> <td>   0.005</td>
</tr>
<tr>
  <th>Model:</th>                   <td>OLS</td>       <th>  Adj. R-squared:    </th> <td>   0.005</td>
</tr>
<tr>
  <th>Method:</th>             <td>Least Squares</td>  <th>  F-statistic:       </th> <td>   76.86</td>
</tr>
<tr>
  <th>Date:</th>             <td>Tue, 25 Jan 2022</td> <th>  Prob (F-statistic):</th> <td>5.03e-34</td>
</tr>
<tr>
  <th>Time:</th>                 <td>17:09:57</td>     <th>  Log-Likelihood:    </th> <td> -825.97</td>
</tr>
<tr>
  <th>No. Observations:</th>      <td> 30847</td>      <th>  AIC:               </th> <td>   1658.</td>
</tr>
<tr>
  <th>Df Residuals:</th>          <td> 30844</td>      <th>  BIC:               </th> <td>   1683.</td>
</tr>
<tr>
  <th>Df Model:</th>              <td>     2</td>      <th>                     </th>     <td> </td>   
</tr>
<tr>
  <th>Covariance Type:</th>      <td>nonrobust</td>    <th>                     </th>     <td> </td>   
</tr>
</table>
<table class="simpletable">
<tr>
                  <td></td>                     <th>coef</th>     <th>std err</th>      <th>t</th>      <th>P>|t|</th>  <th>[0.025</th>    <th>0.975]</th>  
</tr>
<tr>
  <th>Intercept</th>                         <td>    0.3664</td> <td>    0.002</td> <td>  157.011</td> <td> 0.000</td> <td>    0.362</td> <td>    0.371</td>
</tr>
<tr>
  <th>C(leg, Treatment('leg0'))[T.leg1]</th> <td>    0.0427</td> <td>    0.004</td> <td>   11.676</td> <td> 0.000</td> <td>    0.036</td> <td>    0.050</td>
</tr>
<tr>
  <th>C(leg, Treatment('leg0'))[T.leg2]</th> <td>    0.0296</td> <td>    0.003</td> <td>    9.025</td> <td> 0.000</td> <td>    0.023</td> <td>    0.036</td>
</tr>
</table>
<table class="simpletable">
<tr>
  <th>Omnibus:</th>       <td>1495.897</td> <th>  Durbin-Watson:     </th> <td>   0.694</td>
</tr>
<tr>
  <th>Prob(Omnibus):</th>  <td> 0.000</td>  <th>  Jarque-Bera (JB):  </th> <td>1720.585</td>
</tr>
<tr>
  <th>Skew:</th>           <td>-0.570</td>  <th>  Prob(JB):          </th> <td>    0.00</td>
</tr>
<tr>
  <th>Kurtosis:</th>       <td> 3.201</td>  <th>  Cond. No.          </th> <td>    3.63</td>
</tr>
</table><br/><br/>Notes:<br/>[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.</details>


Click on "full output" to see... well, full output. :-)  
Note that the values reported in full output are not the same as in the summary table above. Given a linear coefficient for the logged halfwidth $a$ (in full output), a percent change effect is $(10^a-1)*100\%$. You can use this formula if you wish to report confidence intervals from the full output.

### Controlling for the effect of roi

A devil might say it is possible that a few misbehaving rois with large halfwidths were, for some reason, not firing during the first leg, but only in the later ones, and that the whole effect is only due to presence of such rois. In fact, the devil might suspect all other rois shortened and the effect  we see is solely due to the presence of the misbehaing rois. To address that, we need to asses the effect for each roi individually and then "average" over them. For that we use mixed linear models, but I created a function that takes care of this.  
We also might not care about rois that have less than, say, 10 events in each roi.


```python
result = mf.get_hw_effect(events_short, ref_leg="leg0", control_for_roi=True, minEvents=10)
```

    There are 131 rois, which have at least 10 events in all legs.



<style type="text/css">
            /* Fix details summary arrow
               not shown in Firefox
               due to bootstrap
               display: block;
             */
            summary {
                display: list-item;
            }
            </style><div>reference leg: leg0</div><table><tr><td style="text-align:left">Obtained effect for leg '<span style='font-family:monospace'>leg1</span>':</td><td>+13%,</td><td>with the p-value of 4.3e-45</td></tr><tr><td style="text-align:left">Obtained effect for leg '<span style='font-family:monospace'>leg2</span>':</td><td>+2.7%,</td><td>with the p-value of 0.00052</td></tr></table><span style='font-size:70%'>~ 8 warning(s) encountered when fitting you might want to check full output ~</span><details><summary style='color:navy'>full output</summary><table class="simpletable">
<tr>
       <td>Model:</td>       <td>MixedLM</td> <td>Dependent Variable:</td> <td>log10_hw</td>
</tr>
<tr>
  <td>No. Observations:</td>  <td>26892</td>        <td>Method:</td>         <td>REML</td>  
</tr>
<tr>
     <td>No. Groups:</td>      <td>131</td>         <td>Scale:</td>         <td>0.0538</td> 
</tr>
<tr>
  <td>Min. group size:</td>    <td>55</td>      <td>Log-Likelihood:</td>   <td>898.7931</td>
</tr>
<tr>
  <td>Max. group size:</td>    <td>473</td>       <td>Converged:</td>         <td>No</td>   
</tr>
<tr>
  <td>Mean group size:</td>   <td>205.3</td>           <td></td>               <td></td>    
</tr>
</table>
<table class="simpletable">
<tr>
                  <td></td>                  <th>Coef.</th> <th>Std.Err.</th>    <th>z</th>   <th>P>|z|</th> <th>[0.025</th> <th>0.975]</th>
</tr>
<tr>
  <th>Intercept</th>                         <td>0.383</td>   <td>0.009</td>  <td>44.441</td> <td>0.000</td>  <td>0.366</td>  <td>0.400</td>
</tr>
<tr>
  <th>C(leg, Treatment('leg0'))[T.leg1]</th> <td>0.052</td>   <td>0.004</td>  <td>14.091</td> <td>0.000</td>  <td>0.045</td>  <td>0.059</td>
</tr>
<tr>
  <th>C(leg, Treatment('leg0'))[T.leg2]</th> <td>0.011</td>   <td>0.003</td>   <td>3.471</td> <td>0.001</td>  <td>0.005</td>  <td>0.018</td>
</tr>
<tr>
  <th>roi Var</th>                           <td>0.009</td>   <td>0.006</td>     <td></td>      <td></td>       <td></td>       <td></td>   
</tr>
</table></details>


The effects stay almost the same even if we include all rois.


```python
result = mf.get_hw_effect(events_short, ref_leg="leg0", control_for_roi=True, minEvents=0)
```

    There are 310 rois, which have at least 0 events in all legs.



<style type="text/css">
            /* Fix details summary arrow
               not shown in Firefox
               due to bootstrap
               display: block;
             */
            summary {
                display: list-item;
            }
            </style><div>reference leg: leg0</div><table><tr><td style="text-align:left">Obtained effect for leg '<span style='font-family:monospace'>leg1</span>':</td><td>+13%,</td><td>with the p-value of 4e-54</td></tr><tr><td style="text-align:left">Obtained effect for leg '<span style='font-family:monospace'>leg2</span>':</td><td>+3.3%,</td><td>with the p-value of 5.4e-06</td></tr></table><details><summary style='color:navy'>full output</summary><table class="simpletable">
<tr>
       <td>Model:</td>       <td>MixedLM</td> <td>Dependent Variable:</td> <td>log10_hw</td> 
</tr>
<tr>
  <td>No. Observations:</td>  <td>30847</td>        <td>Method:</td>         <td>REML</td>   
</tr>
<tr>
     <td>No. Groups:</td>      <td>289</td>         <td>Scale:</td>         <td>0.0526</td>  
</tr>
<tr>
  <td>Min. group size:</td>     <td>1</td>      <td>Log-Likelihood:</td>   <td>1266.2863</td>
</tr>
<tr>
  <td>Max. group size:</td>    <td>473</td>       <td>Converged:</td>         <td>Yes</td>   
</tr>
<tr>
  <td>Mean group size:</td>   <td>106.7</td>           <td></td>               <td></td>     
</tr>
</table>
<table class="simpletable">
<tr>
                  <td></td>                  <th>Coef.</th> <th>Std.Err.</th>    <th>z</th>   <th>P>|z|</th> <th>[0.025</th> <th>0.975]</th>
</tr>
<tr>
  <th>Intercept</th>                         <td>0.421</td>   <td>0.008</td>  <td>51.361</td> <td>0.000</td>  <td>0.405</td>  <td>0.438</td>
</tr>
<tr>
  <th>C(leg, Treatment('leg0'))[T.leg1]</th> <td>0.054</td>   <td>0.003</td>  <td>15.491</td> <td>0.000</td>  <td>0.047</td>  <td>0.060</td>
</tr>
<tr>
  <th>C(leg, Treatment('leg0'))[T.leg2]</th> <td>0.014</td>   <td>0.003</td>   <td>4.549</td> <td>0.000</td>  <td>0.008</td>  <td>0.020</td>
</tr>
<tr>
  <th>roi Var</th>                           <td>0.015</td>   <td>0.007</td>     <td></td>      <td></td>       <td></td>       <td></td>   
</tr>
</table></details>


It is enouraging that the effects stay almost the same even if we include all rois and all Events. (Until now, we have been doing this only on events shorter that 10s)


```python
result = mf.get_hw_effect(Events, ref_leg="leg0", control_for_roi=True, minEvents=0)
```

    There are 339 rois, which have at least 0 events in all legs.



<style type="text/css">
            /* Fix details summary arrow
               not shown in Firefox
               due to bootstrap
               display: block;
             */
            summary {
                display: list-item;
            }
            </style><div>reference leg: leg0</div><table><tr><td style="text-align:left">Obtained effect for leg '<span style='font-family:monospace'>leg1</span>':</td><td>+9.8%,</td><td>with the p-value of 2.4e-23</td></tr><tr><td style="text-align:left">Obtained effect for leg '<span style='font-family:monospace'>leg2</span>':</td><td>+1.1%,</td><td>with the p-value of 0.21</td></tr></table><details><summary style='color:navy'>full output</summary><table class="simpletable">
<tr>
       <td>Model:</td>       <td>MixedLM</td> <td>Dependent Variable:</td>  <td>log10_hw</td> 
</tr>
<tr>
  <td>No. Observations:</td>  <td>31736</td>        <td>Method:</td>          <td>REML</td>   
</tr>
<tr>
     <td>No. Groups:</td>      <td>318</td>         <td>Scale:</td>          <td>0.0745</td>  
</tr>
<tr>
  <td>Min. group size:</td>     <td>1</td>      <td>Log-Likelihood:</td>   <td>-4452.2958</td>
</tr>
<tr>
  <td>Max. group size:</td>    <td>478</td>       <td>Converged:</td>          <td>Yes</td>   
</tr>
<tr>
  <td>Mean group size:</td>   <td>99.8</td>            <td></td>                <td></td>     
</tr>
</table>
<table class="simpletable">
<tr>
                  <td></td>                  <th>Coef.</th> <th>Std.Err.</th>    <th>z</th>   <th>P>|z|</th> <th>[0.025</th> <th>0.975]</th>
</tr>
<tr>
  <th>Intercept</th>                         <td>0.578</td>   <td>0.020</td>  <td>29.246</td> <td>0.000</td>  <td>0.539</td>  <td>0.617</td>
</tr>
<tr>
  <th>C(leg, Treatment('leg0'))[T.leg1]</th> <td>0.041</td>   <td>0.004</td>   <td>9.956</td> <td>0.000</td>  <td>0.033</td>  <td>0.049</td>
</tr>
<tr>
  <th>C(leg, Treatment('leg0'))[T.leg2]</th> <td>0.005</td>   <td>0.004</td>   <td>1.245</td> <td>0.213</td> <td>-0.003</td>  <td>0.012</td>
</tr>
<tr>
  <th>roi Var</th>                           <td>0.110</td>   <td>0.043</td>     <td></td>      <td></td>       <td></td>       <td></td>   
</tr>
</table></details>


#### Poor man's mixed linear model

The above statistical analysis is a bit detached from intuition, so let's convince ourselves it makes sense.  
First, let's find out which rois are the ones it is based on.


```python
minEvents = 10
legs = legDict.keys()
chooseRois = [roi for roi, ev in events_short.groupby("roi") if all([len(ev.query(f"leg=='{leg}'")) >= minEvents for leg in legs])]
```


```python
len(chooseRois)
```




    131




```python
regions.plotEdges(ix=chooseRois)
```


![png](temp_files/temp_57_0.png)



```python
# Let's check the leg0 vs leg1
leg = "leg1"
ref_leg = "leg0"
x = np.array([events_short.query(f"roi=={roi} and leg=='{ref_leg}'")["halfwidth"].median() for roi in chooseRois])
y = np.array([events_short.query(f"roi=={roi} and leg=='{leg}'"    )["halfwidth"].median() for roi in chooseRois])
```


```python
ax = plt.subplot(111)
ax.plot(x,y,".")
xr = np.linspace(x.min(),x.max())
ax.plot(xr,xr,"k--", label=r"$y = x$ (no effect)")
ax.plot(xr,1.13*xr,"C3", label=r"$y = 1.13x$ (effect is +13%, as inferred above)")

# See linear_regression example from the 
# documentation for explanation of this
myres = mf.fitLine(x,y)
yfit = myres['params'][0] + myres['params'][1]*xr
yefit = myres["details"]["se_fit"](xr)
upperCI = yfit + yefit
lowerCI = yfit - yefit
ax.plot(xr,yfit, "C0",label="inference from medians (95%CI)")
ax.plot(xr,upperCI, "C0",lw=.5)
ax.plot(xr,lowerCI, "C0",lw=.5)
ax.fill_between(xr,upperCI,lowerCI,color="C0",alpha = .5)

# cosmetics of the figure
ax.set_aspect("equal")
ax.set_xlabel(ref_leg)
ax.set_ylabel(leg)
ax.set_title("median halfwidth per leg")
ax.legend(loc=(1.01,.01))
```




    <matplotlib.legend.Legend at 0x7fe7a2e39908>




![png](temp_files/temp_59_1.png)



```python
np.mean(y/x), np.median(y/x)
```




    (1.2114176307930447, 1.1919868681893064)



So, naively, one arrives to the effect size of approximately +20%. The fact that the values are not the same should not surprise us, since this analysis treats all rois equally and is based only on their medians. But, it is encouraging that the inferred effects are in the same ballpark.

The same game for leg2 vs leg0 is left as an excercise for the careful reader :-D

## Events rate

Let us compare the rate with which the events occur in the three legs.  
And, let's do it in two (in this case arbitrarily chosen) halfwith domains: from 1-8s ("bursts") and above 10s ("long").


```python
hwRegions = {
    "bursts": ( 1,8),
    "long":   (10,200),
}
```

_The function below does the following:_
First, it will discard rois that have less than `minEvents` events in the whole set.  
Then, for each of the legs, and for each domain of durations, it calculates the rate of firing as number / (duration in minutes) / (number of active rois).  
For each leg and each halfwidth region, this is done many times (`Nbootstrap`), with randomly chosen rois (some may repeat, yes), in order to obtain also an estimate of the variance associated with this rate. In addition, we know our segmentation algorithm results in more rois than cells. So if we took all rois, we may obtain a more precise result than it actually is. To reduce this artefact, we can use not all rois, but a subset of them in each subsampling. This is regulated by `reduceRoi` argument. We typically see automatic segmentation without manual curation result in approx 2 times more rois than cells (so I put `reduceRoi=2` below), but if you have carefully used merging and deleting of rois in examiner (or are otherwise confident your rois represent actual independent cells), you can put it to `1`.


```python
evpmpar = mf.get_events_per_min_per_nrois(
    Events,
    hwRegions,
    minEvents=10,
    reduceRois=2,
    Nbootstrap=30
)
```

    There are 244 rois with more than 10 events (72%).
    Out of them, 122 are sampled 30 times to estimate mean and std of the firing rate.


Let have a look at what we got:


```python
evpmpar.sort_values("hw_region")
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>hw_region</th>
      <th>leg</th>
      <th>epmpr</th>
      <th>epmpr_std</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>bursts</td>
      <td>leg0</td>
      <td>3.257027</td>
      <td>0.225021</td>
    </tr>
    <tr>
      <th>2</th>
      <td>bursts</td>
      <td>leg1</td>
      <td>3.369956</td>
      <td>0.198765</td>
    </tr>
    <tr>
      <th>4</th>
      <td>bursts</td>
      <td>leg2</td>
      <td>2.877914</td>
      <td>0.102780</td>
    </tr>
    <tr>
      <th>1</th>
      <td>long</td>
      <td>leg0</td>
      <td>0.099410</td>
      <td>0.007028</td>
    </tr>
    <tr>
      <th>3</th>
      <td>long</td>
      <td>leg1</td>
      <td>0.087642</td>
      <td>0.009161</td>
    </tr>
    <tr>
      <th>5</th>
      <td>long</td>
      <td>leg2</td>
      <td>0.083378</td>
      <td>0.007488</td>
    </tr>
  </tbody>
</table>
</div>



and more visually


```python
fig, axs = plt.subplots(1,len(hwRegions),figsize=(len(hwRegions)*(len(legs)+2)/2,3))
fig.suptitle("Events per minute per active roi*")
for hwr,ax in zip(hwRegions,axs):
    y = [evpmpar.query(f"leg=='{leg}' and hw_region=='{hwr}'")['epmpr'].values for leg in legs]
    ye= [evpmpar.query(f"leg=='{leg}' and hw_region=='{hwr}'")['epmpr_std'].values for leg in legs]
    y = np.array(y).flatten()
    ye= np.array(ye).flatten()
    x = np.arange(len(y))
    ax.errorbar(x,y,ye, color="k", marker="o",ls="none")
    ax.set_xlim(-.5,x[-1]+.5)
    ax.set_xticks(x)
    ax.set_xticklabels(legs)
    ax.set_title(hwr + r" ($%g < \tau_{1/2} \leq %g$)"%hwRegions[hwr])
plt.tight_layout()
ax = axs[0]
ax.text(0,-.3,f"* active rois are those that have at least {minEvents} events recorded",transform=ax.transAxes)
```




    Text(0, -0.3, '* active rois are those that have at least 10 events recorded')




![png](temp_files/temp_71_1.png)


Note that the error bars here represent the uncertainty due to finite number of rois only. A potentially much larger source of uncertainty is the finite duration of observation. Namely, for a system with large synchronicity, the number of events within some time period depends greatly on how many bursts happened, which depends essentially on their frequency. If the frequency is high and the experiment captured for many bursts in a leg, I might have bootstrapped over time slices too to capture that. Unfortunatelly, this is typically not the case.  

But, definitely _the_ largest source of uncertainty is biology itself: mouse, preparation, islet, etc. For this reason, it is much more important to average over different experiments, than rely on a single one. 


```python

```


```python

```
