#' % Event detection and processing
#' % Srdjan Sarikas
#' % 29 Jan 2021
#'
#' # Basics
#' ## Single trace
#' So, we have a trace with events to detect:

#+echo=False
import pandas as pd
pd.set_option("display.precision", 3)
import islets
regions = islets.load_regions("small_example_rois.pkl", baremin=True, calcInterest=False )
#'
#+results='hidden'
roi = 226
regions.plot_trace(roi)

#' We estimate the events to be a few seconds long. So, filtering at a timescale of say 10s could produce a nice separation between the uninteresting component with features slower than 10s and a faster, interesting one.

#+results='hidden'
ts = 10  # timescale
regions.fast_filter_traces(ts)
regions.plot_trace(roi, cols=["trace","slower_10","faster_10"])

#' Filtering also provides a z-score, aka [_standard score_](https://en.wikipedia.org/wiki/Standard_score), of the signal. This measures the significance of a particular event in the units of standard deviation $\sigma$. For example, $z=2$ corresponds to a p-value of 0.025, while for $z=3$, $p_v\sim$0.015.

#+results='hidden'
regions.fast_filter_traces(ts)
regions.plot_trace(roi, cols=["zScore_10"])

#' To detect the excursion which we see in the plot, we need to choose a threshold value for the `zScore`. This trace is particularly prominent, so we can choose quite a high value, say $z_{th} = 5$.

regions.detect_events(ts, z_th=5)

#' The result of this command is saved within the object itself and is just a [pandas DataFrame](https://pythonbasics.org/pandas-dataframe/):


#+results='raw'
print(regions.events["10"].head())

#' Each event is characterized by it `height` (in z-score), the width at the half of this heigth (`halfwidth`), and the timepoint when this is happens. There is also the time of the peak itself `peakpoint` (often, but not necessarily at the center of the peak), and of course what `roi` it refers to.

#' For illustration, we restrict to the one we are interested in

my_events = regions.events["10"].query("roi==%i"%roi)
print(my_events)

#" ## Basic plotting of events

#' The basic command for plotting the events is

#+results="hidden"
# I am presuming you have imported islets previously
islets.EventDistillery.plot_events(my_events)

#' We may want to plot the trace together with the events, but then we need pass the [axes](https://realpython.com/python-matplotlib-guide/#the-matplotlib-object-hierarchy) as an argument. Also, the default colormap above might not look great on white background, so we can choose [another one](https://matplotlib.org/tutorials/colors/colormaps.html).
ax_trace = regions.plot_trace(roi)
islets.EventDistillery.plot_events(my_events, ax=ax_trace, cmap="Greys_r")

#' The same function can be used to plot all detected events
islets.EventDistillery.plot_events((regions.events["10"]))

#' You may want to examine the events more closely using the interactive app
#+evaluate=False
# works only in jupyter notebook
regions.examine_events(my_events, x="halfwidth", y="height")
# like this
regions.examine_events(my_events, x="t0", y="halfwidth")
# you can of course all pass it all detected events
regions.examine_events(regions.events["10"], x="t0", y="halfwidth")

#' An interesting static graph can also be obtained using native plotting library
regions.events["10"].plot.hexbin(
    x="t0",
    y="halfwidth",
    figsize=(14,4),
    )
#' For other useful options of the hexbin plot see (documentation)[https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.hexbin.html].
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'