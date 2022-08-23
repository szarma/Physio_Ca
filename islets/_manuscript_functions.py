import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import warnings
from IPython.display import display, HTML
from matplotlib.lines import Line2D
from scipy.stats import distributions as dst
from scipy.stats import ttest_1samp
from typing import Dict, List, Optional, Tuple, Iterable
import islets
from sklearn.mixture import GaussianMixture
from warnings import warn

myGrey = (.2,)*3

def cleanup_minor_ticks(ax, which):
    ax.get_figure().canvas.draw()
    for z in which:
        minors = eval(f"plt.Axes.get_{z}ticks")(ax,minor=True)
        majors = eval(f"plt.Axes.get_{z}ticks")(ax)
        min_, max_ = eval(f"plt.Axes.get_{z}lim")(ax)
        min_ = majors[np.searchsorted(majors, min_)]
        max_ = majors[np.searchsorted(majors, max_-1)]
        m = [v for v in minors if v>min_ and v<max_]
        eval(f"plt.Axes.set_{z}ticks")(ax, m, minor=True)

def gm_fit(x, ncs=None, crit="bic", debug=False, mindist = 0):
    if ncs is None:
        ncs = [1,2,3]
    gms = [GaussianMixture(nc,
                           means_init=np.percentile(x,np.linspace(1,99,nc)).reshape(-1,1)
                          ).fit(x.reshape(-1,1)) for nc in ncs]
    iopt = np.argmin([getattr(gm,crit)(x.reshape(-1,1)) for gm in gms])
    gm = gms[iopt]
    if mindist>0:
        while True:
            gm = gms[iopt]
            if gm.n_components==1:
                break
            dists = np.diff(sorted(gm.means_.flat))
            if dists.min()<mindist:
                iopt -= 1
                if iopt<0:
                    warn("Could not satisfy minimal distance between components.")
                    break
            else:
                break
    if debug:
        return gm, gms
    else:
        return gm

def expFormat(x, spacer=r"\times"):
    if x<.01:
        v,e = ("%.1e"%x).split("e")
        e = int(e)
        return r"%s{%s}10^{%i}"%(v,spacer, e)
    else:
        return "%.3f"%x

def showLineScan(linescan, axs, pixel_limits, Tind = 1, slice_ = slice(None), offset = 1, nrebin = None):
    data = linescan.data[slice_]
    if nrebin is None:
        nrebin = data.shape[1]//100
    x = islets.numeric.rebin(data,nrebin,1, func=np.sum)
    for j in range(x.shape[0]):
        x[j] -= np.percentile(x[j],10)
    ax = axs[0]
    distance = 10
    physicalSize = linescan.metadata["pxSize"]*data.shape[0]
    txtOffset = physicalSize-1.5*distance
    tmax = data.shape[1]/linescan.Freq#linescan.time[-1]
    ax.imshow(x, cmap="turbo", vmin = 0, vmax = x.max(), aspect="auto", extent = (0,tmax, 0, physicalSize), origin = "lower")

    ax.plot([tmax*.05]*2,[txtOffset,txtOffset+distance],"lightgrey")
    ax.text(tmax*.05,txtOffset+distance/2," "+str(distance)+"µm", va="center", color="w")
    # ax.imshow(x, cmap="Greys", vmin = 0, vmax = x.max(), aspect="auto")
    ax.plot([tmax*.05, tmax*.03+Tind],[.03*physicalSize]*2,"lightgrey")
    ax.text(tmax*.05+Tind/2,.03*physicalSize,str(Tind)+"s\n", ha="center", va="center", color="w", )
    ax1 = axs[1]
    off = 0
    for px0,px1 in pixel_limits:
        emphasize_region(ax,ax.get_xlim(), [px0*physicalSize/x.shape[0],px1*physicalSize/x.shape[0]], color="lightgrey", extend = (-.01,0),lw=1)
        ax1.plot(np.linspace(0,tmax,x.shape[1]), x[px0:px1].mean(0) + off, lw=.5, color=myGrey)
        # ar = ax.arrow(tmax*1.05, np.mean([px0, px1])*physicalSize/x.shape[0], tmax*.2, 0, color=myGrey, lw=1, clip_on=False, head_width = tmax, head_length = .2)
        # ar.set_clip_on(False)
        ax.set_xlim(0,tmax)
        off += offset
    
    mystyle_axes(ax)
    mystyle_axes(ax1)
#     ln = Line2D([tmax*1.05, tmax*1.3],[np.mean([px0, px1])*physicalSize/x.shape[0],]*2, color=myGrey, lw=1)
#     ln.set_clip_on(False)
#     ax.add_line(ln)
    #place_indicator(ax1,Tind,position=(tmax-1.5*Tind,1))
    return x
#
def fitLine(x, y, alpha=0.05, newx=None, plotFlag=False, show0=False, summary=False, ax=None, prediction=True):
    """ Fit a curve to the data using a least squares 1st order polynomial fit
    mostly copied from: central.scipy.org/item/50/1/line-fit-with-confidence-intervals
    """
    from scipy.stats import t
    # Summary data
    n = len(x)  # number of samples

    Sxx = np.sum(x ** 2) - np.sum(x) ** 2 / n
    Sxy = np.sum(x * y) - np.sum(x) * np.sum(y) / n
    mean_x = np.mean(x)
    mean_y = np.mean(y)

    # Linefit
    b = Sxy / Sxx
    a = mean_y - b * mean_x

    # Residuals
    fit = lambda xx: a + b * xx
    residuals = y - fit(x)

    var_res = np.sum(residuals ** 2) / (n - 2)
    sd_res = np.sqrt(var_res)

    # Confidence intervals
    se_b = sd_res / np.sqrt(Sxx)
    se_a = sd_res * np.sqrt(np.sum(x ** 2) / (n * Sxx))

    df = n - 2  # degrees of freedom
    tval = t.isf(alpha / 2., df)  # appropriate t value

    ci_a = a + tval * se_a * np.array([-1, 1])
    ci_b = b + tval * se_b * np.array([-1, 1])

    # create series of new test x-values to predict for
    npts = 100
    if show0:
        px = np.linspace(0., np.max(x), num = npts)
    else:
        px = np.linspace(np.min(x), np.max(x), num = npts)

    se_fit = lambda xi: sd_res * np.sqrt(1. / n + (xi - mean_x) ** 2 / Sxx)
    se_predict = lambda xi: sd_res * np.sqrt(1 + 1. / n + (xi - mean_x) ** 2 / Sxx)

    if summary:
        print('Fitting y = a + b*x')
        print('a={0:5.4f}+/-{1:5.4f}, b={2:5.4f}+/-{3:5.4f}'.format(a, tval * se_a, b, tval * se_b))
        print('Confidence intervals: ci_a=({0:5.4f} - {1:5.4f}), ci_b=({2:5.4f} - {3:5.4f})'.format(ci_a[0], ci_a[1],
                                                                                                    ci_b[0], ci_b[1]))
        print('Residuals: variance = {0:5.4f}, standard deviation = {1:5.4f}'.format(var_res, sd_res))
        print('alpha = {0:.3f}, tval = {1:5.4f}, df={2:d}'.format(alpha, tval, df))

    # Return info
    ri = {'residuals': residuals,
          'var_res': var_res,
          'sd_res': sd_res,
          'alpha': alpha,
          'tval': tval,
          'df': df,
          'se_fit': se_fit,
          }

    if plotFlag:
        # Plot the data
        if ax is None:
            plt.figure()
            ax = plt.subplot(111)

        ax.plot(px, fit(px), 'k', label = 'Regression line')
        ax.plot(x, y, 'r.', label = 'Sample observations')

        x.sort()
        limit = (1 - alpha) * 100
        ax.plot(px, fit(px) + tval * se_fit(px), 'r--', label = 'Confidence limit ({0:.1f}%)'.format(limit))
        ax.plot(px, fit(px) - tval * se_fit(px), 'r--')
        if prediction:
            ax.plot(px, fit(px) + tval * se_predict(px), 'c--', label = 'Prediction limit ({0:.1f}%)'.format(limit))
            ax.plot(px, fit(px) - tval * se_predict(px), 'c--')

        ax.set_xlabel('X values')
        ax.set_ylabel('Y values')
        # ax.set_xlim(0,plt.xlim()[1]*1.1)
        ax.set_title('Linear regression and confidence limits')

        # configure legend
        leg = ax.legend(loc = 0)
        # leg = plt.gca().get_legend()
        ltext = leg.get_texts()
        plt.setp(ltext, fontsize = 10)

    output = {
        "func": "y = a + b*x",
        "params": (a, b),
        "confidence interval": (ci_a, ci_b),
        "t value": tval,
        "standard error": (se_a, se_b),
        "details": ri,
    }
    if newx is not None and len(newx):
        newx = np.array(newx)

        print('Example: x = {0}+/-{1} =&gt; se_fit = {2:5.4f}, se_predict = {3:6.5f}'
              .format(newx[0], tval * se_predict(newx[0]), se_fit(newx[0]), se_predict(newx[0])))

        newy = (fit(newx), fit(newx) - se_predict(newx), fit(newx) + se_predict(newx))
        output["newy"] = newy
    return output


protocolColorScheme = {
    ("glucose", "6mM"): (0.98,) * 3,
    ("glucose", "8mM"): (0.8,) * 3,
    ("Ca", "2mM"): tuple(1 - (1 - np.array(plt.cm.Oranges(0))) / 2),
    ("Ca", "0.4mM"): plt.cm.Oranges(.1),
    ("Ca", "0mM"): plt.cm.Oranges(.3),
    ("isradipine",  "5uM"): "lightblue",
    ("isradipine",  "5µM"): "lightblue",
    ("diazoxide","100uM"): "xkcd:pale rose",
    ("ryanodine","100nM"): "#f6cefc",
    ("ryanodine","100uM"): "xkcd:light purple",
    ("ryanodine","10uM"): "xkcd:light purple",
    ("acetylcholine","10uM"): "xkcd:light pastel green",
    ("caffeine","10mM"): "xkcd:light brown",
    ("xestospongin","3uM"): "goldenrod",
}


def beautify_protocol(protocol):
    from .protocol import Protocol
    protocol = protocol.copy()
    colors = []
    for i in protocol.index:
        colors += [protocolColorScheme.get((protocol.loc[i, "compound"], protocol.loc[i, "concentration"]), "white")]
    protocol["color"] = colors
    for comp, df in protocol.groupby("compound"):
        units = df["unit"].unique()
        if len(units) == 1:
            unit = units[0]
            # protocol.loc[df.index, "compound"] = [f"{comp} [{unit}]"] * len(df)
            protocol.loc[df.index, "concentration"] = [df["concentration"].iloc[0]] + [f"{c}" for c in df["conc"][1:]]
    # protocol = protocol.replace("Ca [mM]", r"Ca${}^{2\!+}$[mM] ")
    protocol = protocol.replace("Ca", r"Ca${}^{2\!+}$ ")
    # protocol = Protocol(protocol)
    return protocol


def ruler(fig, margin=.5, grid=True):
    from matplotlib.ticker import MultipleLocator
    figwidth = fig.get_figwidth()
    figheight = fig.get_figheight()
    dax = fig.add_axes([0, 0, 1, 1],
                       zorder = -2,
                       facecolor = (.98,) * 3  # "whitesmoke"
                       )
    dax.axvline(margin, color = "salmon", lw = .5, )
    dax.axvline(figwidth - margin, color = "salmon", lw = .5, )
    dax.axhline(margin, color = "salmon", lw = .5, )
    dax.axhline(figheight - margin, color = "salmon", lw = .5, )
    dax.set_xlim(0, figwidth)
    dax.set_ylim(0, figheight)
    dax.set_xticks(np.arange(0, figwidth + 1e-10))
    dax.set_yticks(np.arange(0, figheight + 1e-10))
    mystyle_axes(dax, retain = ["bottom","right","top", "left"], bounded = [False] * 4)
    dax.xaxis.set_minor_locator(MultipleLocator(.2))
    dax.yaxis.set_minor_locator(MultipleLocator(.2))
    if grid:
        dax.grid(color = (.9,) * 3, which = "major")
        dax.grid(color = (.94,) * 3, which = "minor")


def silent_fit(model):
    with warnings.catch_warnings(record = True) as w:
        result = model.fit()
        return result, w


def get_rate_effect(events, quantity="epmpr", ref_leg=None, legs=None, groups = "experiment"):
    events = events.copy()
    if legs is None:
        legs = events["leg"].dropna().unique()
    if ref_leg is None:
        ref_leg = events[~events["leg"].isna()].iloc[0]['leg']
    log_quantity = f"log10_{quantity}"
    events[log_quantity] = np.log10(events[quantity])
    model = sm.MixedLM.from_formula(f"{log_quantity} ~ C(leg, Treatment('{ref_leg}')) + 1", data = events,
                                    groups = groups)
    result, warns = silent_fit(model)
    legs = [leg for leg in legs if leg != ref_leg]
    # summary = get_summary(model, legs, ref_leg)
    html_summary = \
        """<style type="text/css">
            /* Fix details summary arrow
               not shown in Firefox
               due to bootstrap
               display: block;
             */
            summary {
                display: list-item;
            }
            </style>"""
    html_summary += f"<div>reference leg: {ref_leg}</div>"
    html_summary += "<table>"
    for leg in legs:
        coef = result.params[f"C(leg, Treatment('{ref_leg}'))[T.{leg}]"]
        pv = result.pvalues[f"C(leg, Treatment('{ref_leg}'))[T.{leg}]"]
        effect_in_pc = (10 ** coef - 1) * 100
        html_summary += f"""<tr><td style="text-align:left">Obtained effect for leg '<span style='font-family:monospace'>{leg}</span>':</td><td>%+.2g%%,</td><td>with the p-value of %.2g</td></tr>""" % (
            effect_in_pc, pv)
    #         html_summary += "<br>"
    html_summary += "</table>"
    if len(warns):
        html_summary += f"<span style='font-size:70%'>~ {len(warns)} warning(s) encountered when fitting you might want to check full output ~</span>"
    summary = result.summary()
    html_summary = f"{html_summary}<details><summary style='color:navy'>full output</summary>{summary._repr_html_()}</details>"
    display(HTML(html_summary))
    return result


def get_effect(events, quantity, ref_leg=None, legs=None, control_for_roi=True, minEvents=10, groups="roi"):
    events = events.copy()
    if legs is None:
        legs = events["leg"].dropna().unique()
    if ref_leg is None:
        ref_leg = events[~events["leg"].isna()].iloc[0]['leg']
    depvar = f"log10_{quantity}"
    events[depvar] = np.log10(events[quantity])
    if control_for_roi:
        events.dropna(inplace = True)
        chooseRois = [roi for roi, ev in events.groupby("roi") if
                      all([len(ev.query(f"leg=='{leg}'")) >= minEvents for leg in legs])]
        print(f"There are {len(chooseRois)} rois, which have at least {minEvents} events in all legs.")
        events = events[events.roi.isin(chooseRois)]
        model = sm.MixedLM.from_formula(
            f"{depvar} ~ C(leg, Treatment('{ref_leg}')) + 1",
            data = events,
            groups = groups,
            vc_formula = None if groups == "roi" else {"roi": "0 + C(roi)"},

        )
    else:
        model = sm.OLS.from_formula(f"{depvar} ~ C(leg, Treatment('{ref_leg}')) + 1", data = events)
    result, warns = silent_fit(model)
    legs = [leg for leg in legs if leg != ref_leg]
    # summary = get_summary(model, legs, ref_leg)
    html_summary = \
        """<style type="text/css">
            /* Fix details summary arrow
               not shown in Firefox
               due to bootstrap
               display: block;
             */
            summary {
                display: list-item;
            }
            </style>"""
    html_summary += f"<div>reference leg: {ref_leg}</div>"
    html_summary += "<table>"
    for leg in legs:
        coef = result.params[f"C(leg, Treatment('{ref_leg}'))[T.{leg}]"]
        pv = result.pvalues[f"C(leg, Treatment('{ref_leg}'))[T.{leg}]"]
        effect_in_pc = (10 ** coef - 1) * 100
        html_summary += f"""<tr><td style="text-align:left">Obtained effect for leg '<span style='font-family:monospace'>{leg}</span>':</td><td>%+.2g%%,</td><td>with the p-value of %.2g</td></tr>""" % (
            effect_in_pc, pv)
    #         html_summary += "<br>"
    html_summary += "</table>"
    if len(warns):
        html_summary += f"<span style='font-size:70%'>~ {len(warns)} warning(s) encountered when fitting you might want to check full output ~</span>"
    summary = result.summary()
    html_summary = f"{html_summary}<details><summary style='color:navy'>full output</summary>{summary._repr_html_()}</details>"
    display(HTML(html_summary))
    return result

def get_hw_effect(events, ref_leg=None, legs=None, control_for_roi=True, minEvents=10, groups="roi"):
    print ("This function is deprecated, please use get_effect(..., quantity = 'halfiwdth'), which then also accepts other quantities like 'height' and 'auc'.")
    events = events.copy()
    if legs is None:
        legs = events["leg"].dropna().unique()
    if ref_leg is None:
        ref_leg = events[~events["leg"].isna()].iloc[0]['leg']
    events["log10_hw"] = np.log10(events["halfwidth"])
    if control_for_roi:
        events.dropna(inplace = True)
        chooseRois = [roi for roi, ev in events.groupby("roi") if
                      all([len(ev.query(f"leg=='{leg}'")) >= minEvents for leg in legs])]
        print(f"There are {len(chooseRois)} rois, which have at least {minEvents} events in all legs.")
        events = events[events.roi.isin(chooseRois)]
        model = sm.MixedLM.from_formula(
            f"log10_hw ~ C(leg, Treatment('{ref_leg}')) + 1",
            data = events,
            groups = groups,
            vc_formula = None if groups == "roi" else {"roi": "0 + C(roi)"},

        )
    else:
        model = sm.OLS.from_formula(f"log10_hw ~ C(leg, Treatment('{ref_leg}')) + 1", data = events)
    result, warns = silent_fit(model)
    legs = [leg for leg in legs if leg != ref_leg]
    # summary = get_summary(model, legs, ref_leg)
    html_summary = \
        """<style type="text/css">
            /* Fix details summary arrow
               not shown in Firefox
               due to bootstrap
               display: block;
             */
            summary {
                display: list-item;
            }
            </style>"""
    html_summary += f"<div>reference leg: {ref_leg}</div>"
    html_summary += "<table>"
    for leg in legs:
        coef = result.params[f"C(leg, Treatment('{ref_leg}'))[T.{leg}]"]
        pv = result.pvalues[f"C(leg, Treatment('{ref_leg}'))[T.{leg}]"]
        effect_in_pc = (10 ** coef - 1) * 100
        html_summary += f"""<tr><td style="text-align:left">Obtained effect for leg '<span style='font-family:monospace'>{leg}</span>':</td><td>%+.2g%%,</td><td>with the p-value of %.2g</td></tr>""" % (
            effect_in_pc, pv)
    #         html_summary += "<br>"
    html_summary += "</table>"
    if len(warns):
        html_summary += f"<span style='font-size:70%'>~ {len(warns)} warning(s) encountered when fitting you might want to check full output ~</span>"
    summary = result.summary()
    html_summary = f"{html_summary}<details><summary style='color:navy'>full output</summary>{summary._repr_html_()}</details>"
    display(HTML(html_summary))
    return result

# def get_summary(model, legs, ref_leg=None):
#     if ref_leg is None:
#         ref_leg = legs[0]
#         legs = legs[1:]
#     result, warns = silent_fit(model)
#     html_summary = \
#          """<style type="text/css">
#             /* Fix details summary arrow
#                not shown in Firefox
#                due to bootstrap
#                display: block;
#              */
#             summary {
#                 display: list-item;
#             }
#             </style>"""
#     html_summary += "<table>"
#     for leg in legs:
#         coef = result.params[f"C(leg, Treatment('{ref_leg}'))[T.{leg}]"]
#         pv = result.pvalues[f"C(leg, Treatment('{ref_leg}'))[T.{leg}]"]
#         effect_in_pc = (10 ** coef - 1) * 100
#         html_summary += f"""<tr><td style="text-align:left">Obtained effect for leg '<span style='font-family:monospace'>{leg}</span>':</td><td>%+.2g%%,</td><td>with the p-value of %.2g</td></tr>""" % (
#         effect_in_pc, pv)
#     #         html_summary += "<br>"
#     html_summary += "</table>"
#     if len(warns):
#         html_summary += f"~ {len(warns)} warning(s) encountered when fitting you might want to check full output ~"
#     summary = result.summary()
#     return f"{html_summary}<details><summary>full output</summary>{summary._repr_html_()}</details>"


def get_tidy(Events, hwmin, hwmax, legKeys=None):
    if legKeys is None:
        legKeys = sorted([leg for leg in Events['leg'].unique() if isinstance(leg, str)])
    rois = Events.roi.unique()
    Nevents = pd.DataFrame(
        index = rois,
        columns = legKeys
    )
    for leg in legKeys:
        Nevents[leg] = [len(Events.query(f"roi=={roi} and leg=='{leg}' and halfwidth>={hwmin} and halfwidth<{hwmax}"))
                        for roi in rois]

    activeRois = Nevents.index[(Nevents > 1).all(1)]
    tmp = []
    for leg in legKeys:
        ev = Events.query(f"leg=='{leg}'")
        duration = ev["peakpoint"].max() - ev["peakpoint"].min()
        for roi in activeRois:
            ev = Events.query(f"roi=={roi} and leg=='{leg}' and halfwidth>={hwmin} and halfwidth<{hwmax}")
            assert len(ev) == Nevents.loc[roi, leg]
            tmp += [{
                "leg": leg,
                "roi": roi,
                "nevents": len(ev),
                "halfwidth": ev["halfwidth"].median(),  # if len(ev)>3 else np.nan,
                "height": ev["height"].median(),  # if len(ev)>3 else np.nan,
                "duration": duration
            }]
    return pd.DataFrame(tmp)


def histogram_of_hw(Events, hwbinEdges, legDict, legColorDict, ax=None, orientation="bottom", hist_scale="left",
                    control_line=False, n_sigma=False):
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    if "leg" not in Events.columns:
        islets.EventDistillery.define_legs(Events, legDict)
    for leg in legDict:
        ev = Events.query(f"leg=='{leg}'")
        h = np.histogram(ev["halfwidth"], hwbinEdges)[0]
        herr = np.sqrt(h)
        hmin = h / (legDict[leg][1] - legDict[leg][0]) * 60  # events per minute
        hmin_err = herr / (legDict[leg][1] - legDict[leg][0]) * 60
        if orientation in ["bottom", "top"]:
            ax.plot(
                np.repeat(hwbinEdges, 2),
                [0] + list(np.repeat(hmin, 2)) + [0], "-", c = legColorDict[leg],
                label = leg)
            if n_sigma:
                ax.fill_between(
                    np.repeat(hwbinEdges, 2),
                    [0] + list(np.repeat(hmin + n_sigma * hmin_err, 2)) + [0],
                    [0] + list(np.repeat(hmin - n_sigma * hmin_err, 2)) + [0],
                    color = legColorDict[leg],
                    alpha = .1,
                )

            ax.set_xscale("log")
            if hist_scale == "top":
                ax.set_ylim(h.max() * 1.05, 0)
            ax.set_xlim(hwbinEdges[0], hwbinEdges[-1])
            ax.set_ylabel("events/min")
            if control_line:
                x = np.sqrt(hwbinEdges[:-1] * hwbinEdges[1:])
                ax.step(x, 60 * dst.norm.sf(Events.z_max.min()) / x,  # *np.diff(hwbinEdges),
                        c = legColorDict[leg], ls = "--", where = "mid")
        else:
            ax.plot([0] + list(np.repeat(hmin, 2)) + [0], np.repeat(hwbinEdges, 2), "-", c = legColorDict[leg],
                    label = leg)
            ax.set_yscale("log")
            if hist_scale == "right":
                ax.set_xlim(h.max() * 1.05, 0)
            ax.set_ylim(hwbinEdges[0], hwbinEdges[-1])
            ax.set_xlabel("events/min")
    # mystyle_axes(ax, retain = [orientation, hist_scale], bounded = [False, True])

    legend = ax.legend(frameon = False)


def boxplots_per_leg(roi_stats, legDict, legColorDict, pvalue=None, ax=None, **bx_kwargs):
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    for ileg, leg in enumerate(legDict):
        box_kwargs = dict(
            notch = True,
            showfliers = False,
            widths = .4,
            whis = (5, 95)
        )
        if bx_kwargs is not None:
            box_kwargs.update(bx_kwargs)
        bxs = ax.boxplot(roi_stats[leg], positions = [ileg], **box_kwargs)
        for el in bxs:
            for ln in bxs[el]:
                ln.set_color(legColorDict[leg])
    if pvalue is not None:
        for kp in pvalue:
            i, j = kp
            ki = list(legDict)[i]
            kj = list(legDict)[j]
            x = roi_stats[kj] - roi_stats[ki]
            # std = x.std(ddof=1)
            # mu = x.mean()
            # sem = std/np.sqrt(len(x))
            # assert sem==pd.Series(x).sem()
            # p = dst.norm.sf(mu/sem)
            ypos = pvalue[kp].get("ypos", np.percentile(roi_stats[[ki, kj]], 99.7))
            alternative = pvalue[kp].get("alternative", "two-sided")
            p = ttest_1samp(x, popmean = 0, alternative = alternative).pvalue
            ax.plot([i, j], [ypos] * 2, color = "k")
            if p < .005:
                plabel = r"   $p < 10^{%i}$" % int(np.log10(p))
            else:
                plabel = r"   $p = %.1g$" % p
            ax.text((i + j) / 2, ypos, plabel, ha = "center", va = "bottom")
    ax.set_xticks(range(len(legDict)))
    ax.set_xticklabels(list(legDict))


def fancy_barchart(ax, df, c=(.4,) * 3, orientation="horizontal", plabel_offset=.1, value="mean", error=None):
    c = df.get("color", c)
    y = df[value].values
    if error is None:
        yerr = np.zeros_like(y)
    else:
        yerr = df[error].values
    nticks = len(df.index)
    if orientation == "horizontal":
        ax.barh(range(nticks), y, color = c)
        ax.errorbar(y, range(nticks), xerr = yerr, ls = 'none', c = c)
        ax.set_yticks([])
        for jt, tk in enumerate(df.index):
            ax.text(0, jt, tk + " ", ha = "right", va = "center")
    elif orientation == "vertical":
        ax.bar(range(nticks), y, color = c)
        ax.errorbar(range(nticks), y, yerr = yerr, ls = 'none', c = c)
        ax.set_xticks([])
    else:
        raise ValueError("orientation can only be 'horizonal' or 'vertical'")
    if error is None:
        return None
    comparisons = [(i, j) for i in range(len(df.index)) for j in range(i + 1, len(df.index))]
    comparisons = sorted(comparisons, key = lambda ij: ((y + yerr)[list(ij)].max(), np.diff(ij)[0]))

    offset = max(y + 2 * yerr) * plabel_offset
    xpos0 = 0

    for i, j in comparisons:
        dy = y[j] - y[i]
        sigma = (yerr[i] ** 2 + yerr[j] ** 2) ** .5
        z = np.abs(dy) / sigma
        p = dst.norm.sf(z)
        xpos = max([y[k] + 2 * yerr[k] for k in [i, j]]) + offset / 2
        while xpos <= xpos0 + offset:
            xpos += offset
        # print(i, j, xpos)
        xpos0 = xpos
        ln = Line2D([xpos - offset / 5, xpos, xpos, xpos - offset / 5], [i, i, j, j], c = "k", lw = .7)
        ln.set_clip_on(False)
        ax.add_line(ln)
        if p < .01:
            txt = "%.0e" % p
            txt, exponent = txt.split("e")
            txt = txt.split(".")[0]
            if txt == "1":
                plabel = r" $p = 10^{%i}$" % (int(exponent))
            else:
                plabel = r" $p = %s{\cdot}10^{%i}$" % (txt, int(exponent))
            # plabel = r" $p < 10^{%i}$" % int(np.log10(p))
            # plabel = "*"*max(int(np.log10(p)),3)
            # plabel = r" $%i\sigma$"%int(z)
        else:
            plabel = r" $p = %.1g$" % p
        if orientation == "horizontal":
            ax.text(xpos, (i + j) / 2, plabel,
                    # rotation="90",
                    va = "center"
                    )


def mystyle_axes(ax, retain=None, xlim=None, ylim=None, bounded=None):
    if retain is None:
        retain = []
    if bounded is None:
        bounded = [True] * len(retain)
    boundedDict = dict(zip(retain, bounded))
    # renderer = ax.get_figure().canvas.get_renderer()
    # ax.draw(renderer)
    for sp in ax.spines:
        ax.spines[sp].set_visible(sp in retain)
    if xlim is None:
        xlim = ax.get_xlim()
    xt = ax.get_xticks()
    if ylim is None:
        ylim = ax.get_ylim()
    yt = ax.get_yticks()
    for sp in ["top", "bottom"]:
        if len(xt) and sp in retain and boundedDict[sp]:
            xbounds = xt[np.searchsorted(xt, xlim[0])], xt[np.searchsorted(xt, xlim[1] + 1e-10) - 1]
            ax.spines[sp].set_bounds(xbounds)
    for sp in ["left", "right"]:
        if len(yt) and sp in retain and boundedDict[sp]:
            ybounds = yt[np.searchsorted(yt, ylim[0])], yt[np.searchsorted(yt, ylim[1] + 1e-10) - 1]
            ax.spines[sp].set_bounds(ybounds)
    if "top" in retain:
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')
    if "right" in retain:
        ax.yaxis.set_ticks_position('right')
        ax.yaxis.set_label_position('right')
    if "bottom" not in retain and "top" not in retain:
        ax.set_xticks([])
    if "right" not in retain and "left" not in retain:
        ax.set_yticks([])


def emphasize_region(ax, x, y, extend=(0, 0), clip_on = False, **line_kwargs):
    xb, xe = x
    yb, ye = y
    assert xe > xb
    assert ye > yb

    dx = xe - xb
    dx = extend[0] * dx
    xe = xe + dx
    xb = xb - dx

    dy = ye - yb
    dy = extend[1] * dy
    ye = ye + dy
    yb = yb - dy
    if clip_on:
        ln = ax.plot(
            [xb, xe, xe, xb, xb],
            [yb, yb, ye, ye, yb],
            **line_kwargs
        )[0]
    else:
        ln = Line2D(
            [xb, xe, xe, xb, xb],
            [yb, yb, ye, ye, yb],
            **line_kwargs
        )
        ln.set_clip_on(False)
        ax.add_line(ln)
    return ln

def plot_events(Events,
                regions,
                timeUnits="s",
                plottype="hexbin",
                only_legs=False,
                protocol=None,
                legColorDict=None,
                axheight=2,
                offset=2,
                ySpineOffset=100,
                default_color = myGrey,
                **kwargs):
    if protocol is None:
        protocol = regions.protocol.copy()

    duration = regions.time[-1] - regions.time[0]

    if only_legs and "leg" in Events.columns:
        # take only events that happen during one of the specified legs
        Events = Events[~Events.leg.isna()].copy()

    ### creation and formatting of figure
    axwidth = (duration + ySpineOffset) / 500

    figwidth, figheight = axwidth + 2 * offset, axheight + 2 * offset
    fig = plt.figure(figsize = (figwidth, figheight))
    ax = fig.add_axes([offset / figwidth, offset / figheight, axwidth / figwidth, axheight / figheight])
    # axp = fig.add_axes([offset / figwidth, (axheight + 1.1 * offset) / figheight, axwidth / figwidth,
    #                     offset / 10 / figheight * len(protocol.compound.unique())])
    # axp.get_shared_x_axes().join(ax, axp)
    # mystyle_axes(axp)
    # from .utils import add_protocol
    # add_protocol(axp, protocol)
    if timeUnits[0] == "m":
        # duration = duration/60
        ax.set_xlabel("time [min]")
        for col in ["t_begin", "t_end"]:
            protocol[col] = protocol[col] / 60
    elif timeUnits[0] == "s":
        ax.set_xlabel("time [s]")
    ax.set_ylabel(r"$\tau_{1/2}$, halfwidth [s]")
    axp = protocol.plot_protocol(ax, only_number = False)
    output = [fig, (ax, axp)]
    ## plotting
    if plottype == "hexbin":
        cax = fig.add_axes(
            [(1.1 * offset + axwidth) / figwidth, (offset) / figheight, .2 / figwidth, axheight / figheight])
        x, y = Events["peakpoint"].values, Events["halfwidth"].values
        if timeUnits == "min":
            x = x.copy() / 60
        if "bins" not in kwargs:
            kwargs.update({"vmin": 0})
        if "gridsize" not in kwargs:
            kwargs.update({"gridsize": (1 + int(duration / 20), 50)})
        if "mincnt" not in kwargs:
            kwargs.update({"mincnt": 1})
        if "cmap" not in kwargs:
            kwargs.update({"cmap": "hot"})

        hx = ax.hexbin(x, y, yscale = "log", **kwargs)
        cax = plt.colorbar(hx, cax = cax)
        cax.set_label("bin count")
        output[1] = output[1] + (cax,)

    elif plottype == "scatter":
        if "markersize" not in kwargs and "ms" not in kwargs:
            kwargs.update({"markersize": 2})
        ax.set_yscale("log")
        if "leg" in Events.columns:
            iterlist = list(Events.groupby("leg")) + [(None, Events[Events.leg.isna()])]
        else:
            iterlist = [(None, Events)]
        for leg, ev in iterlist:
            x, y = ev["peakpoint"].values.copy(), ev["halfwidth"]
            if timeUnits == "min":
                x = x / 60
            ax.plot(x, y, ".", c = default_color if leg is None else legColorDict[leg], **kwargs)
    else:
        raise ValueError("plottype can only be 'scatter' or 'hexbin")
    ax.set_xlim(regions.time[0] - ySpineOffset, regions.time[-1])
    axp.set_xlim(regions.time[0] - ySpineOffset, regions.time[-1])
    mystyle_axes(ax, ["bottom", "left"], bounded = [True, False])
    mystyle_axes(axp)
    return output


def get_events_per_min_per_nrois(Events, hwRegions, minEvents, reduceRois, Nbootstrap=50, seed = 0):
    ### active rois are those that have at least 10 events
    activeRois = [roi for roi, evroi in Events.groupby("roi") if len(evroi) > minEvents]
    pc = len(activeRois) / len(Events['roi'].unique()) * 100
    print(f"There are {len(activeRois)} rois with more than {minEvents} events ({pc:.0f}%).")
    nrois = len(activeRois) // reduceRois
    print(f"Out of them, {nrois} are sampled {Nbootstrap} times to estimate mean and std of the firing rate.")
    ev_pm_par = []  # number of events per min per active roi
    np.random.seed(seed)
    for leg, ev in Events.groupby("leg"):
        duration = ev["peakpoint"].max() - ev["peakpoint"].min()
        for hwr in hwRegions:
            hwb = hwRegions[hwr]
            tmp = []
            for j in range(Nbootstrap):
                acts = np.random.choice(activeRois, nrois, replace = False)
                ## take only events that belong to active rois
                evs = Events[Events.roi.isin(acts)]
                evs = evs.query(f"leg=='{leg}' and halfwidth>={hwb[0]} and halfwidth<{hwb[1]}")
                tmp += [len(evs)]

            ev_pm_par += [{
                "hw_region": hwr,
                "leg": leg,
                "epmpr": np.mean(tmp) / duration / nrois * 60,
                "epmpr_std": np.std(tmp) / duration / nrois * 60,
            }]
    ev_pm_par = pd.DataFrame(ev_pm_par)
    return ev_pm_par

def get_activity(Events, minEvents, reduceRois, Nbootstrap=50, seed = 0):
    ### active rois are those that have at least 10 events
    activeRois = [roi for roi, evroi in Events.groupby("roi") if len(evroi) > minEvents]
    pc = len(activeRois) / len(Events['roi'].unique()) * 100
    print(f"There are {len(activeRois)} rois with more than {minEvents} events ({pc:.0f}%).")
    nrois = len(activeRois) // reduceRois
    print(f"Out of them, {nrois} are sampled {Nbootstrap} times to estimate mean and std of the activity.")
    out = []  # number of events per min per active roi
    np.random.seed(seed)
    for leg, ev in Events.groupby("leg"):
        duration = ev["peakpoint"].max() - ev["peakpoint"].min()
        tmp = []
        for j in range(Nbootstrap):
            acts = np.random.choice(activeRois, nrois, replace = False)
            ## take only events that belong to active rois
            evs = Events[Events.roi.isin(acts)]
            evs = evs.query(f"leg=='{leg}'")
            tmp += [evs['auc'].sum()]

        out += [{
            "leg": leg,
            "activity": np.mean(tmp) / duration / nrois * 60,
            "activity_std": np.std(tmp) / duration / nrois * 60,
        }]
    out = pd.DataFrame(out)
    return out

def get_trace_activity(regions: islets.Regions,
                       legs: dict,
                       ts: float = 300,
                       reduceRois: int = 2,
                       Nbootstrap: int = 30,
                       seed: int = 0,
                       verbose: bool = True
                       ) -> pd.DataFrame:
    regions.fast_filter_traces(ts,Npoints=int(3*ts), verbose=verbose, write=True)
    ratio = np.vstack(regions.df['faster_%g'%ts]) / np.vstack(regions.df['slower_%g'%ts])
    t = regions.showTime.get('%g' % ts, regions.time)
    dt = np.diff(t).mean()
    nrois = len(regions.df)
    AUCs = []
    np.random.seed(seed)
    nrois_choose = int(np.ceil(nrois / reduceRois))
    for leg in legs:
        begin, end = legs[leg]
        duration = end-begin
        fltr = (t > begin) & (t < end)
        integrals = ratio[:,fltr].sum(1)*dt/duration
        if Nbootstrap==1 and reduceRois==1:
            AUCs += [{"leg": leg, "activity": integrals.sum(), "activity_std": np.nan}]
        else:
            x = np.ones(nrois)
            for j in range(Nbootstrap):
                subset = np.random.choice(nrois, size=nrois_choose)
                x[j] = integrals[subset].sum()
            AUCs += [{"leg":leg, "activity": x.mean(), "activity_std": x.std()}]
    return pd.DataFrame(AUCs)




def big_plot(regions: islets.Regions,
             Events: pd.DataFrame,
             rois: Optional[Iterable] = None,
             plot_sum: bool = False,
             frame=False,
             labels=False,
             legColors=None,
             cdf=False,
             ):
    figwidth, figheight = 8.27 * 1.5, 11.69 / 2 * 1.5
    fig, axs = plt.subplots(3, 2, figsize = (figwidth, figheight), gridspec_kw = {"width_ratios": [1, 3]}, )
    if frame:
        ax = fig.add_axes([0, 0, 1, 1], zorder = -1, facecolor = "none")
        ax.plot([0, 0, 1, 1, 0], [0, 1, 1, 0, 0], "k-")
        ax.set_xlim(-.001, 1.001)
        ax.set_ylim(-.001, 1.001)
        mystyle_axes(ax)
    # plt.subplots_adjust(left = .5/figwidth, right=1-.5/figwidth, bottom = .5/figheight, top = 1-.5/figheight)
    # mf.ruler(fig)
    if rois is None:
        ccs = islets.utils.get_ccs(regions,
                                   timescales = np.array([10]),
                                   col = "faster",
                                   skip_sequential_filtering = "faster_10" in regions.df.columns
                                   )
        ccs = pd.DataFrame(ccs["all"]).T
        ccs.set_index(regions.df.index, inplace = True)
        ranks = ccs.rank(ascending = False)
        mostranked = ranks.sum(1).sort_values().index[:3]
        rois = mostranked
    regions.plotTraces(rois, axs = axs[0], protocol = False, labels = labels)
    #     axs[0,1].set_xlim()
    axs[0, 1].set_xlabel("time [s]")
    mystyle_axes(axs[0, 1], retain = ["top"], bounded = [False])
    ax = axs.flat[0]
    ax.set_xticks([])
    ax.set_yticks([])

    ax = axs.flat[3]
    hxb = ax.hexbin(Events["peakpoint"], Events["halfwidth"], yscale = "log", cmap = "hot", gridsize = (300, 100),
                    mincnt = 1)
    ax.set_ylabel("halfwidth [s]")
    if not plot_sum:
        ax.set_xlabel("time [s]")
    pos = ax.get_position()
    cax = fig.add_axes([pos.x0 + pos.width, pos.y0 + .1 * pos.height, pos.width / 40, pos.height * .8])
    plt.colorbar(hxb, cax = cax, label = "bin counts")
    if "leg" in Events.columns:
        legs = sorted(Events.dropna()["leg"].unique())
        if legColors is None:
            legColors = dict(zip(legs, ["C%i" % (j % 10) for j in range(len(legs))]))
        for leg, ev in Events.groupby("leg"):
            ys = ev['halfwidth'].min() * .95, ev['halfwidth'].max() * 1.05
            xs = ev['peakpoint'].min(), ev['peakpoint'].max()
            emphasize_region(ax, xs, ys, color = legColors[leg])

    ax.set_xlim(axs.flat[1].get_xlim())
    mystyle_axes(ax, retain = ["left", "bottom"][:2 - int(plot_sum)], bounded = [False, False])
    axp = regions.protocol.plot_protocol(ax = ax)
    mystyle_axes(axp)
    axp.set_xlim(axs.flat[1].get_xlim())

    axtraces = axs[0, 1]
    pos_traces = axtraces.get_position()
    pos_axp = axp.get_position()
    posy0tr = pos_axp.y0 + pos_axp.height + .01
    dy = pos_traces.y0 - posy0tr
    axtraces.set_position([pos_traces.x0, pos_traces.y0 - dy, pos_traces.width, pos_traces.height + dy])

    tmin = regions.time[0]
    for t in set(regions.protocol.dropna()[["t_begin", "t_end"]].values.flatten()):
        if t <= tmin + 5: continue
        for ax in axs[:, 1]:
            ax.axvline(t, lw = .5, color = (.5,) * 3, ls = "dotted")
    epmpr = get_events_per_min_per_nrois(Events, {"all": (Events["halfwidth"].min(), Events["halfwidth"].max())},
                                         minEvents = 10, reduceRois = 2)
    pos = axs[1, 0].get_position()
    axev = fig.add_axes([pos.x0, pos.y0, pos.width * .4, pos.height])
    axev.set_ylabel("events / min / active roi")
    axbx = fig.add_axes([pos.x0 + pos.width * .6, pos.y0, pos.width * .4, pos.height])
    for i, row in epmpr.iterrows():
        axev.bar([i], [row['epmpr']], color = legColors[row["leg"]])
        axev.errorbar([i], [row['epmpr']], [row['epmpr_std']], color = legColors[row["leg"]])
        axev.text(i + .5, 0, row["leg"] + " ", rotation = 45, va = "top", ha = "right")
        axbx.text(i + .5, 0, row["leg"] + " ", rotation = 45, va = "top", ha = "right")
    for il, leg in enumerate(legs):
        ev = Events.query(f"leg=='{leg}'")
        bxs = axbx.boxplot(ev["halfwidth"], positions = [il], showfliers = False, widths = .7, notch = True)
        for el in bxs:
            for ln in bxs[el]:
                ln.set_color(legColors[leg])
    hwmin, hwmax = axbx.get_ylim()
    dhw = hwmax - hwmin
    hwmin -= dhw / 5
    hwmin = max(0, hwmin)
    axbx.set_ylim(hwmin, hwmax)
    mystyle_axes(axev, ["left"])
    mystyle_axes(axbx, ["right"], bounded = [False])
    axs[1, 0].remove()

    pos = axs[2, 0].get_position()
    axs[2, 0].remove()
    if cdf:
        axcdf = fig.add_axes([pos.x0, pos.y0, pos.width, pos.height * .8])
        for leg, ev in Events.groupby("leg"):
            x = ev["halfwidth"].values.copy()
            x.sort()
            axcdf.plot(x, np.linspace(0, 1, len(x)), color = legColors[leg])
        axcdf.set_xscale("log")
        axcdf.set_xlabel("halfwidth [s]")
        axcdf.set_ylabel("CDF")
        mystyle_axes(axcdf, ["left", "bottom"], bounded = [True, False])
    ax = axs.flat[-1]
    if plot_sum:
        nrebin = int(regions.Freq / 1)
        t = regions.time
        y = np.sum([regions.df.loc[i, "trace"] * regions.df.loc[i, "size"] for i in regions.df.index],
                   0)  # /regions.df["size"].sum()

        if nrebin > 1:
            t = islets.numeric.rebin(t, nrebin)
            y = islets.numeric.rebin(y, nrebin)
            ax.plot(t, y, color = "k")
        ax.set_xlim(axs.flat[1].get_xlim())
        mystyle_axes(ax, retain = ["bottom"], bounded = [False])
        ax.set_xlabel("time [s]")
    else:
        ax.remove()


def place_indicator(ax, duration, position=(0, 1), unit="s", linekwargs={}, fontkwargs={}, text_vert_offset=-.12):
    axx = ax.get_figure().add_axes(ax.get_position(), facecolor = "none")
    axx.set_xlim(ax.get_xlim())
    if position[0] == "left":
        x = axx.get_xlim()[0]
    elif position[0] == "right":
        x = axx.get_xlim()[1] - duration
    else:
        x = position[0]

    # if position[1]=="top":
    #    y = axx.get_ylim()[1]
    # else:
    #    y = axx.get_ylim()[0]
    y = position[1]
    if "lw" not in linekwargs:
        linekwargs["lw"] = 2
    if "color" not in linekwargs:
        linekwargs["color"] = "black"
    #     print ([x,x+duration],[y]*2,**linekwargs)
    ln = Line2D([x, x + duration], [y] * 2, solid_capstyle="butt",**linekwargs)
    ln.set_clip_on(False)
    axx.add_artist(ln)
    if unit == "s":
        txt = "%gs"%duration
    if unit == "min":
        txt = " %gmin"%(duration / 60)
    fontkwargs.update({"ha": "center"})
    axx.text(x + duration / 2, y + text_vert_offset, txt, **fontkwargs)
    axx.set_ylim(-.1, 1.1)
    mystyle_axes(axx)

def show_line_across_axes(list_of_coordinates_with_axes, fig, **linekwargs):
    xy = []
    for x,y,ax in list_of_coordinates_with_axes:
        xy += [
            fig.transFigure.inverted().transform(
                ax.transData.transform((x,y))
            )
        ]
    x,y = np.array(xy).T
    fig.add_artist(Line2D(x,y,**linekwargs))

import matplotlib.collections as mcoll
import matplotlib.path as mpath

def colorline(
    x, y, z=None, cmap=plt.get_cmap('copper'), norm=plt.Normalize(0.0, 1.0),
        linewidth=3, alpha=1.0):
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    """

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])

    z = np.asarray(z)

    segments = make_segments(x, y)
    lc = mcoll.LineCollection(segments, array=z, cmap=cmap, norm=norm,
                              linewidth=linewidth, alpha=alpha)

    ax = plt.gca()
    ax.add_collection(lc)

    return lc


def make_segments(x, y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

def plot_closeup_traces(ax, regions, rois, closeup, **tracekwargs):
    regions.plotTraces(rois,axs=ax,labels=False,protocol=False, **tracekwargs)
    for l in ax.lines:
        x,y = l.get_data()
        fltr = (x>closeup[0]) & (x<closeup[1])
        l.set_data((x[fltr], y[fltr]))
    ax.relim()

def paper_plot(Events, regions,
               roi_ax_data = dict(),
               events_ax_data = dict(),
               protocol = None,
               format = "A4",
               draft = False,
               axroidim = 2,
               axhexHeight = 2.6,
               trace_ax_data = None,
               caxratio = .8,
               dpi = 72
               ):
    from matplotlib.ticker import AutoMinorLocator
    Events = Events.sort_values("z_max").copy()
    if protocol is None:
        protocol = regions.protocol
    figwidth, figheight = {
        "A4": (2*8.27,11.69*2),
        "A5": (11.69, 2*8.27),
        "Letter": (8.5*2,11.*2),
    }[format]
    fig_discard = plt.figure(figsize = (figwidth, figheight), dpi = dpi)
    fig = plt.figure(figsize=(figwidth, figheight), dpi = dpi)
    if draft:
        ruler(fig, margin = 1)
        # ruler(fig_discard, margin = 1)
    axtrWidth = figwidth - 5-axroidim

    if "discr_qty" in roi_ax_data and "threshold" in roi_ax_data:
        validRois = regions.df.query(f"{roi_ax_data['discr_qty']}>={roi_ax_data['threshold']}").index
        discard_Events = Events[~Events.roi.isin(validRois)].copy()
        valid_Events = Events[Events.roi.isin(validRois)].copy()
    else:
        roi_ax_data = {"discr_qty": "activity", "threshold":0}
        valid_Events = Events.copy()
        discard_Events = Events[:0].copy()

    ## Rois plot
    axrois = fig.add_axes([(1.2)/figwidth, 1-(axroidim+1.2)/figheight, axroidim/figwidth, axroidim/figheight], facecolor="none", label="rois")
    axcolor = fig.add_axes([(1.2) / figwidth, 1-(axroidim*1.57+1.2)/figheight, axroidim / figwidth, axroidim * .4 / figheight], label="roi colors")
    axc1, h = plot_rois_colored_acc_quantity(regions, [axrois, axcolor], **roi_ax_data)
    axc1.set_label("roi colorbar")



    # Events plot
    axhex = fig.add_axes([(3+axroidim)/figwidth, 1-(axhexHeight+1.2)/figheight, axtrWidth/figwidth, axhexHeight/figheight], facecolor="none", label="events")
    caxhex = fig.add_axes([(3+axroidim+axtrWidth+.2)/figwidth, 1-(axhexHeight*(1-(1-caxratio)/2)+1.2)/figheight, .2/figwidth, axhexHeight*caxratio/figheight], facecolor="none",label="events cax")
    axhexes = [axhex]
    caxhexes = [caxhex]
    eventsList = [valid_Events]
    if len(discard_Events):
        axhexes = [fig_discard.add_axes(axhex.get_position())] + axhexes
        caxhexes = [fig_discard.add_axes(caxhex.get_position())] + caxhexes
        eventsList = [discard_Events] + eventsList
    if "style" in events_ax_data:
        events_plot = events_ax_data.pop("style")
    else:
        events_plot = "scatter"

    for ax_, cax_, Events in zip(axhexes, caxhexes, eventsList):
        if events_plot == "scatter":
            defaults = dict(s = 10, linewidth = .3, edgecolor = 'k', vmax = 10, cmap = "magma_r")
            defaults.update(events_ax_data)
            pth = ax_.scatter(
                Events["peakpoint"],
                Events["halfwidth"],
                c = Events["z_max"],
                **defaults
            )
            ax_.set_yscale("log")
            ax_.set_xlim(-10,regions.time[-1]+10)
            clb = plt.colorbar(pth, cax=cax_)
            cax_.text(0,1.1,r"$z$-score", transform = cax_.transAxes)

        elif events_plot == 'hexbin':
            if "granularity" in events_ax_data:
                granularity = events_ax_data.pop("granularity")
            else:
                granularity = 40
            if "cmap" not in events_ax_data:
                events_ax_data["cmap"] = "hot"
            pth = ax_.hexbin(
                Events["peakpoint"],
                Events["halfwidth"],
                yscale="log",
                mincnt = 1,
                gridsize=(np.array([axtrWidth*1.2, axhexHeight])*granularity).astype(int),
                extent = (regions.time[0]-10, regions.time[-1]+10, -1.05,2.4),
                **events_ax_data
            )
            ax_.set_xlim(-10, regions.time[-1] + 10)
            clb = plt.colorbar(pth,cax=cax_)
            cax_.text(0, 1.1, "bin count", transform = cax_.transAxes)
        else:
            raise ValueError(f"Only 'scatter' and 'hexbin' allowed as styles in `events_ax_data`")
        ax_.spines['left'].set_position(("outward",5))
        ax_.spines['top'].set_position(("outward", 5))
        ax_.set_yticks([.1, 1, 10, 100])
        ax_.set_yticklabels([.1, 1, 10, 100])
        ax_.set_ylabel(r"halfwidth, $\tau_{1/2}$ [s]")
        ax_.text(0, 1.1, "time [s]  ", transform = ax_.transAxes, ha = "right", va="center")
        ax_.xaxis.set_minor_locator(AutoMinorLocator())
        mystyle_axes(ax_, retain = ["left","top"], bounded=[False,False])

        axp = protocol.plot_protocol( ax=ax_, position = "bottom", )

    ## Trace plot
    if trace_ax_data is not None and len(trace_ax_data):
        axtrHeight = trace_ax_data.pop("height",1)/figheight
        axtrOffset = trace_ax_data.pop("spacing",.1)/figheight
        roi_label_position = trace_ax_data.pop("roi_label_position",-.01)
        showRois = trace_ax_data.pop('show_rois')
        axppos = axp.get_position()
        axtr = fig.add_axes([axppos.x0,axppos.y0-axtrOffset-axtrHeight, axtrWidth/figwidth, axtrHeight], label="traces", facecolor='none')
        regions.plotTraces(showRois, axs = axtr, protocol = False, labels = True, **trace_ax_data)
        axtr.set_xticks(axhex.get_xticks())
        axtr.set_xticklabels(axhex.get_xticklabels())
        axtr.set_xticks(axhex.get_xticks(minor=True), minor=True)
        axtr.set_xlim(axhex.get_xlim())
        mystyle_axes(axtr, retain=['bottom'], bounded=[False])

        regions.plotEdges(showRois, color="k", ax=axrois, image=False)
        showRoiDF = pd.DataFrame(regions.df.loc[showRois, 'pixels'].apply(lambda xi: np.mean(xi, 0)))
        showRoiDF.columns = ["center"]
        showRoiDF["x0"] = [xy[1] for xy in showRoiDF['center']]
        showRoiDF["y0"] = [xy[0] for xy in showRoiDF['center']]
        showRoiDF["tan"] = (showRoiDF["y0"] - regions.image.shape[0] / 2) / showRoiDF["x0"]
        rpos = np.arange(0, len(showRois)) * regions.image.shape[0]/10
        rpos += -rpos.mean() + max(axrois.get_ylim()) / 2
        for yr, roi in zip(rpos, showRoiDF.sort_values("tan").index):
            x = regions.image.shape[1] * roi_label_position
            # axrois.plot([x, showRoiDF.loc[roi, "x0"]], [yr, showRoiDF.loc[roi, "y0"]], "k-", lw = .7)
            # axrois.plot(x, yr, showRoiDF.loc[roi, "marker"], ms = 5, mfc = "w", mec = "k", mew = .7)
            axrois.annotate(str(roi),
                            xy = showRoiDF.loc[roi, ["x0","y0"]],
                            xytext = (x, yr),
                            arrowprops = dict(arrowstyle = "-"),
                            fontsize=plt.rcParams['font.size']*.8,
                            horizontalalignment = "right",
                            verticalalignment = "center",
                            )
    if len(discard_Events):
        ax = fig_discard.get_axes()[-1]
        ax.text(1,0,"\n\n above: discarded events (left of the threshold in the lower plot)",transform=ax.transAxes, ha='right', va='top')
        ax.text(1,0,"\n\n\n\n This message will need to cropped manually for embedding this figure into document.",transform=ax.transAxes, ha='right', va='top', size="x-small")
        pos = ax.get_position()
        fig_discard.add_artist(Line2D([1-pos.x0-pos.width,pos.x0+pos.width],[pos.y0-.05]*2,color='k'),)
        fig_discard.show()
        fig.show()
    else:
        plt.close(fig_discard)
        fig_discard = None
    axs = {ax._label: ax for ax in fig.get_axes() if ax._label is not None}
    return (fig, axs), fig_discard

def dummy():
    fig = plt.figure()
    plt.close(fig)

def plot_rois_colored_acc_quantity(regions, axs, discr_qty, cmap = "RdBu", vmin = None, vmax = None, hist_kwargs = None, rel_height = .2, rel_pos = 1.1, scaleFontSize=None, ceiling=None, threshold=None):
    axrois, axcolor = axs
    if vmin is None:
        vmin = np.nanpercentile(regions.df[discr_qty],1)
    if vmax is None:
        vmax = np.nanpercentile(regions.df[discr_qty],95)
    if hist_kwargs is None:
        hist_kwargs = {"bins":50, "color":myGrey}
    if scaleFontSize is None:
        scaleFontSize = plt.rcParams['font.size']*.8
    regions.color_according_to(discr_qty, cmap=cmap, vmin = vmin, vmax = vmax)
    regions.plotEdges(ax=axrois, scaleFontSize=scaleFontSize, separate=True, fill=True, alpha =1, lw=0, smoothness = 2)
    regions.plotEdges(ax=axrois,color="grey", image=False,lw=.3, separate=True, smoothness = 2)
    axrois.set_xticks([])
    axrois.set_yticks([])
    if ceiling is None:
        ceiling = 1.1 * vmax
    h = axcolor.hist(np.minimum(regions.df[discr_qty], ceiling*.98),**hist_kwargs)[0]
    axc_pos = axcolor.get_position()
    fig = axrois.get_figure()
    axc1 = fig.add_axes([axc_pos.x0, axc_pos.y0+axc_pos.height*rel_pos, axc_pos.width, axc_pos.height*rel_height])
    min_, max_ = axcolor.get_xlim()
    axrange = max_ - min_
    min_ -= axrange*.3
    max_ += axrange*.3
    x = np.linspace(min_, max_, 200)
    x = (x-vmin)/(vmax-vmin)
    axc1.imshow(plt.cm.get_cmap(cmap)([x]), aspect="auto", extent = (min_, max_, 0, 1))
    axc1.set_xlim(axcolor.get_xlim())
    axc1.set_yticks([])
    axc1.set_xticks([])
    if threshold is not None:
        axcolor.axvline(threshold,color="r", ls="dotted", lw=1)
    return axc1, h


def check_roi_heterogeneity(Events, regions, quantity = "halfwidth", min_events = 10, transform = None,
                            boxplot_kwargs = dict(flierprops = dict(markersize = 1))
                            ):
    fig, axs = plt.subplots(2, 2, figsize = (20, 8), gridspec_kw={"width_ratios":[1,4]})
    if transform is None: transform = lambda xi: xi
    X = [transform(dfr[quantity].values) for roi, dfr in Events.groupby("roi") if len(dfr) > min_events]
    X = sorted(X, key = np.median)
    bx = axs[0,1].boxplot(X, positions = np.arange(len(X)), **boxplot_kwargs)
    axs[1,1].plot(list(map(len, X)))
    axs[1,1].set_xlim(axs[0,1].get_xlim())
    hetvar = "heterogeneity_var"
    x = []
    for roi in regions.df.index:
        ev = Events.query(f"roi=={roi}")
        if len(ev)>min_events:
            x += [np.median(transform(ev[quantity]))]
        else:
            x += [np.nan]
    regions.df[hetvar] = x
    plot_rois_colored_acc_quantity(regions, axs[::-1,0], hetvar, )

    return (fig, axs, bx)

def get_event_shapes(events, regions, n_traces_min = 10, plot = False):
    from scipy import signal
    timescaleDict = {it: t for (it, t), _ in events.groupby(["its", "ts"])}
    events = events.copy()
    events["trace"] = [islets.utils.get_evshape(row, regions) for ir, row in events.iterrows()]
    events["trace_points"] = events["trace"].apply(len)

    its = sorted(events.its.unique())
    nps = events["trace_points"].value_counts().sort_index()
    M = np.array([[len(events.query(f"trace_points=={n} and its=={its_}")) for n in nps.index] for its_ in its])
    if plot:
        plt.imshow(M,aspect="auto")
        plt.colorbar()
        plt.yticks(np.arange(len(its)), [timescaleDict[i] for i in its]);
        plt.xticks(np.arange(len(nps)), nps.index);
        plt.xlabel("trace points")
        plt.ylabel("filtering timescale index")
        plt.show()

    largest_values = np.argsort((-M).flat)
    nr = 4
    nc = 10
    Shapes = []
    if plot:
        fig, axs = plt.subplots(
            nr, nc, figsize = (nc * 2, nr * 2),
            sharex = False, sharey = True
        )
        ia = 0

    for ij in largest_values:
        i, j = ij // M.shape[1], ij % M.shape[1]
        value = M[i, j]
        if value < n_traces_min:
            continue
        if plot: ax = axs.flat[ia]
        n = nps.index[j]
        it = its[i]
        curTs = timescaleDict[it]
        traces = np.vstack(events.query(f"trace_points=={n} and its=={it}")["trace"].values)
        t = regions.showTime.get("%g" % curTs, regions.time)[:traces.shape[1]]
        X = np.zeros_like(traces)
        for j, tr in enumerate(traces):
            x = tr - tr.mean()
            x /= x.std()
            if j < len(X):
                X[j, :] = x
        x = X.mean(0)
        x -= np.percentile(x, 5)
        if plot:
            ax.plot(t, x, "C2", )
        dt = t[1] - t[0]
        xsmooth = islets.numeric.runningAverage(x,len(x)//10//2*2+1)
        pks, pks_d = signal.find_peaks(xsmooth, prominence = .2)
        w, h, t0, te = np.vstack(signal.peak_widths(xsmooth, pks, )) * dt
        h /= dt
        ix = np.argmax(pks_d["prominences"])
        midline = [t0[ix] + dt / 2, te[ix] + dt / 2], [h[ix]] * 2
        if plot:
            ax.plot(*midline, "C1.-")
            txt = f"ts: {curTs}Hz\nL: {n}\nN: {len(traces)}"
            txt += "\nhw: %.3g" % w[ix]
            ax.text(1, 1, txt, transform = ax.transAxes, va = "top", ha = "right")
        Shapes += [{
            "filtering timescale": curTs,
            "halfwidth": w[ix],
            "shape": (t, x, xsmooth),
            "midline": midline,
            "N traces": len(traces),
            "N points": n
        }]
        if plot:
            ia += 1
            if ia >= axs.size // 1: break
    if plot:
        for ax in axs.flat[ia:]:
            ax.remove()
    Shapes = pd.DataFrame(Shapes)
    if plot:
        return (Shapes, fig, axs)
    else:
        return (Shapes, )

def get_axsrow(ax, naxs, height_in_inches,  spacing = .1, rel_offset = -1.5):
    fig = ax.get_figure()
    figwidth, figheight = fig.get_size_inches()
    axpos = ax.get_position()
    axclspWidth = axpos.width/(naxs+(naxs-1)*spacing)
    axclspHeight = height_in_inches/figheight
    axclsp = [
        fig.add_axes([axpos.x0+(spacing+1)*axclspWidth*jx, axpos.y0+axclspHeight*rel_offset, axclspWidth, axclspHeight])
        for jx in range(naxs)
    ]
    return axclsp