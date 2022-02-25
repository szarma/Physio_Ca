import matplotlib.pyplot as plt
from scipy.stats import ttest_1samp
import numpy as np
import islets
import pandas as pd
from scipy.stats import distributions as dst
from matplotlib.lines import Line2D
import warnings

from IPython.display import display, HTML
import statsmodels.api as sm

myGrey = (.2,)*3
default_fontsize = 7
figlabeldict = dict(fontsize=8,fontweight="bold")

def showLineScan(linescan, axs, pixel_limits, Tind = 1, slice_ = None, offset = 1):
    data = linescan.data[slice_]
    x = islets.numeric.rebin(data,100,1, func=np.sum)
    for j in range(x.shape[0]):
        x[j] -= np.percentile(x[j],10)
    ax = axs[0]
    distance = 10
    physicalSize = linescan.metadata["pxSize"]*data.shape[0]
    txtOffset = physicalSize-1.5*distance
    tmax = linescan.time[-1]
    ax.imshow(x, cmap="turbo", vmin = 0, vmax = x.max(), aspect="auto", extent = (0,tmax, 0, physicalSize), origin = "lower")

    ax.plot([tmax*.05]*2,[txtOffset,txtOffset+distance],"lightgrey")
    ax.text(tmax*.05,txtOffset+distance/2," "+str(distance)+"µm", va="center", color="w")
    # ax.imshow(x, cmap="Greys", vmin = 0, vmax = x.max(), aspect="auto")
    ax.plot([tmax*.05, tmax*.03+Tind],[.03*physicalSize]*2,"lightgrey")
    ax.text(tmax*.05+Tind/2,.03*physicalSize,str(Tind)+"s\n", ha="center", va="center", color="w", fontsize=7)
    ax1 = axs[1]
    off = 0
    for px0,px1 in pixel_limits:
        emphasize_region(ax,ax.get_xlim(), [px0*physicalSize/x.shape[0],px1*physicalSize/x.shape[0]], color="lightgrey", extend = (-.01,0),lw=1)
        ax1.plot(np.linspace(0,tmax,x.shape[1]), x[px0:px1].mean(0) + off, lw=.5, color=myGrey)
        ar = ax.arrow(tmax*1.05, np.mean([px0, px1])*physicalSize/x.shape[0], tmax*.2, 0, color=myGrey, lw=1, clip_on=False, head_width = tmax, head_length = .2)
        ar.set_clip_on(False)
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
    ("diazoxide","100uM"): "xkcd:pale rose",
    ("ryanodine","100nM"): "xkcd:very light purple",
    ("ryanodine","100uM"): "xkcd:light purple",
    ("acetylcholine","10uM"): "xkcd:light pastel green",
    ("caffeine","10mM"): "xkcd:light brown",
    ("xestospongin","3uM"): "goldenrod",
}


def beautify_protocol(protocol):
    from .protocol import Protocol
    colors = []
    for i in protocol.index:
        colors += [protocolColorScheme[(protocol.loc[i, "compound"], protocol.loc[i, "concentration"])]]
    protocol["color"] = colors
    for comp, df in protocol.groupby("compound"):
        units = df["unit"].unique()
        if len(units) == 1:
            unit = units[0]
            # protocol.loc[df.index, "compound"] = [f"{comp} [{unit}]"] * len(df)
            protocol.loc[df.index, "concentration"] = [df["concentration"].iloc[0]] + [f"{c}" for c in df["conc"][1:]]
    # protocol = protocol.replace("Ca [mM]", r"Ca${}^{2\!+}$[mM] ")
    protocol = protocol.replace("Ca", r"Ca${}^{2\!+}$ ")
    protocol = Protocol(protocol)
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


def get_rate_effect(events, ref_leg=None, legs=None, groups = "experiment"):
    events = events.copy()
    if legs is None:
        legs = events["leg"].dropna().unique()
    if ref_leg is None:
        ref_leg = events[~events["leg"].isna()].iloc[0]['leg']
    events["log10_epmpr"] = np.log10(events["epmpr"])
    model = sm.MixedLM.from_formula(f"log10_epmpr ~ C(leg, Treatment('{ref_leg}')) + 1", data = events,
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

    dy = xe - xb
    dy = extend[1] * dy
    ye = ye + dy
    yb = yb - dy
    if clip_on:
        ax.plot(
            [xb, xe, xe, xb, xb],
            [yb, yb, ye, ye, yb],
            **line_kwargs
        )
    else:
        ln = Line2D(
            [xb, xe, xe, xb, xb],
            [yb, yb, ye, ye, yb],
            **line_kwargs
        )
        ln.set_clip_on(False)
        ax.add_line(ln)


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
        for leg, ev in list(Events.groupby("leg")) + [(None, Events[Events.leg.isna()])]:
            x, y = ev["peakpoint"].values.copy(), ev["halfwidth"]
            if timeUnits == "min":
                x = x / 60
            ax.plot(x, y, ".", c = "lightgrey" if leg is None else legColorDict[leg], **kwargs)
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


from typing import Dict, List, Optional, Tuple, Iterable


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
    axx.plot([x, x + duration], [y] * 2, **linekwargs)
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

def plot_closeup_traces(ax, regions, rois, closeup, **tracekwargs):
    regions.plotTraces(rois,axs=ax,labels=False,protocol=False, **tracekwargs)
    for l in ax.lines:
        x,y = l.get_data()
        fltr = (x>closeup[0]) & (x<closeup[1])
        l.set_data((x[fltr], y[fltr]))
    ax.set_xlim(closeup)
    ax.relim()

def plot_rois_colored_acc_quantity(regions, axs, discr_qty, cmap = "turbo", vmin = None, vmax = None, hist_kwargs = None, rel_height = .1, rel_pos = 1.05):
    axrois, axcolor = axs
    if vmin is None:
        vmin = np.nanpercentile(regions.df[discr_qty],1)
    if vmax is None:
        vmax = np.nanpercentile(regions.df[discr_qty],95)
    if hist_kwargs is None:
        hist_kwargs = {"bins":50, "color":myGrey}
    regions.color_according_to(discr_qty, cmap=cmap, vmin = vmin, vmax = vmax)
    regions.plotEdges(ax=axrois, scaleFontSize=6, separate=True, fill=True, alpha =1)
    regions.plotEdges(ax=axrois, separate=True, fill=False, image=False)
    axrois.set_xticks([])
    axrois.set_yticks([])
    h = axcolor.hist(regions.df[discr_qty],**hist_kwargs)[0]
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