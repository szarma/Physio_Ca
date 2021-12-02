import matplotlib.pyplot as plt
from scipy.stats import ttest_1samp
import numpy as np
import islets
import pandas as pd
from scipy.stats import distributions as dst
from matplotlib.lines import Line2D


def histogram_of_hw(Events, hwbinEdges, legDict, legColorDict, ax=None, orientation="bottom", hist_scale="left", control_line=False, n_sigma=False):
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
                    [0] + list(np.repeat(hmin + n_sigma*hmin_err, 2)) + [0],
                    [0] + list(np.repeat(hmin - n_sigma*hmin_err, 2)) + [0],
                    color = legColorDict[leg],
                    alpha = .1,
                )

            ax.set_xscale("log")
            if hist_scale == "top":
                ax.set_ylim(h.max() * 1.05, 0)
            ax.set_xlim(hwbinEdges[0], hwbinEdges[-1])
            ax.set_ylabel("events/min")
            if control_line:
                x = np.sqrt(hwbinEdges[:-1]*hwbinEdges[1:])
                ax.step(x, 60* dst.norm.sf(Events.z_max.min()) / x, #*np.diff(hwbinEdges),
                        c = legColorDict[leg], ls = "--", where="mid")
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


def boxplots_per_leg(roi_stats,legDict,legColorDict, pvalue=None, ax=None, **bx_kwargs):
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    for ileg,leg in enumerate(legDict):
        box_kwargs = dict(
            notch = True,
            showfliers = False,
            widths = .4,
            whis = (5, 95)
        )
        if bx_kwargs is not None:
            box_kwargs.update(bx_kwargs)
        bxs = ax.boxplot(roi_stats[leg], positions=[ileg], **box_kwargs)
        for el in bxs:
            for ln in bxs[el]:
                ln.set_color(legColorDict[leg])
    if pvalue is not None:
        for kp in pvalue:
            i,j = kp
            ki = list(legDict)[i]
            kj = list(legDict)[j]
            x = roi_stats[kj]-roi_stats[ki]
            #std = x.std(ddof=1)
            #mu = x.mean()
            #sem = std/np.sqrt(len(x))
            #assert sem==pd.Series(x).sem()
            #p = dst.norm.sf(mu/sem)
            ypos = pvalue[kp].get("ypos",np.percentile(roi_stats[[ki,kj]],99.7))
            alternative = pvalue[kp].get("alternative","two-sided")
            p = ttest_1samp(x, popmean=0, alternative=alternative).pvalue
            ax.plot([i,j],[ypos]*2,color="k")
            if p<.005:
                plabel = r"   $p < 10^{%i}$"%int(np.log10(p))
            else:
                plabel = r"   $p = %.1g$"%p
            ax.text((i+j)/2,ypos,plabel, ha="center",va="bottom")
    ax.set_xticks(range(len(legDict)))
    ax.set_xticklabels(list(legDict))

def fancy_barchart(ax, df, c=(.4,) * 3, orientation="horizontal", plabel_offset=.1):
    c = df.get("color", c)
    y = df["mean"].values
    yerr = df["std"].values
    nticks = len(df.index)
    if orientation == "horizontal":
        ax.barh(range(nticks), y, color = c)
        ax.errorbar(y, range(nticks), xerr = yerr, ls = 'none', c = c)
        ax.set_yticks([])
        for jt,tk in enumerate(df.index):
            ax.text(0, jt, tk+" ", ha="right",va="center")
    elif orientation == "vertical":
        ax.bar(range(nticks), y, color = c)
        ax.errorbar(range(nticks), y, yerr = yerr, ls = 'none', c = c)
        ax.set_xticks([])
    else:
        raise ValueError("orientation can only be 'horizonal' or 'vertical'")

    comparisons = [(i, j) for i in range(len(df.index)) for j in range(i + 1, len(df.index))]
    comparisons = sorted(comparisons, key = lambda ij: ((y + yerr)[list(ij)].max(), np.diff(ij)[0]))

    offset = max(y + 2*yerr) * plabel_offset
    xpos0 = 0

    for i, j in comparisons:
        dy = y[j] - y[i]
        sigma = (yerr[i] ** 2 + yerr[j] ** 2) ** .5
        z = np.abs(dy) / sigma
        p = dst.norm.sf(z)
        xpos = max([y[k] + 2 * yerr[k] for k in [i, j]])+offset/2
        while xpos <= xpos0 + offset:
            xpos += offset
        # print(i, j, xpos)
        xpos0 = xpos
        ln = Line2D([xpos - offset / 5, xpos, xpos, xpos - offset / 5], [i, i, j, j], c = "k", lw=.7)
        ln.set_clip_on(False)
        ax.add_line(ln)
        if p < .01:
            txt = "%.0e"%p
            txt,exponent = txt.split("e")
            txt = txt.split(".")[0]
            if txt=="1":
                plabel = r" $p = 10^{%i}$"%(int(exponent))
            else:
                plabel = r" $p = %s{\cdot}10^{%i}$"%(txt,int(exponent))
            # plabel = r" $p < 10^{%i}$" % int(np.log10(p))
            # plabel = "*"*max(int(np.log10(p)),3)
            # plabel = r" $%i\sigma$"%int(z)
        else:
            plabel = r" $p = %.1g$" % p
        if orientation=="horizontal":
            ax.text(xpos, (i + j) / 2, plabel,
                    # rotation="90",
                    va="center"
                    )

def emphasize_region(ax,x,y,extend=(0,0),**line_kwargs):
    xb, xe = x
    yb, ye = y
    assert xe>xb
    assert ye>yb

    dx = xe - xb
    dx = extend[0] * dx
    xe = xe + dx
    xb = xb - dx

    dy = xe - xb
    dy = extend[1] * dy
    ye = ye + dy
    yb = yb - dy
    ln = Line2D(
        [xb,xe,xe,xb,xb],
        [yb,yb,ye,ye,yb],
        **line_kwargs
    )
    ln.set_clip_on(False)
    ax.add_line(ln)

def plot_events(Events, regions, timeUnits="s", plottype="scatter", only_legs=False, protocol=None, legColorDict=None, axheight = 2, offset = 2,**kwargs):
    from .general_functions import  mystyle_axes
    from .utils import add_protocol

    if protocol is None:
        protocol = regions.protocol.copy()

    duration = regions.time[-1] - regions.time[0]

    if only_legs and "leg" in Events.columns:
        # take only events that happen during one of the specified legs
        Events = Events[~Events.leg.isna()].copy()

    ### creation and formatting of figure
    axwidth = duration / 500

    figwidth, figheight = axwidth + 2 * offset, axheight + 2 * offset
    fig = plt.figure(figsize = (figwidth, figheight))
    ax = fig.add_axes([offset / figwidth, offset / figheight, axwidth / figwidth, axheight / figheight])
    axp = fig.add_axes([offset / figwidth, (axheight + 1.1 * offset) / figheight, axwidth / figwidth,
                        offset / 10 / figheight * len(protocol.compound.unique())])
    axp.get_shared_x_axes().join(ax, axp)
    mystyle_axes(axp)
    output = [fig, (ax, axp)]
    if timeUnits[0]=="m":
        # duration = duration/60
        ax.set_xlabel("time [min]")
        for col in ["t_begin", "t_end"]:
            protocol[col] = protocol[col]/60
    elif timeUnits[0]=="s":
        ax.set_xlabel("time [s]")
    ax.set_ylabel(r"$\tau_{1/2}$, halfwidth [s]")
    add_protocol(axp, protocol)
    ## plotting
    if plottype=="hexbin":
        cax = fig.add_axes([(1.1*offset + axwidth) / figwidth, (offset) / figheight, .2 / figwidth, axheight/figheight])
        x, y = Events["peakpoint"].values, Events["halfwidth"].values
        if timeUnits == "min":
            x = x.copy() / 60
        if "bins" not in kwargs:
            kwargs.update({"vmin":0})
        hx = ax.hexbin(x,y,yscale="log", mincnt = 1, gridsize=(1+int(duration/20), 50),**kwargs)
        cax = plt.colorbar(hx, cax=cax)
        cax.set_label("bin count")
        output[1] = output[1] + (cax,)

    elif plottype=="scatter":
        ax.set_yscale("log")
        for leg,ev in list(Events.groupby("leg"))+[(None, Events[Events.leg.isna()])]:
            x,y = ev["peakpoint"].values.copy(), ev["halfwidth"]
            if timeUnits=="min":
                x = x/60
            ax.plot(x,y,".",c="lightgrey" if leg is None else legColorDict[leg],**kwargs)
    else:
        raise ValueError("plottype can only be 'scatter' or 'hexbin")

    mystyle_axes(ax,["bottom","left"],bounded = [True,False])
    return output

def get_events_per_min_per_nrois(Events, hwRegions):
    nrois = len(Events.roi.unique())
    ev_pm_par = []  # number of events per min per active roi
    for leg,ev in Events.groupby("leg"):
        duration = ev["peakpoint"].max() - ev["peakpoint"].min()
        for hwr in hwRegions:
            hwb = hwRegions[hwr]
            ev = Events.query(f"leg=='{leg}' and halfwidth>{hwb[0]} and halfwidth<{hwb[1]}")

            ev_pm_par += [{
                "hw_region": hwr,
                "leg": leg,
                "mean": len(ev) / duration / nrois * 60,
                "std": len(ev) ** .5 / duration / nrois * 60,
            }]
    ev_pm_par = pd.DataFrame(ev_pm_par)
    return ev_pm_par