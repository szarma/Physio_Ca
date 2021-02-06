import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import median_abs_deviation


def define_legs(events, legs):
    if "tend" not in events.columns:
        events["tend"] = events['t0']+events["halfwidth"]
    for leg in legs:
        t0,tend = legs[leg]
        indices = events.query(f"t0>{t0} and tend<{tend}").index
        events.loc[indices,"leg"] = leg
    events['leg'] = events['leg'].astype("category")

def event_diff(events_0, events_1, hw_toll=.2):
    indexmap = {}
    remainder = []
    for i0, row0 in events_0.iterrows():
        df = events_1.query(f"roi=={row0.roi}")
        df = df.query(f"t0>{row0.t0-hw_toll*row0.halfwidth}")
        df = df.query(f"t0<{row0.t0+hw_toll*row0.halfwidth}")
        df = df.query(f"tend>{row0.tend-hw_toll*row0.halfwidth}")
        df = df.query(f"tend<{row0.tend+hw_toll*row0.halfwidth}")
        if len(df):
            indexmap[i0] = set(df.index)
        else:
            remainder += [i0]
    return indexmap, remainder


def sequential_filtering(
        regions,
        timescales=None,
        verbose=True,
        smooth=None,
        filt_cutoff=2.,
        z_th=3,
        npass=3,
        debug=False,
        rescale_z=True
        ):
    if timescales is None:
        timescales = 2. ** np.arange(-1, 20, .25)
    timescales = np.unique([float("%.3g" % ts) for ts in timescales])
    timescales = timescales[timescales < regions.time[-1] / 3]
    timescales = timescales[timescales > 5/regions.Freq]
    if verbose:
        print("defined:", "  ".join(["%g" % ts for ts in timescales]))
    regions.timescales = timescales
    allEvents = []
    for its in range(len(timescales)):
        ts = timescales[its]
        k = "%g"%ts
        if verbose:
            print ("#"*30+f"\t {its}:  ts = {ts}s")
        # regions.infer_TwoParFit(ts=ts, verbose=verbose)
        regions.fast_filter_traces(ts,
                                   filt_cutoff=filt_cutoff,
                                   verbose=verbose,
                                   npass=npass
                                  )
        from .numeric import mad2std
        if "correctZ" in regions.df.columns:
            regions.df["zScore_"+k] = [regions.df.loc[ix,"zScore_"+k]/regions.df.loc[ix,"correctZ"] for ix in regions.df.index]
        if rescale_z:
            corr = regions.df["zScore_" + k].apply(median_abs_deviation).median()*mad2std
            if verbose:
                print("correcting z by", corr)
            regions.df["zScore_" + k] = [z / corr for z in regions.df["zScore_" + k]]
        regions.detect_events(ts, z_th=z_th, smooth=smooth, verbose=verbose, debug=debug)
        spikeDF = regions.events[k]
        spikeDF["ts"] = ts
        spikeDF["its"] = its
        for ttmp in allEvents:
            spikeDF.index += len(ttmp.index)
        allEvents += [spikeDF]
    allEvents = pd.concat(allEvents)
    if verbose: print ("#"*30+"\n"+"#"*10+"   DONE   "+"#"*10+"\n"+"#"*30)
    return allEvents


def distill_events(
        candidateEvents,
        regions,
        **kwargs):
#     from tqdm.notebook import tqdm
    from tqdm import tqdm
    # tqdm().pandas()
    DF_filt = []
    for roi,dfroi in tqdm(candidateEvents.groupby("roi")):
        df_filt = distill_events_per_roi(dfroi, regions, **kwargs)[0]
        DF_filt += [df_filt]
    return pd.concat(DF_filt)


def distill_events_per_roi(roiEvents,
                           regions,
                           plot=False,
                           halfwidth_toll=.2,
                           freqShow=5,
                           take_best=False,
                           plotSlows=None,
                           minconfidence=4,
                           small_halfwidth=2,
                           #movement_pvalue_threshold=.0,
                           require_contingency=False,
                           figsize=(13, 10),
                           ):
    # from scipy.stats import ttest_1samp
    roi = roiEvents["roi"].unique()
    assert len(roi) <= 1
    if len(roi)==0:
        return pd.DataFrame()
    roi = roi[0]
    if "tend" not in roiEvents.columns:
        roiEvents["tend"] = roiEvents["t0"] + roiEvents["halfwidth"]
    origtrace = regions.df.loc[roi, "trace"]
    tsslows = {float(col.split("_")[1]): regions.df.loc[roi, col] for col in regions.df.columns if "slower" in col}
    roiEvents = pd.DataFrame(roiEvents)
    roiEvents["status"] = "ok"
    if minconfidence==1:
        take_best=True
    if take_best:
        require_contingency = False
        small_halfwidth = np.inf
    if plot:
        fig, axs = plt.subplots(3, 1, figsize=figsize, sharex=True, sharey=True)
        from matplotlib.ticker import MultipleLocator
        def draw_spike(row, ax):
            ax.fill(
                row.t0 + row.halfwidth * np.array([0, 0, 1, 1, 0]),
                row.its + np.array([0, 1, 1, 0, 0]) * .8 - .5,
                color=row.color,
                zorder=-1
                # linewidth=.7,
            )
        for ix, row in roiEvents.iterrows():
            draw_spike(row, axs[1])
    colorDict = {}
    for ix, row in roiEvents.iterrows():
        t = regions.showTime.get("%g" % row.ts, regions.time)
        it0 = np.searchsorted(t, row.t0) - 1
        ite = np.searchsorted(t, row.t0 + row.halfwidth) + 1
        spikeLength = ite - it0
        if spikeLength < 3 and not take_best:
            status = "too short"
            roiEvents.loc[ix, "status"] = status
            if plot:
                if status not in colorDict:
                    colorDict[status] = axs[0].plot([], label=status)[0].get_color()
                row = pd.Series(row)
                row.color = colorDict[status]
                draw_spike(row, axs[0])
                continue
        if (row.t0-row.halfwidth/2 < regions.time[0]) or (row.tend+row.halfwidth/2 > regions.time[-1]):
            status = "too close to begging/end"
            roiEvents.loc[ix, "status"] = status
            if plot:
                if status not in colorDict:
                    colorDict[status] = axs[0].plot([], label=status)[0].get_color()
                row = pd.Series(row)
                row.color = colorDict[status]
                draw_spike(row, axs[0])
                continue
    g = nx.DiGraph()
    for isp, spike in roiEvents.query("status=='ok'").sort_values("halfwidth").iterrows():
        toll = spike.halfwidth * halfwidth_toll
        conflicts = roiEvents.query(f"t0>{spike.t0 - toll}")
        conflicts = conflicts.query(f"t0<{spike.t0 + toll}")
        if require_contingency and row.halfwidth>small_halfwidth:
            conflicts = conflicts.query(f"its>={spike.its-1}")
            conflicts = conflicts.query(f"its<={spike.its+1}")

        conflicts = conflicts.query(f"tend>{spike.tend - toll}")
        conflicts = conflicts.query(f"tend<{spike.tend + toll}")
        # conflicts = conflicts.query(f"index!={spike.name}")
        if len(conflicts.index)==1 and not take_best:
            status = "unique"
            roiEvents.loc[isp, "status"] = status
            if plot:
                if status not in colorDict:
                    colorDict[status] = axs[0].plot([], label=status)[0].get_color()
                row = pd.Series(spike)
                row.color = colorDict[status]
                draw_spike(row, axs[0])
        for other in list(conflicts.index):
            g.add_edge(isp, other)
    if plot:
        ax = axs[1]
        nx.draw_networkx(g,
                         ax=ax,
                         with_labels=False,
                         pos={node: (
                         roiEvents.loc[node, "t0"] + .5 * roiEvents.loc[node, "halfwidth"], roiEvents.loc[node, "its"])
                              for node in g.nodes},
                         node_size=0,
                         edge_color="grey",
                         width=.8,
                         arrows=False,
                         )
        # ax.set_yticks(np.arange(len(timescales)));
        # ax.set_yticklabels(timescales);
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    ######################################################
    df_filt = []
    for ixs in nx.connected_components(g.to_undirected()):
        if len(ixs)==1 and not take_best: continue
        row = pd.Series(roiEvents.loc[list(ixs)].sort_values("height").iloc[-1])
        if not take_best:
            for col in ["t0","tend","height","color","coltrans"]:
                if col in row.index:
                    row[col] = np.median(np.vstack(roiEvents.loc[list(ixs), col]),axis=0)
                    if len(row[col])==1:
                        row[col] = row[col][0]
                    #row[col] = np.mean(roiEvents.loc[list(ixs), col],axis=0)
        row.halfwidth = row.tend-row.t0
        for ix in ixs:
            if ix != row.name:
                roiEvents.loc[ix,"status"]="merged"
                roiEvents.loc[ix,"attractor"]=row.name
        # TODO: movement detection
        # dt = roiEvents.loc[ixs].sort_values("ts")[["t0","tend"]].diff(axis=0).dropna().values.flatten()
        # if row.halfwidth>small_halfwidth and p:
        #     tstat = ttest_1samp(dt,popmean=0)
        #     if tstat.pvalue<movement_pvalue_threshold:
        #         status = "movement"
        #         roiEvents.loc[row.name,"status"] = status
        #         if plot:
        #             if status not in colorDict:
        #                 colorDict[status] = axs[0].plot([], label=status)[0].get_color()
        #             row = pd.Series(row)
        #             row.color = colorDict[status]
        #             draw_spike(row, axs[0])
        #         continue

        if len(ixs)<minconfidence and row.halfwidth>small_halfwidth:
            status = "too few calls at large timescale"
            roiEvents.loc[row.name, "status"] = status
            if plot:
                if status not in colorDict:
                    colorDict[status] = axs[0].plot([], label=status)[0].get_color()
                row = pd.Series(row)
                row.color = colorDict[status]
                draw_spike(row, axs[0])
            continue
        df_filt += [row]
    df_filt = pd.DataFrame(df_filt).sort_index()
    for col in ["loghalfwidth", "score", "status","attractor"]:
        if col in df_filt.columns:
            del df_filt[col]
    if plot:
        for ix, row in df_filt.iterrows():
            draw_spike(row, axs[2])
        nr = np.round(regions.Freq / freqShow)
        if nr > 1:
            nr = int(nr)
            from .numeric import rebin
            x = rebin(regions.time, nr)
            y = rebin(origtrace, nr)
        else:
            x, y = regions.time, origtrace
        axxs = [ax.twinx() for ax in axs]
        axxs[0].get_shared_y_axes().join(*axxs)
        for axx in axxs:
            axx.set_ylabel("light intensity")
            axx.plot(x, y, c="navy", lw=.7)
            if plotSlows is not None:
                for ts in plotSlows:
                    axx.plot(regions.showTime.get("%g" % ts, regions.time), tsslows[ts])

        yl = axs[1].get_ylim()[::-1]
        for ax in axs:
            # ax.set_facecolor("k")
            ax.set_ylabel("filtering timescale [s]")
            if hasattr(regions, "timescales"):
                timescales = regions.timescales
                # ax.set_yticks(np.arange(len(timescales)))
                # ax.set_yticklabels(["%g"%ts if int(np.log2(ts))==np.log2(ts) else "" for ts in timescales]);
                ax.set_yticks([jt for jt,ts in enumerate(timescales) if int(np.log2(ts))==np.log2(ts)])
                ax.set_yticklabels(["%g"%ts for ts in timescales if int(np.log2(ts))==np.log2(ts)]);
            else:
                ax.set_yticks(roiEvents.its.unique());
                ax.set_yticklabels(["%g"%ts for ts in roiEvents.ts.unique()]);
            ax.set_xlim(roiEvents.t0.min() - 1, roiEvents.tend.max() + 1)
            ax.set_ylim(yl)
        axs[2].set_xlabel("time [s]")
        if len(colorDict):
            lg = axs[0].legend(loc=1)
            fr = lg.get_frame()
            fr.set_facecolor("w")
            fr.set_alpha(.7)
        for ax, txt in zip(axs, ["discarded","all","distilled"]):
            ax.text(.01,.9,txt, fontsize=14, transform=ax.transAxes, bbox=dict(facecolor="w",alpha=.7,linewidth=0))
            ax.yaxis.set_minor_locator(MultipleLocator(1))
        plt.subplots_adjust(hspace=0.01)
    else:
        axs = None
    if "attractor" in roiEvents.columns:
        roiEvents['attractor']= roiEvents['attractor'].astype(pd.Int64Dtype())
    return df_filt, None, roiEvents


def plot_events(events,
                ax=None,
                min_height=.3,
                allrois=None,
                cmap="turbo",
                modify=False,
                **hillkwargs):
    if not modify:
        events = events.copy()
    if "tend" not in events.columns:
        events["tend"] = events["t0"]+events['halfwidth']
    if "coltrans" not in events.columns:
        from .numeric import hillCurve
        events["coltrans"] = hillCurve(events.halfwidth,**hillkwargs)
    if ax is None:
        h = .15 * len(events.roi.unique())+1
        fig = plt.figure(figsize=(15,h))
        ax = fig.add_axes([0,1/h,1,1-1/h])
        ax.set_facecolor("k")
        ax.set_yticks([])
        ax.set_xlabel("time [s]")
        for sp in ["left","right","top"]: ax.spines[sp].set_visible(False)
    if "color" not in events.columns:
        events["color"] = list(plt.cm.get_cmap(cmap)(events.coltrans))
    if allrois is None:
        allrois = events.roi.unique()
    ir=0
    for roi in allrois:
        evroi = events.query(f"roi=={roi}")
        for ix,row in evroi.sort_values("halfwidth", ascending=False).iterrows():
            h = (row.coltrans+min_height*(1+min_height))/(1+min_height)
            h = min(1,h)
            ax.fill(
                row.t0+row.halfwidth*np.array([0,0,1,1,0]),
                -(-.5+np.array([0,1,1,0,0]))*h+ir,
                color=row.color,
                alpha = row.alpha if "alpha" in row.index else 1
            )
        ir+=1
    return ax