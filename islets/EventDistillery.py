import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import median_abs_deviation
from tqdm import tqdm


def define_legs(events, legs):
    if "leg" in events.columns:
        del events["leg"]
#     if "tend" not in events.columns:
#         events["tend"] = events['t0']+events["halfwidth"]
    for leg in legs:
        t0,tend = legs[leg]
        indices = events.query(f"peakpoint>{t0} and peakpoint<{tend}").index
        events.loc[indices,"leg"] = [leg]*len(indices)
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
        verbose=False,
        smooth=None,
        filt_cutoff=2.,
        z_th=3,
        npass=3,
        debug=False,
        rescale_z="soft"
        ):
    if timescales is None:
        timescales = 2. ** np.arange(2, 8+1e-10, .25)
    timescales = np.unique([float("%.3g" % ts) for ts in timescales])
    timescales = timescales[timescales < regions.time[-1] / 3]
    timescales = timescales[timescales > 5/regions.Freq]
    if verbose:
        print("defined:", "  ".join(["%g" % ts for ts in timescales]))
    regions.timescales = timescales
    allEvents = []
    trange = range(len(timescales))
    if not verbose:
        trange = tqdm(trange)
    for its in trange:
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
        if rescale_z=="hard":
            rescale_factor = [ median_abs_deviation(z) * mad2std for z in regions.df["zScore_" + k] ]
            if verbose:
                print("correcting each z: mean %.3g, std %.3g"%(np.mean(rescale_factor), np.std(rescale_factor)))
            regions.df["zScore_" + k] = [ z/rescale_factor[i] for i,z in enumerate(regions.df["zScore_" + k]) ]
        if rescale_z=="soft":
            # corr = regions.df["zScore_" + k].apply(median_abs_deviation).median()*mad2std
            corr = regions.df["zScore_" + k].apply(np.std).median()
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
    if len(allEvents):
        allEvents = pd.concat(allEvents)
    else:
        allEvents = pd.DataFrame()
    if verbose:
        print ("#"*30+"\n"+"#"*10+"   DONE   "+"#"*10+"\n"+"#"*30)
    return allEvents

def distill_events(
        candidateEvents,
        regions,
        **kwargs):
    global iterf
    def iterf(roi_dfroi_):
        return distill_events_per_roi(roi_dfroi_[1], regions, **kwargs)[0]

    from multiprocessing import Pool
    from tqdm import tqdm

    groupby = candidateEvents.groupby("roi")

    DF_filt = []
    pool = Pool(10)
    try:
        with tqdm(total=len(groupby)) as pbar:
            for res in pool.imap_unordered(iterf, groupby):
                DF_filt += [res]
                pbar.update()
    finally:
        pool.close()
        pool.join()
    if len(DF_filt):
        DF_filt = pd.concat(DF_filt)
    else:
        DF_filt = pd.DataFrame()
    return DF_filt


def draw_event(row, ax, **kwargs):
    ax.barh(row.its, row.halfwidth, left=row.t0, align="center", color=row["color"], **kwargs)

# def draw_event(row, ax, edge=False,zorder=3, height=.8):
    # if edge:
    #     ax.plot(
    #         row['t0'] + row["halfwidth"] * np.array([0, 0, 1, 1, 0]),
    #         row["its"] + (np.array([0, 1, 1, 0, 0]) - .5) * height,
    #         color=row['color'],
    #         zorder=zorder,
    #         linewidth=.5,
    #     )
    # else:
    #     ax.fill(
    #         row['t0'] + row["halfwidth"] * np.array([0, 0, 1, 1, 0]),
    #         row["its"] + (np.array([0, 1, 1, 0, 0]) - .5) * height,
    #         color=row['color'],
    #         zorder=zorder,
    #         linewidth=0,
    #     )

def distill_events_per_roi(roiEvents,
                           regions,
                           plot=False,
                           opt_col = "z_max",
                           halfwidth_toll=.2,
                           discardCloseToEdge=True,
                           freqShow=5,
                           take_best=False,
                           plotSlows=None,
                           minconfidence=4,
                           small_halfwidth=2,
                           #movement_pvalue_threshold=.0,
                           require_contingency=False,
                           figsize=(13, 10),
                           plotTrace=slice(3),
                           ):
    # from scipy.stats import ttest_1samp
    roi = roiEvents["roi"].unique()
    assert len(roi) <= 1
    if len(roi)==0:
        return pd.DataFrame()
    roi = roi[0]
    if "tend" not in roiEvents.columns:
        roiEvents["tend"] = roiEvents["t0"] + roiEvents["halfwidth"]
    # hatches = ['//','\\',"++","**"]
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
        if not hasattr(regions, "showTime"): regions.showTime = {}
        fig, axs = plt.subplots(3, 1, figsize=figsize, sharex=True, sharey=True)
        from matplotlib.ticker import MultipleLocator
        for ix, row in roiEvents.iterrows():
            draw_event(row, axs[1])
    colorDict = {
        # "too short": plt.cm.Greys(.9),
        # "too close to beginning/end": plt.cm.Greys(.6),
        # "too few calls": plt.cm.Greys(.5),
        # "solitary": plt.cm.Greys(.3),
    }
    for ix, row in roiEvents.iterrows():
        t = regions.showTime.get("%g" % row.ts, regions.time)
        it0 = np.searchsorted(t, row.t0) - 1
        ite = np.searchsorted(t, row.t0 + row.halfwidth) + 1
        eventLength = ite - it0
        if eventLength < 3 and not take_best:
            status = "too short"
            roiEvents.loc[ix, "status"] = status
            if plot:
                if status not in colorDict:
                    colorDict[status] = axs[0].plot([], label=status)[0].get_color()
                # else:
                #     axs[0].plot([], label=status,c=colorDict[status])
                row = pd.Series(row)
                row.color = colorDict[status]
                # print (row.color)
                draw_event(row, axs[0])
                continue
        if discardCloseToEdge and ((row.t0-row.halfwidth/2 < regions.time[0]) or (row.tend+row.halfwidth/2 > regions.time[-1])):
            status = "too close to beginning/end"
            roiEvents.loc[ix, "status"] = status
            if plot:
                if status not in colorDict:
                    colorDict[status] = axs[0].plot([], label=status)[0].get_color()
                row = pd.Series(row)
                row.color = colorDict[status]
                draw_event(row, axs[0])
                continue
    g = nx.DiGraph()
    for isp, event in roiEvents.query("status=='ok'").sort_values("halfwidth").iterrows():
        toll = event.halfwidth * halfwidth_toll
        conflicts = roiEvents.query(f"t0>{event.t0 - toll}")
        conflicts = conflicts.query(f"t0<{event.t0 + toll}")
        if require_contingency and row.halfwidth>small_halfwidth:
            conflicts = conflicts.query(f"its>={event.its-1}")
            conflicts = conflicts.query(f"its<={event.its+1}")

        conflicts = conflicts.query(f"tend>{event.tend - toll}")
        conflicts = conflicts.query(f"tend<{event.tend + toll}")
        # conflicts = conflicts.query(f"index!={event.name}")
        for other in list(conflicts.index):
            g.add_edge(isp, other)
    if plot:
        ax = axs[1]
        lines = nx.draw_networkx_edges(g,
                         ax=ax,
                         pos={node: (
                         roiEvents.loc[node, "t0"] + .5 * roiEvents.loc[node, "halfwidth"], roiEvents.loc[node, "its"])
                              for node in g.nodes},
                         edge_color="grey",
                         width=.8,
                         arrows=False,
                         )
        lines.set_zorder(5)
        # ax.set_yticks(np.arange(len(timescales)));
        # ax.set_yticklabels(timescales);
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    ######################################################
    df_filt = []
    for ixs in nx.connected_components(g.to_undirected()):
        if len(ixs)==1 and not take_best: 
            status = "solitary"
            isp = list(ixs)[0]
            roiEvents.loc[isp, "status"] = status
            roiEvents.loc[isp, "attractor"] = isp
            if plot:
                if status not in colorDict:
                    colorDict[status] = axs[0].plot([], label=status)[0].get_color()
                # else:
                #     axs[0].plot([], label=status,c=colorDict[status])
                #     axs[0].plot([], label=status,c=colorDict[status])
                row = pd.Series(roiEvents.loc[isp])
                row.color = colorDict[status]
                draw_event(row, axs[0])
            continue
        row = pd.Series(roiEvents.loc[list(ixs)].sort_values(opt_col).iloc[-1])
        if not take_best:
            for col in ["t0","tend","z_max","color","coltrans"]:
                if col in row.index:
                    row[col] = np.median(np.vstack(roiEvents.loc[list(ixs), col]),axis=0)
                    if len(row[col])==1:
                        row[col] = row[col][0]
                    #row[col] = np.mean(roiEvents.loc[list(ixs), col],axis=0)
        row.halfwidth = row.tend-row.t0
        row["clique_size"] = len(ixs)
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
        #             draw_event(row, axs[0])
        #         continue

        if len(ixs)<minconfidence and row.halfwidth>small_halfwidth:
            status = "too few calls"
            roiEvents.loc[row.name, "status"] = status
            if plot:
                if status not in colorDict:
                    colorDict[status] = axs[0].plot([], label=status)[0].get_color()
                # else:
                #     axs[0].plot([], label=status,c=colorDict[status])
                row = pd.Series(row)
                row.color = colorDict[status]
                draw_event(row, axs[0])
            continue
        df_filt += [row]
    df_filt = pd.DataFrame(df_filt).sort_index()
    for col in ["loghalfwidth", "score", "status","attractor"]:
        if col in df_filt.columns:
            del df_filt[col]
    if plot:
        for ix, row in df_filt.iterrows():
            draw_event(row, axs[2])

        yl = axs[1].get_ylim()[::-1]
        for ax in axs:
            # ax.set_facecolor("k")
            ax.set_ylabel("filtering\ntimescale [s]")
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
            ax.text(.01,.95,txt, fontsize=12, transform=ax.transAxes, bbox=dict(facecolor="w",alpha=.7,linewidth=0),va="top")
            ax.yaxis.set_minor_locator(MultipleLocator(1))
        plt.subplots_adjust(hspace=0.01)
    else:
        axs = None
    if "attractor" in roiEvents.columns:
        roiEvents['attractor']= roiEvents['attractor'].astype(pd.Int64Dtype())


    if plot and plotTrace:
        nr = np.round(regions.Freq / freqShow)
        if nr > 1:
            nr = int(nr)
            from .numeric import rebin
            x = rebin(regions.time, nr)
            y = rebin(origtrace, nr)
        else:
            x, y = regions.time, origtrace
        axxs = [ax.twinx() for ax in axs[plotTrace]]
        axxs[0].get_shared_y_axes().join(*axxs)
        for axx in axxs:
            # axx.set_ylabel("light intensity")
            axx.plot(x, y, c=plt.cm.Greys(.9), lw=.7, zorder=-1)
            axx.set_yticks([])
            if plotSlows is not None:
                for ts in plotSlows:
                    axx.plot(regions.showTime.get("%g" % ts, regions.time), tsslows[ts])
    return df_filt, axs, roiEvents


def plot_events(events,
                ax=None,
                min_height=.3,
                allrois=None,
                cmap="turbo_r",
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
        h = .05 * len(events.roi.unique())+1
        fig = plt.figure(figsize=(15,h))
        ax = fig.add_axes([0,1/h,1,1-1/h])
        ax.set_facecolor("k")
        ax.set_yticks([])
        ax.set_xlabel("time [s]")
        for sp in ["left","right","top"]: ax.spines[sp].set_visible(False)
    if "color" not in events.columns:
        events["color"] = list(plt.cm.get_cmap(cmap)(events["coltrans"]))
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
    ax.set_xlim(0,None)
    return ax