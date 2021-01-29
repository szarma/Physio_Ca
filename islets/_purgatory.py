import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


def filter_events_per_roi(roiEvents,
                          regions,
                          plot=False,
                          halfwidth_toll=.33,
                          freqShow=5,
                          plotSlows=None,
                          ):
    roi = roiEvents["roi"].unique()
    assert len(roi)==1
    roi = roi[0]
    origtrace = regions.df.loc[roi,"trace"]
    tsslows = {float(col.split("_")[1]): regions.df.loc[roi,col] for col in regions.df.columns if "slower" in col}
    roiEvents = pd.DataFrame(roiEvents)
    roiEvents["status"] = "ok"

    if plot:
        fig, axs = plt.subplots(3,1,figsize=(13,14), sharex=True, sharey=True)

        def draw_spike(row, ax):
            ax.fill(
                row.t0+row.halfwidth*np.array([0,0,1,1,0]),
                row.its+np.array([0,1,1,0,0])*.8-.5,
                color=row.color,
                #linewidth=.7,
            )
        for ix, row in roiEvents.iterrows():
            draw_spike(row, axs[1])
    colorDict = {}
    for ix, row in roiEvents.iterrows():
        t = regions.showTime.get("%g"%row.ts, regions.time)
        it0 = np.searchsorted(t,row.t0)-1
        ite = np.searchsorted(t,row.t0+row.halfwidth)+1
        spikeLength = ite-it0
        if spikeLength<3:
            status = "too short"
            roiEvents.loc[ix,"status"] = status
            if status not in colorDict:
                colorDict[status] = axs[0].plot([],label=status)[0].get_color()
            row = pd.Series(row)
            row.color = colorDict[status]
            if plot: draw_spike(row,axs[0])
            continue
        slowtrace = tsslows[row.ts][it0:ite].copy()
        # if row.ts*4 in tsslows:
        #     slowtrace-=tsslows[row.ts*4][it0:ite]
        # else:
        #     if row.ts*2 in tsslows:
        #         slowtrace-=tsslows[row.ts*2][it0:ite]
        ddslow = np.diff(np.diff(slowtrace))
        ddorig = np.diff(np.diff(origtrace[it0:ite]))
        # if ddorig.mean()>0:
        if spikeLength>20 and ddslow.mean()>0.1 and ddslow.mean()>-ddorig.mean()*.5:
            print (ddslow.mean(), ddorig.mean())
            status = "spurious"
            roiEvents.loc[ix,"status"] = status
            if status not in colorDict:
                colorDict[status] = axs[0].plot([],label=status)[0].get_color()
            row = pd.Series(row)
            row.color = colorDict[status]
            if plot: draw_spike(row, axs[0])
            if plot: axs[0].text(row.t0+row.halfwidth/2, row.its+.3, "x", va="bottom", ha="center")
            continue

        tcur = t[it0:ite]
        bkg = tsslows[row.ts][it0]+np.arange(spikeLength)*(tsslows[row.ts][ite-1]-tsslows[row.ts][it0])/(spikeLength-1)
        purePeak = origtrace[it0:ite]-bkg
        amplitude = np.percentile(purePeak, 95)
        if amplitude/bkg.mean()<.1:
#         if amplitude<1.*np.abs(bkg[-1]-bkg[0]):
            status = "too steep"
            roiEvents.loc[ix,"status"] = status
            if plot:
                if status not in colorDict:
                    colorDict[status] = axs[0].plot([],label=status)[0].get_color()
                row = pd.Series(row)
                row.color = colorDict[status]
                draw_spike(row,axs[0])
                axs[0].text(row.t0+row.halfwidth/2, row.its+.3, ("/" if bkg[-1]>bkg[0] else "\\"), va="bottom", ha="center")

    if len(colorDict):
        axs[0].legend()
    g = nx.DiGraph()
    for isp, spike in roiEvents.query("status=='ok'").sort_values("halfwidth").iterrows():
        conflicts = roiEvents.query(  f"t0>{spike.t0   - spike.halfwidth*halfwidth_toll}")
        conflicts = conflicts.query(  f"t0<{spike.t0   + spike.halfwidth*halfwidth_toll}")
        conflicts = conflicts.query(f"tend>{spike.tend - spike.halfwidth*halfwidth_toll}")
        conflicts = conflicts.query(f"tend<{spike.tend + spike.halfwidth*halfwidth_toll}")
        conflicts = conflicts.query(f"index!={spike.name}")
        g.add_node(isp)
        for other in list(conflicts.index):
            g.add_edge(isp, other)
    if plot:
        ax = axs[1]
        nx.draw_networkx(g,
                         ax=ax,
                         with_labels=False,
                         pos={node: (roiEvents.loc[node, "t0"]+.5*roiEvents.loc[node, "halfwidth"], roiEvents.loc[node, "its"]) for node in g.nodes},
                         node_size=0,
                         edge_color="grey"
                        )
        #ax.set_yticks(np.arange(len(timescales)));
        #ax.set_yticklabels(timescales);
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    ######################################################
    df_filt = []
    for ixs in nx.connected_components(g.to_undirected()):
        row = pd.Series(roiEvents.loc[list(ixs)].sort_values("height").iloc[-1])
        for col in ["t0","tend","height"]:
            row.col = roiEvents.loc[list(ixs),col].mean()
        row.halfwidth = row.tend-row.t0
        df_filt += [row]
    df_filt = pd.DataFrame(df_filt)
    for col in ["loghalfwidth","score","status"]:
        if col in df_filt.columns:
            del df_filt[col]
    if plot:
#         if len(df_filt):
#             df_filt["color"] = list(plt.cm.hot(coltrans(df_filt.height, vmax=10)))
        for ix, row in df_filt.iterrows():
            draw_spike(row, axs[2])
        nr = np.round(regions.Freq/freqShow)
        if nr>1:
            nr = int(nr)
            from .numeric import rebin
            x = rebin(regions.time, nr)
            y = rebin(origtrace, nr)
        else:
            x,y = regions.time, origtrace
#         y = y - np.percentile(y,1)
#         y = y / np.percentile(y,99)
#         y = y * roiEvents.its.max() + roiEvents.its.max()/3
        for ax in axs:
            ax.set_facecolor("k")
            axx = ax.twinx()
            axx.set_ylabel("light intensity")
            axx.plot(x, y, c="grey", lw=.7)
            if plotSlows is not None:
                for ts in plotSlows:
                    axx.plot(regions.showTime.get("%g"%ts,regions.time), tsslows[ts])
            ax.set_ylabel("filtering timescale [s]")
            ax.set_xlabel("time [s]]")
            ax.set_ylim(ax.get_ylim()[::-1])
        axs[0].set_title("discarded")
        axs[1].set_title("all")
        axs[2].set_title("filtered")
        if hasattr(regions, "timescales"):
            timescales=regions.timescales
            ax.set_yticks(np.arange(len(timescales)))
            ax.set_yticklabels(["%g"%ts for ts in timescales]);
        else:
            ax.set_yticks(roiEvents.its.unique());
            ax.set_yticklabels(["%g"%ts for ts in roiEvents.ts.unique()]);
        ax.set_xlim(roiEvents.t0.min()-1, roiEvents.tend.max()+1)
        plt.tight_layout()
        df_filt, axs
    return df_filt, None