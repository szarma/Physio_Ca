from collections import OrderedDict
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from islets.numeric import climb
from islets.utils import multi_map


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


def showEventEdges(ks, evIndices, dims, resc=0):
    imm = putRois(ks,evIndices,dims)
    return getEdges(imm)

def getEdges(image,rescale=0):
    from scipy.ndimage import binary_fill_holes
    from skimage.feature import canny
    from skimage.morphology import remove_small_objects
    edges = canny(image)
    fill_coins = binary_fill_holes(edges)
    roiEdges = ndi.label(remove_small_objects(fill_coins, 21))[0]
    return mark_boundaries(image*rescale, roiEdges)


def getEvents(tfMatrix_,DistMatrix_):
    import networkx as nx
    G = nx.from_scipy_sparse_matrix(getEvents_(tfMatrix_,DistMatrix_))
    events = [cmp for cmp in nx.connected_components(G)
              if len(cmp)>1 or allGraph[tuple(cmp)*2]]
    return events

@jit
def runningMode(x_, filterSize, boundary="reflective"):
    if filterSize%2==0:
        raise ValueError("filter size needs to be odd number")
    delta = filterSize//2
    out = np.zeros_like(x_)
    allowedBoundaries = ["reflective","equal"]
    if boundary not in allowedBoundaries:
        raise ValueError(f"boundary {boundary} not recognized. Allowed values are: "+repr(allowedBoundaries))
    if boundary=="reflective":
        x_ = np.concatenate((x_[:delta][::-1], x_, x_[-delta:][::-1]))
    if boundary=="equal":
        x_ = np.concatenate(([x_[0]]*delta, x_, [x_[-1]]*delta))
    for i in range(len(out)):
        c,v = np.histogram(x_[i:i+filterSize], nbins)
        out[i] = v[np.argmax(c)]
    return out


def crawlDict_restr(image, pixels, diag=False, processes=10, verbose=False):
    global iterfcd
    # noinspection PyRedeclaration
    def iterfcd(ij):
        return climb((ij[0],ij[1]), image, diag=diag,)
    R_ = multi_map(iterfcd,pixels, processes=processes)
    A_ = [ij+r for ij,r in zip(pixels,R_) if r in pixels]
    A_ = [el for el in A_ if el[-1] is not None]
    B_ = OrderedDict()
    for (i0,j0,i1,j1) in A_:
        if (i1,j1) not in B_:
            B_[(i1,j1)] = []
        B_[(i1,j1)] += [(i0,j0)]
    return B_


def get_crawl_dict(image, pixels, diag=False, offset=(0,0)):
    # noinspection PyRedeclaration
    R_ = [climb((i,j), image, diag=diag,) for i,j in pixels]
    A_ = [ij+r for ij,r in zip(pixels,R_) if r in pixels]
    A_ = [el for el in A_ if el[-1] is not None]
    B_ = OrderedDict()
    for (i0,j0,i1,j1) in A_:
        i0 += offset[0]
        j0 += offset[1]
        i1 += offset[0]
        j1 += offset[1]
        if (i1,j1) not in B_:
            B_[(i1,j1)] = []
        B_[(i1,j1)] += [(i0,j0)]
    return B_


    def infer_gain(self, ts=None, plot=False, verbose=False, ret_points=False):
        if ts is None:
            minDt = np.diff(self.time).mean()
            freq = 1 / minDt
            ts = 30 / freq
        absSlow, absFast, _ = self.fast_filter_traces(ts, write=False, normalize=False, filt_cutoff=0)
        di = 30
        slow_est, fast_vars = [], []
        for i in range(absFast.shape[0]):
            for j in range(di, absFast.shape[1] - di, absFast.shape[1] // 30):
                slow_est += [absSlow[i, j]]
                fast_vars += [absFast[i, j - di:j + di].var()]
        fast_vars = np.array(fast_vars)
        slow_est = np.array(slow_est)
        logbs = np.log(np.logspace(np.log10(np.percentile(slow_est, 2)), np.log10(np.percentile(slow_est, 98))))
        d = np.digitize(np.log(slow_est), logbs)
        x = np.array([slow_est[d == i].mean() for i in np.unique(d)])
        y = np.array([np.median(fast_vars[d == i]) if (d == i).sum() > 10 else np.nan for i in np.unique(d)])
        x = x[np.isfinite(y)]
        y = y[np.isfinite(y)]
        if ret_points:
            return x, y
        gain = np.mean(y / x)
        slow_est[slow_est <= 0] = np.nan
        if plot:
            ax = plt.subplot(111)
            ax.hexbin(slow_est, fast_vars, bins="log",
                      xscale="log",
                      yscale="log",
                      cmap="Greys",
                      mincnt=1
                      )
            c = ax.plot(x, y, "o", mfc="none")[0].get_color()
            ax.plot(x, x * gain, c=c)

        if verbose: print("initial estimate of the gain is", gain)

        for _ in range(5):
            fast_vars[fast_vars > 5 * gain * slow_est] = np.nan
            if np.isnan(fast_vars).any():
                x = np.array([slow_est[d == i].mean() for i in np.unique(d)])
                y = np.array([np.nanmedian(fast_vars[d == i]) if (d == i).sum() > 10 else np.nan for i in np.unique(d)])
                x = x[np.isfinite(y)]
                y = y[np.isfinite(y)]
                y[y <= 0] = y[y > 0].min()
                gain = np.nanmean(y / x)
                if verbose: print("revised estimate of the gain", gain)
                if plot:
                    c = ax.plot(x, y, "o", mfc="none")[0].get_color()
                    ax.plot(x, x * gain, c=c)
        if plot:
            ax.set_title("gain inference")
            ax.set_xlabel("window means")
            ax.set_ylabel("window variances")
            ax.grid()
        self.gain = gain