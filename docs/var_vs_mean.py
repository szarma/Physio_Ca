import islets
import numpy as np
import matplotlib.pyplot as plt
from copy import copy, deepcopy
from matplotlib.ticker import MultipleLocator, LogLocator
from matplotlib.colors import LogNorm
from matplotlib import patches

def get_means_and_variances(regions, ts, Nt = 100):
    s,f,z = regions.fast_filter_traces(ts,
                                       verbose=True,
                                       dilate=True,
                                       write=False,
                                       normalize=False,
                                       npass=0
                                      )
    dn = int(np.ceil(regions.Freq*ts))
    m,v = np.zeros((s.shape[0],Nt)),np.zeros((s.shape[0],Nt))
    for i in range(s.shape[0]):
        its = np.random.randint(dn,s.shape[1]-dn,Nt)
        for jt,it in enumerate(its):
            m[i,jt] = s[i,it]
            v[i,jt] = f[i,it-dn:it+dn].var()
    return m,v


plt.rcParams["font.size"] = 8
cmap = copy(plt.cm.Greys)
cmap.set_bad("lime")

regions = {
    "photon counting mode": islets.load_regions(
        "/home/jupyter-srdjan/notebooks/receptors/bull_rois.pkl",
        baremin=True, calcInterest=False),
}
regions["standard mode"] = islets.load_regions(
#     "/data/Sandra/2019/2019_07_01/Experiment22.lif_analysis/Series025/6_rois.pkl",
    "/data/Sandra/2019/2019_10_15/Experiment38a.lif_analysis/Series018/7_rois.pkl", # the pear 18
#     "/data/Sandra/2019/2019_08_09/Experiment27.lif_analysis/Series025/8_rois.pkl", # the burn
#         "/data/Sandra/2019/2019_10_15/Experiment38a.lif_analysis/Series020/2021_01_11_7_rois.pkl", # the pear
#         "/data/Sandra/2019/2019_10_15/Experiment38a.lif_analysis/Series020/7_rois.pkl", # the pear
        baremin=True, calcInterest=False)



r=regions['standard mode']
r.merge_closest(
    mergeSizeTh=30,
    # plot=True,
    mergeDist=2
)
r.reassign_peaks(r.statImages["mean"])
indices = r.df.query("peakValue>15").index
# r.plotEdges(indices)
r.df = r.df.loc[indices]
r.update()
r.infer_gain()



Npoints = 100
hexbincmap = "YlGnBu_r"
hexvmax = Npoints*1
bbox = {"alpha":.8,"facecolor":"w","edgecolor":"none","pad":0}
figsize = np.array((7,4.))*1.2
fig, axs = plt.subplots(2,4,figsize=figsize,  gridspec_kw=dict(width_ratios=[1,.3,1,1]))
for ax in axs[:,1]: ax.remove()
axs = axs[:,[0,2,3]]
for ax in axs.flat:
    ax.set_aspect("equal")
alpha = .15
ar = figsize[1]/figsize[0]
plt.subplots_adjust(top=1-alpha*ar, bottom=alpha*ar, right=1-alpha, left=alpha,wspace=.1, hspace=.0)

for k,ax in zip(regions, [axs[0,0], axs[1,0]]):
    r = regions[k]
    im = r.statImages["mean"].copy()
    hp = r.statImages["highperc"].copy()
    vmin = np.percentile(im,10)
    im[hp==hp.max()] = np.nan
    r.df["Nsatur"] = [sum([np.isnan(im[px]) for px in pxs]) for pxs in r.df.pixels]
    ax.imshow(im, cmap=cmap, norm=LogNorm(vmin=vmin))
    ax.set_xticks([])
    ax.set_yticks([])
    r.plotEdges(ax=ax,image=False, lw=.5, scaleFontSize=0, color="navy")
    ax.text(.01,-.01,"%.0fmin @ %.fHz"%(r.time[-1]/60,r.Freq),va="top", transform=ax.transAxes)
    ax.set_ylabel(k.replace(" mode","\n"))
r.df["Nsatur"] = [sum([np.isnan(im[px]) for px in pxs]) for pxs in r.df.pixels]
r.plotEdges(r.df.query("Nsatur>1").index, lw=.8, color="r",ax=ax, image=False, scaleFontSize=0)

axs[-1,-1].get_shared_x_axes().join(*axs[:,1:].flat)
axs[-1,-1].get_shared_y_axes().join(*axs[:,1:].flat)
axs[-1,-1].set_ylim(6,3e6)
axs[-1,-1].set_xlim(6,3e6)
for ik,k in enumerate(regions):
    r = regions[k]
    for ts,ax in zip([1,10],axs[ik,1:3]):
        m,v = get_means_and_variances(r,ts,Npoints)
        pc = ax.hexbin(m,v,xscale="log", yscale="log", mincnt=1, cmap=hexbincmap, vmax=hexvmax,zorder=3)
        xl = ax.get_xlim()
        ax.plot(xl,xl,c="C1",lw=.7,zorder=5)
        if hasattr(r,"gain"):
            xl = m.min()/2, m.max()*2
            ax.plot(xl,np.array(xl)*r.gain,c="C1",ls="--",lw=1,zorder=5)
        if ts>1:
            ff = r.df.Nsatur<30
            mrest, vrest = m[~ff], v[~ff]
            m,v = m[ff], v[ff]
            ax.plot(mrest.flat,vrest.flat,".", color="red", ms=1,zorder=3)
        pc = ax.hexbin(m,v,xscale="log", yscale="log", mincnt=1, cmap=hexbincmap, vmax=hexvmax,zorder=4)
        filtPars = r.prep_filtering(ts,write=False)
        if filtPars["Nrebin"]>1:
            txt = r"$N_{\rm rebin} = %i$"%filtPars["Nrebin"]
            ax.text(.95,.03, txt, ha="right", transform=ax.transAxes,zorder=4,bbox=bbox)
        
axs[0,2].set_xticklabels([])
axs[0,2].set_yticklabels([])
axs[1,2].set_yticklabels([])
axs[0,1].set_xticklabels([])
fig.set_facecolor("y")
axs[0,1].set_ylabel(r"variance $\sigma^2$")
axs[1,1].set_ylabel(r"variance $\sigma^2$")
axs[1,1].set_xlabel(r"mean @ 1s")
axs[1,2].set_xlabel(r"mean @ 10s")

for ax in axs[:,1:].flat:
    ax.xaxis.set_major_locator(LogLocator(numticks=10))
    ax.yaxis.set_major_locator(LogLocator(numticks=10))
    ax.minorticks_off()
    ax.grid(lw=.5)

ax = axs[0,-1]
ellipse = patches.Ellipse([.45,.7], .6, .25,angle=45,transform=ax.transAxes, fill=False,edgecolor="darkred",zorder=10)
ax.add_patch(ellipse)
ax.text(.1,.8,"signal",
        transform=ax.transAxes,
        fontdict=dict(color="darkred"),
        zorder=4,
        bbox=bbox
       )


ax = axs[1,2]
ellipse = patches.Ellipse([.77,.75], .25, .15, angle=45, transform=ax.transAxes, fill=False, edgecolor="darkred",zorder=20)
ax.add_patch(ellipse)
ax.text(.95,.55,"saturated",
        transform=ax.transAxes,
        fontdict=dict(color="darkred"),
        bbox={"alpha":.8,"facecolor":"w","edgecolor":"none","pad":0},
        zorder=10,
        ha='right',
       )
txt = f"Npoints={Npoints}"
for k in regions:
    txt += f", Nrois({k.split()[0]})=%i"%regions[k].df.shape[0]
fig.text(0,1,txt,va="top", fontsize=6)
fig.savefig("figures/mean_vs_var_tmp.pdf")
fig.savefig("figures/mean_vs_var_tmp.png",dpi=300)

