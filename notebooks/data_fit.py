import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize_scalar
from scipy.stats import distributions as dst
import logging
from matplotlib import ticker

def inhibition(c, top=1, ic50=1, h=1, bottom=0):
    return bottom + (top-bottom)/(1+(c/ic50)**h)

def stimulation(c, top=1, ec50=1, h=1, bottom=0):
    return bottom + (top-bottom)*(c/ec50)**h/(1+(c/ec50)**h)


class Fit():
    def __init__(self, std_tidy, function=None, bounds=None, bootstrap=100):
        
        self.data = std_tidy
        self.minNonZeroConc = self.data.query("concentration>0")['concentration'].min()
        
        if function is None:
            function = stimulation
        self.function = function
        
        if bounds is None:
            bounds = [
                np.zeros(len(self.function.__defaults__)),
                [
                    self.data['measurement'].median()*1000,
                    self.data['concentration'].max()*100,
                    10,
                    self.data['measurement'].median()*10,
                ]
            ]
        
        pars, pars_cov = curve_fit(
            self.function,
            self.data['concentration'].values,
            self.data['measurement'],
            bounds = bounds
        )
        self.pars = pars
        self.pars_cov = pars_cov
        self.bootstrap = bootstrap
        
        if self.bootstrap:
            multi_pars = []
            for n in range(bootstrap):
                np.random.seed(n)
                idxs = self.data.sample(len(self.data), replace=True).index
                multi_pars += [
                    curve_fit(
                    self.function,
                    self.data.loc[idxs, 'concentration'].values,
                    self.data.loc[idxs, 'measurement'].values,
                    bounds = bounds
                    )[0]
                ]
            self.pars_bootstrapped = multi_pars
        
    def infer_value(self, y, pars=None, debug=False):
        if pars is None:
            pars = self.pars
        def toMin(x, y):
            return np.abs(y-self.function(x, *pars))
        ret = minimize_scalar(toMin, args=(y,), bounds=(0,self.data["concentration"].max()*10), method="bounded")
        if debug:
            print(ret)
        if ret.success and ret.x>self.minNonZeroConc/2:
            return ret.x
        else:
            return np.nan
        
    def infer_values(self, ys, pars=None):
        return np.array([self.infer_value(y, pars) for y in ys])
    
    def infer_confidence_interval(self, y, alpha=0.05):
        assert self.bootstrap>0
        val = self.infer_value(y,self.pars)
        if np.isnan(val):
            return np.nan, np.nan
        vals = [self.infer_value(y, pars) for pars in self.pars_bootstrapped]
        Nnans = np.isnan(vals).sum()
        if Nnans:
            logging.warning(
                "Concentration cannot be inferred for some of the values."
                "The inferred confidence interval shold be considered a lower bound."
            )
        if self.bootstrap-Nnans<30:
            mean =  np.nanmean(vals)
            std = np.nanstd(vals)
            return dst.norm(mean, std).ppf(alpha/2), dst.norm(mean, std).ppf(1-alpha/2)
        else:
            return tuple(np.nanpercentile(vals, [100*(alpha/2), 100*(1-alpha/2)] ))
    
    def infer_confidence_intervals(self, ys, alpha=0.05):
        vals = np.array([self.infer_confidence_interval(y, alpha=alpha) for y in ys])
        return vals.T
        
    def plot(self, x_conc=None, infer_data = None, alpha=0.05, show_confidence_intervals=False, show_relative_error=False):
        minNonZeroConc = self.minNonZeroConc
        maxConc = self.data['concentration'].max()
        if infer_data is not None:
            nans = infer_data["inferred_concentration"].isna()
            if sum(~nans):
                minNonZeroConc = min(infer_data[~nans]['inferred_concentration'].min(), minNonZeroConc)
                maxConc = max(maxConc, infer_data[~nans]["inferred_concentration"].max())
        if x_conc is None:
            x_conc = np.array([0]+list(np.geomspace(minNonZeroConc*.7, maxConc*2)))
        fig, axs = plt.subplots(1,2,figsize=(7,4), gridspec_kw={"width_ratios":[1,15]}, dpi = 200)
        for ax in axs:
            ax.plot(x_conc, self.function(x_conc, *self.pars), "k--",
                    label='fit' if ax is axs[0] else None)
        ylim = ax.get_ylim()[1]
        if self.bootstrap:
            yl = np.geomspace(self.pars[-1], ylim, 100)
            l,u = self.infer_confidence_intervals(yl, alpha = alpha)
            # l =
            # u =
            ax.fill_betweenx(yl,
                             np.clip(u, x_conc[1], x_conc[-1]),
                             np.clip(l, x_conc[1], x_conc[-1]),
                             color='grey',
                             alpha = .3,
                             zorder=0,
                             )
            ax.plot(u, yl, color='k', lw=.5)
            ax.plot(l, yl, color='k', lw=.5)
            if show_relative_error:
                axx = ax.twinx()
                v = self.infer_values(yl)
                axx.plot(v, (u-l)/v/2*100, ls="dotted", color='grey')
                axx.set_ylabel("relative error of\ninferred concentration [%]")
                axx.set_yscale("symlog")
                axx.spines["left"].set_visible(False)
                axx.yaxis.set_major_formatter(ticker.ScalarFormatter())
                y_minor = ticker.LogLocator(base = 10.0, subs = np.arange(.1,1,.1), numticks = 10)
                axx.yaxis.set_minor_locator(y_minor)

        for rep, df in self.data.groupby("replicate"):
            axs[1].plot(df['concentration'], df["measurement"], "o", mfc="none", label=rep, mew=2, zorder=-2)
            axs[0].plot(df['concentration'], df["measurement"], "o", mfc="none", mew=2, zorder=-2)
        if show_relative_error:
            axs[1].legend(
                title = "std curve\nreplicate",
                loc = "upper center",
                frameon=False
            )
        else:
            axs[1].legend(
                title = "std curve\nreplicate",
                loc = (1.05, 0)
            )

        
        plt.setp(axs[0],zorder=10)
        axs[0].set_xlim(-minNonZeroConc/10,minNonZeroConc/10)
        axs[0].set_xticks([0])
        axs[0].set_xticklabels([0])
        if infer_data is not None:
            if nans.sum():
                axs[0].plot(
                    [-minNonZeroConc/10]*nans.sum(),
                    infer_data[nans]["measurement"],
                    "+",
                    color='darkred',
                    label='values outside of\nregion of validity'
                )
            if (~nans).sum():
                for ax in axs:
                    ax.plot(
                        infer_data[~nans]["inferred_concentration"],
                        infer_data[~nans]["measurement"],
                        "*",
                        color='darkred',
                        label='successfully inferred\ndata points',
                        mfc="w",
                        ms = 10
                    )
                if show_confidence_intervals and 'inferred_lowerCI' in infer_data.columns and 'inferred_upperCI' in infer_data.columns:
                    ax.hlines(
                        infer_data['measurement'],
                        infer_data['inferred_lowerCI'],
                        infer_data['inferred_upperCI'],
                        color='darkred'
                    )
                    ax.plot(
                        list(infer_data['inferred_lowerCI']) + list(infer_data['inferred_upperCI']),
                        list(infer_data['measurement'])*2,
                        "|",
                        color='darkred'
                    )
                    
            axs[0].set_xlim(-2*minNonZeroConc/10,2*minNonZeroConc/10)
        

        axs[0].text(1,0,"/", ha="center", va="center", transform = axs[0].transAxes)
        axs[0].text(1,1,"/", ha="center", va="center", transform = axs[0].transAxes)
        axs[1].text(0,0,"/", ha="center", va="center", transform = axs[1].transAxes)
        axs[1].text(0,1,"/", ha="center", va="center", transform = axs[1].transAxes)
        axs[1].set_xscale("log")
        axs[1].set_xlabel("concentration")
        axs[0].set_ylabel("measurement")
        plt.subplots_adjust(wspace=.05, right = .85)
        axs[0].spines['right'].set_visible(False)
        axs[1].spines['left'].set_visible(False)
        axs[1].set_yticks([])
        axs[0].legend( loc=2, frameon=False)
        for ax in axs:
            ax.set_ylim(0,ylim)
        xt = axs[1].get_xticks()
        axs[1].set_xticks(xt)
        axs[1].set_xticklabels(["%g"%t for t in xt])
        axs[1].set_xlim(x_conc[1], x_conc[-1]*.95)
        return fig