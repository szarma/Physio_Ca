import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .numeric import rebin


class LineScan:
    def __init__(self, data,
                 metadata=None,
                 freq=None,
                 verbose=False,
                 name=None,
                 detrend=True
                ):
        self.data = data.astype("float32")
        if metadata is None:
            if freq is None:
                raise ValueError("If you do not pass metadata, you need to at least provide frequency.")
            else:
                self.metadata = pd.Series({
                    "SizeX":data.shape[0],
                    "pxSize":data.shape[0],
                    "pxUnit":"px",
                    "Frequency":freq
                })
        else:
            assert metadata["SizeX"]== data.shape[0]
            self.metadata = metadata
        self.Freq = self.metadata.Frequency
        self.time = np.arange(0,data.shape[1])/self.Freq
        self.name = name
        if detrend:
            self.detrend()
        
    def plot(self,
        tmax=None,
        timeScales = [None,3],
        distances = [5,10,20,50,100],
        times = [.1,.5,1,2,5,10,30],
        save=False,
        Npoints = 1000,
        verbose=False,
            axs = None
        ):
        if tmax is None:
            tmax = self.time[-1]
        physicalSize = self.metadata[["SizeX","pxSize"]].prod()
        distance = distances[np.searchsorted(distances, physicalSize*.15)-1]
        Tind = times[np.searchsorted(times, tmax*.2)-1]
        txtOffset = physicalSize-1.2*distance
        nr = len(timeScales)
        if axs is None:
            fig, axs = plt.subplots(nr,1,
                                    figsize=(4,3*nr),
                                    sharey=True,
                                    sharex=True,
                                    squeeze = False
                                   )
            axs = axs.flatten()
            fig.patch.set_facecolor('k')
            plt.subplots_adjust(wspace=.02, hspace=.02)
        nRebin = int(np.round(len(self.time)/Npoints))
        if nRebin>1:
            if verbose: print ("rebinning by", nRebin)
            newdata = rebin(self.data, nRebin, axis=1)
            ls = LineScan(newdata, freq=self.Freq/nRebin)
        else:
            ls = self
        for it in range(len(timeScales)):
            ts = timeScales[it]
            if ts is None:
                txt = "raw"
                x = ls.data
            else:
                txt = "detrended"
                if ls.time[-1]>ts:
                    s,x,z = ls.fast_filter_traces(ts,write=False)
                else:
                    if not hasattr(ls, "detrended"):
                        ls.detrend()
                    x = ls.detrended
            ax = axs[it]
            ax.set_yticks([])
            ax.set_xticks([])
            ax.imshow(x,extent=(0, x.shape[1]/ls.Freq,0,physicalSize,),
                      aspect="auto",
                      cmap="hot",
                      vmin=0
                     )
            ax.set_ylabel(txt,color="lightgrey")
            if it==len(timeScales)-1:
                ax.plot([tmax*.03]*2,[txtOffset,txtOffset+distance],"lightgrey")
                ax.text(tmax*.03,txtOffset+distance/2,"  "+str(distance)+"µm", va="center", color="w")
                ax.plot([tmax*.03, tmax*.03+Tind],[.03*physicalSize]*2,"lightgrey")
                ax.text(tmax*.03+Tind/2,.03*physicalSize,str(Tind)+"s\n", ha="center", va="center", color="w")
            if it==0 and self.name is not None:
                ax.set_title(self.name,color='w')
        ax.set_xlim(0,tmax)
        if save:
            fig.tight_layout()
            fig.savefig(save, dpi=150, facecolor='k')
            plt.close(fig)

    def detrend(self, fast=True, n=None, processes=1, func=None, points=None):
        from .numeric import fit_bleaching
        from .utils import multi_map
        if n is None:
            n = self.data.shape[0]
        self.detrended = np.zeros_like(self.data)
        if fast:
            trend = np.mean(self.data[:n],1)
        else:
            global iterf
            if points is None:
                points = 100
            def iterf(xi):
                return fit_bleaching(xi, func=func, points=points)
            trend = multi_map(iterf, self.data[:n], processes = processes)
            if processes==1:
                trend = list(trend)
        self.trend = np.array(trend)
        for i in range(n):
            self.detrended[i] = self.data[i] - trend[i]
    
    def examine(self, debug=False):
        from .linescanner import examine
        return examine(self, debug=debug)
    
    def fast_filter_traces(self,
                           ironTimeScale,
                           z_sp=0,
                           order=5,
                           Npoints=None,
                           write=True,
                           verbose=False,
                           dilation=False
                          ):
        from .numeric import fast_filter
        if hasattr(self,"gain"):
            gain = self.gain
        else:
            if z_sp==0:
                gain = 1
            else:
                raise ValueError("z_sp allowed to be non zero only if gain is provided. If the data represent raw photon counts, put gain to 1.")
            
        x_slow, x_fast, zscore = fast_filter(
            self.data,
            ironTimeScale= ironTimeScale,
            gain=gain,
            freq=self.Freq,
            z_sp=z_sp,
            order=5,
            verbose=verbose,
            npass=2,
            nreflect="auto",
            dilation=dilation
              )
        if write:
            setattr(self, "slower_%g"%ironTimeScale, x_slow)
            setattr(self, "faster_%g"%ironTimeScale, x_fast)
            setattr(self, "zScore_%g"%ironTimeScale, zscore)
        else:
            return x_slow, x_fast, zscore
        
    def infer_gain(self, plot=False, verbose=False):
        ts = 100/self.Freq
        absSlow, absFast, _ = self.fast_filter_traces(ts, write=False, z_sp=0)
        di = 100
        slow_est, fast_vars = [],[]
        for i in range(absFast.shape[0]):
            for j in range(di, absFast.shape[1]-di, absFast.shape[1]//30):
                slow_est  += [absSlow[i,j]]
                fast_vars += [absFast[i,j-di:j+di].var()]
        fast_vars = np.array(fast_vars)
        slow_est = np.array(slow_est)
        slow_est[slow_est<=0] = slow_est[slow_est>0].min()
        logbs = np.log(np.logspace(np.log10(np.percentile(slow_est,2)),np.log10(min(np.percentile(slow_est,98), self.data.max()*.4))))
        d = np.digitize(np.log(slow_est), logbs)
        x = np.array([slow_est[d==i].mean() for i in np.unique(d)])
        y = np.array([np.median(fast_vars[d==i]) for i in np.unique(d)])
        gain = np.mean(y/x)
#         gain = np.exp(np.mean(np.log(y)-np.log(x)))
        
        if plot:
            ax = plt.subplot(111)
            ax.hexbin(slow_est, fast_vars, bins="log",
                      xscale="log",
                      yscale="log",
                      cmap="Greys",
                      mincnt=1
                     )
            c = ax.plot(x,y,"o",mfc="none")[0].get_color()
            ax.plot(x,x*gain,c=c)
            
        if verbose: print ("initial estimate of the gain is", gain)
        
        for _ in range(5):    
            fast_vars[fast_vars>10*gain*slow_est] = np.nan
            if np.isnan(fast_vars).any():
                y = np.array([np.nanmedian(fast_vars[d==i]) for i in np.unique(d)])
                newgain = np.nanmean(y/x)
                if verbose: print ("revised estimate of the gain", gain)
                if plot:
                    c = ax.plot(x,y,"o",mfc="none")[0].get_color()
                    ax.plot(x,x*gain,c=c)
                if np.abs(newgain/gain-1)<.03:
                    if verbose:
                        print("gain converged.")
                    break
                gain = newgain
        
        self.gain = gain

    