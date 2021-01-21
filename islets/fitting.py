from collections import namedtuple
import numpy as np
from scipy.stats import distributions as dst
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd

def generic_function(fname,pars,t,**kwargs):
    return eval(fname)(*(t,)+pars,**kwargs)

def myShape(t,loc=0,amplitude=1,decay=1,halfwidth=1,offset=0,end=None,fillnan=True,bkgsum=True,spike_cut=0.01):
    if end is None:
        end=offset
    alpha = (halfwidth*decay/2)**2+1
    x0 = alpha-1-halfwidth*decay/2
    x1 = alpha-1+halfwidth*decay/2
    dloc = x0
    dloc /= decay
    yd = dst.gamma.pdf(t,alpha,scale=1./decay,loc=loc-dloc)
    ydmax = dst.gamma.pdf((alpha-1)/decay,alpha,scale=1./decay,)
    Offset = np.zeros_like(yd)
    good = np.abs(yd)>np.abs(ydmax)*spike_cut
    if good.any():
        ibegin = np.where(np.abs(yd)>np.abs(ydmax)*spike_cut)[0][0]
        ibegin = max(0,ibegin-2)
        iend = np.where(good)[0][-1]+1
        iend = min(iend, len(yd))
    else:
        ibegin, iend = 0, len(yd)
    Offset[:ibegin] = np.nan if fillnan else offset
    Offset[iend:] = np.nan if fillnan else end
    Offset[ibegin:iend] = offset + (end-offset)*np.arange(iend-ibegin)/(iend-ibegin)
    yd = amplitude*yd/ydmax
    ysum = yd + Offset
    offset, end = ysum[0], ysum[-1]
    Offset[ibegin:iend] = offset + (end-offset)*np.arange(iend-ibegin)/(iend-ibegin)
    yd = ysum-Offset
    
#     if np.all(np.diff(yd[imax:iend])<0): break
    try:
        imax = np.nanargmax(yd)+1
        iend = imax+np.where(np.diff(yd[imax:])>0)[0][0]
        yd[iend:] = np.nan if fillnan else 0
        ysum = Offset+yd
        ibegin = np.where(np.diff(yd[:imax-2])<0)[0][-1]
        yd[:ibegin] = np.nan if fillnan else 0
        yd[yd<0] = 0
        ysum = Offset+yd
    except:
        pass
    if bkgsum:
        return ysum
    else:
        return yd, Offset

def myLogNorm(t,loc=0,amplitude=1,scale=1,s=.2,offset=0,end=None, reshift=True,fillnan=True,bkgsum=True,spike_cut=0.01):
    if end is None:
        end=offset
#     loc = loc+(np.exp(-s**2)-np.e)*scale
    if reshift:
        inflection_points = scale*np.exp(-5/8*s**2+np.array([-1,1])*s/2*np.sqrt(57/16*s**2+2))
        loc -= inflection_points[0]
#     loc *= scale
    ds = dst.lognorm(loc=loc,scale=scale,s=s)
    yd = ds.pdf(t)
    ydmax = ds.pdf(loc+np.exp(-s**2)*scale)
    Offset = np.zeros_like(yd)
    good = np.abs(yd)>np.abs(ydmax)*spike_cut
    if good.any():
        ibegin = np.where(np.abs(yd)>np.abs(ydmax)*spike_cut)[0][0]
        ibegin = max(0,ibegin-2)
        iend = np.where(good)[0][-1]+1
        iend = min(iend, len(yd))
    else:
        ibegin, iend = 0, len(yd)
    Offset[:ibegin] = np.nan if fillnan else offset
    Offset[iend:] = np.nan if fillnan else end
    Offset[ibegin:iend] = offset + (end-offset)*np.arange(iend-ibegin)/(iend-ibegin)
    yd = amplitude*yd/ydmax
    ysum = yd + Offset
    if fillnan:
        offset, end = ysum[np.isfinite(ysum)][[0,-1]]
    else:
        offset, end = ysum[0], ysum[-1]
    Offset[ibegin:iend] = offset + (end-offset)*np.arange(iend-ibegin)/(iend-ibegin)
    try:
        imax = np.nanargmax(yd)+1
        iend = imax+np.where(np.diff(yd[imax:])>0)[0][0]
        yd[iend:] = np.nan if fillnan else 0
        ysum = Offset+yd
    except:
        pass
    if bkgsum:
        return ysum
    else:
        yd = ysum-Offset
        return yd, Offset

def getInitialPars(t_,x_,p_,fshape, offset=None, end=None, x_slow=None):
    if x_slow is None and (offset is None or end is None):
        pc10 = np.percentile(x_,10)
        if offset is None:
            offset = pc10
        if end is None:
            end = pc10
    if x_slow is not None:
        offset = x_slow[0]
        end = x_slow[-1]
    from .numeric import guessDecayPars
    imax_ = np.argmax(x_)
    if len(x_)>5 and imax_>0:
        amplitude = x_[imax_-1:imax_+1].mean()-offset
    else:
        amplitude = x_[imax_]-offset
    if fshape.__name__ == "myShape":
        r_ = guessDecayPars(x_)[-1]
        myParameters = namedtuple("myParameters", ['loc', 'amplitude', 'decay', 'halfwidth', 'offset', 'end'])
        pars0 = myParameters(decay      = (t_[1]-t_[0])/r_/10,
                             halfwidth = p_["halfwidth"]*.7,
                             amplitude  = amplitude,
                             loc        = p_["t0"] + p_["halfwidth"]*.1,
                             offset     = offset,
                             end        = end
                          )
    if fshape.__name__ == "myLogNorm":
        LNParameters = namedtuple("LNParameters", ['loc', 'amplitude', 'scale', 's', 'offset', 'end'])
        pars0 = LNParameters(s         =.2,
                             scale     = p_["halfwidth"]*1.5,
                             amplitude = amplitude,
                             loc       = p_["t0"] + p_["halfwidth"]*.4,
                             offset     = offset,
                             end        = end
                          )
    return pars0

def toMin(x, t_, y_, shape=myShape):
    yf_ = shape(t_,*x,fillnan=False)
    return ((yf_-y_)**2).sum()

def fit_spikes(t,y,spikeDf,
               y_slow=None,
               ifit=None,
               plot=False,
               colorCode = None,
               z_score_th=None,
               rel_amplitude_threshold = None,
               halfwidth_threshold = None,
               ax=None,
               verbose=False,
               allowjump=True,
              ):
    from scipy.signal import find_peaks, peak_widths
    from scipy.optimize import minimize, basinhopping
    if colorCode is None:
        colorCode = dict(zip(["myShape", "myLogNorm"],["C1","C0"]))
    if plot:
        if ax is None:
            fig, ax = plt.subplots(1,1,figsize=(9,4), sharex=True)
        ax.plot(t, y, lw=.6, color="grey")
        plt.tight_layout()
    if ifit is None:
        ifit = np.inf
    iSpike = 0
    parDf = []
    for ip,p in spikeDf.iterrows():
        i0  = np.searchsorted(t , p.t0-p["halfwidth"]*.7,)
        ie  = np.searchsorted(t , p.t0+p["halfwidth"]*3)
        ycur = y[i0:ie]
        tcur = t[i0:ie]
        if y_slow is None:
            y_slow_curr = None
        else:
            y_slow_curr = y_slow[i0:ie]
            if z_score_th is not None:
                dy = ycur-y_slow_curr
                excursions = ycur>z_score_th*y.std()
                if not any(excursions):
                    if verbose: print(f"Rejected @ ix={ip}: the neighbourhood is too variable")
                    continue
        imax = np.argmax(ycur[:len(ycur)//2+1])
        if plot:
            ax.plot( tcur[imax], ycur[imax], "ko", mfc="w", mew=.5)
#         res = pd.DataFrame(columns=["fshape","pars","score","background","halfwidth"])
        res = pd.DataFrame()

        for fname in colorCode:
            fshape = eval(fname)
            c = colorCode[fname]
            pars0 = getInitialPars(tcur,ycur,p,fshape,x_slow=y_slow_curr)
            if plot:
                y0 = fshape(tcur, *pars0, fillnan=False)
                ax.plot(tcur, y0, c, lw=.7, ls="--", label=fshape.__name__ if iSpike==0 else None)
            lowerB = {"offset":.8}
            parsdict = pars0._asdict()
            bounds = [(parsdict[j]*lowerB.get(j,0.05), parsdict[j]*10) for j in parsdict.keys()]
            bounds = [b if b[0]<b[1] else (b[1]/20, b[1]) for b in bounds]
            x0 = np.array(pars0)
            if not allowjump:
                x0 = x0[:-1]
                bounds = bounds[:-1]
            try:
                opt = minimize(toMin, x0, args=(tcur,ycur,fshape),
                               bounds = bounds
                              )
    #             opt = basinhopping(toMin, np.array(pars0), minimizer_kwargs=dict(
    #                 args=(tcur,ycur,fshape),
    #                 bounds = [(parsdict[j]*lowerB.get(j,0.05), parsdict[j]*10) for j in parsdict.keys()])
    #                               )
    #             except:
    #                 opt = minimize(toMin, np.array(pars0), args=(tcur,ycur,fshape),
    #                                # bounds = [(pars0[j]*.05, pars0[j]*10) for j in range(len(pars0))],
    #                                # method="Nelder-Mead"
    #                               )
                pars = pars0._replace(**dict(zip(parsdict.keys(), opt.x)))
            except ValueError:
                pars = pars0
            if not allowjump:
                pars = pars._replace(end=pars.offset)
            yf,bkg = fshape(tcur,*pars, bkgsum=False, fillnan=True)
            if np.isfinite(yf).sum()<3: continue
            ycur = ycur[np.isfinite(yf)]
            tcur = tcur[np.isfinite(yf)]
            bkg  = bkg [np.isfinite(yf)]
            yf   = yf  [np.isfinite(yf)]
#             yf_th = yf.max()/10
#             if bkg[0]>yf_th or bkg[-1]>yf_th:
#                 if verbose: print(f"Rejected @ ix={ip}: the spike is not prominent enough")
#                 continue
            if not np.isfinite(yf).any(): 
                if verbose: print(f"Rejected @ ix={ip}: all nans(!?)")
                continue
#             if (yf<0).any():
#                 if verbose: print(f"Rejected @ ix={ip}: as negative peak")
#                 continue
            pks = find_peaks(yf)[0]
            if len(pks)==0: 
                if verbose: print(f"Rejected @ ix={ip}: Could not infer the peak. yf={np.round(yf,3)}")
                continue
                hw = p['halfwidth']
                edges = None
            else:
                hw,_,ledge,redge = peak_widths(yf,pks)
                dt = t[1]-t[0]
                hw = hw[0]*dt
                redge = redge[0]*dt+tcur[0]
                ledge = ledge[0]*dt+tcur[0]
                edges = (ledge, redge)
            if np.abs(pars0.loc-pars.loc)>p["halfwidth"]+hw: 
                if verbose: print(f"Rejected @ ix={ip}: inferred peak too far from the candidate")
                continue
            yf += bkg
            res = res.append({
                "fshape": fshape.__name__,
                "pars": pars,
                "chi2_raw": np.mean((ycur-yf)**2/(yf+yf.min()/100)),
                "ycur": ycur,
                "yf": yf,
                "edges": edges,
                "time": tcur,
                "background": bkg,
                "halfwidth": hw,
            }, ignore_index=True)
            if plot:
                ax.plot(tcur, yf+bkg, lw=1.5, c=c, label=fshape.__name__ if iSpike==0 else None)
                ax.plot(tcur, bkg, lw=1.5, c=c)
                ax.plot(tcur, ycur-yf, lw=.5, c=c)
#             except:
#                 continue
        if len(res):
            res = res.sort_values("chi2_raw").iloc[0]
            if plot:
                ax.plot(tcur[imax],ycur[imax],"o",c=colorCode[res.fshape],ms=3.5)
            save=True
            if halfwidth_threshold is not None:
                if res.halfwidth<halfwidth_threshold:
                    save=False
            if save:
                res.name = ip
                parDf += [res]
            else:
                if plot:
                    ax.plot(tcur[imax],ycur[imax],"kx",ms=5)
            
    #     y[i0:ie][np.isfinite(yf)] -= yf[np.isfinite(yf)]
    #     plt.plot(tcur, y[i0:ie], "C1",lw=.5)

        if iSpike==0:
            tmin = t[i0]-10
        iSpike += 1
        if iSpike>=ifit: break
    if plot:
#         ax.set_xlim(tmin, tcur[-1]+10)
        ax.legend()
    parDF = pd.DataFrame(parDf)
    if len(parDF):
        parDF["t0"]       = [par.loc       for par in parDF.pars]
#         parDF["amplitude"] = [par.amplitude for par in parDF.pars]
    return parDF


def multiShape(t_,df_,bkgsum=False):
    out = np.zeros_like(t_)
    for _,row in df_.iterrows():
        f = eval(row.fshape)
        tmp = f(t_,*row.pars, fillnan=False, bkgsum=bkgsum)
        if not bkgsum:
            tmp = tmp[0]
        out += tmp
    return out


def guessPars(t,x):
    x0 = np.percentile(x,5)
    crit = x>x0+4
    if np.any(crit):
        iloc = np.where(crit)[0][0]
    else:
        iloc = np.argmax(x)
#     if iloc<0:iloc=0
#     x0 = np.mean(x[:iloc])
    loc = t[iloc]
    imax = iloc+np.argmax(x[iloc:min(len(x),iloc+10)])
    xmaxAvg = x[max(imax-1,0):min(len(x),imax+2)].mean()
    if xmaxAvg<np.percentile(x,10):
        return False
    Ampl = xmaxAvg-x0
    scale = .6
    s = 1
    return Ampl,x0,loc,scale,s


def myGamma(t,Ampl,loc,scale=.1,s=2,offset=0,end=None):
    if end is None:
        end=offset
    yd = dst.gamma.pdf(t,s,loc=loc,scale=scale)
#     ishoot = np.argmax(np.diff(yd))
#     ibegin = np.where(t<loc)[0][-1]
# #     loc -= (t[ishoot]-t[ibegin])/2
#     yd = dst.lognorm.pdf(t,loc=loc,scale=scale,s=s)
    Offset = np.zeros_like(yd)
    ibegin = np.where(t<loc)[0][-1]
    iend = np.where(yd>yd.max()/100)[0][-1]
#     print (ibegin, iend, offset, end)
    Offset[:ibegin] = np.nan
    Offset[iend:] = np.nan
    Offset[ibegin:iend] = offset + (end-offset)*np.arange(iend-ibegin)/(iend-ibegin)
    yd = Ampl*yd/yd.max()+Offset
    return yd

def multiFun(t,ps):
    if len(ps)==0:
        return np.zeros_like(t)
    else:
        return np.sum([myLogNorm(t,*p) for p in ps], axis=0)

def firstOccurConseq(a_,conseq=1):
    a_=a_.copy()
    icc = np.where(a_)[0]
    if conseq==1:
        return icc[0]
    else:
        for i in range(conseq-1):
            try:
                icc = icc[:-1][np.diff(icc)==1]
            except:
                return -1
        return icc[0]

def EvalModel(ps,t,x):
    xf = multiFun(t,ps)
    return np.sum((x-xf)**2)

def Fit(t,y_,ax=None,nPeaks=30,verbose=False, showFailed=True):
    yoffset = np.percentile(y_,5)
    y = y_-yoffset
    yscale = y.std()*2
    y = y/yscale
    if ax:
        ax.plot(t,y,lw=.5)
    iBegin = 0
    dt = np.diff(t).mean()
#     print (len(t),len(y))
    ps = []
    for ip in range(nPeaks):
        if iBegin>len(y)-10:
            break
        y0 = multiFun(t,ps)
#         if y0 = 0: y0 = y.mean()
        tr,yr = t[iBegin:],(y-y0)[iBegin:]
        p0 = guessPars(tr,yr)
        if not p0:
            if verbose: print ("could not guess, p0=",p0)
            break
        p0 = p0[:-1]
        crit = yr-myLogNorm(tr,*p0)>3
#         iStop = firstOccurConseq(crit,3)
        try:
            iStop = firstOccurConseq(crit,3)
            if iStop<0:
                iStop = 10
            if iStop>len(yr):
                iStop = len(yr)
        except: iStop = len(yr)
        iStop += iBegin
#         if iStop<0: 
#             iStop = min(len(yr),iBegin+1./dt)
#         else:
#         if iStop==iBegin:
#             iStop = iBegin+int(.3/dt)
#         iStop = min(len(y),iStop)
        if verbose: print (ip,iBegin,iStop,p0,t[iBegin],t[iStop-1])
        tr = tr[:iStop-iBegin]
        yr = yr[:iStop-iBegin]
        if ax: 
            ax.plot(tr,multiFun(tr,ps+[p0]),"lightgrey",lw=.8)
        try:
            t0 = p0[2]
            p = curve_fit(myLogNorm,tr,yr,p0=p0,
                          bounds=np.array([
                              (1,np.inf),
                              (-np.inf,np.inf),
                              (t0-.15,t0+.05),
                              (.05,.2),
                              (0.8,2),
                          ])[:len(p0)].T 
                         )[0]
        except:
            if verbose: print ("something was off with fitting")
            iBegin += int(.1/dt)
            continue
#         print (p)
        score    = EvalModel([ ],t[:iStop], y[:iStop]-multiFun(t[:iStop],ps)-y.mean()*int(len(ps)==0))
        newScore = EvalModel([p],t[:iStop], y[:iStop]-multiFun(t[:iStop],ps))
#         score    = EvalModel([ ],tr,yr-yr.mean())
#         newScore = EvalModel([p],tr,yr)

        if newScore>score*.99:
            if verbose: print ("refused since score0 = %.1f and suggested score = %.1f"%(score,newScore))
            if ax and showFailed: ax.plot(tr,multiFun(tr,ps+[p]),color="grey")
            iBegin += int(.1/dt)
            continue
        if verbose: print ("accepted since score0 = %.1f and suggested score = %.1f"%(score,newScore))
        ps += [p]
        score=newScore
        if ax: ax.plot(t,y0+myLogNorm(t,*p))
        iBegin = max(0,iStop - int(.2/dt))
#         iBegin += int(.2/dt)
    return ps

def manyFit(t,x,verbose=False,toll=.1):
    var0 = x.var()
    xf = x.copy()
    PS_ = []
    for i in range(10):
        ps = Fit(t,xf,nPeaks=3)
        if not len(ps):
            print ("stopped because there was nothing else to find")
            break
        xf = xf-multiFun(t,ps)
        var1 = xf.var()
        dr2 = 1-var1/var0
        if dr2<toll:
            print ("stopped because the fve difference was too small: %f --> %f (%f)"%(var0,var1,dr2))
            break
        PS_ += ps
        print ("accepted because the fve difference was large enough: %f --> %f (%f)"%(var0,var1,dr2))
        var0 = var1
#         print (PR_)
    return PS_