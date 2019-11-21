import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import distributions as dst
from scipy.optimize import curve_fit#,minimize,basinhopping
from numba import jit

def decay(time,top,bottom,rate):
#     top,bottom,rate = p
    return bottom+(top-bottom)*np.exp(-time*rate)
from scipy.optimize import curve_fit

def guessDecayPars(y):
    b0 = np.nanmin(y)
    y  = np.log(y-b0)
    r0 = np.diff(y)
    r0 = r0[np.isfinite(r0)]
    r0 = -np.mean(r0)
    t0 = np.nanpercentile(y,99)
    p0 = (np.exp(t0)+b0,b0,r0)
    return p0

def decayfit(x,Ntrials=None):
    if Ntrials is None:
        lx = 10
    else:
        lx = Ntrials
    nx = len(x)//10
    TT = np.arange(len(x))
    tt = TT.copy()
    for j_ in range(lx):
        try:
            p0 = guessDecayPars(x)
            ff = np.isfinite(x)
            popt = curve_fit(decay,tt[ff],x[ff],p0=p0)[0]
            expDecay = decay(TT,*popt)
            return expDecay
        except:
            x = x[:-nx]
            tt = tt[:-nx]
    return p0

@jit 
def percFilter(x_,perc,filterSize):
    if filterSize%2==0:
        raise ValueError("filter size needs to be odd number")
    delta = filterSize//2
    out = np.zeros_like(x_)
    x_ = np.hstack([[x_[0]]*delta,x_,[x_[-1]]*delta])
#     print (x_.shape, out.shape)
    for i in range(len(out)):
        out[i] = np.percentile(x_[i:i+filterSize],perc)
    return out
@jit 
def runningAverage(x_,filterSize):
    if filterSize%2==0:
        raise ValueError("filter size needs to be odd number")
    delta = filterSize//2
    out = np.zeros_like(x_)
    x_ = np.hstack([[x_[0]]*delta,x_,[x_[-1]]*delta])
#     print (x_.shape, out.shape)
    for i in range(len(out)):
        out[i] = np.nanmean(x_[i:i+filterSize])
    return out
@jit 
def runningStd(x_,filterSize):
    if filterSize%2==0:
        raise ValueError("filter size needs to be odd number")
    delta = filterSize//2
    out = np.zeros_like(x_)
    x_ = np.hstack([[x_[0]]*delta,x_,[x_[-1]]*delta])
#     print (x_.shape, out.shape)
    for i in range(len(out)):
        out[i] = np.nanstd(x_[i:i+filterSize])
    return out

def lowPass(x, wIron, wAvg=None, perc=5,npass = 1):
    ironedx = percFilter(x, perc, wIron)
    if wAvg is not None:
        ironedx = runningAverage(ironedx,wAvg)
    for i in range(npass-1):
        ironedx += lowPass(x-ironedx,wIron,wAvg,perc,)
    return ironedx

def getInterestingROIcenters(image_, Nsample=10000, exponent=3, bandwidth=10, n_jobs=4):
    from sklearn.cluster import MeanShift
    p = image_.flatten()**exponent
    p = p/p.sum()
    yd,xd = np.array([(x//image_.shape[0], x%image_.shape[0]) for x in np.random.choice(np.arange(image_.size),p=p,size=Nsample)]).T
    meanshift = MeanShift(bandwidth,n_jobs=n_jobs)
    meanshift.fit(np.array([xd,yd]).T, )
    return meanshift.cluster_centers_.T