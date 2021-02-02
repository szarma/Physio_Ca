#!/usr/bin/env python
# coding: utf-8

# In[1]:


# %load_ext autoreload
# %autoreload 2


import numpy as np
import matplotlib.pyplot as plt
import islets
from scipy.stats import distributions as dst

get_ipython().run_line_magic('config', "InlineBackend.figure_format='retina'")

plt.rcParams["font.size"] = 16

from scipy.spatial import distance_matrix

import pandas as pd

k = 11
th = k*.45

res = 64
np.random.seed(0)
points = (np.random.rand(100,2)*res).astype(int)

while True:
    dd = distance_matrix(points, points)
    indices = np.array(np.triu_indices_from(dd,1)).T
    closest_neighbors = indices[dd[tuple(indices.T)]<th]
    if len(closest_neighbors)==0: break
    todel = pd.Series(closest_neighbors.flatten()).value_counts().index[0]
    points = points[[j!=todel for j in range(len(points))]]

image = np.zeros((res,res))
image[points.T[0], points.T[1]] = 500
m = islets.cmovie(image.reshape(1,res,res))
image = m.gaussian_blur_2D(k,k,k/6,k/6)[0].copy()
image = np.maximum(image, np.percentile(image[image>0],1))
real  = dst.poisson(mu=image+1+m.gaussian_blur_2D(k*2+1,k*2+1,-1,-1)[0]).rvs().astype(float)
vmax = islets.numeric.robust_max(real)//10*10

fig, axs = plt.subplots(1,3,figsize=(6.3,3), gridspec_kw={"width_ratios":[1,1,.06]})

for ax in axs[:-1]:
    ax.set_xticks([])
    ax.set_yticks([])
ax = axs[0]
ax.plot(points.T[1],points.T[0],"wo", mfc="none", ms=15, alpha = .7)
ax.text(0.01,.03,"ideal", transform=ax.transAxes, color="w")
im = ax.imshow(image, cmap="hot", vmax=vmax, vmin=0)
cticks = range(0,41,5)
cbar = plt.colorbar(im, cax=axs[-1], 
#                     ticks=cticks
                   )
ax = axs[1]
ax.text(0.01,.03,"real", transform=ax.transAxes, color="w")
ax.imshow(real, cmap="hot",  vmax=vmax, vmin=0)
cbar.ax.set_yticklabels([str(i) for i in cticks])  # vertically oriented colorbar
plt.tight_layout()
plt.savefig("figures/regions_ideal_real.pdf")


# In[2]:


k0 = int(k/2)+1
k0


# In[3]:


ks = [int(k0*.6//2*2+1),int(k0*1.7/2)*2+1,k0//2*2+1]
ks


# In[4]:


blur = {}
bkg  = {}
regions = {}
for ir,kk in enumerate(ks):
    blur[kk] = islets.cmovie(real.copy().reshape(1,res,res)).gaussian_blur_2D(kk,kk,-1,-1)[0]
    bkg[kk]  = islets.cmovie(real.copy().reshape(1,res,res)).gaussian_blur_2D(kk*2+1,kk*2+1,-1,-1)[0]
    dimage = blur[kk]-bkg[kk]
    dimage = dimage/islets.numeric.robust_max(dimage)
    regions[kk] = islets.Regions(dimage,mode="custom", verbose=True, img_th=.1)


# In[5]:


fig, axs = plt.subplots(3,4,figsize=(9,7.8), gridspec_kw={"width_ratios":[1,1,1,.06]})
for ir,kk in enumerate(blur):
    axr = axs[ir]
    axr[0].text(1.05,.5,r"$-$",transform=axr[0].transAxes, va="center",)
#     axr[0].set_ylabel(r"kernel = %i"%kk)
    txt = dict(zip(ks, ["too small","too large","adequate"]))[kk]
    axr[0].set_ylabel(txt)
    axr[1].text(1.05,.5,r"$=$",transform=axr[1].transAxes, va="center",)
    axr[0].imshow(blur[kk], cmap="hot", vmin=0, vmax=vmax)
    axr[1].imshow( bkg[kk], cmap="hot", vmin=0, vmax=vmax)
    dimage = blur[kk]-bkg[kk]
    dimage = dimage/islets.numeric.robust_max(dimage)
    im = axr[2].imshow(dimage, cmap="bwr",  vmax = 2,  vmin =-2)
    if ir==0:
        cbar = plt.colorbar(im, cax=axr[-1], ticks=[-1.5,0,1.5])
        cbar.ax.set_yticklabels(["-","bkg","+"])  # vertically oriented colorbar
    else:
        axr[-1].remove()
    for ax in axr[:-1]:
        ax.set_xticks([])
        ax.set_yticks([])
    
    regions[kk].plotEdges(ax=axr[2], image=False, color="k")
for ax,txt in zip(axs[0,:-1],["blurred","inferred\nbackground", "difference"]):
    ax.text(.5,1.15,txt, transform=ax.transAxes, va="center",ha="center")
axs[0,0].text(1.05,1.15,r"$-$",transform=axs[0,0].transAxes, va="center",)
axs[0,1].text(1.05,1.15,r"$=$",transform=axs[0,1].transAxes, va="center",)
# for ax in axs[:,2]:
#     ax.plot(points.T[1],points.T[0],"g+",ms=4)
ax = axs[1,0]
ax.text(-.2,.5,"kernel size",transform=ax.transAxes,va="center",ha="center", rotation=90)
plt.tight_layout()
# fig.set_facecolor("y")
plt.savefig("figures/regions_kernels.pdf")


# In[ ]:




