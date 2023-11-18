
import numpy as np

from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as pe
import cmasher as cmr
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.offsetbox import AnchoredText
from latex_form import latex_float
from termcolor import colored
import matplotlib.gridspec as gridspec
interval = np.linspace(0.12,0.9)
colors = cmr.pride(interval)
cmap = LinearSegmentedColormap.from_list('name', colors)


def line_background(lw,col):
    return [pe.Stroke(linewidth=lw, foreground=col), pe.Normal()]


###########################################################
Rate_e1 = np.loadtxt('/home/parracho/Desktop/Project/Ultralight_Dark_Matter/Real_pseudo/Rates/Real_pseudo/Rates(mh2<mH)3.dat',unpack=True,skiprows=1)
Rate_e2 = np.loadtxt('/home/parracho/Desktop/Project/Ultralight_Dark_Matter/Real_pseudo/Rates/Real_pseudo/Rates(mh2<mH)1.dat',unpack=True,skiprows=1)
Rate_e3 = np.loadtxt('/home/parracho/Desktop/Project/Ultralight_Dark_Matter/Real_pseudo/Rates/Real_pseudo/Rates(mh2<mH)2.dat',unpack=True,skiprows=1)


fig, axes= plt.subplots(nrows=1,ncols=1, figsize=(25, 15), dpi=100)
fig.set_figheight(10)
fig.set_figwidth(15)
plt.subplots_adjust(left=0.3,right=0.7,bottom=0.2,top=0.8,hspace = 1)
plt.style.use('classic')
plt.rcParams["font.family"] = "lmodern"
# plt.rcParams["font.weight"] = "ultralight"
plt.rcParams["font.size"] = "9"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rc('text', usetex=True)
plt.rcParams['axes.linewidth'] = 2.5
mpl.rcParams['hatch.linewidth'] = 3  # previous pdf hatch linewidth
mpl.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
ax = plt.gca()

axes.loglog(1 / Rate_e1[0, :], Rate_e1[1, :], lw=4, label='\\boldmath$\\alpha_1={}$'.format(latex_float(format(Rate_e1[6, 1],".1e"))),rasterized=True,color="b")
axes.loglog(1 / Rate_e2[0, :], Rate_e2[1, :], lw=4, label='\\boldmath$\\alpha_2={}$'.format(latex_float(format(Rate_e2[6, 1],".1e"))),rasterized=True,color="r")
axes.loglog(1 / Rate_e3[0, :], Rate_e3[1, :], lw=4, label='\\boldmath$\\alpha_3={}$'.format(latex_float(format(Rate_e3[6, 1],".1e"))),rasterized=True,color="g")

axes.axhline(y=1, color='k', ls='dashed', lw=4)
axes.axvline(x=2.7*2*(Rate_e1[9,1])**-1, color='k', alpha=0.4, ls='dashed', lw=4)
axes.axvline(x=2.7*2*(Rate_e1[8,1])**-1, color='k', alpha=0.4, ls='dashed', lw=4)



axes.text(3*1e-4, 2, r'$\boldsymbol{\Gamma/H = 1}$', color='k', fontsize=25, zorder=11, rasterized=True)
axes.text(0.02, 1e-29, r'$M_h={}$ GeV'.format(Rate_e1[9, 1]), color='k', fontsize=25, zorder=11, rasterized=True,rotation=90)
axes.text(8.4, 1e-33, r'$M_{{h_2}}={:.2f}$ GeV'.format(Rate_e1[8, 1]), color='k', fontsize=25, zorder=11, rasterized=True,rotation=90)


axes.set_ylim(1.0e-37, 1.0e14)
axes.set_xlim(1.0e-4, 1.0e3)

anchored_text_1 = AnchoredText(r'$\mathbf{M_\theta ={%s}}$ GeV' % latex_float(Rate_e1[7, 1]), loc=3,prop=dict(size=25),frameon=True)



axes.add_artist(anchored_text_1)
#axes.add_artist(anchored_text_2)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
axes.set_xlabel('$T^{-1}$ [GeV]$^{-1}$',fontsize=25)
axes.set_ylabel('$\\Gamma/H$',fontsize=25)
axes.grid(ls='-', alpha=0.1)
axes.legend(loc=2)

axes.legend(loc="upper right",prop={"size":20})

plt.figure(1)
plt.tight_layout()
plt.show()

fig.savefig("ratedm_lesser.pdf")
#############################################################
