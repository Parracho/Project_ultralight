
import numpy as np

from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as pe
import cmasher as cmr
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.offsetbox import AnchoredText
from latex_form import latex_float

interval = np.linspace(0.12,0.9)
colors = cmr.pride(interval)
cmap = LinearSegmentedColormap.from_list('name', colors)


def line_background(lw,col):
    return [pe.Stroke(linewidth=lw, foreground=col), pe.Normal()]


###########################################################
Rate_e1 = np.loadtxt('/home/parracho/Desktop/Project/Ultralight_Dark_Matter/Real_pseudo/Rates/New_higgs/Rates(mh2>mH)1.dat',unpack=True,skiprows=1)
Rate_e2 = np.loadtxt('/home/parracho/Desktop/Project/Ultralight_Dark_Matter/Real_pseudo/Rates/New_higgs/Rates(mh2>mH)2.dat',unpack=True,skiprows=1)


[fig, axes] = plt.subplots(nrows=1, ncols=2, figsize=(40, 15), dpi=100)
plt.subplots_adjust(hspace = 1)
plt.style.use('classic')
plt.rcParams["font.family"] = "serif"
# plt.rcParams["font.weight"] = "ultralight"
plt.rcParams["font.size"] = "9"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rc('text', usetex=True)
plt.rcParams['axes.linewidth'] = 2.5
mpl.rcParams['hatch.linewidth'] = 3  # previous pdf hatch linewidth
mpl.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
ax = plt.gca()

axes[0].loglog(1 / Rate_e1[0, :], Rate_e1[1, :], lw=4, label='$\\alpha_1={}$'.format(latex_float(Rate_e1[6, 1])),rasterized=True)
axes[1].loglog(1 / Rate_e2[0, :], Rate_e2[1, :], lw=4, label='$\\alpha_2={}$'.format(latex_float(Rate_e2[6, 1])),rasterized=True)

axes[0].axhline(y=1, color='k', alpha=0.4, ls='dashed', lw=4)
axes[1].axhline(y=1, color='k', alpha=0.4, ls='dashed', lw=4)


axes[0].text(10, 2, r'$\boldsymbol{\Gamma/H = 1}$', color='k', fontsize=9, zorder=11, rasterized=True)
axes[1].text(10, 2, r'$\boldsymbol{\Gamma/H = 1}$', color='k', fontsize=9, zorder=11, rasterized=True)


axes[0].set_ylim(1.0e-100, 1.0e-20)
axes[0].set_xlim(1.0e-4, 1.0e2)

axes[1].set_ylim(1.0e-100, 1.0e-20)
axes[1].set_xlim(1.0e-4, 1.0e2)

anchored_text_1 = AnchoredText(r'$\mathbf{M_\theta ={%s}}$ GeV' % latex_float(Rate_e1[7, 1]) + '\n'\
                            r'$\mathbf{M_{h_2} ={%s}}$ GeV' %
                            Rate_e1[8, 1] + '\n'\
                            r'$\mathbf{\nu_s = {%s}}$ GeV' % Rate_e1[5, 1] + '\n'\
                            r'$\mathbf{\lambda_H = {%s}}$' %
                            Rate_e1[4, 1] + '\n'\
                            r'$\mathbf{\lambda_\phi = {%s}}$' % Rate_e1[3, 1] + '\n'\
                            r'$\mathbf{\lambda_{H\phi} = {%s}}$' %
                            latex_float(Rate_e1[2, 1]), loc=3)

anchored_text_2 = AnchoredText(r'$\mathbf{M_\theta ={%s}}$ GeV' % latex_float(Rate_e2[7, 1]) + '\n' \
                            r'$\mathbf{M_{h_2} ={%s}}$ GeV' %
                            Rate_e2[8, 1] + '\n' \
                            r'$\mathbf{\nu_s = {%s}}$ GeV' % Rate_e2[5, 1] + '\n' \
                            r'$\mathbf{\lambda_H = {%s}}$' %
                            Rate_e2[4, 1] + '\n' \
                            r'$\mathbf{\lambda_\phi = {%s}}$' % Rate_e2[3, 1] + '\n' \
                            r'$\mathbf{\lambda_{H\phi} = {%s}}$' %
                            latex_float(Rate_e2[2, 1]), loc=3)



axes[0].add_artist(anchored_text_1)
axes[1].add_artist(anchored_text_2)


axes[0].set_xlabel('$T^{-1}$ [GeV]$^{-1}$')
axes[0].set_ylabel('$\Gamma/H$')
axes[0].grid(ls='-', alpha=0.1)
axes[0].legend(loc=4)

axes[1].set_xlabel('$T^{-1}$ [GeV]$^{-1}$')
axes[1].set_ylabel('$\Gamma/H$')
axes[1].grid(ls='-', alpha=0.1)
axes[1].legend(loc=4)



axes[0].legend(loc="best")
axes[1].legend(loc="best")
plt.figure(1)

plt.show()

fig.savefig("rateh2.pdf")