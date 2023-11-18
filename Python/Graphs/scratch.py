
import numpy as np

from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as pe
import cmasher as cmr
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.offsetbox import AnchoredText
from latex_form import latex_float
import matplotlib.gridspec as gridspec
interval = np.linspace(0.12,0.9)
colors = cmr.pride(interval)
cmap = LinearSegmentedColormap.from_list('name', colors)


def line_background(lw,col):
    return [pe.Stroke(linewidth=lw, foreground=col), pe.Normal()]


###########################################################
Rate_e1 = np.loadtxt(r'/home/parracho/Desktop/Project/Ultralight_Dark_Matter/Real_pseudo/Cuba-4.2.2/Ultralight (new higgs)/build/bin/Rates-0.dat',unpack=True,skiprows=1)



fig, axes= plt.subplots(nrows=1,ncols=1, figsize=(40, 15), dpi=100)
plt.subplots_adjust(left=0.3,right=0.7,hspace = 1)
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

axes.loglog(1 / Rate_e1[0, :], Rate_e1[1, :], lw=4, label='$\\alpha_1={}$'.format(latex_float(format(Rate_e1[6, 1],".1e"))),rasterized=True)


axes.axhline(y=1, color='k', alpha=0.4, ls='dashed', lw=4)


axes.text(10, 2, r'$\boldsymbol{\Gamma/H = 1}$', color='k', fontsize=9, zorder=11, rasterized=True)


axes.set_ylim(1.0e-60, 1.0e-35)
axes.set_xlim(1.0e-4, 1.0e2)


anchored_text_1 = AnchoredText(r'$\mathbf{M_\theta ={%s}}$ GeV' % latex_float(Rate_e1[7, 1]) + '\n'\
                            r'$\mathbf{M_{h_2} ={%s}}$ GeV' %
                            Rate_e1[8, 1] + '\n'\
                            r'$\mathbf{\nu_s = {%s}}$ GeV' % Rate_e1[5, 1] + '\n'\
                            r'$\mathbf{\lambda_H = {%s}}$' %
                            Rate_e1[4, 1] + '\n'\
                            r'$\mathbf{\lambda_\phi = {%s}}$' % Rate_e1[3, 1] + '\n'\
                            r'$\mathbf{\lambda_{H\phi}= {%s}}$' %
                            latex_float(Rate_e1[2, 1]), loc=3)



axes.add_artist(anchored_text_1)


axes.set_xlabel('$T^{-1}$ [GeV]$^{-1}$')
axes.set_ylabel('$\Gamma/H$')
axes.grid(ls='-', alpha=0.1)
axes.legend(loc=2)

axes.legend(loc="upper left",prop={"size":10})

plt.figure(1)

plt.show()

fig.savefig("rate_scratch.pdf")
#############################################################
