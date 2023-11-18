
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
Rate_e1 = np.loadtxt('/home/parracho/Desktop/Project/Ultralight_Dark_Matter/Real_pseudo/Rates/Real_pseudo/Rates(mh2=mH)1.dat',unpack=True,skiprows=1)
Rate_e2 = np.loadtxt('/home/parracho/Desktop/Project/Ultralight_Dark_Matter/Real_pseudo/Rates/Real_pseudo/Rates(mh2=mH)2.dat',unpack=True,skiprows=1)
Rate_e3 = np.loadtxt('/home/parracho/Desktop/Project/Ultralight_Dark_Matter/Real_pseudo/Rates/Real_pseudo/Rates(mh2=mH)3.dat',unpack=True,skiprows=1)
Rate_g1 = np.loadtxt('/home/parracho/Desktop/Project/Ultralight_Dark_Matter/Real_pseudo/Rates/Real_pseudo/Rates(mh2>mH)1.dat',unpack=True,skiprows=1)
Rate_g2 = np.loadtxt('/home/parracho/Desktop/Project/Ultralight_Dark_Matter/Real_pseudo/Rates/Real_pseudo/Rates(mh2>mH)2.dat',unpack=True,skiprows=1)
Rate_g3 = np.loadtxt('/home/parracho/Desktop/Project/Ultralight_Dark_Matter/Real_pseudo/Rates/Real_pseudo/Rates(mh2>mH)3.dat',unpack=True,skiprows=1)
Rate_f1 = np.loadtxt('/home/parracho/Desktop/Project/Ultralight_Dark_Matter/Real_pseudo/Rates/Real_pseudo/Rates(mh2<mH)1.dat',unpack=True,skiprows=1)
Rate_f2 = np.loadtxt('/home/parracho/Desktop/Project/Ultralight_Dark_Matter/Real_pseudo/Rates/Real_pseudo/Rates(mh2<mH)2.dat',unpack=True,skiprows=1)
Rate_f3 = np.loadtxt('/home/parracho/Desktop/Project/Ultralight_Dark_Matter/Real_pseudo/Rates/Real_pseudo/Rates(mh2<mH)3.dat',unpack=True,skiprows=1)

axes=[]

gs = gridspec.GridSpec(2, 4)
gs.update(left=0.2,right=0.8,top=0.9,bottom=0.1,hspace=0.2,wspace=0.5)
axes.append(plt.subplot(gs[0, :2], ))
axes.append(plt.subplot(gs[0, 2:]))
axes.append(plt.subplot(gs[1, 1:3]))
#fig1, axes= plt.subplots(nrows=2,ncols=2, figsize=(40, 15), dpi=180)
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

axes[0].loglog(1 / Rate_e1[0, :], Rate_e1[1, :], lw=4, label='$\\alpha_1={:.3e}$'.format(latex_float(Rate_e1[6, 1])),rasterized=True)
axes[0].loglog(1 / Rate_e2[0, :], Rate_e2[1, :], lw=4, label='$\\alpha_2={:.3e}$'.format(latex_float(Rate_e2[6, 1])),rasterized=True)
axes[0].loglog(1 / Rate_e3[0, :], Rate_e3[1, :], lw=4, label='$\\alpha_3={:.3e}$'.format(latex_float(Rate_e3[6, 1])),rasterized=True)
axes[1].loglog(1 / Rate_g1[0, :], Rate_g1[1, :], lw=4, label='$\\alpha_1={:.3e}$'.format(latex_float(Rate_g1[6, 1])),rasterized=True)
axes[1].loglog(1 / Rate_g2[0, :], Rate_g2[1, :], lw=4, label='$\\alpha_2={:.3e}$'.format(latex_float(Rate_g2[6, 1])),rasterized=True)
axes[1].loglog(1 / Rate_g3[0, :], Rate_g3[1, :], lw=4, label='$\\alpha_3={:.3e}$'.format(latex_float(Rate_g3[6, 1])),rasterized=True)
axes[2].loglog(1 / Rate_f1[0, :], Rate_f1[1, :], lw=4, label='$\\alpha_1={:.3e}$'.format(latex_float(Rate_f1[6, 1])),rasterized=True)
axes[2].loglog(1 / Rate_f2[0, :], Rate_f2[1, :], lw=4, label='$\\alpha_2={:.3e}$'.format(latex_float(Rate_f2[6, 1])),rasterized=True)
axes[2].loglog(1 / Rate_f3[0, :], Rate_f3[1, :], lw=4, label='$\\alpha_3={:.3e}$'.format(latex_float(Rate_f3[6, 1])),rasterized=True)

axes[0].axhline(y=1, color='k', alpha=0.4, ls='dashed', lw=4)
axes[1].axhline(y=1, color='k', alpha=0.4, ls='dashed', lw=4)
axes[2].axhline(y=1, color='k', alpha=0.4, ls='dashed', lw=4)

axes[0].text(10, 2, r'$\boldsymbol{\Gamma/H = 1}$', color='k', fontsize=9, zorder=11, rasterized=True)
axes[1].text(10, 2, r'$\boldsymbol{\Gamma/H = 1}$', color='k', fontsize=9, zorder=11, rasterized=True)
axes[2].text(10, 2, r'$\boldsymbol{\Gamma/H = 1}$', color='k', fontsize=9, zorder=11, rasterized=True)

axes[0].set_ylim(1.0e-25, 1.0e7)
axes[0].set_xlim(1.0e-4, 1.0e2)

axes[1].set_ylim(1.0e-25, 1.0e7)
axes[1].set_xlim(1.0e-9, 1.0e2)

axes[2].set_ylim(1.0e-10, 1.0e14)
axes[2].set_xlim(1.0e-4, 1.0e2)

anchored_text_1 = AnchoredText(r'$\mathbf{M_\theta ={%s}}$ GeV' % latex_float(Rate_e1[7, 1]) + '\n'\
                            r'$\mathbf{M_{h_2} ={%s}}$ GeV' %
                            Rate_e1[8, 1] + '\n'\
                            r'$\mathbf{\nu_s = {%s}}$ GeV' % Rate_e1[5, 1] + '\n'\
                            r'$\mathbf{\lambda_H = {%s}}$' %
                            Rate_e1[4, 1] + '\n'\
                            r'$\mathbf{\lambda_\phi = {%s}}$' % Rate_e1[3, 1] + '\n'\
                            r'$\mathbf{\lambda_{H\phi}^{(1)} = {%s}}$' %
                            latex_float(Rate_e1[2, 1])+ '\n'\
                            r'$\mathbf{\lambda_{H\phi}^{(2)}= {%s}}$' %
                            latex_float(Rate_e2[2, 1])+ '\n'\
                            r'$\mathbf{\lambda_{H\phi}^{(3)}= {%s}}$' %
                            latex_float(Rate_e3[2, 1]), loc=3)

anchored_text_2 = AnchoredText(r'$\mathbf{M_\theta ={%s}}$ GeV' % latex_float(Rate_g1[7, 1]) + '\n' \
                            r'$\mathbf{M_{h_2} ={%s}}$ GeV' %
                            Rate_g1[8, 1] + '\n' \
                            r'$\mathbf{\nu_s = {%s}}$ GeV' % Rate_g1[5, 1] + '\n' \
                            r'$\mathbf{\lambda_H = {%s}}$' %
                            Rate_g1[4, 1] + '\n' \
                            r'$\mathbf{\lambda_\phi = {%s}}$' % Rate_g1[3, 1] + '\n' \
                            r'$\mathbf{\lambda_{H\phi}^{(1)} = {%s}}$' %
                            latex_float(Rate_g1[2, 1])+ '\n'\
                            r'$\mathbf{\lambda_{H\phi}^{(2)}= {%s}}$' %
                            latex_float(Rate_g2[2, 1])+ '\n'\
                            r'$\mathbf{\lambda_{H\phi}^{(3)}= {%s}}$' %
                            latex_float(Rate_g3[2, 1]), loc=3)

anchored_text_3 = AnchoredText(r'$\mathbf{M_\theta ={%s}}$ GeV' % latex_float(Rate_f1[7, 1]) + '\n'\
                            r'$\mathbf{M_{h_2} ={%s}}$ GeV' %
                            Rate_f1[8, 1] + '\n'\
                            r'$\mathbf{\nu_s = {%s}}$ GeV' % Rate_f1[5, 1] + '\n'\
                            r'$\mathbf{\lambda_H = {%s}}$' %
                            Rate_f1[4, 1] + '\n'\
                            r'$\mathbf{\lambda_\phi = {%s}}$' % Rate_f1[3, 1] + '\n'\
                            r'$\mathbf{\lambda_{H\phi}^{(1)} = {%s}}$' %
                            latex_float(Rate_f1[2, 1])+ '\n'\
                            r'$\mathbf{\lambda_{H\phi}^{(2)}= {%s}}$' %
                            latex_float(Rate_f2[2, 1])+ '\n'\
                            r'$\mathbf{\lambda_{H\phi}^{(3)}= {%s}}$' %
                            latex_float(Rate_f3[2, 1]), loc=3)

axes[0].add_artist(anchored_text_1)
axes[1].add_artist(anchored_text_2)
axes[2].add_artist(anchored_text_3)

axes[0].set_xlabel('$T^{-1}$ [GeV]$^{-1}$')
axes[0].set_ylabel('$\Gamma/H$')
axes[0].grid(ls='-', alpha=0.1)
axes[0].legend(loc=2)

axes[1].set_xlabel('$T^{-1}$ [GeV]$^{-1}$')
axes[1].set_ylabel('$\Gamma/H$')
axes[1].grid(ls='-', alpha=0.1)
axes[1].legend(loc=2)

axes[2].set_xlabel('$T^{-1}$ [GeV]$^{-1}$')
axes[2].set_ylabel('$\Gamma/H$')
axes[2].grid(ls='-', alpha=0.1)
axes[2].legend(loc=2)

axes[0].legend(loc="upper left",prop={"size":8})
axes[1].legend(loc="upper left",prop={"size":8})
axes[2].legend(loc="best")
plt.figure(1)

plt.show()

fig.savefig("ratedm.pdf")
#############################################################
