
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

mtheta=[1e4,1e-1,1e-6]
vsigma=3e1/np.sqrt(50)

xmin=0
xmax=0.20
xn=1000
x = np.linspace(xmin,xmax,xn)
fermi=[]
bose=[]
maxw=[]
for i in x:
    try:
        ferm =np.sqrt((i/0.11)*mtheta[0]**(-1/2))*1e-3/vsigma
        bos =np.sqrt((i/0.11)*mtheta[1]**(-1/2))*1e-3/vsigma
        max = np.sqrt((i / 0.11) * mtheta[2] ** (-1 / 2)) * 1e-3 / vsigma
        fermi.append(ferm)
        bose.append(bos)
        maxw.append(max)
    except:
        break

axes.plot(fermi,x,lw=4, label='\\boldmath$ m_\\theta={}$ eV'.format(latex_float(format(mtheta[0]*1e-14,".1e"))),rasterized=True,color="b")
axes.plot(bose,x,lw=4, label='\\boldmath$ m_\\theta={}$ eV'.format(latex_float(format(mtheta[1]*1e-14,".1e"))),rasterized=True,color="r")
axes.plot(maxw,x,lw=4, label='\\boldmath$ m_\\theta={}$ eV'.format(latex_float(format(mtheta[2]*1e-14,".1e"))),rasterized=True,color="g")
axes.set_xscale("log", basex=10)
axes.set_ylim(0, 0.20)
axes.set_xlim(1.0e-10, 1.0e-1)

axes.axhline(y=0.11, color='k', alpha=0.4, ls='dashed', lw=4)


axes.text(7e-10,0.112, r'$\Omega h^2 \simeq 0.11$', color='k', fontsize=25, zorder=11, rasterized=True)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
axes.set_xlabel('$\Theta_{osc}$',fontsize=25)
axes.set_ylabel('$\Omega h^2$',fontsize=25)
axes.grid(ls='-', alpha=0.1)
axes.legend(loc=2)

axes.legend(loc="upper left",prop={"size":20})

plt.figure(1)
plt.tight_layout()
plt.show()

fig.savefig("relic_density.pdf")
#############################################################



