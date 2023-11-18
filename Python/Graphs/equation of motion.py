import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
from matplotlib import rc
from matplotlib.pyplot import gca,show
from IPython.display import Math, display
import numpy as np
import scipy as sp
import sympy as sym
from scipy.special import kv,zeta, polygamma, factorial, erf
from scipy import integrate
from matplotlib.offsetbox import AnchoredText
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as pe
import cmasher as cmr
import mpmath as mp
from scipy.integrate import odeint, solve_ivp, ode, quad
import math
from mpmath import *
from sympy import Eq, Symbol, solve
from scipy.interpolate import griddata
import constants
from ipynb.fs.full.parameters import gs
from ipynb.fs.full.parameters import ge


k = 0.301/np.sqrt(80)*constants.MP


def Hubble(T):
    return 1.66 *np.sqrt(100) * T**2/constants.MP

def HI(T,H_I):
    T_I = H_I/2/np.pi
    if T > T_I:
        return H_I
    else:
        return (1.66 *np.sqrt(100) * T**2/constants.MP)

def time(T):
    return 0.301*(gs(T))**(-1/2)*constants.MP/T**2

def m_(x):
    return 1.68e-7*(0.4**4/vev**2)*(0.4*x)**(6.68)


def energy_density(m,phidot,phi,x):
    result = np.zeros(len(x))
    for i in range(len(result)):
        result[i] = (1/2/k/x[i])**2*phidot[i]**2/2 + m**2*phi[i]**2/2
    return result

def pressure(m,phidot,phi,x):
    result = np.zeros(len(phidot))
    for i in range(len(result)):
        result[i] = (1/2/k/x[i])**2*phidot[i]**2/2 - m**2*phi[i]**2/2
    return result

def Abundance(m,theta_i,gosc):
    return 7.02e-65*m**2*theta_i**2/gosc/(m/np.sqrt(gosc))**(3/2)/8.15e-47


def dy_(y,t):
    """
    The right-hand side of the damped oscillator ODE
    """
    x, p = y[0], y[1]
    dx = p
    dp = (1/t - 6/t)*p - 4*k**2*t**2*m**2*x
    return [dx, dp]


def dy_Quartico(y,t):
    """
    The right-hand side of the damped oscillator ODE
    """
    x, p = y[0], y[1]
    dx = p
    dp = (1/t - 6/t)*p - 4*k**2*t**2*m**2*x + 4*k**2*m**2/6/vev**2 * x**3
    return [dx, dp]


def Number_of_particles(theta_dot, theta, t):
    result = list(np.zeros(len(theta)))
    for i in range(len(result)):
        result[i] = (1/2*(theta_dot[i]/2/k/t[i])**2 + m_(t[i])**2*(1 - np.cos(theta[i])))*(k*t[i]**2)**(3/2)/m_(t[i])
    return result


m = 1.0e-29
t = np.logspace(-16,7,10000)

y0 = [1.0e0, 0.0]

y_ULA = odeint(dy_, y0,t)

#y_ULA = solve_ivp(dy_,t_eval = t,y0 = [1.0e1,0.0], t_span=[1.0e-15,1.0e10])

t_ = np.logspace(-16,7, 10000)


y0 = [1.0e0, 0.0]

vev = 1.0e-2

y_ULA_quartico = odeint(dy_Quartico, y0,t_)
#y_ULA_quartico = solve_ivp(dy_Quartico,t_span = [1.0e-17,1.0e15],y0 = [1.0e10,0.0],method = 'DOP853')


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

#plt.axvline(x = 1/np.sqrt(m*constants.MP/1.66/2.7/np.sqrt(10)), color = 'gray', lw = 3, ls = 'dashed')



axes.plot(t,y_ULA[:,0] , lw = 5)
#plt.plot(y_ULA.t, y_ULA.y[0,:], lw = 3)



axes.plot(t_,y_ULA_quartico[:,0],color = 'r' ,lw = 2)

#plt.plot(y_ULA_quartico.t, y_ULA_quartico.y[0,:], lw = 3)


anchored_text = AnchoredText(r"$\mathbf{m_\theta = 1.0 \times 10^{-15}}$ eV" + '\n'\
                             r"$\mathbf{\Theta_i = 1}$", loc=3,prop=dict(size=25),frameon=True)
#ax.add_artist(anchored_text)

axes.add_artist(anchored_text)
#plt.text(1.4e3, 0.0e0, r'$\mathbf{T_{osc} \simeq 0.4}$ MeV', fontsize=14, color = 'k', alpha = 1, rotation = 90, weight = 'light')

#axes.set_xscale("log", basex=10)
axes.axvline(x=5e5, color='k', alpha=0.4, ls='dashed', lw=4)
axes.axhline(y=0, color='k')
axes.text(6e5, 0.1, r'$\boldsymbol{T_{\textrm{osc}}}$', color='k', fontsize=25, zorder=11, rasterized=True,rotation=90,alpha=0.5)


axes.set_ylabel('$\\Theta(T)$',fontsize=20)
axes.set_xlabel('$T^{-1}$ [GeV]$^{-1}$',fontsize=20)


axes.set_xlim(1.0e-20,8.0e6)

axes.set_ylim(-1.1,1.1e0)


plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
axes.grid(ls='-', alpha=0.1)
axes.legend(loc=2)

plt.figure(1)
plt.tight_layout()
plt.show()

#fig.savefig("relic_density.pdf")
#fig.savefig("theta.pdf", bbox_inches='tight')



