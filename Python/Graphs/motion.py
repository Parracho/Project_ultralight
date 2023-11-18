import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt
def f(u,t):
    mtheta=10**(-20)
    return (u[1],-3/2*t*u[1]+mtheta**2/4*u[0])



thetai=10
# initial condition
y0 = [thetai,0]

# time points
t = np.linspace(0.01,1e18,200000)

# solve ODE
us = odeint(f,y0,t)
ys=us[:,0]
# plot
plt.ylim(-1, 20)
plt.xlim(0.01, 1.0e18)
plt.plot(t,ys)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()