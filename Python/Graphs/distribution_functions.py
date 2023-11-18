import math

import matplotlib.pyplot as plt
import numpy as np

xmin=0.001
xmax=10
xn=100000
x = np.linspace(xmin,xmax,xn, dtype=np.float128)
fermi=[]
bose=[]
maxw=[]
for i in x:
    try:
        ferm =(1/(math.exp(i)+1))
        bos = (1 / (math.exp(i) - 1))
        max = (math.exp(-i))
        fermi.append(1-ferm)
        bose.append(1+bos)
        maxw.append(1)
    except:
        break

plt.plot(x[0:len(fermi)],fermi,x[0:len(fermi)],bose,x[0:len(fermi)],maxw)
#plt.plot([xmin,xmax],[1,1],"--")
plt.xlim([10e-3,10e0])
plt.ylim([10e-3,10e1])
plt.xlabel("$E_i/T$")
plt.ylabel("$1 \pm f_i^{eq}$")
plt.grid()
plt.legend(["Fermi-Dirac","Bose-Einstein","Maxwell-Boltzmann"])
plt.yscale('log')
plt.xscale('log')
plt.show()