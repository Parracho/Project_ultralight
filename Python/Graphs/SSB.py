import math

import matplotlib.pyplot as plt
import numpy as np


phi = np.linspace(-10,10,1000, dtype=np.float128)
lam=0.5
miu=0.1
miub=-10
V=[]
Vb=[]

for x in phi:
    Vi=lam/4*(x**4)+miu**2*x**2
    Vib = lam / 4 * (x ** 4) + miub ** 2 * x ** 2 +
    V.append(Vi)
    Vb.append(Vib)


plt.plot(phi,V,phi,Vb)
plt.xlim([-1,1])
plt.ylim([-1,5])
plt.xlabel("$\phi$")
plt.ylabel("$V(\phi)$")
plt.legend(["Stable ground state","Unstable ground state"])
plt.grid()
plt.show()

