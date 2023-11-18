import numpy as np

MP     = 2.43*10**18 # Plack Mass GeV
Gf     = 1.16 *10**(-5) # GeV-2
alpha1 = 127.9
MZ     = 91.19 # Z Boson Mass
alpha_ = 1.0/alpha1 
ee     = np.sqrt(4.0*np.pi*alpha_) # Eletric charge
vev    = 1/np.sqrt(np.sqrt(2.0)*Gf) #Vev of Standard 
A      = np.sqrt(np.pi*alpha_)*vev
thetaW = np.arcsin(2*A/MZ)/2.0 #Weinberg angle
MW     = A/np.sin(thetaW) #W boson mass
SW     = np.sin(thetaW) # Sin of Weinberg angle
g      = ee/np.sin(thetaW) #Coupling of SM