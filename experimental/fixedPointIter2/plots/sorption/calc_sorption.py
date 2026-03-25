import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from scipy import interpolate
import pandas as pd
import sys
from matplotlib.pyplot import figure

font = {'size'   : 20}
plt.rc('font', **font)
mpl.rcParams['mathtext.default'] = 'regular'
cols = ['r', 'b', 'g']

#1) get kd value from Freundlich sorption parameters for Alanine and Acetate 

#Cs = kf *Cl**m 
#CS in [micromol C/kg soil] --> [mol C/ cm**3]
#CL in [micromol C/dm**3 W] --> [mol C/ cm**3 W]
labels = ['Alanine', 'Acetate', 'Glucose']
m_ = [0.8, 0.14, 0.52] 
kf_ = [1.13, 10.09, 6.47] #micromol**(1-m) (dm**3)**m kg**-1
rhob = 1.4 #[g/cm**3]
g_per_molC = 12.01
#Cs  = Cs_ *  10**-9 * rhob 
#Cl = Cl_ * 10**-9

#From Freundlich sorption parameters to kd value 
#kd = kf*Cl**(m-1)
Cl = np.linspace(1e-8, 5e-5, 1000)*10**9 #[mol C / cm** W] --> [micromol C/dm**3 W]
kd_lim = [0.25, 0.005, 0.11]

figure(figsize=(8, 6), dpi=80)
for i in range(0, len(m_)): 
    m = m_[i]
    kf = kf_[i]
    x = Cl*10**-9
    y = kf*Cl**(m-1)*rhob
    plt.plot(x,y, color = cols[i], label = labels[i])
    plt.plot([0,1],[kd_lim[i], kd_lim[i]], color = cols[i], linestyle = '--', label = labels[i])
    print(y[-10:], labels[i])

plt.xlim([1e-8, 5e-5])
plt.xscale('log')
plt.xlabel('Dissolved exudate concentration (mol C $cm^{-3}$ water)') 
plt.ylabel('kd value ($cm^3$ W $cm^{-3}$ soil)') 
plt.legend()
plt.show()


#2) define kads and kdes from kd values 
#kd = kads/kdes
kd_l = np.asarray([0.25, 0.005, 0.11]) #cm^3 C / cm^3 soil 


#Alanine (Amino acids) - high
#Desorption is faster than adsorption → weak net sorption
#Equilibrium is reached quickly (hours to <1 day)
#Alanine remains mobile

#Acetate (Carboxylates) - low
#Very weak retention
#Rapid equilibration
#Nearly conservative transport
#Retardation factor very close to 1
#Transport dominated by advection + biodegradation rather than sorption

#Glucose (Sugars) - medium
#Fast reversible sorption
#Weak retention (consistent with low Kd)
#Rapid desorption → high mobility
#Transport dominated more by biodegradation than sorption

kads_l = np.asarray([0.9*g_per_molC/rhob, 0.072*g_per_molC/rhob, 0.8*g_per_molC/rhob]) #[cm^3/mol/d]
kdes_l = np.asarray([5, 20, 10]) #[1/d]

print(labels, ' loam')
print(kd_l)
print(kads_l/(g_per_molC/rhob))
print(kdes_l)


#3) scaling the parameters for sand 
#Sand vs loam: 
#much less sorption capacity, slower adsorption, faster release --> much higher mobility an leaching risk 
silt_clay_loam = 67.4 #%
silt_clay_sand = 8.2 #%
f_scale = silt_clay_sand / silt_clay_loam 
kd_s = kd_l * f_scale 
kads_s = kads_l*f_scale 
kdes_s = kdes_l/f_scale 
print(labels, ' sand')
print(kd_s) 
print(kads_s/(g_per_molC/rhob))
print(kdes_s) 