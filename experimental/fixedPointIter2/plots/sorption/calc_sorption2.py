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

font = {'size'   : 25}
plt.rc('font', **font)
mpl.rcParams['mathtext.default'] = 'regular'
cols = ['r', 'b', 'g']

#1) get kd value from Freundlich sorption parameters for Alanine, Glucose and Acetate 

#Cs = kf *Cl**m 
#CS in [micromol C/kg soil] --> [mol C/ cm**3]
#CL in [micromol C/dm**3 W] --> [mol C/ cm**3 W]
labels = ['Amino acids', 'Carboxylates', 'Sugars']
m_ = [0.8, 0.14, 0.52] 
kf_ = [1.13, 10.09, 6.47] #micromol**(1-m) (dm**3)**m kg**-1
t_eq = np.asarray([60,28,400])/60/24 #day
rhob = 1.4 #[g/cm**3]
g_per_molC = 12.01
#Cs  = Cs_ *  10**-9 * rhob 
#Cl = Cl_ * 10**-9

#From Freundlich sorption parameters to kd value 
#kd = kf*Cl**(m-1)
Cl = np.linspace(0, 1e-5, 1000)*10**9 #[mol C / cm** W] --> [micromol C/dm**3 W]
kd = np.zeros((len(m_)))

figure(figsize=(8, 6), dpi=80)
for i in range(0, len(m_)): 
    m = m_[i]
    kf = kf_[i]
    x = Cl*10**-9
    y = (kf*Cl**(m))*10**-9
    plt.plot(x,y, color = cols[i], label = labels[i])
    
    idx = np.argmax(x>3e-6)
    x1 = x[idx:]
    y1 = y[idx:]
    z = np.polyfit(x1, y1, 1)
    slope = z[0]
    p = np.poly1d(z)
    plt.plot(x1, p(x1), color = cols[i], linestyle = '--', linewidth = 2)
    t = plt.text(x1[int(len(x1)/2)], y1[int(len(x1)/2)], '$K_d$ = '+str('{:.2E}'.format(slope)), color = cols[i]) #'{:.2E}'.format(slope)
    t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))
    kd[i] = slope #cm^3 C / cm^3 soil 

# plt.xlim([1e-8, 5e-5])
plt.xlabel('$C_l$ (mol C $cm^{-3}$ water)') 
plt.ylabel('$C_s$ (mol C $cm^{-3}$ soil)') 
plt.legend()
plt.show()


#2) define kads and kdes from kd values 
#kd = kads/kdes
kdes = 4/t_eq #[1/d]
kads = kd*kdes #[1/d]
np.set_printoptions(formatter={'float': lambda x: format(x, '.2E')})
print(labels, ' sand')
print(kd)
print(kads) 
print(kdes) 


#3) scaling the parameters for sand 
f_scale = np.mean([2.6/19.5,0.15/0.85, 0.25/1.32])
kd_s = kd * f_scale 
kads_s = kads*f_scale 
kdes_s = kdes
print(labels, ' sand')
print(kd_s) 
print(kads_s)
print(kdes_s) 