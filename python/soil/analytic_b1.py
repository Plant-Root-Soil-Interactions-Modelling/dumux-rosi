"""
produces the analytical solution Figure 2abc
from Vanderborght et al. (2005)

D. Leitner, 2018
"""
import sys; sys.path.append("../modules"); sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import functional.van_genuchten as vg

import numpy as np
from scipy import optimize
from scipy import integrate
import matplotlib.pyplot as plt

sand = vg.Parameters([0.045, 0.43, 0.15, 3, 1.1574e-04 * 100 * 3600 * 24])
loam = vg.Parameters([0.08, 0.43, 0.04, 1.6, 5.7870e-06 * 100 * 3600 * 24])
clay = vg.Parameters([0.1, 0.4, 0.01, 1.1, 1.1574e-06 * 100 * 3600 * 24])

Jw = 0.5;  # constant downward flow rate [cm/d]

#
# The lower part (constant pressure)
#
Ks = lambda psi: vg.hydraulic_conductivity(psi, sand) - Jw
Kl = lambda psi: vg.hydraulic_conductivity(psi, loam) - Jw
Kc = lambda psi: vg.hydraulic_conductivity(psi, clay) - Jw
psi_s = optimize.brentq(Ks, -100, 0)
psi_l = optimize.brentq(Kl, -100, 0)
psi_c = optimize.brentq(Kc, -100, 0)

 #
# The upper part
#
Ks = lambda psi: vg.hydraulic_conductivity(psi, sand)
Kl = lambda psi: vg.hydraulic_conductivity(psi, loam)
Kc = lambda psi: vg.hydraulic_conductivity(psi, clay)
# integrand Eqn [14]
Fs = lambda psi: 1. / (Jw / Ks(psi) - 1.)
Fl = lambda psi: 1. / (Jw / Kl(psi) - 1.)
Fc = lambda psi: 1. / (Jw / Kc(psi) - 1.)

N = 100  # resolution
dz = np.ones(N,)

psiA = np.linspace(-45, psi_s, N)  # Loam (on sand)
for i in range(0, N):
    dz[i], err = integrate.quad(Fl, psi_s, psiA[i])
zA = dz - 50.

psiB = np.linspace(psi_l, psi_s - 1e-10, N)  # Sand (on loam)
for  i in range(0, N):
    dz[i], err = integrate.quad(Fs, psi_l, psiB[i]);
zB = dz - 50;

psiC = np.linspace(psi_s, psi_c - 1e-3, 100);  # Clay (on sand)
for  i in range(0, N):
    dz[i] , err = integrate.quad(Fc, psi_s, psiC[i]);
zC = dz - 50;

#
# prepare plot
#
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

z_ = np.linspace(-50, -100, 2)
ax1.plot([psi_s, psi_s], z_, 'b')  # sand
ax1.plot(psiA, zA, 'b');
ax1.set_xlabel('$\psi$ (cm)')
ax1.set_ylabel('Depth (cm)')
ax1.set_xlim(-50, -10)
ax1.set_ylim(-100, 0)

z_ = np.linspace(-50, -100, 2)
ax2.plot([psi_l, psi_l], z_, 'b')
ax2.plot(psiB, zB, 'b');
ax2.set_xlabel('$\psi$ (cm)')  # loam
ax2.set_xlim(-50, -10)
ax2.set_ylim(-60, -40)

z_ = np.linspace(-50, -100, 2)
ax3.plot([psi_s, psi_s], z_, 'b')  # sand
ax3.plot(psiC, zC, 'b');
ax3.set_xlabel('$\psi$ (cm)')
ax3.set_xlim(-28, -8)
ax3.set_ylim(-100, 0)

if __name__ == "__main__":
    plt.show()

