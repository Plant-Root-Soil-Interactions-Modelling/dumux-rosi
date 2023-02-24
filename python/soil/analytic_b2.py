"""
produces the analytical solution Figure 3
from Vanderborght et al. (2005)

D. Leitner, 2018
"""
import sys; sys.path.append("../modules"); sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import functional.van_genuchten as vg

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

loam = vg.Parameters([0.08, 0.43, 0.04, 1.6, 50])

Jw = -0.5  # cm/day

K = lambda psi: vg.hydraulic_conductivity(psi, loam)
F = lambda psi: 1. / (Jw / K(psi) - 1.)

psi_ = np.linspace(0, -300, 300)  # psi(-54) = 0
dz = np.zeros(len(psi_),)
for i in range(0, len(psi_)):
    dz[i], err = integrate.quad(F, 0, psi_[i])
z1 = dz - 53  # this value is not clear to me any more

plt.plot(psi_, z1)
plt.xlabel('$\psi$ (cm)')
plt.ylabel('Depth (cm)')
plt.xlim(-300, 0)
plt.ylim(-60, 0)

if __name__ == "__main__":
    plt.show()
