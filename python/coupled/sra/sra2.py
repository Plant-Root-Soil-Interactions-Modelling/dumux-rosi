import sys; sys.path.append("../../modules/"); sys.path.append("../../../../CPlantBox/");  sys.path.append("../../../../CPlantBox/src/python_modules")
sys.path.append("../")

import numpy as np
import matplotlib.pyplot as plt

import van_genuchten as vg
from root_conductivities import *


def mfp(h, soil):
#     return vg.matric_flux_potential(h, soil)
    return vg.fast_mfp[soil](h)


def imfp(mfp, soil):
#     return vg.matric_potential_mfp(h, soil)
    return vg.fast_imfp[soil](mfp)


def sra_nostress(r, p, q_root, q_out, r_in, r_out, soil):
    rho = r_out / r_in
    mfp = vg.fast_mfp[soil](p) + (q_root * r_in - q_out * r_out) * (r ** 2 / r_in ** 2 / (2 * (1 - rho ** 2)) + rho ** 2 / (1 - rho ** 2) * (np.log(r_out / r) - 0.5)) + q_out * r_out * np.log(r / r_out)
    return vg.fast_imfp[soil](mfp)


def sra_stress(r, p, q_out, r_in, r_out, soil):
    rho = r_out / r_in
    mfp = (vg.fast_mfp[soil](p) + q_out * r_out * np.log(1 / rho)) * ((r ** 2 / r_in ** 2 - 1 + 2 * rho ** 2 * np.log(r_in / r)) / (rho ** 2 - 1 + 2 * rho ** 2 * np.log(1 / rho))) + q_out * r_out * np.log(r / r_in)
    return vg.fast_imfp[soil](mfp)


def sra_flux(p, q_root, q_out, r_in, r_out, soil):
    r = r_in
    dx = 1.e-5  # [cm]
    h0 = sra_nostress(r, p, q_root, q_out, r_in, r_out, soil)
    h1 = sra_nostress(r + dx, p, q_root, q_out, r_in, r_out, soil)
    hc = vg.hydraulic_conductivity(h0, soil)
    f = hc * (h1 - h0) / dx
    return f


def stressed_flux(p, q_out, r_in, r_out, soil):
    r = r_in
    dx = 1.e-5  # [cm]
    h0 = -15000
    h1 = sra_stress(r + dx, p, q_out, r_in, r_out, soil)
    hc = vg.hydraulic_conductivity(h0, soil)
    f = hc * (h1 - h0) / dx
    return f


sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soils = [sand, loam, clay]
soil_names = ["sand", "loam", "clay"]

fig, (ax) = plt.subplots(1, 3)

a = 0.04
L = 0.9
r_ = np.linspace(a, L, 200)

q_out = 0.  # something wrong with positive q_out

for j, s in enumerate(soils):
    sp = vg.Parameters(s)
    vg.create_mfp_lookup(sp)
    h0 = np.zeros(r_.shape)
    h0s = np.zeros(r_.shape)
    q_root = stressed_flux(-700, q_out, a, L, sp)
    print(soil_names[j], "flux", q_root)
    for i, r in enumerate(r_):
         h0[i] = sra_nostress(r, -700, q_root, q_out, a, L, sp)
         h0s[i] = sra_stress(r, -700, q_out, a, L, sp)

    h0s[0] = -15000
    ax[j].plot(r_, h0, "b", label="no stress")
    ax[j].plot(r_, h0s, "r", label="stress")
    ax[j].set_title(soil_names[j])
    print(soil_names[j], "flux", q_root, h0[0], h0s[1])
    ax[j].legend()

plt.show()

