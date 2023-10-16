import sys; sys.path.append("../../../CPlantBox/src");

import functional.van_genuchten as vg

import numpy as np
from scipy.optimize import fsolve

"""
    creates a 4d look up table, keeps soil and kr constant varies rx, sx, radial conductivity, and geometry factor rho
"""


def soil_root_interface(rx, sx, inner_kr, rho, sp):
    """
    rx             xylem matric potential [cm]
    sx             bulk soil matric potential [cm]
    inner_kr       root radius times hydraulic conductivity [cm/day] 
    rho            geometry factor [1]
    sp             soil van Genuchten parameters (type vg.Parameters)
    """

    k_soilfun = lambda hsoil, hint: (vg.fast_mfp[sp](hsoil) - vg.fast_mfp[sp](hint)) / (hsoil - hint)

    # rho = outer_r / inner_r  # Eqn [5]
    rho2 = rho * rho  # rho squared
    # b = 2 * (rho2 - 1) / (1 + 2 * rho2 * (np.log(rho) - 0.5))  # Eqn [4]
    b = 2 * (rho2 - 1) / (1 - 0.53 * 0.53 * rho2 + 2 * rho2 * (np.log(rho) + np.log(0.53)))  # Eqn [7]

    fun = lambda x: (inner_kr * rx + b * sx * k_soilfun(sx, x)) / (b * k_soilfun(sx, x) + inner_kr) - x
    rsx = fsolve(fun, rx)

    return rsx


# sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
# clay = [0.1, 0.4, 0.01, 1.1, 10]
hydrus_loam = [0.078, 0.43, 0.036, 1.56, 24.96]
hydrus_clay = [0.068, 0.38, 0.008, 1.09, 4.8]
hydrus_sand = [0.045, 0.43, 0.145, 2.68, 712.8]
hydrus_sandyloam = [0.065, 0.41, 0.075, 1.89, 106.1]

filename = "table_loam" 

# # for single root comparison
# alpha = 0.0383  # (cm-1) soil
# n = 1.3774
# Ks = 60.  # (cm d-1)
# soil = [0.025, 0.403, alpha, n, Ks]

# soil = [0.08, 0.43, 0.018, 1.8, 28.46]  # standard loam
# soil = [0.08, 0.43, 0.04, 1.6, 50]  # c12 loam
soil = loam 

sp = vg.Parameters(soil)
vg.create_mfp_lookup(sp, -1.e5, 1000)  # creates the matrix flux potential look up table

rxn = 150
rx_ = -np.logspace(np.log10(1.), np.log10(16000), rxn)
rx_ = rx_ + np.ones((rxn,))
rx_ = rx_[::-1]

sxn = 150
sx_ = -np.logspace(np.log10(1.), np.log10(16000), sxn)
sx_ = sx_ + np.ones((sxn,))
sx_ = sx_[::-1]

akrn = 100
akrn_ = np.logspace(np.log10(1.e-7), np.log10(1.e-4), akrn)

rhon = 30
rho_ = np.logspace(np.log10(1.), np.log10(200.), rhon)

print("calculating", rxn * sxn * rhon * akrn, "supporting points")

interface = np.zeros((rxn, sxn, akrn, rhon))
for i, rx in enumerate(rx_):
    print(i)
    for j, sx in enumerate(sx_):
        for k, akr in enumerate(akrn_):
            for l, rho in enumerate(rho_):
                interface[i, j, k, l] = soil_root_interface(rx, sx, akr, rho, sp)

np.save(filename, interface)  # data
np.save(filename + "_", [rx_, sx_, akrn_, rho_, soil])  # meta data
