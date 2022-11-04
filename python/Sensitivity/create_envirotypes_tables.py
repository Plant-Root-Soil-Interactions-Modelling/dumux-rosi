import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/"); sys.path.append("../modules/fv/");

import van_genuchten as vg
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


def create_table(filename, soil):

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

    print(filename, "calculating", rxn * sxn * rhon * akrn, "supporting points")

    interface = np.zeros((rxn, sxn, akrn, rhon))
    for i, rx in enumerate(rx_):
        print(i)
        for j, sx in enumerate(sx_):
            for k, akr in enumerate(akrn_):
                for l, rho in enumerate(rho_):
                    interface[i, j, k, l] = soil_root_interface(rx, sx, akr, rho, sp)

    np.save(filename, interface)  # data
    np.save(filename + "_", [rx_, sx_, akrn_, rho_, soil])  # meta data


if __name__ == "__main__":

    # from file "Van Genuchten.xls"

    # soil0 = [0.0809, 0.52, 0.0071, 1.5734, 99.49]
    # create_table("envirotype0", soil0)

    # soil1 = [0.0874, 0.5359, 0.0087, 1.5231, 93]
    # create_table("envirotype1", soil1)

    # soil36 = [0.0942, 0.5569, 0.0089, 1.4974, 87.79]
    # create_table("envirotype36", soil36)

    soil5 = [0.0539, 0.5193, 0.024, 1.4046, 208.78]
    create_table("data/envirotype5", soil5)

    # soil59 = [0.0675, 0.5109, 0.0111, 1.4756, 107.63]
    # create_table("envirotype59", soil59)

