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
    soil = {}
    soil[0] = [0.0639, 0.3698, 0.0096, 1.4646, 4.47]
    soil[1] = [0.0619, 0.3417, 0.0132, 1.3258, 2.03]
    soil[36] = [0.0760, 0.3863, 0.0091, 1.4430, 2.99]
    soil[5] = [ 0.0451, 0.3721, 0.0325, 1.4393, 34.03]
    soil[59] = [0.0534, 0.3744, 0.0171, 1.4138, 13.09]

    create_table("envirotype0", soil[0])
    # create_table("envirotype1", soil[1])
    # create_table("envirotype36", soil[36])
    # create_table("data/envirotype5", soil[5])
    # create_table("envirotype59", soil[59])

