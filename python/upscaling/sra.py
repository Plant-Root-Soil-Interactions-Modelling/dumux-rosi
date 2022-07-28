"""
functions for the steady rate approach
"""

import numpy as np

from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import fsolve


def open_sra_lookup(filename):
    """ opens the look-up table from a file, to quickly find soil root interface potential """
    sra_table = np.load(filename + ".npy")
    x = np.load(filename + "_.npy", allow_pickle = True)
    kx_ = x[0]
    sx_ = x[1]
    inner_ = x[2]
    outer_ = x[3]
    # kr_ = x[4]
    # soil = x[5]
    # print(np.min(sx_), np.max(sx_))
    # print("shape", sra_table.shape)
    # print(sra_table.shape)
    # print(inner_.shape)
    # print(outer_.shape)
    return RegularGridInterpolator((kx_, sx_, inner_, outer_), sra_table)  # method = "nearest" fill_value = None , bounds_error=False


def soil_root_interface(rx, sx, inner_kr, rho, sp):
    """
    finds potential at the soil root interface
    
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


def soil_root_interface_table(rx, sx, inner_kr_, rho_, f):
    """
    finds potential at the soil root interface
        
    rx             xylem matric potential [cm]
    sx             bulk soil matric potential [cm]
    inner_kr       root radius times hydraulic conductivity [cm/day] 
    rho            geometry factor [1]
    sp             soil van Genuchten parameters (type vg.Parameters)
    """
    rsx = f((rx, sx, inner_kr_ , rho_))
    return rsx
