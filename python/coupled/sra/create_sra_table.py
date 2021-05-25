import sys; sys.path.append("../../modules/"); sys.path.append("../../../../CPlantBox/");  sys.path.append("../../../../CPlantBox/src/python_modules")
from gevent.resolver._addresses import AddressSyntaxError
sys.path.append("../../../build-cmake/cpp/python_binding/"); sys.path.append("../../modules/fv/");

import van_genuchten as vg
import numpy as np
from scipy.optimize import fsolve

"""
    creates a 4d look up table, keeps soil and kr constant varies rx, sx, innter radius, and outer radius
"""


def soil_root_interface(rx, sx, inner_r, outer_r, kr, z, sp):
  
    hintmin = -1.e5
    hintmax = -2.
       
    k_soilfun = lambda hsoil, hint: (vg.fast_mfp[sp](hsoil) - vg.fast_mfp[sp](hint)) / (hsoil - hint)
                              
    if (sx - z) < hintmax:
              
        if (sx - z) > hintmin:
            
            rho = outer_r / inner_r  # Eqn [5]            
            rho2 = rho * rho  # rho squared
            # b = 2 * (rho2 - 1) / (1 + 2 * rho2 * (np.log(rho) - 0.5))  # Eqn [4]
            b = 2 * (rho2 - 1) / (1 - 0.53 * 0.53 * rho2 + 2 * rho2 * (np.log(rho) + np.log(0.53)))  # Eqn [7]
                  
            # fun = lambda x: (inner_r * kr * rx + b * sx * k_soilfun(sx, x)) / (b * k_soilfun(sx, x) + inner_r * kr) - x                                                 
            fun = lambda x: (inner_r * kr * rx + b * sx * k_soilfun(sx, x)) / (b * k_soilfun(sx, x) + inner_r * kr) - x 
            rsx = fsolve(fun, rx)
                                          
        else: 
            rsx = rx
    else:
        rsx = sx
    
    return rsx


kr = 1.73e-4 
soil = [0.08, 0.43, 0.04, 1.6, 50]  # loam
sp = vg.Parameters(soil)  
vg.create_mfp_lookup(sp, -1.e5, 1000)  # creates the matrix flux potential look up table

rxn = 150
rx_ = -np.logspace(np.log10(1.), np.log10(16000), rxn) 
rx_ = rx_ + np.ones((rxn,))
rx_ = rx_[::-1]

sxn = 150
sx_ = -np.logspace(np.log10(1.), np.log10(32000), sxn) 
sx_ = sx_ + np.ones((sxn,))
sx_ = sx_[::-1]

innern = 30
inner_ = np.logspace(np.log10(0.01), np.log10(0.2), innern)  # np.linspace(0.01, 0.2, innern)
print(inner_)

outern = 60
outer_ = np.logspace(np.log10(0.1), np.log10(20.), outern)  # np.linspace(0.1, 1., outern)

np.save("table_", [rx_, sx_, inner_, outer_, kr, soil])  # additional data

print("calculating", rxn * sxn * innern * outern, "supporting points")

interface = np.zeros((rxn, sxn, innern, outern))
for i, rx in enumerate(rx_): 
    print(i)
    for j, sx in enumerate(sx_): 
        for k, inner in enumerate(inner_): 
            for l, outer in enumerate(outer_): 
                interface[i, j, k, l] = soil_root_interface(rx, sx, inner, outer, kr, 0., sp)

np.save("table", interface)  # data

