import sys; sys.path.append("../../modules/"); sys.path.append("../../../../CPlantBox/");  sys.path.append("../../../../CPlantBox/src/python_modules")
sys.path.append("../../../build-cmake/cpp/python_binding/"); sys.path.append("../../modules/fv/");

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def open_sra_lookup(filename):
    """ opens the look from a file """
    sra_table = np.load(filename + ".npy")
    x = np.load(filename + "_.npy", allow_pickle=True)
    kx_ = x[0]
    sx_ = x[1]
    inner_ = x[2]
    outer_ = x[3]
    print(np.min(sx_), np.max(sx_))
    # kr_ = x[4]
    # soil = x[5] 
    print(sra_table.shape)
    print(inner_.shape)
    print(outer_.shape)
    return RegularGridInterpolator((kx_, sx_, outer_, inner_), sra_table)  # method = "nearest" fill_value = None , bounds_error=False

# TODO create table 2 again, witch outer_, and inner_ !!!!
