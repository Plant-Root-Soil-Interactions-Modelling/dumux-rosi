"""
    checks the transpiration curves
"""

import sys; sys.path.append("../../../modules/"); sys.path.append("../../../../../CPlantBox/");  sys.path.append("../../../../../CPlantBox/src/python_modules")
sys.path.append("../../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../modules/fv/");
sys.path.append("../"); sys.path.append("../scenarios/");

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

Tpavg = 0.6


def jan_transpiration(ctin, dtin = 60. / (24.*60 * 60)):
    return np.maximum(0., np.pi * Tpavg * (np.cos(2 * np.pi * (ctin - 0.5)) + np.cos(2 * np.pi * ((ctin + dtin) - 0.5))) / 2)


def jan_transpiration2(ctin):
    return np.maximum(0., np.pi * Tpavg * np.cos(2 * np.pi * (ctin - 0.5)))


t = np.linspace(0, 20, 1000)
dt = 60. / (24.*60 * 60) * np.ones(t.shape)

f = lambda x: jan_transpiration2(x)
print("daily", integrate.quad(f, 0, 1))

plt.plot(t, jan_transpiration(t, dt))
plt.show()
