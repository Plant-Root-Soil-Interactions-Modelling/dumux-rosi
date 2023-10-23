import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

"""
plots a soil retention curve
"""

import functional.van_genuchten as vg
import numpy as np
import matplotlib.pyplot as plt

sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil = clay

wilting_point = -15000

sp = vg.Parameters(soil)

print("matric potential at residual water content", vg.pressure_head(0.1 * 1.1, sp))
print("water content at -1664.67 cm", vg.water_content(-1664.67, sp))
print("water content at wilting point", vg.water_content(wilting_point, sp))

k = vg.hydraulic_conductivity(wilting_point, sp)
print("hydraulic conductivity at wilting point", k, "cm2/day")

h_ = np.linspace(-1e3, 0., 1000)
theta = np.zeros(h_.shape)
for j, h in enumerate(h_):
    theta[j] = vg.water_content(h, sp)

theta_ = np.linspace(soil[0], soil[1], 1000)
k = np.zeros(theta_.shape)
for j, t in enumerate(theta_):
    k[j] = vg.hydraulic_conductivity(t, sp)

plt.semilogy(theta, -h_)
plt.xlabel("water content $\Theta$ (1)")
plt.ylabel("soil water potential $-\psi$ (cm)")
plt.show()

plt.plot(theta_, k)
plt.xlabel("water content $\Theta$ (1)")
plt.ylabel(" hydraulic conductivity $k$ (cm/day)")
plt.show()

