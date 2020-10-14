import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")
sys.path.append("../solvers/")  # for pure python solvers

import van_genuchten as vg

sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil = clay

wilting_point = -15000

sp = vg.Parameters(soil)

print("pressure at residual water content", vg.pressure_head(0.1 * 1.1, sp))
print("water content at -1664.67", vg.water_content(-1664.67, sp))
print("pressure at wilting point", vg.water_content(wilting_point, sp))

k = vg.hydraulic_conductivity(wilting_point, sp)
print(k)
