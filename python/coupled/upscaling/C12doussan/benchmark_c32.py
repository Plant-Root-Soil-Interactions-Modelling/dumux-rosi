""" 
Jan's new scenarios with the upscaled sink (without SRA)

1d soil, 2 cm thick layers, dynamic conductivities (see root_conductivities.py)
"""
import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

from scenario_setup import *

import timeit
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy import sparse
from scipy.linalg import norm
import scipy.sparse.linalg as LA
import visualisation.vtk_plot as vp

""" 
Benchmark M3.2 Root system: steady state small root system solved with the Python/cpp Hybrid solver
(does not work in parallel)
"""
fig, (ax1) = plt.subplots(1, 1)

""" Parameters """
kz = 4.32e-2  # [cm^3/day]
kr = 1.728e-4  # [1/day]
p_s = -200  # static soil pressure [cm]
p0 = -500  # dircichlet bc at top
sim_time = 14  # [day] for task b

""" root problem """
r0 = XylemFluxPython("../../../grids/RootSystem.rsml")
r0.setKr([kr])
r0.setKx([kz])

soil_index = lambda x, y, z: 0
r0.rs.setSoilGrid(soil_index)

r = HydraulicsDoussan("../../../grids/RootSystem.rsml")
r.setKr([kr])
r.setKx([kz])
r.rs.setSoilGrid(soil_index)

nodes = r.get_nodes()
segments = r.rs.segments
n = len(nodes)
ns = len(segments)

""" Numerical solution (a) """

Id = sparse.identity(ns).tocsc()  # identity matrix

print()
kx_ = np.divide(r.getKx(sim_time), r.rs.segLength())  # / dl
Kx = sparse.diags(kx_).tocsc()
print("Kx", Kx.shape, Kx[0, 0], Kx[1, 1])

kr_ = np.array(r.getEffKr(sim_time))  # times surface (2 a pi length)
Kr = sparse.diags(kr_).tocsc()
print("Kr", Kr.shape, Kr[0, 0], Kr[1, 1])

IM = r.get_incidence_matrix().tocsc()
IMt = IM.transpose().tocsc()

# sum_columns = np.sum(IM, axis = 0)
# tips = np.argwhere(sum_columns == 1)  # tip node indices

L = IMt @ Kx @ IM  # Laplacian, Hess Eqn (10)
L = L[1:, 1:]  # == L_{N-1}
print("L", L.shape)

print("inv start")

A_dirichlet = (L + Kr).tocsc()
Ainv_dirichlet = sparse.linalg.inv(A_dirichlet).todense()  # dense

A_neumann = A_dirichlet
A_neumann[0, 0] -= kx_[0]
Ainv_neumann = sparse.linalg.inv(A_neumann).todense()  # dense

# C_comp_dirichlet = Kr @ (Id - Ainv_dirichlet @ Kr)  # Neumann, Hess, Eqn (24)
# c_dirichlet = (Kr @ Ainv_dirichlet)[:, 0] * (-kx_[0])  # # Hess (25)
# print("C_comp_dirichlet", type(C_comp_dirichlet), C_comp_dirichlet.shape)
# print("c_dirichlet", type(c_dirichlet), c_dirichlet.shape)
#
# C_comp_neumann = Kr @ (Id - Ainv_neumann @ Kr)  # Neumann, Hess, Eqn (32)
# c_neumann = (Kr @ Ainv_neumann)[:, 0]  # Hess (33)
# print("C_comp_neumann", type(C_comp_neumann), C_comp_neumann.shape)
# print("c_neumann", type(c_neumann), c_neumann.shape)

print("inv stop")

hs = np.ones((nodes.shape[0] - 1, 1)) * p_s
rx_a0 = r0.solve_dirichlet(0., p0, p_s, [p_s], True)

rx_a1 = (Ainv_dirichlet).dot(Kr.dot(hs)) + Ainv_dirichlet[:, 0] * kx_[0] * p0  # wrong, we need total potentials in hs

for i in range(0, len(nodes) - 1):  # from matric to total
    hs[i, 0] += nodes[i + 1][2]

rx_a3 = Ainv_dirichlet.dot(Kr.dot(hs)) + Ainv_dirichlet[:, 0] * kx_[0] * p0  #   # Hess Eqn (29)
rx_a4 = Ainv_neumann.dot(Kr.dot(hs)) + Ainv_neumann[:, 0] * -1.3636

for i in range(0, len(nodes) - 1):  # from total to matric
    rx_a3[i, 0] -= nodes[i + 1][2]
    rx_a4[i, 0] -= nodes[i + 1][2]

print("Transpiration (meunier)", r0.collar_flux(0., rx_a0, [p_s]), "cm3/day")
print("Transpiration (tot)", r.collar_flux(0., rx_a3), "cm3/day")
print("Transpiration (neumann)", r.collar_flux(0., rx_a4), "cm3/day")
# np.savetxt("results_m32a", np.vstack((nodes[1:, 2], rx_a1)), delimiter = ',')

ax1.plot(rx_a0, nodes[:, 2] , "r.", label = "Meunier")
ax1.plot(rx_a3, nodes[1:, 2] , "g.", label = "Doussan (tot)")
ax1.plot(rx_a4, nodes[1:, 2] , "b.", label = "Neumann (tot)")
# ax1.plot(rx_a2, nodes[1:, 2] , "b.", label = "Doussan (grav)")

# ax1.plot(rx_a2, nodes[1:, 2] , "b*", label = "Doussan (Neumann)")
ax1.set_xlabel("Xylem pressure (cm)")
ax1.set_ylabel("Depth (m)")
ax1.set_title("Constant conductivities")
ax1.legend()

plt.show()
