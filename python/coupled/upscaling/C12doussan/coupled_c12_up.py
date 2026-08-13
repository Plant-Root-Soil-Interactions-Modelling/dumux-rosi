""" 
Benchmark C1.2 for a static root system in soil (1D or 3D)
with upscaled C_comp matrix (for calculating the fluxes) 
using Doussan approach in HESS paper notation

similar speed for 3d, great speedup for 1d
"""
import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import timeit
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse

from rhizo_models import plot_transpiration
from scenario_setup import *

r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("3D")
name = "results/c12_up"

# r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("3D", age_dependent = True)
# dt = 36 / (24 * 3600)  # unstable otherwise
# name = "results/c12b_up"

# r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("1D")
# name = "results/c12_up1d"

nodes = r.get_nodes()
ns = len(r.rs.segments)

""" Doussan """
A_dirichlet, Kr, kx0 = r.doussan_system_matrix(rs_age)
Id = sparse.identity(ns).tocsc()  # identity matrix

print("invert matrix start")
Ainv_dirichlet = sparse.linalg.inv(A_dirichlet).todense()  # dense
A_neumann = A_dirichlet.copy()
A_neumann[0, 0] -= kx0
Ainv_neumann = sparse.linalg.inv(A_neumann).todense()  # dense

C_comp_dirichlet = Kr @ (Id - Ainv_dirichlet @ Kr)  # Neumann, Hess, Eqn (24)
c_dirichlet = (Kr @ Ainv_dirichlet)[:, 0] * (-kx0)  # # Hess (25)
# print("C_comp_dirichlet", type(C_comp_dirichlet), C_comp_dirichlet.shape)
# print("c_dirichlet", type(c_dirichlet), c_dirichlet.shape)
C_comp_neumann = Kr @ (Id - Ainv_neumann @ Kr)  # Neumann, Hess, Eqn (32)
c_neumann = (Kr @ Ainv_neumann)[:, 0]  # Hess (33)
# print("C_comp_neumann", type(C_comp_neumann), C_comp_neumann.shape)
# print("c_neumann", type(c_neumann), c_neumann.shape)

print("invert matrix stop")

print("upscaling start")

B, soil2matrix, matrix2soil = r.get_soil_matrix()
nmax = len(matrix2soil)  # dimension of the upsaled problem

Bt = B.transpose()
# print(Bt.shape, C_comp_neumann.shape, B.shape)
# print(np.sum(Bt @ B), (Bt @ B).shape)
BBt_inv = sparse.linalg.inv(B @ Bt)  # sparse

C_comp_neumann_up = B @ C_comp_neumann @ Bt
c_neumann_up = B @ c_neumann
C_comp_dirichlet_up = B @ C_comp_dirichlet @ Bt
c_dirichlet_up = B @ c_dirichlet
# print(C_comp_neumann_up.shape, type(C_comp_neumann_up))

print("upscaling end")

""" Numerical solution"""
start_time = timeit.default_timer()
rs_age = np.max(r.get_ages())
x_, y_, z_ = [], [], []
sink1d = []

N = round(sim_time / dt)
t = 0.

rx = [0]

sx = s.getSolutionHead()  # inital condition, solverbase.py
centers = s.getCellCenters()

for i in range(0, N):

    t_pot = -trans * sinusoidal(t)  # potential transpiration ...
    print("t_pot", t_pot)

    # hs = np.transpose(np.array([[sx[mapping[j]][0] for j in range(0, ns)]]))
    # for j in range(0, len(nodes) - 1):  # from matric to total
    #     hs[j, 0] += nodes[j + 1][2]
    #
    # q_dirichlet0 = -(C_comp_dirichlet.dot(hs) + c_dirichlet * wilting_point)
    #
    # if np.sum(q_dirichlet0) > t_pot:
    #     print("dirichlet", np.sum(q_dirichlet0), t_pot)
    #     fluxes = r.sumSegFluxes(q_dirichlet0[:, 0])
    # else:
    #     q_neumann0 = -(C_comp_neumann.dot(hs) - c_neumann * t_pot)
    #     print("neumann", np.sum(q_neumann0), t_pot)
    #     fluxes = r.sumSegFluxes(q_neumann0[:, 0])

    hs_ = np.zeros((nmax, 1))
    for j in soil2matrix.keys():
            hs_[soil2matrix[j]] += sx[j] + centers[j, 2]

    # hs = np.transpose(np.array([[sx[mapping[j]][0] for j in range(0, ns)]]))
    # for j in range(0, len(nodes) - 1):  # from matric to total
    #     hs[j, 0] += nodes[j + 1][2]
    # hs_ = BBt_inv.dot(B.dot(hs))  # can be proven using pseudo inverse

    q_dirichlet0_up = -(C_comp_dirichlet_up.dot(hs_) + c_dirichlet_up * wilting_point)

    if np.sum(q_dirichlet0_up) > t_pot:
        print("q_dirichlet0_up", np.sum(q_dirichlet0_up), t_pot)
        fluxes = {}
        for j in range(0, nmax):
            fluxes[matrix2soil[j]] = q_dirichlet0_up[j, 0]
    else:
        # print("C_comp_neumann_up", C_comp_neumann_up.shape)
        # print("c_neumann_up", c_neumann_up.shape)
        # print("hs_", hs_.shape)
        q_neumann0_up = -(C_comp_neumann_up.dot(hs_) - c_neumann_up * t_pot)
        print("q_neumann0_up", q_neumann0_up.shape, np.sum(q_neumann0_up), t_pot)
        fluxes = {}
        for j in range(0, nmax):
            fluxes[matrix2soil[j]] = q_neumann0_up[j, 0]

    # print()
    # hx = Ainv_neumann.dot(Kr.dot(hs)) + Ainv_neumann[:, 0] * t_pot  #   # Hess Eqn (29)
    # hx0 = hx[0, 0] + t_pot / kx_[0]  # /kx*l
    # print("hx0", hx0, hx[0, 0], hx[1, 0], hx[2, 0])
    # hxd = Ainv_dirichlet.dot(Kr.dot(hs)) + Ainv_dirichlet[:, 0] * kx_[0] * wilting_point
    # print("hxd", hxd[0, 0], hxd[1, 0], hxd[2, 0])
    # print()
    #
    # print()
    # q_neumann = -Kr.dot(hs - hx)
    # print("q_neumann", q_neumann.shape, np.min(q_neumann), np.max(q_neumann), np.argmin(q_neumann), np.sum(q_neumann))
    # print("q_neumann0", q_neumann.shape, np.min(q_neumann0), np.max(q_neumann0), np.argmin(q_neumann0), np.sum(q_neumann0))
    # q_dirichlet = -Kr.dot(hs - hxd)
    # print("q_dirichlet", q_dirichlet.shape, np.min(q_dirichlet), np.max(q_dirichlet), np.argmin(q_dirichlet), np.sum(q_dirichlet))
    # print("q_dirichlet", q_dirichlet0.shape, np.min(q_dirichlet0), np.max(q_dirichlet0), np.argmin(q_dirichlet0), np.sum(q_dirichlet0))
    # # print("sum_q", np.sum(q), q[0], q[1])
    # # hx = [0]
    # print()

    water = s.getWaterVolume()
    s.setSource(fluxes.copy())  # richards.py
    s.solve(dt)
    soil_water = (s.getWaterVolume() - water) / dt

    old_sx = sx.copy()
    sx = s.getSolutionHead()  # richards.py

    if  i % skip == 0:
        x_.append(t)
        sum_flux = 0.
        for f in fluxes.values():
            sum_flux += f
        # print("sum_flux", sum_flux.shape, sum_flux)
        y_.append(soil_water)  # cm3/day
        z_.append(sum_flux)  # cm3/day
        n = round(float(i) / float(N) * 100.)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days {:g}\n"
              .format(np.min(sx), np.max(sx), np.min(rx), np.max(rx), s.simTime, rx[0]))

        # """ Additional sink plot """
        # if i % 60 == 0:  # every 6h
        #     ana = pb.SegmentAnalyser(r.rs)
        #     fluxes = r.segFluxes(rs_age + t, rx, old_sx, False, cells = True)  # cm3/day WRONG sx, should be old sx
        #     ana.addData("fluxes", fluxes)  # cut off for vizualisation
        #     flux1d = ana.distribution("fluxes", max_b[2], min_b[2], 15, False)
        #     sink1d.append(np.array(flux1d))
        #     print("\nSink integral = ", np.sum(np.array(flux1d)), "\n")

    t += dt

s.writeDumuxVTK(name)

""" Plot """
print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")
plot_transpiration(x_, y_, z_, lambda t: trans * sinusoidal(t))
np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter = ';')

# vp.plot_roots_and_soil(r.rs, "pressure head", rx, s, False, min_b, max_b, cell_number, name)  # VTK vizualisation
# sink1d = np.array(sink1d)
# np.save("sink1d", sink1d)
# print(sink1d.shape)
# print(sink1d[-1])

