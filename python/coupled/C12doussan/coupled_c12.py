""" 
Benchmark C1.2 for a static root system in soil (1D or 3D) 
with the classic sink using Doussan approach in HESS paper notation, solves for q 

same results as version with Meunier approach (coupled/coupled_c12.py)
"""
import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

from rhizo_models import plot_transpiration
from scenario_setup import *

import timeit
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import sparse

r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("3D")
name = "results/c12"

# r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("3D", age_dependent = True)
# dt = 36 / (24 * 3600)  # unstable
# name = "results/c12b"

# r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("1D")
# name = "results/c12_1d"

nodes = r.get_nodes()
ns = len(r.rs.segments)

""" Doussan """
A_dirichlet, Kr, kx0 = r.doussan_system_matrix(rs_age)
Id = sparse.identity(ns).tocsc()  # identity matrix

print("invert matrix start")

# A_dirichlet = A.tocsc()
Ainv_dirichlet = sparse.linalg.inv(A_dirichlet)  # dense
# Ainv_dirichlet = scipy.linalg.inv(A_dirichlet) # this does not work for dense A, no idea why
A_neumann = A_dirichlet.copy()
A_neumann[0, 0] -= kx0
Ainv_neumann = sparse.linalg.inv(A_neumann).todense()  # dense
# Ainv_neumann = scipy.linalg.inv(A_neumann)

C_comp_dirichlet = Kr @ (Id - Ainv_dirichlet @ Kr)  # Neumann, Hess, Eqn (24)
c_dirichlet = (Kr @ Ainv_dirichlet)[:, 0] * (-kx0)  # # Hess (25)
# print("C_comp_dirichlet", type(C_comp_dirichlet), C_comp_dirichlet.shape)
# print("c_dirichlet", type(c_dirichlet), c_dirichlet.shape)

C_comp_neumann = Kr @ (Id - Ainv_neumann @ Kr)  # Neumann, Hess, Eqn (32)
c_neumann = (Kr @ Ainv_neumann)[:, 0]  # Hess (33)
# print("C_comp_neumann", type(C_comp_neumann), C_comp_neumann.shape)
# print("c_neumann", type(c_neumann), c_neumann.shape)

print("invert matrix stop")

""" Numerical solution """
start_time = timeit.default_timer()
x_, y_, z_ = [], [], []
sink1d = []
sx = s.getSolutionHead()  # inital condition, solverbase.py

N = round(sim_time / dt)
t = 0.

rx = [0]

for i in range(0, N):

    t_pot = -trans * sinusoidal(t)  # potential transpiration ...
    print("t_pot", t_pot)

    hs = np.transpose(np.array([[sx[mapping[j]] for j in range(0, ns)]]))
    for j in range(0, len(nodes) - 1):  # from matric to total
        hs[j, 0] += nodes[j + 1][2]
    # print("hs", hs.shape, np.min(hs), np.max(hs))

    q_dirichlet0 = -(C_comp_dirichlet.dot(hs) + c_dirichlet * wilting_point)

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

    if np.sum(q_dirichlet0) > t_pot:
        print("dirichlet", np.sum(q_dirichlet0), t_pot)
        fluxes = r.sumSegFluxes(q_dirichlet0[:, 0])
    else:
        q_neumann0 = -(C_comp_neumann.dot(hs) - c_neumann * t_pot)
        print("neumann", np.sum(q_neumann0), t_pot)
        fluxes = r.sumSegFluxes(q_neumann0[:, 0])

    water = s.getWaterVolume()
    s.setSource(fluxes.copy())  # richards.py
    s.solve(dt)
    soil_water = (s.getWaterVolume() - water) / dt

    old_sx = sx.copy()
    sx = s.getSolutionHead()  # richards.py
    water = s.getWaterVolume()

    if  i % skip == 0:
        min_sx = np.min(sx)
        min_rx = np.min(rx)
        max_sx = np.max(sx)
        max_rx = np.max(rx)
        x_.append(t)
        sum_flux = 0.
        for f in fluxes.values():
            sum_flux += f
        y_.append(soil_water)  # cm3/day (soil uptake)
        z_.append(sum_flux)  # cm3/day (root system uptake)
        n = round(float(i) / float(N) * 100.)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days {:g}\n"
              .format(min_sx, max_sx, min_rx, max_rx, s.simTime, rx[0]))

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

