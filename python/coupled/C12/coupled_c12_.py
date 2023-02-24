""" 
Benchmark C1.2 for a static root system in soil (1D or 3D) 
with the classic sink using Doussan approach in HESS (experimental) paper notation with no matrix inversion 
"""
import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import timeit
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import sparse

# import vtk_plot as vp
from rhizo_models import plot_transpiration
from scenario_setup import *

r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("3D")
name = "results/c12_"

# r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("3D", age_dependent = True)
# dt = 36 / (24 * 3600)  # unstable
# name = "results/c12b"

# r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("1D")
# name = "results/c12_1d"

nodes = r.get_nodes()
ns = len(r.rs.segments)

""" Doussan """
Id = sparse.identity(ns).tocsc()  # identity matrix
A, Kr, Kx = r.doussan_system_matrix(rs_age)
kx_ = np.divide(r.getKx(rs_age), r.rs.segLength())  # / dl
A_n = A.copy()
A_n[0, 0] -= kx_[0]
Kr_inv = sparse.linalg.inv(Kr)

""" see HESS (experimental) """
A_q = A @ Kr_inv
A_nq = A_n @ Kr_inv
B = A - Kr
Bn = A_n - Kr
# print("A_q", A_q.shape, "A_nq", A_nq.shape, "B", B.shape)

""" Numerical solution """
start_time = timeit.default_timer()
x_, y_, w_, cf = [], [], [], []
sink1d = []
sx = s.getSolutionHead()  # inital condition, solverbase.py

N = round(sim_time / dt)
t = 0.

rx = [0]

for i in range(0, N):

    t_pot = -trans * sinusoidal(t)  # potential transpiration ...
    print("t_pot", t_pot)

    hs = np.transpose(np.array([[sx[mapping[j]][0] for j in range(0, ns)]]))
    for j in range(0, len(nodes) - 1):  # from matric to total
        hs[j, 0] += nodes[j + 1][2]

    b = B.dot(hs)
    b[0, 0] -= kx_[0] * wilting_point
    q_dirichlet = -sparse.linalg.spsolve(A_q, b)
    print("q_dirichlet", q_dirichlet.shape, np.sum(q_dirichlet))

    if np.sum(q_dirichlet) > t_pot:
        print("dirichlet", np.sum(q_dirichlet), t_pot)
        fluxes = r.sumSegFluxes(q_dirichlet)
    else:
        b = Bn.dot(hs)
        b[0, 0] -= t_pot
        q_neumann = -sparse.linalg.spsolve(A_nq, b)
        print("neumann", q_neumann.shape, np.sum(q_neumann), t_pot)
        fluxes = r.sumSegFluxes(q_neumann)

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
        cf.append(sum_flux)  # cm3/day (root system uptake)
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
plot_transpiration(x_, y_, cf, lambda t: trans * sinusoidal(t))
np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter = ';')

# vp.plot_roots_and_soil(r.rs, "pressure head", rx, s, False, min_b, max_b, cell_number, name)  # VTK vizualisation
# sink1d = np.array(sink1d)
# np.save("sink1d", sink1d)
# print(sink1d.shape)
# print(sink1d[-1])

