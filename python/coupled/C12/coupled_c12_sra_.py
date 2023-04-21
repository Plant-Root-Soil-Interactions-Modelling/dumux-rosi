""" 
Benchmark C1.2 for a static root system in soil (1D or 3D)

with rhizosphere model using steady rate approach and fixed-point-iteration in HESS paper notation

!without matrix inversion! 
"""
import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import timeit
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse

import visualisation.vtk_plot as vp
from rhizo_models import plot_transpiration
from scenario_setup import *

r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("3D", age_dependent = False)
name = "results/c12_sra"

# r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("3D", age_dependent = True)
# name = "results/c12b_sra"

# r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("1D")
# name = "results/c12_sra_1d"

max_error = 10.
max_iter = 100

ns = len(r.rs.segments)
nodes = r.get_nodes()
inner_ = r.rs.radii
outer_ = r.rs.segOuterRadii()
rho_ = np.divide(outer_, np.array(inner_))
rho_ = np.expand_dims(rho_, axis = 1)
kr_ = r.getKr(rs_age)
inner_kr_ = np.multiply(inner_, kr_)  # multiply for table look up
inner_kr_ = np.expand_dims(inner_kr_, axis = 1)
inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)  ############################################ (too keep within table)
inner_kr_ = np.minimum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-4)  ############################################ (too keep within table)

""" Doussan """
A_d, Kr, kx0 = r.doussan_system_matrix(rs_age)
Id = sparse.identity(ns).tocsc()  # identity matrix

A_n = A_d.copy()
A_n[0, 0] -= kx0
Kr_inv = sparse.linalg.inv(Kr)  # Kr is a diagonal matrix, thus Kr_inv sparse

A_dq = A_d @ Kr_inv
A_nq = A_n @ Kr_inv
Bd = A_d - Kr
Bn = A_n - Kr

""" Numerical solution """
start_time = timeit.default_timer()
rs_age = np.max(r.get_ages())
x_, y_, z_ = [], [], []
sink1d = []
sx = s.getSolutionHead()  # inital condition, solverbase.py

N = round(sim_time / dt)
t = 0.

t_pot = -trans * sinusoidal(0.)
hs = np.transpose(np.array([[sx[mapping[j]][0] for j in range(0, ns)]]))
for j in range(0, len(nodes) - 1):  # from matric to total
    hs[j, 0] += nodes[j + 1][2]

# code can be reduced, always neumann for t_pot = 0, i.e. L79,81, & 82
b = Bd.dot(hs)
b[0, 0] -= kx0 * wilting_point
q_dirichlet = -sparse.linalg.spsolve(A_dq, b)
if np.sum(q_dirichlet) > t_pot:
    rx = hs - Kr_inv.dot(-q_dirichlet[:, np.newaxis])
else:
    b = Bn.dot(hs)
    b[0, 0] -= t_pot
    q_neumann = -sparse.linalg.spsolve(A_nq, b)
    rx = hs - Kr_inv.dot(-q_neumann[:, np.newaxis])

for i in range(0, N):

    t_pot = -trans * sinusoidal(t)  # potential transpiration ...
    print("t_pot", t_pot)

    hs = np.transpose(np.array([[sx[mapping[j]][0] for j in range(0, ns)]]))
    for j in range(0, len(nodes) - 1):  # from matric to total
        hs[j, 0] += nodes[j + 1][2]

    err = 1.e6
    c = 1
    rx_old = rx  # in total potential
    while err > max_error and c < max_iter:  # max_iter = 1000 , err>100 works fine

        """ interpolation """
        wall_interpolation = timeit.default_timer()

        for j in range(0, len(nodes) - 1):  # from total to matric
            rx[j, 0] -= nodes[j + 1][2]

        # print("rx", rx.shape, type(rx), "hs", hs.shape, type(hs), "inner_kr", inner_kr_.shape, type(inner_kr_), "rho_", rho_.shape, type(rho_))
        # print(c, np.min(rx))
        rsx = soil_root_interface_table(rx, hs, inner_kr_, rho_, sra_table_lookup)

        for j in range(0, len(nodes) - 1):  # from matric to total
            rsx[j, 0] += nodes[j + 1][2]
        # print("rsx", rsx.shape, rsx[0], rx[0], hs[0])
        wall_interpolation = timeit.default_timer() - wall_interpolation

        """ xylem matric potential """
        wall_xylem = timeit.default_timer()

        b = Bd.dot(rsx)
        b[0, 0] -= kx0 * wilting_point
        q_dirichlet = -sparse.linalg.spsolve(A_dq, b)
        if np.sum(q_dirichlet) > t_pot:
            rx = rsx - Kr_inv.dot(-q_dirichlet[:, np.newaxis])
        else:
            b = Bn.dot(rsx)
            b[0, 0] -= t_pot
            q_neumann = -sparse.linalg.spsolve(A_nq, b)
            rx = rsx - Kr_inv.dot(-q_neumann[:, np.newaxis])

        err = np.linalg.norm(rx - rx_old)
        wall_xylem = timeit.default_timer() - wall_xylem
        rx_old = rx.copy()
        c += 1

    if np.sum(q_dirichlet) > t_pot:
        print("dirichlet", c, err, np.sum(q_dirichlet), t_pot)
        fluxes = r.sumSegFluxes(q_dirichlet)
    else:
        print("neumann", c, err, np.sum(q_neumann), t_pot)
        fluxes = r.sumSegFluxes(q_neumann)

    water = s.getWaterVolume()
    s.setSource(fluxes.copy())  # richards.py
    s.solve(dt)
    soil_water = (s.getWaterVolume() - water) / dt

    old_sx = sx.copy()
    sx = s.getSolutionHead()  # richards.py
    water = s.getWaterVolume()

    if  i % skip == 0:
        x_.append(t)
        sum_flux = 0.
        for f in fluxes.values():
            sum_flux += f
        y_.append(soil_water)  # cm3/day (soil uptake)
        z_.append(sum_flux)  # cm3/day (root system uptake)
        n = round(float(i) / float(N) * 100.)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days {:g}\n"
              .format(np.min(sx), np.max(sx), np.min(rx), np.max(rx), s.simTime, rx[0, 0]))

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

