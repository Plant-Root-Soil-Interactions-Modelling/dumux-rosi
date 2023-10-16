""" 
Benchmark C1.2 for a static root system in soil (1D or 3D)

with the steady rate approach and fixed-point-iteration in HESS paper notation
"""
import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import timeit
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse

from rhizo_models import plot_transpiration
from scenario_setup import *
import aggregated_rs as agg


def double_(rsx):
    rsx2 = np.array([ 0. if i % 2 == 0 else rsx[int(i / 2)] for i in range(0, ns)])
    return rsx2

# r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("3D", age_dependent = False)
# name = "results/c12_sra"

# r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("3D", age_dependent = True)
# name = "results/c12b_sra"


r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("1D")
min_b, max_b, cell_number = get_domain1D()
name = "results/c12_agg_1d"

""" 
Initialize aggregated hydraulic model 
"""
r = agg.create_aggregated_rs(r, rs_age, min_b, max_b, cell_number)
nodes = r.get_nodes()
picker = lambda x, y, z: s.pick([0., 0., z])  # reset mapper, since it is 1D
r.rs.setSoilGrid(picker)  # maps segment

max_error = 10.
max_iter = 1000

""" for fixed mapping """
types = r.rs.subTypes
outer_r = r.rs.segOuterRadii()
inner_r = r.rs.radii
print("inner_r", np.min(inner_r), np.max(outer_r))

rho_ = np.divide(outer_r, np.array(inner_r))
ns = len(outer_r)
mapping = np.array([r.rs.seg2cell[j] for j in range(0, ns)])  # redo mapping

kr_ = np.array(r.getEffKr(rs_age))  # times surface (2 a pi length)
inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up

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

""" Numerical solution """
start_time = timeit.default_timer()
rs_age = np.max(r.get_ages())
x_, y_, z_ = [], [], []
sink1d = []
sx = s.getSolutionHead()  # inital condition, solverbase.py

N = round(sim_time / dt)
t = 0.

t_pot = -trans * sinusoidal(0)
hs = np.transpose(np.array([[sx[mapping[j]][0] for j in range(0, ns)]]))
for j in range(0, len(nodes) - 1):  # from matric to total
    hs[j, 0] += nodes[j + 1][2]
rx = Ainv_dirichlet.dot(Kr.dot(hs)) + Ainv_dirichlet[:, 0] * kx0 * wilting_point
q_dirichlet = -Kr.dot(hs - rx)
if np.sum(q_dirichlet) < t_pot:
    rx = Ainv_neumann.dot(Kr.dot(hs)) + Ainv_neumann[:, 0] * t_pot  #   # Hess Eqn (29)

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
        # print(rx.shape, type(rx), hs.shape, type(hs), inner_kr_.shape, type(inner_kr_), rho_.shape, type(rho_))
        rsx = soil_root_interface_table(rx, hs, inner_kr_, rho_, sra_table_lookup)

        for j in range(0, len(nodes) - 1):  # from matric to total
            rsx[j, 0] += nodes[j + 1][2]
        # print("rsx", rsx.shape, rsx[0], rx[0], hs[0])
        wall_interpolation = timeit.default_timer() - wall_interpolation

        """ xylem matric potential """
        wall_xylem = timeit.default_timer()

        rx = Ainv_dirichlet.dot(Kr.dot(rsx)) + Ainv_dirichlet[:, 0] * kx_[0] * wilting_point
        q_dirichlet = -Kr.dot(rsx - rx)
        if np.sum(q_dirichlet) <= t_pot:
            rx = Ainv_neumann.dot(Kr.dot(rsx)) + Ainv_neumann[:, 0] * t_pot  #   # Hess Eqn (29)

        err = np.linalg.norm(rx - rx_old)
        wall_xylem = timeit.default_timer() - wall_xylem
        rx_old = rx.copy()
        c += 1

    if np.sum(q_dirichlet) > t_pot:
        print("dirichlet", c, err, np.sum(q_dirichlet), t_pot)
        fluxes = r.sumSegFluxes(q_dirichlet[:, 0])
    else:
        q_neumann = -Kr.dot(rsx - rx)
        print("neumann", c, err, np.sum(q_neumann), t_pot)
        fluxes = r.sumSegFluxes(q_neumann[:, 0])

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
              .format(min_sx, max_sx, min_rx, max_rx, s.simTime, rx[0, 0]))

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

