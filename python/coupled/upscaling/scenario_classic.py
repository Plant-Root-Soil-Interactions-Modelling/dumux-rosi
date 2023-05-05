""" 
New reference scenario for upscaling (DL 29.3.2023)

static root system in soil (1D or 3D) outer radii with Voronoi method, or via densities 

full hydraulic model with perirhizal nonlinear resistances 
(steady rate approach and fixed-point-iteration in new manuscript notation)
"""
import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import argparse
import timeit
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse

from functional.xylem_flux import sinusoidal2
import visualisation.vtk_plot as vp
from rhizo_models import plot_transpiration
from scenario_setup import *


def simulate_classic(sim_time, r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping):

    print("\nInitial root sytstem age", rs_age)
    # print("rs_age", rs_age)
    # rs_age = r.get_ages(rs_age)
    # print("rs_age", np.max(rs_age))

    dt = 360 / (24 * 3600)  # days
    skip = 1

    max_error = 10
    max_iter = 10

    ns = len(r.rs.segments)
    nodes = r.get_nodes()
    # seg_length = r.rs.segLength()
    assert len(nodes) - 1 == ns, "Number of nodes should be equal one less than number of segments"

    """ Doussan """
    A_dirichlet, Kr, kx0 = r.doussan_system_matrix(rs_age)
    Id = sparse.identity(ns).tocsc()  # identity matrix
    col_ind = r.collar_index()
    print("collar_index (segment index)", r.collar_index())
    print("kx0:", kx0)
    print()

    print("invert matrix start")
    Ainv_dirichlet = sparse.linalg.inv(A_dirichlet).todense()  # dense

    A_neumann = A_dirichlet.copy()
    A_neumann[col_ind, col_ind] -= kx0
    Ainv_neumann = sparse.linalg.inv(A_neumann).todense()  # dense

    C_comp_dirichlet = Kr @ (Id - Ainv_dirichlet @ Kr)  # Neumann, Hess, Eqn (24)
    c_dirichlet = (Kr @ Ainv_dirichlet)[:, col_ind] * (-kx0)  # # Hess (25)
    # print("C_comp_dirichlet", type(C_comp_dirichlet), C_comp_dirichlet.shape)
    # print("c_dirichlet", type(c_dirichlet), c_dirichlet.shape)

    C_comp_neumann = Kr @ (Id - Ainv_neumann @ Kr)  # Neumann, Hess, Eqn (32)
    c_neumann = (Kr @ Ainv_neumann)[:, col_ind]  # Hess (33)
    # print("C_comp_neumann", type(C_comp_neumann), C_comp_neumann.shape)
    # print("c_neumann", type(c_neumann), c_neumann.shape)
    print("invert matrix stop")

    """ Numerical solution """
    start_time = timeit.default_timer()
    x_, y_, z_, sink_, hs_, hx_, hsr_ = [], [], [], [], [], [], []
    sx = s.getSolutionHead()  # inital condition, solverbase.py

    N = round(sim_time / dt)
    t = 0.

    t_pot = -trans * sinusoidal2(0., dt)
    hs = np.transpose(np.array([[sx[mapping[j]][0] for j in range(0, ns)]]))
    for j in range(0, len(nodes) - 1):  # from matric to total
        hs[j, 0] += nodes[j + 1][2]
    rx = Ainv_dirichlet.dot(Kr.dot(hs)) + Ainv_dirichlet[:, 0] * kx0 * wilting_point
    q_dirichlet = -Kr.dot(hs - rx)
    if np.sum(q_dirichlet) < t_pot:
        rx = Ainv_neumann.dot(Kr.dot(hs)) + Ainv_neumann[:, 0] * t_pot  #   # Hess Eqn (29)

    for i in range(0, N):

        t_pot = -trans * sinusoidal2(t, dt)  # potential transpiration ...
        print("t_pot", t_pot)

        hs = np.transpose(np.array([[sx[mapping[j]][0] for j in range(0, ns)]]))
        for j in range(0, ns):  # from matric to total potential
            hs[j, 0] += nodes[j + 1][2]

        rx = Ainv_dirichlet.dot(Kr.dot(hs)) + Ainv_dirichlet[:, 0] * kx0 * wilting_point
        q_dirichlet = -Kr.dot(hs - rx)
        if np.sum(q_dirichlet) <= t_pot:
            rx = Ainv_neumann.dot(Kr.dot(hs)) + Ainv_neumann[:, 0] * t_pot  #   # Hess Eqn (29)

        if np.sum(q_dirichlet) > t_pot:
            print("dirichlet:", np.sum(q_dirichlet), t_pot)
            fluxes = r.sumSegFluxes(q_dirichlet[:, 0])
            sum_root_flux = np.sum(q_dirichlet)
        else:
            q_neumann = -Kr.dot(hs - rx)
            print("neumann:", np.sum(q_neumann), t_pot, np.sum(q_dirichlet))
            fluxes = r.sumSegFluxes(q_neumann[:, 0])
            sum_root_flux = np.sum(q_neumann)

        water = s.getWaterVolume()
        s.setSource(fluxes.copy())  # richards.py
        s.solve(dt)
        sum_soil_flux = (s.getWaterVolume() - water) / dt

        old_sx = sx.copy()
        sx = s.getSolutionHead()  # richards.py
        water = s.getWaterVolume()

        # remember results
        x_.append(t)
        y_.append(sum_soil_flux)  # cm3/day (soil uptake)
        z_.append(sum_root_flux)  # cm3/day (root system uptake)

        if  i % skip == 0:
            n = round(float(i) / float(N) * 100.)

            print("t_pot", t_pot, "q_root", sum_root_flux, "soil", sum_soil_flux)
            print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days {:g}\n"
                  .format(np.min(sx), np.max(sx), np.min(rx), np.max(rx), s.simTime, rx[0, 0]))

            """ remember results """
            sink = np.zeros(sx.shape)
            for k, v in fluxes.items():
                sink[k] += v
            sink_.append(sink)  # [cm3/day] per soil cell
            hs_.append(sx.copy())  # [cm] per soil cell
            hx_.append(rx)  # [cm] per segment

        t += dt

    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return hx_, hs_, sink_, x_, y_, z_, hs_, dt


if __name__ == "__main__":

    # TODO parse from args?
    plant = "soybean"  # soybean, maize
    dim = "3D"  # 1D, 3D
    soil = "hydrus_loam"  #  hydrus_loam, hydrus_clay, hydrus_sand
    outer_method = "surface"  # voronoi, length, surface, volume
    initial = -200  # cm total potential

    sim_time = 7.
    name = "classic_" + plant + "_" + dim + "_" + soil + "_" + outer_method
    print(name, "\n")

    r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = set_scenario(plant, dim, initial, soil, outer_method)
    print("trans", trans)

    hx_, hsr_, sink_, x_, y_, z_, hs_, dt = simulate_classic(sim_time, r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping)

    """ write and plot """
    s.writeDumuxVTK("results/" + name)
    write_files(name, hx_, hsr_, sink_, x_, y_, z_, hs_)
    plot_transpiration(x_, y_, z_, lambda t: trans * sinusoidal2(t, dt), name)
