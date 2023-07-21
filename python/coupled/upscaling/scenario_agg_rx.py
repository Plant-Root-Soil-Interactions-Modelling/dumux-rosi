""" 
static root system in soil (1D or 3D) outer radii with Voronoi method, or via densities 

aggregated hydraulic model with aggregated perirhizal nonlinear resistances 
(aggregated  over soil cells, steady rate approach and fixed-point-iteration on aggregated values, 
in new manuscript notation)
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


def simulate_agg(sim_time, r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping):

    print("\nInitial root sytstem age (agg)", rs_age)
    # print("rs_age", rs_age)
    # rs_age = r.get_ages(rs_age)
    # print("rs_age", np.max(rs_age))

    dt = 360 / (24 * 3600)  # days
    skip = 10

    max_error = 10
    max_iter = 100

    ns = len(r.rs.segments)
    nodes = r.get_nodes()
    seg_length = r.rs.segLength()
    assert len(nodes) - 1 == ns, "number of nodes should be equal one less than number of segments"

    """ Fetch rhizosphere model params """
    kr_ = np.array(r.getKr(rs_age))
    inner_ = r.rs.radii
    inner_kr_ = np.multiply(inner_, kr_)  # multiply for table look up
    inner_kr_ = np.expand_dims(inner_kr_, axis = 1)
    # inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)
    inner_kr_ = np.minimum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-4)
    rho_ = np.maximum(rho_, np.ones(rho_.shape) * 1.)
    rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200.)

    """ Doussan """
    A_dirichlet, Kr, kx0 = r.doussan_system_matrix(rs_age)
    Id = sparse.identity(ns).tocsc()  # identity matrix

    print("inv start")
    Ainv_dirichlet = sparse.linalg.inv(A_dirichlet).todense()  # dense
    A_neumann = A_dirichlet.copy()
    A_neumann[0, 0] -= kx0
    Ainv_neumann = sparse.linalg.inv(A_neumann).todense()  # dense
    print("inv stop")

    print("up start")
    B, soil2matrix, matrix2soil = r.get_soil_matrix()
    nmax = len(matrix2soil)

    Bt = B.transpose()
    BBt_inv = sparse.linalg.inv(B @ Bt)  # sparse

    C_comp_dirichlet = Kr @ (Id - Ainv_dirichlet @ Kr)  # Neumann, Hess, Eqn (24) ->[25]
    c_dirichlet = (Kr @ Ainv_dirichlet)[:, 0] * (-kx0)  # # Hess (25) -> [26]
    C_comp_neumann = Kr @ (Id - Ainv_neumann @ Kr)  # Neumann, Hess, Eqn (32)
    c_neumann = (Kr @ Ainv_neumann)[:, 0]  # Hess (33)

    AinvKr_dirichlet_up = (((B @ Ainv_dirichlet) @ Kr) @ Bt)
    Ainv_dirichlet_up = B @ Ainv_dirichlet
    C_comp_dirichlet_up = B @ C_comp_dirichlet @ Bt
    c_dirichlet_up = B @ c_dirichlet
    # print(C_comp_neumann_up.shape, type(C_comp_neumann_up))

    AinvKr_neumann_up = B @ ((Ainv_neumann) @ Kr) @ Bt
    Ainv_neumann_up = B @ Ainv_neumann
    C_comp_neumann_up = B @ C_comp_neumann @ Bt
    c_neumann_up = B @ c_neumann

    Kr_up = B @ Kr @ Bt  # sparse
    # Kr_up_inv = sparse.linalg.inv(Kr_up)

    inner_kr_up = BBt_inv.dot(B.dot(inner_kr_))
    inner_kr_up = np.maximum(inner_kr_up, np.ones(inner_kr_up.shape) * 1.e-7)  ############################################ (too keep within table)
    inner_kr_up = np.minimum(inner_kr_up, np.ones(inner_kr_up.shape) * 1.e-4)  ############################################ (too keep within table)

    rho_up = BBt_inv.dot(B.dot(rho_))

    print("up end")

    """ Numerical solution (a) """
    start_time = timeit.default_timer()
    x_, y_, z_, sink_, sx_, hx_, hsr_ = [], [], [], [], [], [], []
    sx = s.getSolutionHead()  # inital condition, solverbase.py

    N = round(sim_time / dt)
    t = 0.
    rx = [0]

    centers = s.getCellCenters()

    t_pot = -trans * sinusoidal2(t, dt)  # potential transpiration ...
    hs_ = np.zeros((nmax, 1))  # sx -> hs_ # soil cell indices to soil matrix indices
    for j in soil2matrix.keys():
            hs_[soil2matrix[j]] += sx[j] + centers[j, 2]  # from matric to total

    hxd = BBt_inv.dot(AinvKr_dirichlet_up.dot(hs_) + Ainv_dirichlet_up[:, 0] * kx0 * wilting_point)
    q_dirichlet_up = -Kr_up.dot(hs_ - hxd)
    rx = hxd
    # if np.sum(q_dirichlet_up) > t_pot:
    #     rx = hxd
    # else:
    #     rx = BBt_inv.dot(AinvKr_neumann_up.dot(hs_) + Ainv_neumann_up[:, 0] * t_pot)

    for i in range(0, N):

        t_pot = -trans * sinusoidal2(t, dt)  # potential transpiration ..

        hs_ = np.zeros((nmax, 1))  # sx -> hs_ # soil cell indices to soil matrix indices
        for j in soil2matrix.keys():
                hs_[soil2matrix[j]] += sx[j] + centers[j, 2]  # from matric to total

        wall_iteration = timeit.default_timer()
        err = 1.e6
        c = 1
        rx_old = rx
        while err > max_error and c < max_iter:

            # print("rx", np.min(rx), np.max(rx))

            """ interpolation """
            wall_interpolation = timeit.default_timer()

            # print(rx.shape, type(rx), hs_.shape, type(hs_), inner_kr_.shape, type(inner_kr_), rho_.shape, type(rho_))
            for j in soil2matrix.keys():  # from total to matric
                    rx[soil2matrix[j]] -= centers[j, 2]

            rsx = soil_root_interface_table(rx, hs_, inner_kr_up, rho_up, sra_table_lookup)  # in matric potential

            for j in soil2matrix.keys():  # from matric to total
                    rsx[soil2matrix[j]] += centers[j, 2]

            wall_interpolation = timeit.default_timer() - wall_interpolation

            """ xylem matric potential """
            wall_xylem = timeit.default_timer()

            hxd = BBt_inv.dot(AinvKr_dirichlet_up.dot(rsx) + Ainv_dirichlet_up[:, 0] * kx0 * wilting_point)
            q_dirichlet_up = -Kr_up.dot(rsx - hxd)
            rx = hxd
            # if np.sum(q_dirichlet_up) > t_pot:
            #
            # else:
            #     rx = BBt_inv.dot(AinvKr_neumann_up.dot(rsx) + Ainv_neumann_up[:, 0] * t_pot)

            err = np.linalg.norm(rx - rx_old)
            wall_xylem = timeit.default_timer() - wall_xylem
            rx_old = rx.copy()
            c += 1

        wall_soil = timeit.default_timer()

        fluxes = {}
        for j in range(0, nmax):
            fluxes[matrix2soil[j]] = q_dirichlet_up[j, 0]
        sum_root_flux = np.sum(q_dirichlet_up)

        # if np.sum(q_dirichlet_up) > t_pot:
        #     print("q_dirichlet_up", c, err, np.sum(q_dirichlet_up), t_pot)
        #     fluxes = {}
        #     for j in range(0, nmax):
        #         fluxes[matrix2soil[j]] = q_dirichlet_up[j, 0]
        #     sum_root_flux = np.sum(q_dirichlet_up)
        # else:
        #     # print("C_comp_neumann_up", C_comp_neumann_up.shape)
        #     # print("c_neumann_up", c_neumann_up.shape)
        #     # print("hs_", hs_.shape)
        #     q_neumann_up = -Kr_up.dot(rsx - rx)
        #     print("q_neumann0_up", c, err, q_neumann_up.shape, np.sum(q_neumann_up), t_pot)
        #     fluxes = {}
        #     for j in range(0, nmax):
        #         fluxes[matrix2soil[j]] = q_neumann_up[j, 0]
        #     sum_root_flux = np.sum(q_neumann_up)

        water = s.getWaterVolume()
        s.setSource(fluxes.copy())  # richards.py
        s.solve(dt)
        sum_soil_flux = (s.getWaterVolume() - water) / dt

        old_sx = sx.copy()
        sx = s.getSolutionHead()  # richards.py

        wall_iteration = timeit.default_timer() - wall_iteration

        wall_soil = timeit.default_timer() - wall_soil

        x_.append(t)
        y_.append(sum_soil_flux)  # cm3/day (soil uptake)
        z_.append(sum_root_flux)  # cm3/day (root system uptake)

        if  i % skip == 0:
            n = round(float(i) / float(N) * 100.)

            # print("rx", rx.shape)
            # print("rsx", rsx.shape)
            # print("nodes", nodes.shape)
            # rx_ = np.zeros((ns, 1))
            # rsx_ = np.zeros((ns, 1))
            # for j in range(0, ns):  # from total to matric
            #     #   print(nodes[j + 1][2])
            #     rx_[j, 0] = rx[j, 0] - nodes[j + 1][2]
            #     rsx_[j, 0] = rsx[j, 0] - nodes[j + 1][2]

            print("number of iterations", c)
            print("t_pot", t_pot, "q_root", sum_root_flux, "soil", sum_soil_flux)
            print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days {:g}\n"
                  .format(np.min(sx), np.max(sx), np.min(rx), np.max(rx), s.simTime, rx[0, 0]))

            print("iteration (interpolation, xylem) : ", wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem))
            print("iteration, soil", wall_iteration / (wall_iteration + wall_soil), wall_soil / (wall_iteration + wall_soil))
            print()

            """ remember results """
            sink = np.zeros(sx.shape)
            for k, v in fluxes.items():
                sink[k] += v
            sink_.append(sink)  # [cm3/day] per soil cell
            sx_.append(sx.copy())  # [cm] per soil cell
            # hx_.append(rx_)  # [cm] total potential per segment
            # hsr_.append(rsx_)  # [cm] total potential per segment

        t += dt

    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return hx_, hsr_, sink_, x_, y_, z_, sx_, dt


def run_agg(sim_time, method, plant, dim, soil, outer_method):

    # hidden parameters...
    initial = -200  # cm

    name = "agg_" + plant + "_" + dim + "_" + soil + "_" + outer_method
    print(name, "\n")

    r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = set_scenario(plant, dim, initial, soil, outer_method)

    hx_, hsr_, sink_, x_, y_, z_, hs_, dt = simulate_agg(sim_time, r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping)

    s.writeDumuxVTK("results/" + name)  # final soil VTU
    write_files(name, hx_, hsr_, sink_, x_, y_, z_, hs_)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = 'Simulation options')
    parser.add_argument('plant', type = str, help = 'soybean or maize or springbarley')
    parser.add_argument('dim', type = str, help = '1D or 3D')
    parser.add_argument('soil', type = str, help = 'soil type (hydrus_loam, hydrus_clay, hydrus_sand or hydrus_sandyloam)')
    parser.add_argument('outer_method', type = str, help = 'how to determine outer radius (voronoi, length, surface, volume)')

    args = parser.parse_args(['springbarley', "1D", "hydrus_loam", "surface"])
    # args = parser.parse_args()

    name = "_agg_" + args.plant + "_" + args.dim + "_" + args.soil + "_" + args.outer_method

    initial = -200  # cm     plot_transpiration(x_, y_, z_, lambda t: trans * sinusoidal2(t, dt))
    sim_time = 7.5

    print(name, "\n")

    r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = set_scenario(args.plant, args.dim, initial, args.soil, args.outer_method)

    print()
    r.rs.write(name + ".vtp")
    print()

    hx_, hsr_, sink_, x_, y_, z_, hs_, dt = simulate_agg(sim_time, r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping)

    plot_transpiration(x_, y_, z_, lambda t: trans * sinusoidal2(t, dt))

    """ write """
    s.writeDumuxVTK("results/" + name)
    write_files(name, hx_, hsr_, sink_, x_, y_, z_, hs_)

#
