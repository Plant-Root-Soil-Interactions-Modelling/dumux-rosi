""" 
New reference scenario for upscaling (DL 29.3.2023)

static root system in soil (1D or 3D) outer radii with Voronoi method, or via densities 

full hydraulic model with perirhizal nonlinear resistances 
(steady rate approach and fixed-point-iteration in new manuscript notation)

solved by 
*) inversion dense matrix (slowest)
*) spsolve
*) LU factorisation (fastest)
"""
import sys; sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");

import argparse
import timeit
import numpy as np
# import matplotlib.pyplot as plt
from scipy import sparse
import scipy.sparse.linalg as LA
import vtk

from functional.xylem_flux import sinusoidal2
import visualisation.vtk_plot as vp

# from rhizo_models import plot_transpiration
from scenario_setup import *


def simulate_sra(sim_time, r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping):

    print("\nInitial root sytstem age (sra)", rs_age)
    sys.stdout.flush()
    # print("rs_age", rs_age)
    # rs_age = r.get_ages(rs_age)
    # print("rs_age", np.max(rs_age))

    dt = 360 / (24 * 3600)  # days
    skip = 10

    max_error = 10
    max_iter = 100

    ns = len(r.rs.segments)
    nodes = r.get_nodes()
    # seg_length = r.rs.segLength()
    assert len(nodes) - 1 == ns, "Number of nodes should be equal one less than number of segments"

    """fetch rhizosphere model params """
    kr_ = np.array(r.get_kr(rs_age))
    kr_ = np.maximum(kr_, np.ones(kr_.shape) * 1.e-4)
    inner_ = r.rs.radii
    inner_kr_ = np.multiply(inner_, kr_)  # multiply for table look up
    inner_kr_ = np.expand_dims(inner_kr_, axis = 1)
    inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)
    inner_kr_ = np.minimum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-4)
    rho_ = np.maximum(rho_, np.ones(rho_.shape) * 1.)
    rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200.)

    """ initialize hydraulic model """
    r.update(rs_age)

    """ numerical solution """
    start_time = timeit.default_timer()
    x_, y_, z_, sink_, hs_, hx_, hsr_ = [], [], [], [], [], [], []
    sx = s.getSolutionHead_()  # inital condition, solverbase.py

    N = round(sim_time / dt)
    t = 0.

    t_pot = -trans * sinusoidal2(0., dt)

    hs = np.transpose(np.array([[sx[mapping[j]] for j in range(0, ns)]]))  # TODO use getter
    for j in range(0, len(nodes) - 1):  # from matric to total
        hs[j, 0] += nodes[j + 1][2]

    rx = r.solve(hs, 0., wilting_point)

    for i in range(0, N):

        t_pot = -trans * sinusoidal2(t, dt)  # potential transpiration ...

        # hs = np.transpose(np.array([[sx[mapping[j]] for j in range(0, ns)]]))
        hs = np.expand_dims(r.get_hs(sx), axis = 1)

        wall_iteration = timeit.default_timer()
        err = 1.e6
        c = 0
        rx_old = rx  # in total potential
        while err > max_error and c < max_iter:  # max_iter = 1000 , err>100 works fine

            """ interpolation """
            wall_interpolation = timeit.default_timer()

            for j in range(0, ns):  # from total to matric
                rx[j, 0] -= nodes[j + 1][2]

            rx = np.maximum(rx, np.ones(rx.shape) * (-15999))
            rx = np.minimum(rx, np.ones(rho_.shape) * 0.)
            hs = np.maximum(hs, np.ones(hs.shape) * (-15999))

            # print("*")
            # print(rx.shape, hs.shape, ns)
            rsx = soil_root_interface_table(rx, hs, inner_kr_, rho_, sra_table_lookup)
            # print("*")
            for j in range(0, ns):  # from matric to total
                rsx[j, 0] += nodes[j + 1][2]

            wall_interpolation = timeit.default_timer() - wall_interpolation

            """ xylem matric potential """
            wall_xylem = timeit.default_timer()

            # print("types...", type(rsx), type(Kr), type(Ainv_dirichlet))
            # print("shapes...", rsx.shape, Kr.shape, Ainv_dirichlet.shape)
            # rx = Ainv_dirichlet.dot(Kr.dot(rsx)) + Ainv_dirichlet[:, collar_index] * kx0 * wilting_point
            # b = Kr.dot(rsx)
            # b[collar_index, 0] += kx0 * wilting_point
            # rx = sparse.linalg.spsolve(A_d, b)
            # rx = np.expand_dims(rx, axis = 1)
            rx = r.solve(rsx, t_pot, wilting_point)
            # print("dirichlet rx", rx.shape, np.min(rx), np.max(rx))
            # print("dirichlet rx2", rx2.shape, np.min(rx2), np.max(rx2))

            wall_xylem = timeit.default_timer() - wall_xylem

            q = r.radial_fluxes(rx, rsx)  # -Kr.dot(rsx - rx)  # both total potentials

            err = np.linalg.norm(rx - rx_old)
            wall_xylem = timeit.default_timer() - wall_xylem
            rx_old = rx.copy()
            c += 1

        wall_iteration = timeit.default_timer() - wall_iteration

        wall_soil = timeit.default_timer()

        fluxes = r.sumSegFluxes(q[:, 0])
        sum_root_flux = np.sum(q)

        water = s.getWaterVolume()
        s.setSource(fluxes.copy())  # richards.py
        s.solve(dt)
        sum_soil_flux = (s.getWaterVolume() - water) / dt
        old_sx = sx.copy()
        sx = s.getSolutionHead_()  # richards.py

        wall_soil = timeit.default_timer() - wall_soil

        x_.append(t)
        y_.append(sum_soil_flux)  # cm3/day (soil uptake)
        z_.append(sum_root_flux)  # cm3/day (root system uptake)

        if  i % skip == 0:

            n = round(float(i) / float(N) * 100.)

            rx_ = np.zeros((ns, 1))
            rsx_ = np.zeros((ns, 1))
            for j in range(0, ns):  # from total to matric
                #   print(nodes[j + 1][2])
                rx_[j, 0] = rx[j, 0] - nodes[j + 1][2]
                rsx_[j, 0] = rsx[j, 0] - nodes[j + 1][2]

            print("number of iterations", c)
            print("t_pot", t_pot, "q_root", sum_root_flux, "soil", sum_soil_flux)

            print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days {:g}\n"
                  .format(np.min(sx), np.max(sx), np.min(rx_), np.max(rx_), s.simTime, rx[r.ci, 0]))
            print("iteration (interpolation, xylem) : ", wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem))
            print("iteration, soil", wall_iteration / (wall_iteration + wall_soil), wall_soil / (wall_iteration + wall_soil))
            print()
            sys.stdout.flush()

            """ remember results """
            sink = np.zeros(sx.shape)
            for k, v in fluxes.items():
                sink[k] += v
            sink_.append(sink)  # [cm3/day] per soil cell
            hs_.append(sx.copy())  # [cm] per soil cell
            hx_.append(rx_)  # [cm] total potential per segment
            hsr_.append(rsx_)  # [cm] total potential per segment

        t += dt

    wall_time = timeit.default_timer() - start_time
    print ("Coupled benchmark solved in ", wall_time, " s")
    sys.stdout.flush()

    return hx_, hsr_, sink_, x_, y_, z_, hs_, dt, wall_time


def run_sra(sim_time, method, plant, dim, soil, outer_method):

    # hidden parameters...
    initial = -200  # cm

    name = "sra_" + plant + "_" + dim + "_" + soil + "_" + outer_method
    print(name, "\n")

    r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = set_scenario(plant, dim, initial, soil, outer_method)

    hx_, hsr_, sink_, x_, y_, z_, hs_, dt, wall_time = simulate_sra(sim_time, r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping)

    s.writeDumuxVTK("results/" + name)  # final soil VTU
    write_files(name, hx_, hsr_, sink_, x_, y_, z_, hs_, wall_time)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = 'Simulation options')
    parser.add_argument('plant', type = str, help = 'soybean or maize or springbarley')
    parser.add_argument('dim', type = str, help = '1D or 3D')
    parser.add_argument('soil', type = str, help = 'soil type (hydrus_loam, hydrus_clay, hydrus_sand or hydrus_sandyloam)')
    parser.add_argument('outer_method', type = str, help = 'how to determine outer radius (voronoi, length, surface, volume)')

    args = parser.parse_args(['maize', "2D", "hydrus_loam", "length"])
    # args = parser.parse_args()

    name = "sra_" + args.plant + "_" + args.dim + "_" + args.soil + "_" + args.outer_method
    print()
    print(name, "\n")

    initial = -200  # cm     plot_transpiration(x_, y_, z_, lambda t: trans * sinusoidal2(t, dt))
    sim_time = 14.5

    r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = set_scenario(args.plant, args.dim, initial, args.soil, args.outer_method)

    hx_, hsr_, sink_, x_, y_, z_, hs_, dt, wall_time = simulate_sra(sim_time, r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping)

    """ write """
    s.writeDumuxVTK("results/" + name)
    write_files(name, hx_, hsr_, sink_, x_, y_, z_, hs_, wall_time)

    # """ plot results """
    # plot_transpiration(x_, y_, z_, lambda t: trans * sinusoidal2(t, dt), name)
    #
    # hx = hx_[-1]
    # sx = hs_[-1]
    # ana = pb.SegmentAnalyser(r.rs.mappedSegments())
    # ana.addData("pressure head", hx)
    # # vp.plot_roots(ana, "hx")
    #
    # ns = len(r.rs.segments)
    # hs = np.transpose(np.array([[sx[mapping[j]] for j in range(0, ns)]]))
    # # print("hx", hx.shape)
    # # print("sx", sx.shape)
    # # print("hs", hs.shape)
    # # print(len(r.rs.subTypes))
    # # print(len(r.rs.organTypes))
    # # print(len(r.rs.nodeCTs))
    #
    # ana.addConductivities(r, rs_age)
    # ana.addFluxes(r, hx, hs, rs_age)  # adds "axial_flux", "radial_flux"
    # if args.plant == "maize":
    #     min_b, max_b, cell_number = maize_(args.dim)
    # elif args.plant == "soybean":
    #     min_b, max_b, cell_number = soybean_(args.dim)
    #
    # # print(hx.shape)
    # vp.plot_roots_and_soil(ana, "radial_flux", None, s, True, min_b, max_b, cell_number)
    # # vp.plot_roots_and_soil(ana, "axial_flux", None, s, True, min_b, max_b, cell_number)  # NOT WORKING .... (?)
    # vp.plot_roots_and_soil(ana, "pressure head", None, s, True, min_b, max_b, cell_number)
