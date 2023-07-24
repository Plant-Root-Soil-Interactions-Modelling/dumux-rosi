""" 
New reference scenario for upscaling (DL 29.3.2023)

static root system in soil (1D or 3D) outer radii with Voronoi method, or via densities 

full hydraulic model with perirhizal nonlinear resistances 

solved by 
Meunier matrix (sovled with spsolve, but LU factorisation also possible, in Sensitivity branch)
"""
import sys; sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");

import argparse
import timeit
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
import vtk

from functional.xylem_flux import sinusoidal2
import visualisation.vtk_plot as vp

from rhizo_models import plot_transpiration
from scenario_setup import *


def simulate_sraOld(sim_time, r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping):

    print("\nInitial root sytstem age", rs_age)
    # print("rs_age", rs_age)
    # rs_age = r.get_ages(rs_age)
    # print("rs_age", np.max(rs_age))

    dt = 360 / (24 * 3600)  # days
    skip = 10

    max_error = 10
    max_iter = 10

    ns = len(r.rs.segments)
    nodes = r.get_nodes()
    # seg_length = r.rs.segLength()
    assert len(nodes) - 1 == ns, "Number of nodes should be equal one less than number of segments"

    """ Fetch rhizosphere model params """
    kr_ = np.array(r.getKr(rs_age))
    kr_ = np.maximum(kr_, np.ones(kr_.shape) * 1.e-4)
    inner_ = r.rs.radii
    inner_kr_ = np.multiply(inner_, kr_)  # multiply for table look up
    inner_kr_ = np.expand_dims(inner_kr_, axis = 1)
    inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)
    inner_kr_ = np.minimum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-4)
    rho_ = np.maximum(rho_, np.ones(rho_.shape) * 1.)
    rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200.)

    collar_index = r.collar_index()

    """ Numerical solution """
    start_time = timeit.default_timer()
    x_, y_, z_, sink_, hs_, hx_, hsr_ = [], [], [], [], [], [], []
    sx = s.getSolutionHead_()  # inital condition, solverbase.py

    N = round(sim_time / dt)
    t = 0.

    t_pot = -trans * sinusoidal2(0., dt)
    rx = r.solve(rs_age, t_pot, 0., sx, True, wilting_point, [])  # xylem_flux.py, cells = True
    rx = np.expand_dims(rx, axis = 1)
    print(np.min(rx), np.max(rx))

    for i in range(0, N):

        t_pot = -trans * sinusoidal2(t, dt)  # potential transpiration ...
        if  i % skip == 0:
            print("t_pot", t_pot)

        err = 1.e6
        c = 0
        rx_old = rx  # in total potential
        while err > max_error and c < max_iter:  # max_iter = 1000 , err>100 works fine

            """ interpolation """
            wall_interpolation = timeit.default_timer()

            # rx = np.maximum(rx, np.ones(rx.shape) * (-15999))
            # rx = np.minimum(rx, np.ones(rho_.shape) * 0.)

            hs = np.transpose(np.array([[sx[mapping[j]] for j in range(0, ns)]]))
            # print(hs.shape, rx[1:].shape, ns)
            rsx = soil_root_interface_table(rx[1:], hs, inner_kr_, rho_, sra_table_lookup)
            # print("rsx", rsx.shape, np.min(rsx), np.max(rsx))

            wall_interpolation = timeit.default_timer() - wall_interpolation

            """ xylem matric potential """
            wall_xylem = timeit.default_timer()

            rx = r.solve(rs_age, t_pot, 0., rsx, False, wilting_point, [])  # xylem_flux.py, cells = True
            rx = np.expand_dims(rx, axis = 1)
            # print(np.min(rx), np.max(rx))

            err = np.linalg.norm(rx - rx_old)
            wall_xylem = timeit.default_timer() - wall_xylem
            rx_old = rx.copy()
            c += 1

        seg_fluxes = r.segFluxes(rs_age, rx, rsx, False, False)
        sum_root_flux = np.sum(seg_fluxes)
        fluxes = r.sumSegFluxes(seg_fluxes)

        water = s.getWaterVolume()
        s.setSource(fluxes.copy())  # richards.py
        s.solve(dt)
        sum_soil_flux = (s.getWaterVolume() - water) / dt
        old_sx = sx.copy()
        sx = s.getSolutionHead_()  # richards.py

        x_.append(t)
        y_.append(sum_soil_flux)  # cm3/day (soil uptake)
        z_.append(sum_root_flux)  # cm3/day (root system uptake)

        if  i % skip == 0:
            n = round(float(i) / float(N) * 100.)

            rx_ = np.zeros((ns, 1))
            for j in range(0, ns):  # from total to matric
                rx_[j, 0] = rx[j, 0]  # - nodes[j + 1][2]

            print("number of iterations", c)
            print("t_pot", t_pot, "q_root", sum_root_flux, "soil", sum_soil_flux)
            print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days {:g}\n"
                  .format(np.min(sx), np.max(sx), np.min(rx_), np.max(rx_), s.simTime, rx[collar_index, 0]))

            """ remember results """
            sink = np.zeros(sx.shape)
            for k, v in fluxes.items():
                sink[k] += v
            sink_.append(sink)  # [cm3/day] per soil cell
            hs_.append(sx.copy())  # [cm] per soil cell
            hx_.append(rx)  # [cm] total potential per segment
            hsr_.append(rsx)  # [cm] total potential per segment

        t += dt

    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return hx_, hsr_, sink_, x_, y_, z_, hs_, dt


def run_sraOld(sim_time, method, plant, dim, soil, outer_method):

    name = "sraOld_" + plant + "_" + dim + "_" + soil + "_" + outer_method
    print(name)

    # hidden parameters...
    initial = -200  # cm     plot_transpiration(x_, y_, z_, lambda t: trans * sinusoidal2(t, dt))

    print(name, "\n")

    r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = set_scenario(plant, dim, initial, soil, outer_method)
    hx_, hsr_, sink_, x_, y_, z_, hs_, dt = simulate_sraOld(sim_time, r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping)

    s.writeDumuxVTK("results/" + name)
    write_files(name, hx_, hsr_, sink_, x_, y_, z_, hs_)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = 'Simulation options')
    parser.add_argument('plant', type = str, help = 'soybean or maize')
    parser.add_argument('dim', type = str, help = '1D or 3D')
    parser.add_argument('soil', type = str, help = 'soil type (hydrus_loam, hydrus_clay, hydrus_sand or hydrus_sandyloam)')
    parser.add_argument('outer_method', type = str, help = 'how to determine outer radius (voronoi, length, surface, volume)')

    # args = parser.parse_args(['maize', "1D", "hydrus_loam", "surface"])
    args = parser.parse_args()

    name = "old_" + args.plant + "_" + args.dim + "_" + args.soil + "_" + args.outer_method
    print(name, "\n")
    
    initial = -200  # cm     plot_transpiration(x_, y_, z_, lambda t: trans * sinusoidal2(t, dt))
    sim_time = 7.5



    r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = set_scenario(args.plant, args.dim, initial, args.soil, args.outer_method)

    hx_, hsr_, sink_, x_, y_, z_, hs_, dt = simulate_sraOld(sim_time, r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping)

    """ write """
    s.writeDumuxVTK("results/" + name)
    write_files(name, hx_, hsr_, sink_, x_, y_, z_, hs_)
    # plot_transpiration(x_, y_, z_, lambda t: trans * sinusoidal2(t, dt), name)
    #
    # """ plot results """
    # hx = hx_[-1]
    # sx = hs_[-1]
    # ana = pb.SegmentAnalyser(r.rs.mappedSegments())
    # ana.addData("pressure head", hx)
    # # vp.plot_roots(ana, "hx")
    #
    # ns = len(r.rs.segments)
    # hs = np.transpose(np.array([[sx[mapping[j]][0] for j in range(0, ns)]]))
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
    # vp.plot_roots_and_soil(ana, "axial_flux", None, s, True, min_b, max_b, cell_number)  # NOT WORKING .... (?)
    # vp.plot_roots_and_soil(ana, "pressure head", None, s, True, min_b, max_b, cell_number)
