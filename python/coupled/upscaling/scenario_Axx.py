""" 
    Most exact root hydraulic model Ayz
    
    x: A: exact, B: aggregated over soil elemnts, C: parallel root system 
    y: A: Voronoi, B: Uniformly distributed
    z: A: 3D, B: 2D, C: 1D

    Daniel Leitner 1.2024
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
import visualisation.vtk_plot as vp

# from rhizo_models import plot_transpiration
from scenario_setup import *


def simulate_sra(sim_time, r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping):
    """
    assumes a static root architecure
    values are retrieved in set_scenario(), in scenario_setup.py
    
    sim_time             final simulation time
    r                    reference to PlantHydraulicModel (holds reference to MappedSegments)
    rho_                 precomputed outer perirhizal radii (
    rs_age               root system age 
    trans                daily transpiration rate [cm3/day], shape of sinusoidal2 is applied
    wilting_point        plant wilting point (normally -15000 cm)
    soil                 van Genuchten parameter set (of type vg.Parameters)
    s                    Macroscopic soil grid, e.g. of type RichardsWrapper(RichardsSP())
    sra_table_lookup     soil look up table for perirhizal model 
    mapping              mapping segment to cell (np.array)
    """

    print("\nInitial root sytstem age (sra)", rs_age)

    dt = 360 / (24 * 3600)  # days
    skip = 10

    max_error = 10
    max_iter = 100

    ns = len(r.rs.segments)
    nodes = r.get_nodes()
    assert len(nodes) - 1 == ns, "Number of nodes should be equal one less than number of segments"

    """fetch rhizosphere model params """
    kr_ = np.array(r.get_kr(rs_age))
    inner_ = r.rs.radii
    inner_kr_ = np.multiply(inner_, kr_)  # multiply for table look up
    inner_kr_ = np.expand_dims(inner_kr_, axis = 1)
    n_lower = np.sum(inner_kr_ < 1.e-7)
    n_uper = np.sum(inner_kr_ > 1.e-4)
    print("inner_kr truncated ", n_lower, n_uper, "segments")
    inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)  # truncate at look up table boundaries
    inner_kr_ = np.minimum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-4)
    rho_ = np.maximum(rho_, np.ones(rho_.shape) * 1.)
    rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200.)
    n_lower = np.sum(rho_ < 1)
    n_uper = np.sum(rho_ > 200)
    print("rho truncated ", n_lower, n_uper, "segments")

    """ initialize hydraulic model """
    r.update(rs_age)

    """ numerical solution """
    start_time = timeit.default_timer()
    x_, y_, z_, sink_, hs_, hx_, hsr_ = [], [], [], [], [], [], []

    N = round(sim_time / dt)
    t = 0.

    sx = s.getSolutionHead_()  # inital condition, solverbase.py
    hs = np.expand_dims(r.get_hs(sx), axis = 1)
    for j in range(0, len(nodes) - 1):  # from matric to total
        hs[j, 0] += nodes[j + 1][2]
    rx = r.solve(hs, 0., wilting_point)

    for i in range(0, N):

        t_pot = -trans * sinusoidal2(t, dt)  # potential transpiration ...
        hs = np.expand_dims(r.get_hs(sx), axis = 1)

        wall_iteration = timeit.default_timer()
        err = 1.e6
        c = 0
        rx_old = rx  # in total potential
        while err > max_error and c < max_iter:

            """ interpolation """
            wall_interpolation = timeit.default_timer()

            for j in range(0, ns):  # from total to matric
                rx[j, 0] -= nodes[j + 1][2]

            rx = np.maximum(rx, np.ones(rx.shape) * (-15999))
            # rx = np.minimum(rx, np.ones(rho_.shape) * 0.)
            # hs = np.maximum(hs, np.ones(hs.shape) * (-15999))

            rsx = soil_root_interface_table(rx, hs, inner_kr_, rho_, sra_table_lookup)
            for j in range(0, ns):  # from matric to total
                rsx[j, 0] += nodes[j + 1][2]

            wall_interpolation = timeit.default_timer() - wall_interpolation

            """ xylem matric potential """
            wall_xylem = timeit.default_timer()

            rx = r.solve(rsx, t_pot, wilting_point)
            q = r.radial_fluxes(rx, rsx)  # -Kr.dot(rsx - rx)  # both total potentials
            wall_xylem = timeit.default_timer() - wall_xylem

            err = np.linalg.norm(rx - rx_old)
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

            # rx_ = np.zeros((ns, 1))
            # rsx_ = np.zeros((ns, 1))
            # for j in range(0, ns):  # from total to matric
            #     #   print(nodes[j + 1][2])
            #     rx_[j, 0] = rx[j, 0] - nodes[j + 1][2]
            #     rsx_[j, 0] = rsx[j, 0] - nodes[j + 1][2]

            print("number of iterations", c)
            print("t_pot", t_pot, "q_root", sum_root_flux, "soil", sum_soil_flux)

            print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days {:g}\n"
                  .format(np.min(sx), np.max(sx), np.min(rx), np.max(rx), s.simTime, rx[r.ci, 0]))
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
            hx_.append(rx)  # [cm] total potential per segment
            hsr_.append(rsx)  # [cm] total potential per segment

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

    # args = parser.parse_args(['springbarley', "1D", "hydrus_sandyloam", "voronoi"])
    args = parser.parse_args()

    name = "sra_" + args.plant + "_" + args.dim + "_" + args.soil + "_" + args.outer_method
    print()
    print(name, "\n")

    initial = -200  # cm
    sim_time = 14.5  # day

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
