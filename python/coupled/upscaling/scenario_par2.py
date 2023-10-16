""" 
static root system in soil (1D or 3D) outer radii with Voronoi method, or via densities 

parallel root system hydraulic model with perirhizal nonlinear resistances 
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
import scipy.sparse.linalg as LA

from functional.xylem_flux import sinusoidal2
import visualisation.vtk_plot as vp
from rhizo_models import plot_transpiration
from scenario_setup import *
import parallel_rs as par


def double_(rsx, rsx2):
    """ inserts dummy values for the artificial segments """
    rsx2[:, 1] = rsx  # rsx2.shape = (ns, 2)
    return np.array(rsx2.flat)  # 0, rsx[0], 0, rsx[1], ...


def simulate_const(s, r, sra_table_lookup, trans, sim_time, dt, wilting_point, rs_age, outer_r):
    """     
    simulates the coupled scenario       
        root architecture is not gowing  
        conductivities are not changing over time
        
    s
    r
    sra_table_lookup             potentials a root soil interface  
    trans
    sim_time    
    dt
    
    TODO recyle factorisation of left hand side ... 
    """
    skip = 10  # for output and results, skip iteration
    max_iter = 10  # maximum for fix point iteration

    max_error = 10
    max_iter = 100

    start_time = timeit.default_timer()

    nodes = r.rs.nodes
    segs = r.rs.segments
    ns = len(r.rs.segments)
    mapping = np.array([r.rs.seg2cell[j] for j in range(0, ns)])

    inner_r = r.rs.radii
    types = r.rs.subTypes
    rho_ = np.divide(outer_r, np.array(inner_r))
    rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200)  ############################################ (too keep within table)

    """ Doussan """
    A_d, Kr, kx0 = r.doussan_system_matrix(rs_age)
    Id = sparse.identity(ns).tocsc()  # identity matrix
    collar_index = r.collar_index()

    print("collar_index (segment index)", r.collar_index())
    print("kx0:", kx0)
    print()

    A_n = A_d.copy()
    A_n[collar_index, collar_index] -= kx0

    print("invert matrix start", ns)
    sys.stdout.flush()

    A_n_splu = LA.splu(A_n)
    A_d_splu = LA.splu(A_d)
    # Ainv_dirichlet = sparse.linalg.inv(A_d).todense()  # dense
    # Ainv_neumann = sparse.linalg.inv(A_n).todense()  # dense
    print("done inverting", "\n")
    sys.stdout.flush()

    """ Numerical solution """
    start_time = timeit.default_timer()
    x_, y_, z_, sink_, hs_, hx_, hsr_ = [], [], [], [], [], [], []

    sx = s.getSolutionHead()  # inital condition, solverbase.py
    cell_centers = s.getCellCenters()
    cell_centers_z = np.array([cell_centers[mapping[2 * j + 1]][2] for j in range(0, int(ns / 2))])
    seg_centers_z = np.array([0.5 * (nodes[segs[2 * j + 1].x].z + nodes[segs[2 * j + 1].y].z)  for j in range(0, int(ns / 2))])
    hsb = np.array([sx[mapping[2 * j + 1]][0] for j in range(0, int(ns / 2))])  # soil bulk matric potential per segment
    # print(list([mapping[2 * j + 1] for j in range(0, int(ns / 2))]))

    kr_ = np.zeros((ns,))
    rsx = hsb.copy()  # initial values for fix point iteration
    rsx2 = np.zeros((rsx.shape[0], 2))

    # r.init_solve_static(rs_age, double_(rsx, rsx2), False, wilting_point, soil_k = [])  # speed up & and forever static...

    rx = r.solve(rs_age, -trans * sinusoidal2(0., dt), 0., double_(rsx, rsx2), False, wilting_point, soil_k = [])
    rx_old = rx.copy()

    kr_ = np.array([r.kr_f(rs_age, types[j], 2, j) for j in range(0, len(outer_r))])
    inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up
    inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)
    inner_kr_ = np.minimum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-4)

    N = int(np.ceil(sim_time / dt))  # number of iterations
    t = 0

    """ simualtion loop """
    for i in range(0, N):

        t_pot = -trans * sinusoidal2(t, dt)  # potential transpiration ...
        if  i % skip == 0:
            print("t_pot", t_pot)

        t = i * dt  # current simulation time

        wall_iteration = timeit.default_timer()
        wall_fixpoint = timeit.default_timer()

        err = 1.e6
        c = 0
        while err > 1 and c < max_iter:

            """ interpolation """
            wall_interpolation = timeit.default_timer()
            rx_ = rx[1::2] - seg_centers_z  # from total matric potential to matric potential
            hsb_ = hsb - cell_centers_z  # from total potential to matric potential

            rx_ = np.maximum(rx_, np.ones(rx_.shape) * (-15999))
            rx_ = np.minimum(rx_, np.ones(rx_.shape) * 0.)
            hsb_ = np.maximum(hsb_, np.ones(hsb_.shape) * (-15999))
            hsb_ = np.minimum(hsb_, np.ones(hsb_.shape) * 0.)

            rsx = soil_root_interface_table(rx_, hsb_, inner_kr_[1::2], rho_[1::2], sra_table_lookup)  # [1::2] every second entry, starting from 1
            rsx = rsx + seg_centers_z  # from matric potential to total matric potential
            wall_interpolation = timeit.default_timer() - wall_interpolation

            """ xylem matric potential """
            wall_xylem = timeit.default_timer()
            rx = r.solve(rs_age, -trans * sinusoidal2(t, dt), 0., double_(rsx, rsx2), False, wilting_point, soil_k = [])  # xylem_flux.py, cells = False
            err = np.linalg.norm(rx - rx_old)
            wall_xylem = timeit.default_timer() - wall_xylem

            rx_old = rx.copy()
            c += 1

        wall_fixpoint = timeit.default_timer() - wall_fixpoint

        wall_soil = timeit.default_timer()

        fluxes = r.segFluxes(rs_age, rx, double_(rsx, rsx2), approx = False, cells = False)
        err2 = np.linalg.norm(-trans * sinusoidal2(t, dt) - np.sum(fluxes))
        if r.last == "neumann":
            if err2 > 1.e-6:
                print("error: potential transpiration differs summed radial fluxes in Neumann case" , err2, -trans * sinusoidal2(t, dt), np.sum(fluxes))

        soil_fluxes = r.sumSegFluxes(fluxes)
        water = s.getWaterVolume()
        s.setSource(soil_fluxes.copy())  # richards.py
        s.solve(dt)
        sum_soil_flux = (s.getWaterVolume() - water) / dt
        sx = s.getSolutionHead()[:, 0]  # richards.py
        hsb = np.array([sx[mapping[2 * j + 1]] for j in range(0, int(ns / 2))])

        wall_soil = timeit.default_timer() - wall_soil

        wall_iteration = timeit.default_timer() - wall_iteration

        x_.append(t)  # day
        y_.append(sum_soil_flux)  # cm3/day
        sum_root_flux = np.sum(fluxes)
        z_.append(sum_root_flux)

        """ remember results ... """
        if i % skip == 0:
            n = round(float(i) / float(N) * 100.)

            rx_ = rx.copy()
            rsx_ = rsx.copy()
            # rx_ = np.zeros((ns, 1))
            # rsx_ = np.zeros((ns, 1))
            # for j in range(0, ns):  # from total to matric
            #     #   print(nodes[j + 1][2])
            #     rx_[j, 0] = rx[j, 0] - nodes[j + 1][2]
            #     rsx_[j, 0] = rsx[j, 0] - nodes[j + 1][2]

            print("number of iterations", c)
            print("t_pot", t_pot, "q_root", sum_root_flux, "soil", sum_soil_flux)
            print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days {:g}\n"
                  .format(np.min(sx), np.max(sx), np.min(rx_), np.max(rx_), s.simTime, rx[collar_index]))  # , 0

            print("iteration (interpolation, xylem) : ", wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem))
            print("iteration, soil", wall_iteration / (wall_iteration + wall_soil), wall_soil / (wall_iteration + wall_soil))
            print()
            sys.stdout.flush()

            sink = np.zeros(sx.shape)
            for k, v in soil_fluxes.items():
                sink[k] += v
            sink_.append(sink)  # cm3/day (per soil cell)
            hs_.append(sx.copy())  # [cm] per soil cell
            hx_.append(rx.copy()[::2])  # cm (per root node)
            hsr_.append(rsx.copy())  # cm (per root segment)

    wall_time = timeit.default_timer() - start_time
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return hx_, hsr_, sink_, x_, y_, z_, hs_, dt, wall_time


def simulate_par(sim_time, r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping, outer_method):  # def simulate_const(s, r, sra_table_lookup, trans, sim_time, dt):

    dt = 360 / (24 * 3600)  # days

    hx_, hsr_, sink_, x_, y_, z_, sx_, dt, wall_time = simulate_const(s, r, sra_table_lookup, trans, sim_time, dt, wilting_point, rs_age, outer_method)  #  rename to an old script ...

    return hx_, hsr_, sink_, x_, y_, z_, sx_, dt, wall_time


def run_par(sim_time, method, plant, dim, soil, outer_method):

    # hidden parameters...
    initial = -200  # cm

    name = "par_" + plant + "_" + dim + "_" + soil + "_" + outer_method
    print(name, "\n")

    r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = set_scenario(plant, dim, initial, soil, outer_method)

    if plant == "maize":
        min_b, max_b, cell_number = maize_(dim)
    elif plant == "springbarley":
        min_b, max_b, cell_number = springbarley_(dim)
    elif plant == "soybean":
        min_b, max_b, cell_number = soybean_(dim)
    r_par, outer_r = par.create_parallel_rs(r, rs_age, s.getCellCenters(), min_b, max_b, cell_number, outer_method)
    if dim == "1D":
        picker = lambda x, y, z: s.pick([0., 0., z])
    else:
        picker = lambda x, y, z: s.pick([x, y, z])
    r_par.rs.setSoilGrid(picker)

    hx_, hsr_, sink_, x_, y_, z_, hs_, dt, wall_time = simulate_par(sim_time, r_par, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping, outer_r)

    s.writeDumuxVTK("results/" + name)  # final soil VTU
    write_files(name, hx_, hsr_, sink_, x_, y_, z_, hs_, wall_time)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = 'Simulation options')
    parser.add_argument('plant', type = str, help = 'soybean or maize or springbarley')
    parser.add_argument('dim', type = str, help = '1D or 3D')
    parser.add_argument('soil', type = str, help = 'soil type (hydrus_loam, hydrus_clay, hydrus_sand or hydrus_sandyloam)')
    parser.add_argument('outer_method', type = str, help = 'how to determine outer radius (voronoi, length, surface, volume)')

    args = parser.parse_args(['springbarley', "1D", "hydrus_loam", "length"])
    # args = parser.parse_args()

    name = "par_" + args.plant + "_" + args.dim + "_" + args.soil + "_" + args.outer_method
    print()
    print(name, "\n")
    sys.stdout.flush()

    initial = -200  # cm     plot_transpiration(x_, y_, z_, lambda t: trans * sinusoidal2(t, dt))
    sim_time = 14.5

    print("setting scenario")
    sys.stdout.flush()

    r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = set_scenario(args.plant, args.dim, initial, args.soil, args.outer_method)

    if args.plant == "maize":
        min_b, max_b, cell_number = maize_(args.dim)
    elif args.plant == "springbarley":
        min_b, max_b, cell_number = springbarley_(args.dim)
    elif args.plant == "soybean":
        min_b, max_b, cell_number = soybean_(args.dim)
    r_par, outer_r = par.create_parallel_rs(r, rs_age, s.getCellCenters(), min_b, max_b, cell_number, args.outer_method)
    if args.dim == "1D":
        picker = lambda x, y, z: s.pick([0., 0., z])
    else:
        picker = lambda x, y, z: s.pick([x, y, z])
    r_par.rs.setSoilGrid(picker)

    print("set_scenario done.")
    sys.stdout.flush()

    hx_, hsr_, sink_, x_, y_, z_, hs_, dt, wall_time = simulate_par(sim_time, r_par, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping, outer_r)
    # """ write """
    s.writeDumuxVTK("results/" + name)
    write_files(name, hx_, hsr_, sink_, x_, y_, z_, hs_, wall_time)

    sys.stdout.flush()
#
