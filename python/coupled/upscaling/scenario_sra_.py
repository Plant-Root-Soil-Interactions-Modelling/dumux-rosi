""" 
New reference scenario for upscaling (DL 29.3.2023)

static root system in soil (1D or 3D) outer radii with Voronoi method, or via densities 

full hydraulic model with perirhizal nonlinear resistances 
(steady rate approach and fixed-point-iteration in new manuscript notation)

! no matrix inversion, not working :-( !
"""
import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import timeit
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse

from functional.xylem_flux import sinusoidal2
import visualisation.vtk_plot as vp
from rhizo_models import plot_transpiration
from scenario_setup import *


def simulate_sra(sim_time, r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping, name):

    print("\nInitial root sytstem age", rs_age)
    # print("rs_age", rs_age)
    # rs_age = r.get_ages(rs_age)
    # print("rs_age", np.max(rs_age))

    dt = 360 / (24 * 3600)  # days
    skip = 1

    max_error = 10
    max_iter = 1000

    ns = len(r.rs.segments)
    nodes = r.get_nodes()
    # seg_length = r.rs.segLength()
    assert len(nodes) - 1 == ns, "Number of nodes should be equal one less than number of segments"

    """ Fetch rhizosphere model params """
    kr_ = np.array(r.getKr(rs_age))
    inner_ = r.rs.radii
    inner_kr_ = np.multiply(inner_, kr_)  # multiply for table look up
    inner_kr_ = np.expand_dims(inner_kr_, axis = 1)
    # inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)
    # inner_kr_ = np.minimum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-4)
    # rho_ = np.maximum(rho_, np.ones(rho_.shape) * 1.)
    # rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200.)

    """ Initialize root hydraulic model """
    A_d, Kr, kx0 = r.doussan_system_matrix(rs_age)
    Id = sparse.identity(ns).tocsc()  # identity matrix

    A_n = A_d.copy()
    A_n[0, 0] -= kx0

    Kr_inv = sparse.diags(np.divide(np.ones(ns,), Kr.diagonal()))
    # Kr_inv = sparse.linalg.inv(Kr)  # Kr is a diagonal matrix, thus Kr_inv sparse

    A_dq = A_d @ Kr_inv
    A_nq = A_n @ Kr_inv
    Bd = A_d - Kr
    Bn = A_n - Kr

    """ Numerical solution """
    start_time = timeit.default_timer()
    x_, y_, c_, sink_, hs_, hx_, hsr_ = [], [], [], [], [], [], []
    sx = s.getSolutionHead()  # inital condition, solverbase.py

    N = round(sim_time / dt)
    t = 0.

    t_pot = -trans * sinusoidal2(0., dt)
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

        t_pot = -trans * sinusoidal2(t, dt)  # potential transpiration ...
        print("t_pot", t_pot)

        hs = np.transpose(np.array([[sx[mapping[j]][0] for j in range(0, ns)]]))
        for j in range(0, ns):  # from matric to total potential
            hs[j, 0] += nodes[j + 1][2]

        err = 1.e6
        c = 0
        rx_old = rx  # in total potential
        while err > max_error and c < max_iter:  # max_iter = 1000 , err>100 works fine

            """ interpolation """
            wall_interpolation = timeit.default_timer()

            for j in range(0, ns):  # from total to matric
                rx[j, 0] -= nodes[j + 1][2]

            # hx = np.maximum(hx, np.ones(hx.shape) * (-15000))
            # hx = np.minimum(hx, np.ones(rho_.shape) * 0.)
            hsr = soil_root_interface_table(rx, hs, inner_kr_, rho_, sra_table_lookup)
            for j in range(0, ns):  # from matric to total
                hsr[j, 0] += nodes[j + 1][2]

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

        soil_fluxes = r.sumSegFluxes(q[:, 0])
        water = s.getWaterVolume()
        s.setSource(soil_fluxes.copy())  # richards.py
        s.solve(dt)
        sum_soil_flux = (s.getWaterVolume() - water) / dt
        sum_root_flux = np.sum(q)

        old_sx = sx.copy()
        sx = s.getSolutionHead()  # richards.py
        water = s.getWaterVolume()

        # remember results
        x_.append(t)
        y_.append(sum_soil_flux)  # cm3/day (soil uptake)
        c_.append(sum_root_flux)  # cm3/day (root system uptake)

        if  i % skip == 0:
            n = round(float(i) / float(N) * 100.)

            print("number of iterations", c)
            print("t_pot", t_pot, "q_root", sum_root_flux, "soil", sum_soil_flux)
            print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days {:g}\n"
                  .format(np.min(sx), np.max(sx), np.min(hx), np.max(hx), s.simTime, hx[0, 0]))

            """ remember results """
            sink = np.zeros(sx.shape)
            for k, v in soil_fluxes.items():
                sink[k] += v
            sink_.append(sink)  # [cm3/day] per soil cell
            hs_.append(sx.copy())  # [cm] per soil cell
            hx_.append(hx)  # [cm] per segment
            hsr_.append(hsr)

        t += dt

    s.writeDumuxVTK("results/" + name)
    write_files(name, hx_, hsr_, sink_, x_, y_, c_, hs_)

    """ Plot """
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")
    plot_transpiration(x_, y_, c_, lambda t: trans * sinusoidal2(t, dt))


if __name__ == "__main__":

    # TODO parse from args?
    plant = "soybean"  # soybean, maize
    dim = "1D"  # 1D, 3D
    soil = "hydrus_loam"  #  hydrus_loam, hydrus_clay, hydrus_sand
    outer_method = "surface"  # voronoi, length, surface, volume
    initial = -200  # cm total potential

    sim_time = 2.
    name = "sra_" + plant + "_" + dim + "_" + soil + "_" + dim + "_" + outer_method
    print(name, "\n")

    r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = set_scenario(plant, dim, initial, soil, outer_method)
    print("trans", trans)

    simulate_sra(sim_time, r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping, name)
