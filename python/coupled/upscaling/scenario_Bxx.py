"""     
    Upscaled aggregated root hydraulic model Byz
    
    x: A: exact, B: aggregated over soil elemnts, C: parallel root system 
    y: A: Voronoi, B: Uniformly distributed
    z: A: 3D, B: 2D, C: 1D

    Daniel Leitner 1.2024
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

    print("\nInitial root sytstem age (agg)", rs_age)

    dt = 360 / (24 * 3600)  # days
    skip = 10

    max_error = 10
    max_iter = 100

    ns = len(r.rs.segments)
    nodes = r.get_nodes()
    seg_length = r.rs.segLength()
    assert len(nodes) - 1 == ns, "number of nodes should be equal one less than number of segments"

    """ Fetch rhizosphere model params """
    kr_ = np.array(r.get_kr(rs_age))
    inner_ = r.rs.radii
    inner_kr_ = np.multiply(inner_, kr_)  # multiply for table look up
    inner_kr_ = np.expand_dims(inner_kr_, axis = 1)
    inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)
    inner_kr_ = np.minimum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-4)
    rho_ = np.maximum(rho_, np.ones(rho_.shape) * 1.)
    rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200.)

    """ Doussan """
    Ad, Kr, kx0 = r.get_doussan_system(rs_age)
    Id = sparse.identity(ns).tocsc()  # identity matrix
    collar_index = r.collar_index()

    print("collar_index (segment index)", r.collar_index())
    print("kx0:", kx0)
    print()

    """ upscaling """
    M, soil2matrix, matrix2soil = r.get_soil_matrix()
    nmax = len(matrix2soil)
    Mt = M.transpose()
    MMt_inv = sparse.linalg.inv(M @ Mt)  # sparse

    print("inv start")
    sys.stdout.flush()

    Ad_inv = sparse.linalg.inv(Ad).todense()  # dense
    Ad_up = M @ Kr @ (Id - Ad_inv @ Kr) @ Mt  # [34]
    cd_up = M @ (Kr @ Ad_inv)[:, collar_index] * (kx0)
    del Ad_inv  # free the beast

    An = Ad.copy()
    An[collar_index, collar_index] -= kx0
    An_inv = sparse.linalg.inv(An).todense()  # dense
    An_up = M @ Kr @ (Id - An_inv @ Kr) @ Mt  # [34]
    cn_up = M @ (Kr @ An_inv)[:, collar_index]
    del An_inv  # free the beast

    print("inv stop")
    sys.stdout.flush()

    # AinvKr_dirichlet_up = (((M @ Ad_inv) @ Kr) @ Mt)
    # Ainv_dirichlet_up = M @ Ad_inv
    # AinvKr_neumann_up = M @ ((An_inv) @ Kr) @ Mt
    # Ainv_neumann_up = M @ An_inv

    Kr_up = M @ Kr @ Mt  # sparse
    Kr_up_inv = sparse.linalg.inv(Kr_up).todense()

    print(Kr_up.shape, Kr_up_inv.shape)
    print("finished Kr_up, Kr_up_inv")
    sys.stdout.flush()

    inner_kr_up = MMt_inv.dot(M.dot(inner_kr_))
    inner_kr_up = np.maximum(inner_kr_up, np.ones(inner_kr_up.shape) * 1.e-7)  ############################################ (too keep within table)
    inner_kr_up = np.minimum(inner_kr_up, np.ones(inner_kr_up.shape) * 1.e-4)  ############################################ (too keep within table)
    rho_up = MMt_inv.dot(M.dot(rho_))

    print("up end")
    sys.stdout.flush()

    """ Numerical solution (a) """
    start_time = timeit.default_timer()
    x_, y_, z_, sink_, sx_, hx_, hsr_ = [], [], [], [], [], [], []
    sx = s.getSolutionHead_()  # inital condition, solverbase.py

    N = round(sim_time / dt)
    t = 0.
    rx = [0]

    centers = s.getCellCenters()
    print("centers", np.min(centers[:, 2]), np.max(centers[:, 2]))

    t_pot = -trans * sinusoidal2(t, dt)  # potential transpiration ...
    # t_pot = -100
    hs_ = np.zeros((nmax, 1))  # sx -> hs_ # soil cell indices to soil matrix indices
    for j in soil2matrix.keys():
            hs_[soil2matrix[j]] += sx[j] + centers[j, 2]  # from matric to total

    q_dirichlet_up = -(Ad_up.dot(hs_) - cd_up * wilting_point)  # there was some SIGN MISTAKE
    rx = hs_ - Kr_up_inv.dot(-q_dirichlet_up)
    print("1) rx dirichlet", np.min(rx), np.max(rx), "q_dirichlet_up", np.sum(q_dirichlet_up))

    # rx = MMt_inv.dot(AinvKr_dirichlet_up.dot(hs_) + Ainv_dirichlet_up[:, 0] * kx0 * wilting_point)
    # q_dirichlet_up = -Kr_up.dot(hs_ - rx)
    # print("2) rx dirichlet", np.min(rx), np.max(rx), "q_dirichlet_up", np.sum(q_dirichlet_up))

    print("q_dirichlet_up", np.sum(q_dirichlet_up), "t_pot", t_pot)
    if np.sum(q_dirichlet_up) > t_pot:
        rx = hs_ - Kr_up_inv.dot(-q_dirichlet_up)
        print("rx dirichlet", np.min(rx), np.max(rx))
    else:
        q_neumann_up = An_up.dot(hs_) + cn_up * t_pot
        rx = hs_ - Kr_up_inv.dot(-q_neumann_up)
        print("q_neumann_up", np.sum(q_neumann_up), "t_pot", t_pot)
        print("rx neumann", np.min(rx), np.max(rx))

    for i in range(0, N):

        t_pot = -trans * sinusoidal2(t, dt)  # potential transpiration ..

        hs_ = np.zeros((nmax, 1))  # sx -> hs_ # soil cell indices to soil matrix indices
        for j in soil2matrix.keys():
                hs_[soil2matrix[j]] += sx[j] + centers[j, 2]

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

            rx = np.maximum(rx, np.ones(rx.shape) * (-15999))
            rx = np.minimum(rx, np.ones(rx.shape) * 0.)
            hs_ = np.maximum(hs_, np.ones(hs_.shape) * (-15999))

            rsx = soil_root_interface_table(rx, hs_, inner_kr_up, rho_up, sra_table_lookup)  # in matric potential

            for j in soil2matrix.keys():  # from matric to total
                    rsx[soil2matrix[j]] += centers[j, 2]

            wall_interpolation = timeit.default_timer() - wall_interpolation

            """ xylem matric potential """
            wall_xylem = timeit.default_timer()

            q_dirichlet_up = -(Ad_up.dot(rsx) - cd_up * wilting_point)  # there was some SIGN MISTAKE
            q_neumann_up = -(An_up.dot(rsx) - cn_up * t_pot)
            # rx = rsx - Kr_up_inv.dot(-q_dirichlet_up)
            # rx = MMt_inv.dot(AinvKr_dirichlet_up.dot(rsx) + Ainv_dirichlet_up[:, 0] * kx0 * wilting_point)
            # q_dirichlet_up2 = -Kr_up.dot(rsx - rx)
            # print("q_dirichlet_up", np.sum(q_dirichlet_up), np.sum(q_dirichlet_up2), "t_pot", t_pot, "rsx", np.min(rsx), np.max(rsx), "sx", np.min(sx), np.max(sx), "hs_", np.min(hs_), np.max(hs_))
            if np.sum(q_dirichlet_up) > t_pot:
                rx = rsx - Kr_up_inv.dot(-q_dirichlet_up)
                # print("rx dirichlet", np.min(rx), np.max(rx))
                # print("rx dirichlet", np.min(rx), np.max(rx))
            else:
                # rx = MMt_inv.dot(AinvKr_neumann_up.dot(rsx) + Ainv_neumann_up[:, 0] * t_pot)
                # q_neumann_up2 = -Kr_up.dot(rsx - rx)
                rx = rsx - Kr_up_inv.dot(-q_neumann_up)
                # rx = rsx - Kr_up_inv.dot(-q_neumann_up)
                # print("rx neumann", np.min(rx), np.max(rx))
                # print("q_neumann_up", np.sum(q_neumann_up), "t_pot", t_pot)
                # print("rx neumann", np.min(rx), np.max(rx))

            err = np.linalg.norm(rx - rx_old)
            wall_xylem = timeit.default_timer() - wall_xylem
            rx_old = rx.copy()
            c += 1

        wall_soil = timeit.default_timer()

        if np.sum(q_dirichlet_up) > t_pot:

            # print("dirichlet", np.sum(q_dirichlet_up), t_pot, np.sum(q_neumann_up))
            # if s.simTime>3.4:
            #     print(np.min(q_neumann_up), np.max(q_neumann_up), np.min(q_neumann_up2), np.max(q_neumann_up2))
            #     plt.plot(q_dirichlet_up, label = "dirichlet")
            #     plt.plot(q_neumann_up, label = "q")
            #     plt.plot(q_neumann_up2, label = "rsx")
            #     # plt.plot(rx, label = "rsx")
            #     # rx2 = rsx - Kr_up_inv.dot(-q_neumann_up)
            #     # plt.plot(rx2, ':', label = "q")
            #     plt.legend()
            #     plt.show()

            # print("q_dirichlet_up", c, err, np.sum(q_dirichlet_up), t_pot)
            fluxes = {}
            for j in range(0, nmax):
                fluxes[matrix2soil[j]] = q_dirichlet_up[j, 0]
            sum_root_flux = np.sum(q_dirichlet_up)
        else:
            # print("neumann", np.sum(q_dirichlet_up), t_pot, np.sum(q_neumann_up)) # , np.sum(q_neumann_up2)

            # if s.simTime>3.4:
                # print(np.min(q_neumann_up), np.max(q_neumann_up), np.min(q_neumann_up2), np.max(q_neumann_up2))
                # plt.plot(q_dirichlet_up, label = "dirichlet")
                # plt.plot(q_neumann_up, label = "q")
                # plt.plot(q_neumann_up2, label = "rsx")
                # # plt.plot(rx, label = "rsx")
                # # rx2 = rsx - Kr_up_inv.dot(-q_neumann_up)
                # # plt.plot(rx2, ':', label = "q")
                # plt.legend()
                # plt.show()

            # q_neumann_up = -Kr_up.dot(rsx - rx)
#             print("q_neumann0_up", c, err, q_neumann_up.shape, np.sum(q_neumann_up), t_pot)
            fluxes = {}
            for j in range(0, nmax):
                fluxes[matrix2soil[j]] = q_neumann_up[j, 0]
            sum_root_flux = np.sum(q_neumann_up)

        water = s.getWaterVolume()
        s.setSource(fluxes.copy())  # richards.py
        s.solve(dt)
        sum_soil_flux = (s.getWaterVolume() - water) / dt

        old_sx = sx.copy()
        sx = s.getSolutionHead_()  # richards.py

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
            sys.stdout.flush()

            """ remember results """
            sink = np.zeros(sx.shape)
            for k, v in fluxes.items():
                sink[k] += v
            sink_.append(sink)  # [cm3/day] per soil cell
            sx_.append(sx.copy())  # [cm] per soil cell
            # hx_.append(rx_)  # [cm] total potential per segment
            # hsr_.append(rsx_)  # [cm] total potential per segment

        t += dt

    wall_time = timeit.default_timer() - start_time
    print ("Coupled benchmark solved in ", wall_time, " s")
    sys.stdout.flush()

    return hx_, hsr_, sink_, x_, y_, z_, sx_, dt, wall_time


def run_agg(sim_time, method, plant, dim, soil, outer_method):

    # hidden parameters...
    initial = -200  # cm

    name = "agg_" + plant + "_" + dim + "_" + soil + "_" + outer_method
    print(name, "\n")

    r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = set_scenario(plant, dim, initial, soil, outer_method)

    hx_, hsr_, sink_, x_, y_, z_, hs_, dt, wall_time = simulate_agg(sim_time, r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping)

    s.writeDumuxVTK("results/" + name)  # final soil VTU
    write_files(name, hx_, hsr_, sink_, x_, y_, z_, hs_, wall_time)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = 'Simulation options')
    parser.add_argument('plant', type = str, help = 'soybean or maize or springbarley')
    parser.add_argument('dim', type = str, help = '1D or 3D')
    parser.add_argument('soil', type = str, help = 'soil type (hydrus_loam, hydrus_clay, hydrus_sand or hydrus_sandyloam)')
    parser.add_argument('outer_method', type = str, help = 'how to determine outer radius (voronoi, length, surface, volume)')

    # args = parser.parse_args(['maize', "1D", "hydrus_clay", "voronoi"])
    args = parser.parse_args()

    name = "agg_" + args.plant + "_" + args.dim + "_" + args.soil + "_" + args.outer_method
    print()
    print(name, "\n")
    sys.stdout.flush()

    initial = -200  # cm     plot_transpiration(x_, y_, z_, lambda t: trans * sinusoidal2(t, dt))
    sim_time = 14.5

    print("setting scenario")
    sys.stdout.flush()

    r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = set_scenario(args.plant, args.dim, initial, args.soil, args.outer_method)

    print("set_scenario done.")
    sys.stdout.flush()

    hx_, hsr_, sink_, x_, y_, z_, hs_, dt, wall_time = simulate_agg(sim_time, r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping)

    # """ write """
    s.writeDumuxVTK("results/" + name)
    write_files(name, hx_, hsr_, sink_, x_, y_, z_, hs_, wall_time)

    sys.stdout.flush()
#
