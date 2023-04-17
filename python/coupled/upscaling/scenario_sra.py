""" 
New reference scenario for upscaling (DL 29.3.2023)

static root system in soil (1D or 3D) outer radii with Voronoi method, or via densities 

full hydraulic model with perirhizal nonlinear resistances 
(steady rate approach and fixed-point-iteration in new manuscript notation)
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


def simulate_sra(r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping, name):

    print("\nInitial root sytstem age", rs_age)
    # print("rs_age", rs_age)
    # rs_age = r.get_ages(rs_age)
    # print("rs_age", np.max(rs_age))

    sim_time = 14.
    dt = 3600 / (24 * 3600)  # days
    skip = 1

    max_error = 10
    max_iter = 1000

    ns = len(r.rs.segments)
    nodes = r.get_nodes()
    seg_length = r.rs.segLength()
    assert len(nodes) - 1 == ns, "number of nodes should be equal one less than number of segments"

    """ Fetch rhizosphere model params """
    kr_ = np.array(r.getKr(rs_age))
    kr_min = np.ones(kr_.shape) * 1.e-6
    kr_[kr_ == 0] = kr_min[kr_ == 0]
    # print(kr_)
    inner_ = r.rs.radii
    inner_kr_ = np.multiply(inner_, kr_)  # multiply for table look up
    inner_kr_ = np.expand_dims(inner_kr_, axis = 1)
    # inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)
    inner_kr_ = np.minimum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-4)
    rho_ = np.maximum(rho_, np.ones(rho_.shape) * 1.)
    rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200.)

    """ Initialize root hydraulic model """
    Id = sparse.identity(ns).tocsc()  # identity matrix
    kx_ = np.divide(r.getKx(rs_age), seg_length)  # / dl (Eqn 5)
    Kx = sparse.diags(kx_).tocsc()
    # print("Kx", Kx.shape, Kx[0, 0], Kx[1, 1])
    # kr_ = np.array(r.getEffKr(rs_age))  # times surface (2 a pi length), (Eqn 7)
    kr_ = 2 * np.pi * np.multiply(np.multiply(kr_, inner_), seg_length)  # Eqn 7
    Kr = sparse.diags(kr_).tocsc()
    # print("Kr", Kr.shape, Kr[0, 0], Kr[1, 1])

    C = r.get_incidence_matrix().tocsc()
    Ct = C.transpose().tocsc()
    L = Ct @ Kx @ C  # Laplacian (Eqn 4)
    L = L[1:, 1:]  # == L_{N-1} as in (Eqn 10 or 14)
    # print("L", L.shape)

    Ad = (L + Kr).tocsc()  # (Eqn 10)
    An = Ad.copy()
    An[0, 0] -= kx_[0]  # (Eqn 14)
    Kr_inv = sparse.diags(np.divide(np.ones(kr_.shape), kr_)).tocsc()
    # print("Kr_inv", np.min(Kr_inv), np.max(Kr_inv))

    Adq = Ad @ Kr_inv  # (Eqn 27)
    Anq = An @ Kr_inv  # (Eqn 30)
    Ad_Kr = Ad - Kr  # part of b (Eqn 27)
    An_Kr = An - Kr  # part of b (Eqn 30)

    """ Numerical solution """
    start_time = timeit.default_timer()
    x_, y_, c_, sink_, hs_, hx_, hsr_ = [], [], [], [], [], [], []
    sx = s.getSolutionHead()  # inital condition, solverbase.py
    # print("sx", sx.shape, np.min(sx), np.max(sx))
    # print(sx[s.pick([0., 0., 0.])], sx[s.pick([0., 0., -100.])])

    N = round(sim_time / dt)
    t = 0.

    # solve neumann no flux
    hs = np.transpose(np.array([[sx[mapping[j]][0] for j in range(0, ns)]]))
    for j in range(0, ns):  # from matric to total potential
        hs[j, 0] += nodes[j + 1][2]
    # print("hs", hs.shape, np.min(hs), np.max(hs))
    b = An_Kr.dot(hs)
    q = -sparse.linalg.spsolve(Anq, b)
    q = np.expand_dims(q, axis = 1)  # column vector
    # print("q", q.shape, np.sum(q)) ################################### wrong with maize ###########################################
    hx = hs - Kr_inv.dot(-q)
    print("hx", np.min(hx), np.max(hx))

    for i in range(0, N):

        t_pot = -trans * sinusoidal2(t, dt)  # potential transpiration ...
        # print("t_pot", t_pot)

        hs = np.transpose(np.array([[sx[mapping[j]][0] for j in range(0, ns)]]))
        for j in range(0, ns):  # from matric to total potential
            hs[j, 0] += nodes[j + 1][2]

        err = 1.e6
        c = 0
        rx_old = hx  # in total potential
        while err > max_error and c < max_iter:  # max_iter = 1000 , err>100 works fine

            """ interpolation """
            wall_interpolation = timeit.default_timer()

            for j in range(0, ns):  # from total to matric
                hx[j, 0] -= nodes[j + 1][2]
            hx = np.maximum(hx, np.ones(hx.shape) * (-15000))
            hsr = soil_root_interface_table(hx, hs, inner_kr_, rho_, sra_table_lookup)
            for j in range(0, ns):  # from matric to total
                hsr[j, 0] += nodes[j + 1][2]

            wall_interpolation = timeit.default_timer() - wall_interpolation

            """ xylem matric potential """
            wall_xylem = timeit.default_timer()

            # solve dirichlet uptake (Eqn 27)
            bd = Ad_Kr.dot(hsr)
            bd[0, 0] -= kx_[0] * wilting_point
            q = -sparse.linalg.spsolve(Adq, bd)
            # print("q_d", np.sum(q))
            if np.sum(q) < t_pot:  # solve neumann uptake (Eqn 30)
                b = An_Kr.dot(hsr)
                b[0, 0] -= t_pot
                q = -sparse.linalg.spsolve(Anq, b)
                # print("q_n", np.sum(q))

            q = np.expand_dims(q, axis = 1)  # column vector
            hx = hsr - Kr_inv.dot(-q)
            # print("hx", np.min(hx), np.max(hx))

            err = np.linalg.norm(hx - rx_old)
            # print("err", err)
            wall_xylem = timeit.default_timer() - wall_xylem
            rx_old = hx.copy()
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
    outer_method = "voronoi"  # voronoi, length, surface, volume
    initial = -400  # cm total potential

    name = "sra_" + plant + "_" + dim + "_" + soil + "_" + dim + "_" + outer_method
    print(name, "\n")

    r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = set_scenario(plant, dim, initial, soil, outer_method)
    simulate_sra(r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping, name)
