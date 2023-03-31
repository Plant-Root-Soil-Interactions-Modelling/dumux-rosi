""" 
static root system in soil (1D or 3D) outer radii with Voronoi method, or via densities 

aggregated hydraulic model with aggregated perirhizal nonlinear resistances 
(aggregated  over soil cells, steady rate approach and fixed-point-iteration on aggregated values, 
in new manuscript notation)
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


def simulate_agg(plant, dim, soil, outer_method, name):

    r, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = set_scenario(plant, dim, soil, outer_method)

    print("\nInitial root sytstem age", rs_age)
    # print("rs_age", rs_age)
    # rs_age = r.get_ages(rs_age)
    # print("rs_age", np.max(rs_age))

    sim_time = 14.
    dt = 360 / (24 * 3600)  # days
    skip = 10

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

    """ Aggregate over soil cells """
    B, soil2matrix, matrix2soil = r.get_soil_matrix()
    nmax = len(matrix2soil)
    Bt = B.transpose()
    BBt_inv = sparse.linalg.inv(B @ Bt)  # sparse
    print(type(BBt_inv))

    AinvKr_neumann_up = (((B @ Ainv_neumann) @ Kr) @ Bt)
    Ainv_neumann_up = B @ Ainv_neumann
    C_comp_neumann_up = B @ C_comp_neumann @ Bt
    c_neumann_up = B @ c_neumann

    AinvKr_dirichlet_up = (((B @ Ainv_dirichlet) @ Kr) @ Bt)
    Ainv_dirichlet_up = B @ Ainv_dirichlet
    C_comp_dirichlet_up = B @ C_comp_dirichlet @ Bt
    c_dirichlet_up = B @ c_dirichlet
    # print(C_comp_neumann_up.shape, type(C_comp_neumann_up))

    Kr_up = B @ Kr @ Bt  # sparse
    Kr_up_inv = sparse.linalg.inv(Kr_up)

    inner_kr_up = BBt_inv.dot(B.dot(inner_kr_))
    inner_kr_up = np.maximum(inner_kr_up, np.ones(inner_kr_up.shape) * 1.e-7)  ############################################ (too keep within table)
    inner_kr_up = np.minimum(inner_kr_up, np.ones(inner_kr_up.shape) * 1.e-4)  ############################################ (too keep within table)

    rho_up = BBt_inv.dot(B.dot(rho_))

    print("up end")

    """ Numerical solution (a) """
    start_time = timeit.default_timer()
    rs_age = np.max(r.get_ages())
    x_, y_, w_, cpx, cps, cf = [], [], [], [], [], []
    sink1d = []

    N = round(sim_time / dt)
    t = 0.

    rx = [0]
    sx = s.getSolutionHead()  # inital condition, solverbase.py
    centers = s.getCellCenters()

    t_pot = -trans * sinusoidal(t)
    hs_ = np.zeros((nmax, 1))  # sx -> hs_ # soil cell indices to soil matrix indices
    for j in soil2matrix.keys():
            hs_[soil2matrix[j]] += sx[j] + centers[j, 2]
    hxd = BBt_inv.dot(AinvKr_dirichlet_up.dot(hs_) + Ainv_dirichlet_up[:, 0] * kx_[0] * wilting_point)
    q_dirichlet_up = -Kr_up.dot(hs_ - hxd)
    if np.sum(q_dirichlet_up) > t_pot:
        rx = hxd
    else:
        rx = BBt_inv.dot(AinvKr_neumann_up.dot(hs_) + Ainv_neumann_up[:, 0] * t_pot)


if __name__ == "__main__":

    # TODO parse from args?
    plant = "soybean"  # soybean, maize
    dim = "3D"  # 1D, 3D
    soil = "jan_comp"  #  loam, clay sand
    outer_method = "surface"  # voronoi, length, surface, volume

    name = "sra_" + plant + "_" + dim + "_" + soil + "_" + dim + "_" + outer_method
    print(name, "\n")
    simulate_agg(plant, dim, soil, outer_method, name)

