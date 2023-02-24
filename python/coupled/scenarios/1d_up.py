""" 
Jan's new scenarios with the upscaled sink (without SRA)

1d soil, 2 cm thick layers, dynamic conductivities (see root_conductivities.py)
"""
import sys; sys.path.append("../../../modules/"); sys.path.append("../../../../../CPlantBox/");  sys.path.append("../../../../../CPlantBox/src/python_modules")
sys.path.append("../../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../modules/fv/");

from scenario_setup import *

import timeit
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy import sparse
import scipy.sparse.linalg as LA

""" 
Initialize  
"""

name = "small_up0"  # name to export resutls
sstr = "_dry"

r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario1D(sstr, doussan = True)
min_b, max_b, cell_number = get_domain1D()

nodes = r.get_nodes()

print("first segment", r.rs.segments[0], r.rs.nodes[0])
print("mapping of first segment", mapping[0], mapping[1])

s_ = r.rs.segments[23]
print("23 segment", s_, r.rs.nodes[s_.x], r.rs.nodes[s_.y])
print("mapping of 23", mapping[23])
print("pick -1", s.pick([0., 0., -1]))
print("pick -2", s.pick([0., 0., -2]))
print("pick -3", s.pick([0., 0., -2.09]))

# """ for fixed root system """
# types = r.rs.subTypes
# inner_ = r.rs.radii
# outer_ = r.rs.segOuterRadii()
# rho_ = np.divide(outer_, np.array(inner_))
# ns = len(types)

""" experimental stuff """
n = len(r.rs.nodes)
ns = len(r.rs.segments)

Id = sparse.identity(ns).tocsc()  # identity matrix

kx_ = np.divide(r.getKx(sim_time), r.rs.segLength())  # / dl
Kx = sparse.diags(kx_).tocsc()
# print("Kx", Kx.shape, Kx[0, 0], Kx[1, 1])

kr_ = np.array(r.getEffKr(sim_time))  # times surface (2 a pi length)
Kr = sparse.diags(kr_).tocsc()
# print("Kr", Kr.shape, Kr[0, 0], Kr[1, 1])

IM = r.get_incidence_matrix().tocsc()
IMt = IM.transpose().tocsc()

L = IMt @ Kx @ IM  # Laplacian, Hess Eqn (10)
L = L[1:, 1:]  # == L_{N-1}
# print("L", L.shape)

print("inv start")

A_dirichlet = (L + Kr).tocsc()
Ainv_dirichlet = sparse.linalg.inv(A_dirichlet).todense()  # dense

A_neumann = A_dirichlet
A_neumann[0, 0] -= kx_[0]
Ainv_neumann = sparse.linalg.inv(A_neumann).todense()  # dense

C_comp_dirichlet = Kr @ (Id - Ainv_dirichlet @ Kr)  # Neumann, Hess, Eqn (24)
c_dirichlet = (Kr @ Ainv_dirichlet)[:, 0] * (-kx_[0])  # # Hess (25)
print("C_comp_dirichlet", type(C_comp_dirichlet), C_comp_dirichlet.shape)
print("c_dirichlet", type(c_dirichlet), c_dirichlet.shape)

C_comp_neumann = Kr @ (Id - Ainv_neumann @ Kr)  # Neumann, Hess, Eqn (32)
c_neumann = (Kr @ Ainv_neumann)[:, 0]  # Hess (33)
print("C_comp_neumann", type(C_comp_neumann), C_comp_neumann.shape)
print("c_neumann", type(c_neumann), c_neumann.shape)

print("inv stop")

# B, soil2matrix = r.get_soil_matrix()
# Bt = B.transpose()
#
# print("up start")
# print(Bt.shape, C_comp_neumann.shape, B.shape)
# C_comp_neumann_up = B @ C_comp_neumann @ Bt
# print(C_comp_neumann_up.shape)
# print("up end")

""" Numerical solution """
NT = int(np.ceil(sim_time / dt))  # number of iterations

start_time = timeit.default_timer()

sx = s.getSolutionHead()

x_, y_, w_, cf = [], [], [], []
psi_x_, psi_s_, sink_, psi_s2_ = [], [], [], []

for i in range(0, NT):

    t = i * dt  # current simulation time

    """ 1. xylem model """
    t_pot = -trans * sinusoidal(t)  # potential transpiration ...
    print("t_pot", t_pot)

    hs = np.transpose(np.array([[sx[mapping[j]][0] for j in range(0, ns)]]))
    print("hs", hs.shape, np.min(hs), np.max(hs), np.argmin(hs))
    for j in range(0, len(nodes) - 1):  # from matric to total
        hs[j, 0] += nodes[j + 1][2]

    print("sx", sx.shape, np.min(sx), np.max(sx), np.argmin(sx))

    hx = Ainv_neumann.dot(Kr.dot(hs)) + Ainv_neumann[:, 0] * t_pot  #   # Hess Eqn (29)

    hx0 = hx[0, 0] + t_pot / kx_[0]  # /kx*l
    print("hx0", hx0, hx[0, 0])
    # hx = (Ainv_dirichlet @ Kr).dot(hs) + Ainv_dirichlet[:, 0] * kx_[0] * wilting_point  #   # Hess Eqn (29)
    # print("hx", hx.shape, hx[0], hx[-1])
    # q = Kr.dot(hs - hx)
    # print("sum_q", np.sum(q), q[0], q[1])
    # hx = [0]

    q_neumann = C_comp_neumann.dot(hs) + c_neumann * t_pot
    # print("q_neumann", q_neumann.shape)

    q_dirichlet = -(C_comp_dirichlet.dot(hs) + c_dirichlet * wilting_point)
    # print("q_dirichlet", q_dirichlet.shape)
    print()
    # print("q", q.shape)
    print("soil min", np.min(hs), "at", np.argmin(hs))
    print()
    # print("hx", np.min(hx), np.max(hx), hx[0], hx[1], hx.shape)
    print("sum(q_neumann)", np.sum(q_neumann), q_neumann[0], q_neumann[1])
    print("sum(q_dirichlet)", np.sum(q_dirichlet), q_dirichlet[0], q_dirichlet[1])
    print()
    print("mapping of first segment", mapping[0], sx[mapping[0]], "soil", hs[0], "root", hx[0])
    print("mapping of 23", mapping[23], sx[mapping[23]], "soil", hs[23], "root", hx[23])

    # fluxes = r.sumSegFluxes(q_dirichlet)
    if hx0 < wilting_point:
        print("dirichlet")
        fluxes = r.sumSegFluxes(q_dirichlet[:, 0])
    else:
        print("neumann")
        fluxes = r.sumSegFluxes(q_neumann[:, 0])
        # print(fluxes)

    # per cell
    # hs = sx to

    # # for # mapping[j]][0] for j in range(0, ns)
    # # Q, soil2matrix[i]
    # hs = np.zeros((C_comp_up.shape[0],))
    # for i in range(0, sx.shape[0]):
    #     hs[soil2matrix[i]] = sx[i]

    """ 2. soil model """
    s.setSource(fluxes.copy())  # richards.py
    s.solve(dt)
    sx = s.getSolutionHead()  # richards.py

    """ remember results ... """
    if i % skip == 0:
        sx_ = sx[:, 0]
        rx = hx
        psi_x_.append(rx.copy())  # cm (per root node)
        psi_s_.append(np.array([sx_[ci] for ci in mapping]))  # cm (per root segment)
        sink = np.zeros(sx_.shape)
        for k, v in fluxes.items():
            sink[k] += v
        sink_.append(sink)  # cm3/day (per soil cell)
        x_.append(t)  # day
        y_.append(np.sum(sink))  # cm3/day
        psi_s2_.append(sx_)  # cm (per soil cell)
        min_sx, min_rx, max_sx, max_rx = np.min(sx), np.min(rx), np.max(sx), np.max(rx)
        n = round(float(i) / float(NT) * 100.)
        print("\n[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g}, {:g}"
                .format(min_sx, max_sx, min_rx, max_rx, np.sum(sink), -trans * sinusoidal(t)))

file1 = 'results/psix_' + name + sstr  # per segment
np.save(file1, np.array(psi_x_))
file2 = 'results/psiinterface_' + name + sstr  # per segment
np.save(file2, np.array(psi_s_))
file3 = 'results/sink_' + name + sstr
np.save(file3, np.array(-np.array(sink_)))
file4 = 'results/transpiration_' + name + sstr
np.save(file4, np.vstack((x_, -np.array(y_))))
file5 = 'results/soil_' + name + sstr
np.save(file5, np.array(psi_s2_))
