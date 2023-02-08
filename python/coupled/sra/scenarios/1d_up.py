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
sstr = "_wet"

r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario1D(sstr, doussan = True)
min_b, max_b, cell_number = get_domain1D()

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

print()
kx_ = np.divide(r.getKx(sim_time), r.rs.segLength())  # / dl
Kx = sparse.diags(kx_).tocsc()
print("kx0", kx_[0], kx_[1])
print("Kx", Kx.shape)

kr_ = np.array(r.getEffKr(sim_time))  # times surface (2 a pi length)
Kr = sparse.diags(kr_).tocsc()
print("kr0", kr_[0], kr_[1], kr_[-1])
print("Kr", Kr.shape)
print()

IM = r.get_incidence_matrix().tocsc()
IMt = IM.transpose().tocsc()

L = IMt @ Kx @ IM  # Laplacian, Hess Eqn (10)
L = L[1:, 1:]  # == L_{N-1}
print("L", L.shape)

print("inv start")

A_dirichlet = (L + Kr).tocsc()
Ainv_dirichlet = sparse.linalg.inv(A_dirichlet)

A_neumann = A_dirichlet
A_neumann[0, 0] -= kx_[0]
Ainv_neumann = sparse.linalg.inv(A_neumann)

C_comp_dirichlet = Kr @ (Id - Ainv_dirichlet @ Kr)  # Neumann, Hess, Eqn (24)
c_dirichlet = (Kr @ Ainv_dirichlet)[:, 0].todense() * (-kx_[0])  # # Hess (25)

C_comp_neumann = Kr @ (Id - Ainv_neumann @ Kr)  # Neumann, Hess, Eqn (32)
c_neumann = (Kr @ Ainv_neumann)[:, 0].todense()  # Hess (33)

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

print(np.min(mapping), np.max(mapping))
dd

for i in range(0, NT):

    t = i * dt  # current simulation time

    """ 1. xylem model """
    t_act = -trans * sinusoidal(t)  # potential transpiration ...

    # per segment
    hs = np.transpose(np.array([[sx[mapping[j]][0] for j in range(0, ns)]]))

    # hx = (Ainv_neumann @ Kr).dot(hs) + Ainv_neumann[:, 0] * t_act  #   # Hess Eqn (29)
    hx = (Ainv_dirichlet @ Kr).dot(hs) + Ainv_dirichlet[:, 0] * kx_[0] * wilting_point  #   # Hess Eqn (29)
    print("hx", hx[0], hx[-1], hx.shape)
    q = Kr.dot(hs - hx)
    print("sum_q", np.sum(q), q[0], q[1])
    print()

    # print("hs", hs.shape)
    # print("c_neumann", c_neumann.shape)
    # print("C_comp_neumann", C_comp_neumann.shape)
    # print("C_comp_neumann.dot(hs)", C_comp_neumann.dot(hs).shape)
    q_neumann = C_comp_neumann.dot(hs) + c_neumann * t_act

     # print("c_dirichlet", c_dirichlet.shape)
    q_dirichlet = -(C_comp_dirichlet.dot(hs) + c_dirichlet * wilting_point)

    # print("q", q.shape)
    print("hs", np.min(hs), np.max(hs), hs[0], hs[1])
    # print("hx", np.min(hx), np.max(hx), hx[0], hx[1], hx.shape)
    print("t_act", t_act)
    print("sum(q_neumann)", np.sum(q_neumann), q_neumann[0], q_neumann[1])
    print("sum(q_dirichlet)", np.sum(q_dirichlet), q_dirichlet[0], q_dirichlet[1])

    if np.sum(q_dirichlet) > t_act:
        fluxes = r.sumSegFluxes(q_dirichlet)
    else:
        fluxes = r.sumSegFluxes(q_neumann)
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
