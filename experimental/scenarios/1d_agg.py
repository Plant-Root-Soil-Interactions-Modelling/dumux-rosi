""" 
Jan's new scenarios with aggregated (cut) root system and SRA fix point iteration

see aggregated rs
"""

import sys; sys.path.append("../.../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from scenario_setup import *
import aggregated_rs as agg

import timeit
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy import sparse
import scipy.sparse.linalg as LA


def double_(rsx):
    rsx2 = np.array([ 0. if i % 2 == 0 else rsx[int(i / 2)] for i in range(0, ns)])
    return rsx2

""" 
Initialize  
"""

name = "small_agg"  # name to export resutls
sstr = "_dry"

r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario1D(sstr)
min_b, max_b, cell_number = get_domain1D()

dt = 60 / (24 * 3600)  # time step [day], 120 schwankt stark
skip = 6  # for output and results, skip iteration

""" 
Initialize aggregated hydraulic model 
"""
r = agg.create_aggregated_rs(r, rs_age, min_b, max_b, cell_number)
nodes = r.rs.nodes
picker = lambda x, y, z: s.pick([0., 0., z])  # reset mapper, since it is 1D
r.rs.setSoilGrid(picker)  # maps segment

""" for fixed mapping """
types = r.rs.subTypes
outer_r = r.rs.segOuterRadii()
inner_r = r.rs.radii
rho_ = np.divide(outer_r, np.array(inner_r))
ns = len(outer_r)
mapping = np.array([r.rs.seg2cell[j] for j in range(0, ns)])  # redo mapping

""" Numerical solution """
NT = int(np.ceil(sim_time / dt))  # number of iterations

start_time = timeit.default_timer()

sx = s.getSolutionHead()  # inital condition, solverbase.py
cell_centers = s.getCellCenters()
hsb = np.array([sx[mapping[2 * j + 1]] for j in range(0, int(ns / 2))])  # soil bulk matric potential per segment
cell_centers_z = np.array([cell_centers[mapping[2 * j + 1]][2] for j in range(0, int(ns / 2))])
kr_ = np.zeros((ns,))
rsx = hsb.copy()  # initial values for fix point iteration
t = 0.
rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., double_(rsx), False, wilting_point, soil_k = [])
rx_old = rx.copy()

x_, y_, w_, cf_ = [], [], [], []
psi_x_, psi_s_, sink_, psi_s2_ = [], [], [], []  # outputs

for i in range(0, NT):

    t = i * dt  # current simulation time

    wall_iteration = timeit.default_timer()
    wall_fixpoint = timeit.default_timer()

    for j in range(0, len(outer_r)):  # determine kr at this time step
        kr_[j] = r.kr_f(rs_age + t, types[j], 2, 2, j)
    inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up

    err = 1.e6
    c = 1
    while err > 1 and c < 100:

        """ interpolation """
        # wall_interpolation = timeit.default_timer()
        rx_ = rx[1::2] - np.array([nodes[ii].z for ii in range(1, ns, 2)])  # from total matric potential to matric potential
        hsb_ = hsb - cell_centers_z  # from total matric potential to matric potential
        rsx = soil_root_interface_table2(rx_, hsb_, inner_kr_[1::2], rho_[1::2], sra_table_lookup)  # [1::2] every second entry, starting from 1
        rsx = rsx + np.array([nodes[ii].z for ii in range(1, ns, 2)])  # from matric potential to total matric potential
        # wall_interpolation = timeit.default_timer() - wall_interpolation

        """ xylem matric potential """
        # wall_xylem = timeit.default_timer()
        rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., double_(rsx), False, wilting_point, soil_k = [])  # xylem_flux.py, cells = False
        err = np.linalg.norm(rx - rx_old)
        # wall_xylem = timeit.default_timer() - wall_xylem
        rx_old = rx.copy()
        c += 1
#         print(c, ": ", np.sum(rx[1:]), np.sum(hsb), np.sum(inner_kr_), np.sum(rho_))
#        print(c, "iterations", wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem))

#    wall_fixpoint = timeit.default_timer() - wall_fixpoint

    fluxes = r.segFluxes(rs_age + t, rx, double_(rsx), approx = False, cells = False)

#     min_rsx = np.min(rsx)  # for console output
#     max_rsx = np.max(rsx)
#     print("from", min_rsx, "to", max_rsx)

    wall_soil = timeit.default_timer()
    soil_fluxes = r.sumSegFluxes(fluxes)
    s.setSource(soil_fluxes.copy())  # richards.py
    s.solve(dt)
    sx = s.getSolutionHead()  # richards.py
    hsb = np.array([sx[mapping[2 * j + 1]][0] for j in range(0, int(ns / 2))])
    wall_soil = timeit.default_timer() - wall_soil

    wall_iteration = timeit.default_timer() - wall_iteration

    """ remember results ... """
    if i % skip == 0:
        print(i / skip)
        x_.append(t)
        psi_x_.append(rx[1:])
        psi_s_.append(rsx.copy())
        sum_flux = 0.
        sink = np.zeros(sx[:, 0].shape)
        for k, v in soil_fluxes.items():
            sink[k] += v
            sum_flux += v
        y_.append(sum_flux)  # cm3/day
        sink_.append(sink)
        dd = np.array(s.getSolutionHead())
        psi_s2_.append(dd[:, 0])

print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

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

print("fin.")
