""" 
Jan's new scenarios with the new SRA fix point iteration sink

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

name = "small_sra"  # name to export resutls
sstr = "_dry"

r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario1D(sstr)

""" for fixed root system """
types = r.rs.subTypes
inner_ = r.rs.radii
outer_ = r.rs.segOuterRadii()
rho_ = np.divide(outer_, np.array(inner_))
ns = len(types)

""" Numerical solution  """
NT = int(np.ceil(sim_time / dt))  # number of iterations

start_time = timeit.default_timer()

sx = s.getSolutionHead()  # inital condition, solverbase.py
hsb = np.array([sx[mapping[j]][0] for j in range(0, ns)])  # soil bulk matric potential per segment
kr_ = np.zeros((ns,))
rsx = hsb.copy()  # initial values for fix point iteration
t = 0
rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, [])  # xylem_flux.py, cells = False
rx_old = rx.copy()

x_, y_, w_, cf_ = [], [], [], []
psi_x_, psi_s_, sink_, psi_s2_ = [], [], [], []  # outputs

for i in range(0, NT):

    t = i * dt  # current simulation time

    wall_iteration = timeit.default_timer()

    wall_fixpoint = timeit.default_timer()
    err = 1.e6
    c = 1

    while err > 1. and c < 100:

        """ interpolation """
        wall_interpolation = timeit.default_timer()
        for j in range(0, ns):  # determine kr at this time step
            kr_[j] = r.kr_f(rs_age + t, types[j])
        inner_kr_ = np.multiply(inner_, kr_)  # multiply for table look up
        rsx = soil_root_interface_table2(rx[1:], hsb, inner_kr_, rho_, sra_table_lookup)
        wall_interpolation = timeit.default_timer() - wall_interpolation

        """ xylem matric potential """
        wall_xylem = timeit.default_timer()
        rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, [])  # xylem_flux.py, cells = False
        err = np.linalg.norm(rx - rx_old)
        wall_xylem = timeit.default_timer() - wall_xylem
        rx_old = rx.copy()
        c += 1
    # print(c, "iterations", wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem))

    wall_fixpoint = timeit.default_timer() - wall_fixpoint

    wall_soil = timeit.default_timer()
    seg_fluxes = r.segFluxes(rs_age + t, rx, rsx, False)  # class XylemFlux is defined in MappedOrganism.h, approx = True
    fluxes = r.sumSegFluxes(seg_fluxes)  # Soil part runs parallel
    s.setSource(fluxes.copy())  # richards.py
    s.solve(dt)
    sx = s.getSolutionHead()  # richards.py
    hsb = np.array([sx[mapping[j]][0] for j in range(0, ns)])  # soil bulk matric potential per segment
    water = s.getWaterVolume()
    wall_soil = timeit.default_timer() - wall_soil

    wall_iteration = timeit.default_timer() - wall_iteration

    if i % skip == 0:  # remember results
        x_.append(t)
        psi_x_.append(rx[1:])
        psi_s_.append(rsx.copy())
        sum_flux = 0.
        sink = np.zeros(sx[:, 0].shape)
        for k, v in fluxes.items():
            sink[k] += v
            sum_flux += v
        y_.append(sum_flux)  # cm3/day
        sink_.append(sink)  # cm3/day
        psi_s2_.append(sx[:, 0])
        w_.append(water)  # cm3

        cf = r.collar_flux(rs_age + t, rx, rsx, k_soil = [], cells = False)
        print("Summed fluxes ", sum_flux, "= collar flux", cf_, "= prescribed", -trans * sinusoidal(t))
        cf_.append(cf)  # cm3/day
        n = round(float(i) / float(NT) * 100.)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil sx [{:g}, {:g}], interface [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days"
              .format(np.min(sx), np.max(sx), np.min(rsx), np.max(rsx), np.min(rx), np.max(rx), s.simTime))
        print("Iteration {:g} took {:g} seconds [{:g} fixpoint iteration, {:g} soil] \n".format(i, wall_iteration, wall_fixpoint, wall_soil))

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

