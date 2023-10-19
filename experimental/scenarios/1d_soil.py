""" 
Jan's new scenario

with the classical sink

complicated MPI support (a non-mpi version of richards_cyl is needed, see script dumux3_nompi.sh)
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

name = "small_soil"  # name to export resutls
sstr = "_dry"

r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario1D(sstr)

""" Numerical solution  """
NT = int(np.ceil(sim_time / dt))  # number of iterations

start_time = timeit.default_timer()

sx = s.getSolutionHead()

x_, y_, w_, cf_ = [], [], [], []
psi_x_, psi_s_, sink_, psi_s2_ = [], [], [], []  # outputs

for i in range(0, NT):

     t = i * dt  # current simulation time

     """ 1. xylem model """
     rx = r.solve(rs_age, -trans * sinusoidal(t), 0., sx, cells = True, wilting_point = wilting_point)  # xylem_flux.py
     fluxes = r.soilFluxes(rs_age, rx, sx, False)  # class XylemFlux is defined in MappedOrganism.h, approx = False

     """ 2. soil model """
     s.setSource(fluxes.copy())  # richards.py
     s.solve(dt)
     sx = s.getSolutionHead()  # richards.py

     """ remember results ... """
     if i % skip == 0:
         sx_ = sx[:, 0]
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
