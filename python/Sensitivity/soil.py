""" 
    Macroscopic soil model 
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver

# from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

import scenario_setup as scenario
from xylem_flux import sinusoidal2
import evapotranspiration as evap

""" parameters   """
min_b = [-38, -7.5, -200.]  # Domain Mais: 60 cm Reihe, 10 cm Pflanzen
max_b = [38, 7.5, 0.]
cell_number = [1, 1, 200]  # 1 cm3

soil0 = [0.0809, 0.52, 0.0071, 1.5734, 99.49]
soil1 = [0.0874, 0.5359, 0.0087, 1.5231, 93]
soil36 = [0.0942, 0.5569, 0.0089, 1.4974, 87.79]
soil5 = [0.0539, 0.5193, 0.024, 1.4046, 208.78]
soil59 = [0.0675, 0.5109, 0.0111, 1.4756, 107.63]
soil_ = soil0

Kc_maize = 1.2  # book "crop evapotranspiration" Allen, et al 1998

area = 76 * 15  # cm2

sim_time = 95  #  [day]
dt = 360 / (24 * 3600)  # time step [day] 20

# start_date = '1995-03-15 00:00:00'
# x_, y_ = evap.net_infiltration_table_beers_pickle('data/95.pkl', start_date, sim_time, evap.lai_maize, Kc_maize)
# trans_maize = evap.get_transpiration_beers_pickle('data/95.pkl', start_date, sim_time, area, evap.lai_maize, Kc_maize)

start_date = '2020-03-15 00:00:00'  # INARI csv data
x_, y_ = evap.net_infiltration_table_beers_csv(start_date, sim_time, evap.lai_maize, Kc_maize)  # evap.lai_noroots

""" set up simulator """
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = -330, p_bot = -130, type = 2, times = x_, net_inf = y_)  # , times = x_, net_inf = y_

""" simulation loop """
c = []  # resulting solute concentration
N = int(np.ceil(sim_time / dt))

for i in range(0, N):
    t = i * dt  # current simulation time
    print(t)
    soil_sol_fluxes = {}  # empy dict
    evap.add_nitrificatin_source(s, soil_sol_fluxes, nit_flux = 1.e-7 * (75 * 16 * 1))
    s.setSource(soil_sol_fluxes.copy(), eq_idx = 1)  # richards.py
    s.solve(dt)
    c.append(s.getSolution_(1))

z = s.getDofCoordinates()

""" nitrate plot """
c = np.transpose(c)
c = c[::-1,:]
c = c[:150]
# c = np.minimum(c, 5.5e-4)

fig, ax = plt.subplots(1, 1, figsize = (18, 10))
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size = '5%', pad = 0.05)

cmap_ = matplotlib.cm.get_cmap('jet')
im = ax.imshow(c, cmap = cmap_, aspect = 'auto', vmin = 0., vmax = 1.e-3, extent = [0 , sim_time, -150., 0.])  #  interpolation = 'bicubic', interpolation = 'nearest', vmin = 0., vmax = 1.e-3,

cb = fig.colorbar(im, cax = cax, orientation = 'vertical')
cb.ax.get_yaxis().labelpad = 30
cb.set_label('Concentration [kg/m3]', rotation = 270)

ax.set_ylabel("depth [cm]")
ax.set_xlabel("time [days]")

print("range", np.min(c), np.max(c), "cm")

plt.tight_layout()
plt.show()
