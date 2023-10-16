""" 
    Macroscopic soil model 
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver

# from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from datetime import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

import scenario_setup as scenario
from xylem_flux import sinusoidal2
import evapotranspiration as evap

""" parameters   """
min_b = [-38, -7.5, -200.]  # Domain Mais: 60 cm Reihe, 10 cm Pflanzen
max_b = [38, 7.5, 0.]
cell_number = [1, 1, 200]  # 1 cm3

soil[0] = [0.0639, 0.3698, 0.0096, 1.4646, 4.47]
soil[1] = [0.0619, 0.3417, 0.0132, 1.3258, 2.03]
soil[36] = [0.0760, 0.3863, 0.0091, 1.4430, 2.99]
soil[5] = [ 0.0451, 0.3721, 0.0325, 1.4393, 34.03]
soil[59] = [0.0534, 0.3744, 0.0171, 1.4138, 13.09]
soil_ = soil[0]

Kc_maize = 1.2  # book "crop evapotranspiration" Allen, et al 1998

area = 76 * 15  # cm2

sim_time = 10  # 95 + 17  #  [day]
dt = 360 / (24 * 3600)  # time step [day] 20

start_date = '1995-03-23 00:00:00'
start_date = datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
end_date = start_date + timedelta(sim_time, 0)
range_ = [str(start_date), str(end_date)]
print("calculating from", range_[0], "to", range_[1])

x_, y_ = evap.net_infiltration_table_beers('data/95.pkl', range_, sim_time, evap.lai_maize, Kc_maize)
t_ = np.array(x_)
y_ = np.array(y_)
# t_ = t_[::2]
# y_ = y_[::2]

""" set up simulator """
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = -330, p_bot = -130, type = 2, times = x_, net_inf = y_)

""" simulation loop """
c = []  # resulting solute concentration
N = int(np.ceil(sim_time / dt))

for i in range(0, N):
    t = i * dt  # current simulation time
    if i % 24 == 0:
        print(t)
    soil_sol_fluxes = {}  # empy dict
    evap.add_nitrificatin_source(s, soil_sol_fluxes, nit_flux = 1.e-6 * area)
    # print(soil_sol_fluxes)
    s.setSource(soil_sol_fluxes.copy(), eq_idx = 1)  # richards.py
    s.solve(dt)
    c.append(s.getSolutionHead_())

z = s.getDofCoordinates()

""" nitrate plot """
c = np.transpose(c)
c = c[:100:-1,:]
times = np.linspace(0., sim_time, N)

fig, ax = plt.subplots(2, 1, figsize = (18, 10), gridspec_kw = {'height_ratios': [1, 3]})
bar = ax[0].bar(t_, -np.array(y_), 0.035)
ax[0].set_ylabel("net infiltration [cm/day]")
ax[0].set_xlim(times[0], times[-1])
# if ylim_ is not None:
#     ax[0].set_ylim(ylim_, 1.)
divider = make_axes_locatable(ax[0])
cax0 = divider.append_axes('right', size = '5%', pad = 0.05)
cax0.axis('off')
divider = make_axes_locatable(ax[1])
cax = divider.append_axes('right', size = '5%', pad = 0.05)
cmap = matplotlib.cm.get_cmap('jet')
im = ax[1].imshow(c, cmap = cmap, aspect = 'auto', extent = [times[0] , times[-1], -100, 0.])  #  interpolation = 'bicubic', interpolation = 'nearest',

cb = fig.colorbar(im, cax = cax, orientation = 'vertical')
cb.ax.get_yaxis().labelpad = 30
cb.set_label('soil matric potential [cm]', rotation = 270)
ax[1].set_ylabel("depth [cm]")
ax[1].set_xlabel("time [days]")
# ax[1].scatter([0, 1, 17, 18, 53, 54], [-100, -100, -100, -100, -100, -100], np.array([40, 40, 40, 40, 40, 40]), color = 'k')  # sol_times = np.array([0., 1., 1., 17., 17., 18. , 18., 53., 53, 54, 54., 1.e3]) [0, 0, 0, 0, 0, 0]
# ax[1].scatter([17.], [-100], np.array([40]), color = 'r')
# ax[1].plot([17., 17.], [0., -100], 'k:')
print("range", np.min(c), np.max(c), "[cm]")
plt.tight_layout()
plt.show()

