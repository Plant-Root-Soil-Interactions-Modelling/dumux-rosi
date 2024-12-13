"""
2d image of soil matric potential or concentration vs time (of a 1d soil), withg net infiltration
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

import evapotranspiration as evap

from datetime import *

Kc_maize = 1.2  # book "crop evapotranspiration" Allen, et al (1998)
Kc_soybean = 1.15  # book "crop evapotranspiration" Allen, et al (1998)

""" pick... """  # name = "maize"
# str_ = "_sra0"
# Kc = Kc_maize
# lai = evap.lai_maize2
# ylim_ = None  # -10

name = "local_soybean_1"
name = "soybean_testsingle_0"
name = "local_singleroot_conductivities64_27"
name = "singleroot_test"
str_ = ""
Kc = Kc_soybean
lai = evap.lai_soybean2
ylim_ = None

# name = "local_soybean"
# str_ = "1"
# Kc = Kc_soybean
# lai = evap.lai_soybean
# ylim_ = None
# name = "maize"
# str_ = "_sra100"
# Kc = Kc_maize
# lai = evap.lai_maize
# ylim_ = None

# name = "maize"
# str_ = "_sra0d"
# Kc = Kc_maize
# lai = evap.lai_maize2
# ylim_ = None  # -10

fname = name + str_

# start_date = '1995-03-14 00:00:00'  # substract 1 day, since inital rs_age
start_date = '2021-05-10 00:00:00'  # INARI csv data

path = "results/"

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

""" load data """
alldata = np.load(path + name + ".npz")
times = alldata["times"]
pot_trans = alldata["pot_trans"]
act_trans = alldata["act_trans"]
sim_time = np.max(times)

# t_, y_ = evap.net_infiltration_table_beers('data/95.pkl', start_date, times[-1], lai, Kc)
t_, y_ = evap.net_infiltration_table_beers_csvS(start_date, sim_time, lai, Kc, initial_age = times[0])
# t_, y_ = evap.net_infiltration_table_beers('data/95.pkl', start_date, times[-1], lai, Kc)
t_ = np.array(t_)
y_ = np.array(y_)
t_ = t_[::2]
y_ = y_[::2]

fname_ = name + str_

depths = alldata["depth"]

data = alldata["psi_s"]  # np.load(path + fname + ".npy")
data = np.transpose(data)
data = data[::-1,:]
yy = abs(int(depths[-1] - 10))  # depth index for visualization
data = data[:yy,:]

""" soil plot """
fig, ax = plt.subplots(2, 1, figsize = (18, 10), gridspec_kw = {'height_ratios': [1, 3]})

bar = ax[0].bar(t_, 10 * np.array(y_), 1 / 24.)
ax[0].set_ylabel("net inf [mm/day]")
ax[0].set_xlim(times[0], times[-1])
if ylim_ is not None:
    ax[0].set_ylim(ylim_, 1.)
divider = make_axes_locatable(ax[0])
cax0 = divider.append_axes('right', size = '5%', pad = 0.05)
cax0.axis('off')

divider = make_axes_locatable(ax[1])
cax = divider.append_axes('right', size = '5%', pad = 0.05)
cmap_reversed = matplotlib.cm.get_cmap('jet_r')
print("data", data.shape)
im = ax[1].imshow(data, cmap = cmap_reversed, vmin = -10000, aspect = 'auto', extent = [times[0] , times[-1], -yy, 0.])  #  interpolation = 'bicubic', interpolation = 'nearest',
ax[1].plot(times[::10], depths, 'k:')
x = np.linspace(0, times[-1], data.shape[1])
y = np.linspace(0, -yy, data.shape[0])
X, Y = np.meshgrid(x, y)
contours = ax[1].contour(X, Y, data, [0.], colors = 'black')
cb = fig.colorbar(im, cax = cax, orientation = 'vertical')
cb.ax.get_yaxis().labelpad = 30
cb.set_label('soil matric potential [cm]', rotation = 270)
ax[1].set_ylabel("depth [cm]")
ax[1].set_xlabel("time [days]")

print("range", np.min(data), np.max(data), "[cm]")
print("final rs depth", depths[-1])
plt.tight_layout()
plt.show()

