"""
2d image of soil matric potential or concentration vs time (of a 1d soil)
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

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

# name = "soybean"
# str_ = "_sra0w"
# Kc = Kc_soybean
# lai = evap.lai_soybean2
# ylim_ = None
# sim_time = 87.5

# name = "local_soybean"
# str_ = "1"
# Kc = Kc_soybean
# lai = evap.lai_soybean
# ylim_ = None
# sim_time = 87.5
#
# name = "maize"
# str_ = "_sra100"
# Kc = Kc_maize
# lai = evap.lai_maize
# ylim_ = None
# sim_time = 95

name = "maize"
str_ = "_sra0d"
Kc = Kc_maize
lai = evap.lai_maize2
ylim_ = None  # -10
sim_time = 95

fname = "soil_" + name + str_

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
fnames_t = np.array(["transpiration_" + name + str_ ])  # data currently not stored for rhizosphere models
times = [np.load(path + n_ + ".npy")  for n_ in fnames_t]
times = times[0][0,:]

# t_, y_ = evap.net_infiltration_table_beers('data/95.pkl', start_date, times[-1], lai, Kc)
t_, y_ = evap.net_infiltration_table_beers_csvS(start_date, sim_time, lai, Kc)
# t_, y_ = evap.net_infiltration_table_beers('data/95.pkl', start_date, times[-1], lai, Kc)
t_ = np.array(t_)
y_ = np.array(y_)
t_ = t_[::2]
y_ = y_[::2]

fnames_depth = np.array(["depth_" + name + str_ ])
depths = [np.load(path + n_ + ".npy")  for n_ in fnames_depth]
depths = depths[0]

data = np.load(path + fname + ".npy")
data = np.transpose(data)
data = data[::-1,:]
data = data[:abs(int(depths[-1] - 10)),:]
# data = np.minimum(data, 0.001)

""" sink plot """
if fname.startswith("soilc_"):
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
    cmap = matplotlib.cm.get_cmap('jet')
    im = ax[1].imshow(data, vmin = 0., vmax = 1.e-3, cmap = cmap, aspect = 'auto', extent = [times[0] , times[-1], depths[-1] - 10, 0.])  #  interpolation = 'bicubic', interpolation = 'nearest',
    ax[1].plot(times[::10], depths, 'k:')
    cb = fig.colorbar(im, cax = cax, orientation = 'vertical')
    cb.ax.get_yaxis().labelpad = 30
    cb.set_label('nitrate concentration [kg/m3]', rotation = 270)
    ax[1].set_ylabel("depth [cm]")
    ax[1].set_xlabel("time [days]")
    # ax[1].scatter([16, 17, 29, 30], [0, 0, 0, 0], [30, 30, 60, 60], color = 'w')
    # ax[1].scatter([1, 18, 54],
    #               [depths[-1] - 10] * 3, 3 * np.array([40] * 3), color = 'k')  # sol_times = np.array([0., 1., 1., 17., 17., 18. , 18., 53., 53, 54, 54., 1.e3]) [0, 0, 0, 0, 0, 0]
    print("data ranges from", np.min(data), "to ", np.max(data), "[kg/m3]'")
    # print("final rs depth", depths[-1])
    plt.tight_layout()
    plt.show()

else:
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
    im = ax[1].imshow(data, cmap = cmap_reversed, aspect = 'auto', extent = [times[0] , times[-1], depths[-1] - 10, 0.])  #  interpolation = 'bicubic', interpolation = 'nearest',
    ax[1].plot(times[::10], depths, 'k:')
    cb = fig.colorbar(im, cax = cax, orientation = 'vertical')
    cb.ax.get_yaxis().labelpad = 30
    cb.set_label('soil matric potential [cm]', rotation = 270)
    ax[1].set_ylabel("depth [cm]")
    ax[1].set_xlabel("time [days]")

    print("range", np.min(data), np.max(data), "[cm]")
    print("final rs depth", depths[-1])
    plt.tight_layout()
    plt.show()

