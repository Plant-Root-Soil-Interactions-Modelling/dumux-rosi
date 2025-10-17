"""
    Dynamic:

    Plots net infiltration and 2D image of soil matric potential or concentration vs time (of a 1d soil)
    in script pick name from L50-
    manually set corresponding net infiltration by setting Kc and lai (see L113)

    Daniel Leitner, 2025  
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import numpy as np
from datetime import *
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import scenario_setup
import evapotranspiration as evap

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


def plot_soil(path, name, show = True):

    ylim_ = None
    soil_min = -6000

    # load data
    alldata = np.load(path + name + ".npz")
    times = alldata["times"]
    pot_trans = alldata["pot_trans"]
    act_trans = alldata["act_trans"]
    depths = alldata["depth"]
    sim_time = np.max(times)
    data = alldata["psi_s"]
    data = np.transpose(data)
    data = data[::-1,:]
    yy = abs(int(depths[-1] - 10))  # depth index for visualization
    data = data[:yy,:]

    # derive net infiltration for setup
    start_date = '2021-05-10 00:00:00'  # INARI csv data
    Kc = 1.15
    lai = evap.lai_soybean2
    t_, y_ = evap.net_infiltration_table_beers_csvS(start_date, sim_time, lai, Kc, initial_age = times[0])
    t_ = np.array(t_)
    y_ = np.array(y_)
    t_ = t_[::2]
    y_ = y_[::2]

    fig, ax = plt.subplots(2, 1, figsize = (18, 10), gridspec_kw = {'height_ratios': [1, 3]})

    # Net infiltration
    bar = ax[0].bar(t_, 10 * np.array(y_), 1 / 24.)
    ax[0].set_ylabel("net inf [mm/day]")
    ax[0].set_xlim(times[0], times[-1])
    if ylim_ is not None:
        ax[0].set_ylim(ylim_, 1.)
    divider = make_axes_locatable(ax[0])
    cax0 = divider.append_axes('right', size = '5%', pad = 0.05)
    cax0.axis('off')

    # Soil
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes('right', size = '5%', pad = 0.05)
    cmap_reversed = matplotlib.cm.get_cmap('jet_r')
    im = ax[1].imshow(data, cmap = cmap_reversed, vmin = soil_min, aspect = 'auto', extent = [times[0] , times[-1], -yy, 0.])  #  interpolation = 'bicubic', interpolation = 'nearest',
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

    print("Soil ranges from", np.min(data), "to", np.max(data), "[cm] matric potential")
    print("Final rs depth", depths[-1])

    plt.tight_layout()
    if not os.path.exists("figures"):
        os.mkdir("figures")
    plt.savefig('figures/' + name + '.png', dpi = 150, bbox_inches = 'tight')
    if show:
        plt.show()


if __name__ == "__main__":

    # """ pick... """  # name = "maize"
    # pick12 = [
    # "soybean_all14_e5288952625e69534d1987146abf0263c78dea24a4f8ebbcfd54b3ec02b33752",
    # "soybean_all14_8583e300050dcde7b7745e83d952833faac5258877ab5607822a1f0e4c639b85",
    # "soybean_all14_2d85b66aacdc462dee8c271e13c52d281aa78b54c32225c310526a7f9e0ec247",
    # "soybean_all14_fb6876837510555b1f91d490edd2707e5024bccd8b0cc89b8156e299648513f6",
    # "soybean_all14_6dcaf82cca09f5c136ca0b9da5c6f539197a85f733211ece6665d9afb69ec4e2",
    # "soybean_all14_6f0184db286042bf501cf953d9f82cbf2c6801c1ada76b8a0ae2911a57cfc189",
    # "soybean_all14_b95151025ffa954589fb8d863d5be5a9110ecd4538f23b5cc77750cbee997ee9",
    # "soybean_all14_a656e46e3de66c57292b1abc890a27cc0712984db218111d14d3ec513319ea70",
    # "soybean_all14_6d244e8c11a6a20ad6d4985c337944e9a074b7ce204cbf27c423cf9abed7973b",
    # "soybean_all14_a817d2cd48002b337996bd68fff0e653b5c9e677fae41fca4bab8f3d84f1205b",
    # "soybean_all14_b4d63e89a73c45e0ef4f338b53b5b1ea82434d1bd9b7de5b71e774f9ef9d5fd2",
    # "soybean_all14_3a7d79c45a73e419323d343f07ee3bec6c631bac462622ac21a73c3823c740d0"
    # ]
    # envirotypes = ["0", "1", "5", "36", "59"]
    # cluster_10_opt = [
    #     "soybean_all14_3564a0f636e9f600fe68bf96ffca4124135106ae4787b9a3332334b04abcdf1a",  # 213
    #     "soybean_all14_f7319dc5d83c72932fd39e4afbf6e50822c2f7bf13b27fc5749c1128642a95d2"  # 138
    #     ]

    # name = pick12[6] + "_36"  # node + envirotype

    # path = "local_soybean_0/"
    # name = "local_soybean_0_001"
    # plot_soil(path, name)

    # path = "local_soybean_1/"
    # name = "local_soybean_1_001"
    # plot_soil(path, name)

    # path = "local_soybean_free_0/"
    # name = "local_soybean_free_0_001"
    # plot_soil(path, name)

    path = "local_soybean_free_1/"
    name = "local_soybean_free_1_001"
    plot_soil(path, name)

