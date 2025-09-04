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
pick12 = [
"soybean_all14_e5288952625e69534d1987146abf0263c78dea24a4f8ebbcfd54b3ec02b33752",
"soybean_all14_8583e300050dcde7b7745e83d952833faac5258877ab5607822a1f0e4c639b85",
"soybean_all14_2d85b66aacdc462dee8c271e13c52d281aa78b54c32225c310526a7f9e0ec247",
"soybean_all14_fb6876837510555b1f91d490edd2707e5024bccd8b0cc89b8156e299648513f6",
"soybean_all14_6dcaf82cca09f5c136ca0b9da5c6f539197a85f733211ece6665d9afb69ec4e2",
"soybean_all14_6f0184db286042bf501cf953d9f82cbf2c6801c1ada76b8a0ae2911a57cfc189",
"soybean_all14_b95151025ffa954589fb8d863d5be5a9110ecd4538f23b5cc77750cbee997ee9",
"soybean_all14_a656e46e3de66c57292b1abc890a27cc0712984db218111d14d3ec513319ea70",
"soybean_all14_6d244e8c11a6a20ad6d4985c337944e9a074b7ce204cbf27c423cf9abed7973b",
"soybean_all14_a817d2cd48002b337996bd68fff0e653b5c9e677fae41fca4bab8f3d84f1205b",
"soybean_all14_b4d63e89a73c45e0ef4f338b53b5b1ea82434d1bd9b7de5b71e774f9ef9d5fd2",
"soybean_all14_3a7d79c45a73e419323d343f07ee3bec6c631bac462622ac21a73c3823c740d0"
]
envirotypes = ["0", "1", "5", "36", "59"]
cluster_10_opt = [
    "soybean_all14_3564a0f636e9f600fe68bf96ffca4124135106ae4787b9a3332334b04abcdf1a",  # 213
    "soybean_all14_f7319dc5d83c72932fd39e4afbf6e50822c2f7bf13b27fc5749c1128642a95d2"  # 138
    ]

# str_ = "_sra0"
# Kc = Kc_maize
# lai = evap.lai_maize2
# ylim_ = None  # -10

name = "local_soybean_1"
name = "soybean_testsingle_0"
name = "local_singleroot_conductivities64_27"
name = "singleroot_test"
name = "local_soybean_radii_1"
name = "soybean_test_0"

name = pick12[6] + "_36"  # node + envirotype
name = cluster_10_opt[0] + "_1"

str_ = ""
Kc = Kc_soybean
lai = evap.lai_soybean2
ylim_ = None

soil_min = -6000

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

print("range", np.min(data), np.max(data), "[cm]")
print("final rs depth", depths[-1])
plt.tight_layout()
plt.savefig('soil' + name[-5:] + '.png')
plt.show()

