"""
    transpiration plot (one column, number of rows as number of filenames)
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import numpy as np
import matplotlib.pyplot as plt

from functional.xylem_flux import sinusoidal2

import evapotranspiration as evap

Kc_maize = 1.2
Kc_soybean = 1.15

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

name = "soybean_water_0"
name = "local_soybean_noFlux_1"
name = "local_soybean_1"
# name = "local_singleroot_conductivities1_1"
name = "singleroot_test"
name = "soybean_test_1"
name = "local_soybean_new_1"

name = pick12[6] + "_1"
str_ = [""]
area = 76 * 3
# start_date = '2021-05-10 00:00:00'  # INARI csv data
# potential_trans = evap.get_transpiration_beers_csvS(start_date, 87.5, area, evap.lai_soybean2, Kc_soybean)  # 87.5
# evap.net_infiltration_table_beers_csv(start_date, 87.5, evap.lai_soybean2, Kc_soybean)
# trans = 1

# name = "local_soybean"
# str_ = ["1"]
# area = 76 * 3
# # start_date = '1995-03-15 00:00:00'
# # potential_trans = evap.get_transpiration_beers_pickle('data/95.pkl', start_date, 87.5, area, evap.lai_soybean2, Kc_soybean)
# start_date = '2021-05-10 00:00:00'  # INARI csv data
# potential_trans = evap.get_transpiration_beers_csvS(start_date, 87.5, area, evap.lai_soybean2, Kc_soybean)
# trans = 1

# name = "test1.0"
# str_ = [""]
# potential_trans = lambda t, dt: trans * sinusoidal2(t, dt) * t / sim_time  # soybean
# area = 76 * 3
# potential_trans = evap.get_transpiration_beers('data/95.pkl', start_date, 87.5, area, evap.lai_soybean2, Kc_soybean)
# trans = 1

# name = "maize"
# str_ = ["_cyl0"]
# area = 76 * 16  # cm2
# potential_trans = evap.get_transpiration_beers('data/95.pkl', start_date, 95, area, evap.lai_maize2, Kc_maize)
# trans = 1
#
# name = "maize"
# str_ = ["_sra0d"]
# area = 76 * 16  # cm2
# # start_date = '1995-03-15 00:00:00'
# # potential_trans = evap.get_transpiration_beers_pickle('data/95.pkl', start_date, 95, area, evap.lai_maize2,2 Kc_maize)
# start_date = '2021-05-10 00:00:00'  # INARI csv data
# potential_trans = evap.get_transpiration_beers_csvS(start_date, 95, area, evap.lai_maize2, Kc_maize)
# # evap.net_infiltration_table_beers_csv(start_date, 95, evap.lai_maize2, Kc_maize)
# trans = 1

rs_age = 1

# fnames = np.array(["transpiration_" + name + s for s in str_ ])
# fnames2 = np.array(["nitrate_" + name + s for s in str_ ])
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

""" transpiration plot """
dt_ = 360 / (24 * 3600)

# # load data
# n = len(fnames)
# data = [np.load(path + n_ + ".npy")  for n_ in fnames]
# try:
#     data2 = [np.load(path + n_ + ".npy")  for n_ in fnames2]
# except:
#     data2 = None

alldata = np.load(path + name + ".npz")
times = alldata["times"]
pot_trans = alldata["pot_trans"]
act_trans = alldata["act_trans"]
print(times)
print(pot_trans)
print(act_trans)
# plt.plot(times, pot_trans)
# plt.plot(times, act_trans)
# plt.show()

n = 1

fig, ax = plt.subplots(n, 1, figsize = (18, 8))
if n == 1:
    ax = [ax]

for i in range(0, n):
    t = times
    ax[i].plot(t, -10 * pot_trans / area, 'k', label = "potential transpiration")  # potential transpiration
    # y = np.maximum(y, 0)
    ax[i].plot(t, -10 * act_trans / area, 'g', label = "actual transpiration")  # actual transpiration  according to soil model
    # if data2 is not None:
    #     t = data[i][0]
    #     c = data2[i]
    #     ax[i].plot(t, 10 * c, 'r:', label = "nitrate uptake 1.e-2*[g/day]")  # actual transpiration  according to soil model
    #     dt = np.diff(t)
    #     cuc = np.cumsum(np.multiply(c[:-1], dt))
    #     print(str_[i], "cumulative nitrate uptake", cuc[-1] , "g", c[0])

    ax[i].set_ylabel("transpiration [mm/day]")
    ax[i].legend(loc = 'upper left')
    ax2 = ax[i].twinx()
    dt = np.diff(t)
    cup = -10 * np.cumsum(np.multiply(act_trans[:-1], dt)) / area
    cum_pot = -10 * np.cumsum(np.multiply(pot_trans[:-1], dt)) / area
    ax2.plot(t[1:], cup, 'c--', label = "cumulative transpiration")  # cumulative transpiration (neumann)
    # ax2.plot(t[1:], cum_pot, 'c--', label = "cumulative potential")
    print("cumulative pot transp", cum_pot[-1])
    ax2.set_ylabel("cumulative [mm]")
    # if data2 is not None:
    #     ax2.plot(t[1:], 10 * cuc, 'r--', label = "cumulative nitrate uptake 1.e-2*[g]")  # cumulative transpiration (neumann)

    ax2.legend(loc = 'center right')
    print("cumulative water uptake", cup[-1], "cm", -10 * np.sum(np.multiply(act_trans[:-1], dt)) / area)

ax[i].set_xlabel("time [day]")

plt.tight_layout()
plt.show()
