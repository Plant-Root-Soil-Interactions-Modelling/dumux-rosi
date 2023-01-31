"""
transpiration plot (one column, number of rows as number of filenames)
"""
import sys; sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

import numpy as np
import matplotlib.pyplot as plt

from xylem_flux import sinusoidal2
import evapotranspiration as evap

Kc_maize = 1.2
Kc_soybean = 1.15

# name = "soybean"
# str_ = ["_sra0d"]
# area = 76 * 3
# # start_date = '1995-03-15 00:00:00'
# # potential_trans = evap.get_transpiration_beers_pickle('data/95.pkl', start_date, 87.5, area, evap.lai_soybean, Kc_soybean)
# start_date = '2021-05-10 00:00:00'  # INARI csv data
# potential_trans = evap.get_transpiration_beers_csvS(start_date, 87.5, area, evap.lai_soybean2, Kc_soybean)  # 87.5
# evap.net_infiltration_table_beers_csv(start_date, 87.5, evap.lai_soybean2, Kc_soybean)
#
# trans = 1
#
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
name = "maize"
str_ = ["_sra0d"]
area = 76 * 16  # cm2
# start_date = '1995-03-15 00:00:00'
# potential_trans = evap.get_transpiration_beers_pickle('data/95.pkl', start_date, 95, area, evap.lai_maize2,2 Kc_maize)
start_date = '2021-05-10 00:00:00'  # INARI csv data
potential_trans = evap.get_transpiration_beers_csvS(start_date, 95, area, evap.lai_maize2, Kc_maize)
# evap.net_infiltration_table_beers_csv(start_date, 95, evap.lai_maize2, Kc_maize)
trans = 1

rs_age = 1

fnames = np.array(["transpiration_" + name + s for s in str_ ])
fnames2 = np.array(["nitrate_" + name + s for s in str_ ])
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

# load data
n = len(fnames)
data = [np.load(path + n_ + ".npy")  for n_ in fnames]
try:
    data2 = [np.load(path + n_ + ".npy")  for n_ in fnames2]
except:
    data2 = None

fig, ax = plt.subplots(n, 1, figsize = (18, 8))
if n == 1:
    ax = [ax]

for i in range(0, n):
    t = data[i][0]
    y = data[i][1]
    if trans > 0:
        pt = 10 * np.array([ -potential_trans(t[i], dt_) / area for i in range(0, t.shape[0]) ])
        ax[i].plot(t, pt, 'k', label = "potential transpiration")  # potential transpiration
    y = np.maximum(y, 0)
    ax[i].plot(t, 10 * y / area, 'g', label = "actual transpiration")  # actual transpiration  according to soil model
    if data2 is not None:
        t = data[i][0]
        c = data2[i]
        ax[i].plot(t, 10 * c, 'r:', label = "nitrate uptake 1.e-2*[g/day]")  # actual transpiration  according to soil model
        dt = np.diff(t)
        cuc = np.cumsum(np.multiply(c[:-1], dt))
        print(str_[i], "cumulative nitrate uptake", cuc[-1] , "g", c[0])

    ax[i].set_ylabel("transpiration [mm/day]")
    ax[i].legend(loc = 'upper left')
    ax2 = ax[i].twinx()
    dt = np.diff(t)
    so = np.array(y) / area
    cup = np.cumsum(np.multiply(so[:-1], dt))
    cum_pot = np.cumsum(np.multiply(pt[:-1], dt))
    ax2.plot(t[1:], 10 * cup, 'c--', label = "cumulative transpiration")  # cumulative transpiration (neumann)
    print("cumulative pot transp", cum_pot[-1])
    ax2.set_ylabel("cumulative [mm]")
    if data2 is not None:
        ax2.plot(t[1:], 10 * cuc, 'r--', label = "cumulative nitrate uptake 1.e-2*[g]")  # cumulative transpiration (neumann)

    ax2.legend(loc = 'center right')
    print(str_[i], "cumulative water uptake", cup[-1], "cm")

ax[i].set_xlabel("time [day]")

plt.tight_layout()
plt.show()
