"""
transpiration plot (one column, number of rows as number of filenames)
"""
import sys;
sys.path.append("../../../../CPlantBox/");
sys.path.append("../../../../CPlantBox/src/")
sys.path.append("../")
sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../modules/");
sys.path.append("/../data");
sys.path.append("../../../../CPlantBox/src/functional/");
sys.path.append("../../../../CPlantBox/src/rsml/");
sys.path.append("../../../../CPlantBox/src/visualisation/")
sys.path.append("../../../../CPlantBox/src/structural/")
sys.path.append("../../../../CPlantBox/src/external/")
import numpy as np
import matplotlib.pyplot as plt
from functional.xylem_flux import sinusoidal2
import scenario_setup as scenario
import evapotranspiration as evap
import os

"""scenario"""
year = 2019
soil_type = "loam"
genotype = "WT"
name = "maize_exudate_2019"


os.chdir("../")
print(os.getcwd())
soil_, min_b, max_b, cell_number, area, Kc = scenario.maize_SPP(soil_type)
sim_time = 50 #154   #  [day]
x_, y_, lai = evap.net_infiltration(year, soil_type, genotype, sim_time, Kc)
potential_trans = evap.get_transpiration(year, sim_time, area, lai, Kc)
os.chdir("plotting/")

rs_age = 1

fnames = np.array(["transpiration_" + name ])
fnames2 = np.array(["carbon_" + name])
path = "../results/"

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
data = [np.load(path + n_ + ".npy", allow_pickle=True)  for n_ in fnames]
try:
    data2 = [np.load(path + n_ + ".npy", allow_pickle=True)  for n_ in fnames2]
except:
    data2 = None

fig, ax = plt.subplots(n, 1, figsize = (18, 8))
if n == 1:
    ax = [ax]

for i in range(0, n):
    t = data[i][0]
    y = data[i][1]
    #if trans > 0:
    pt = 10 * np.array([ -potential_trans(t[i], dt_) / area for i in range(0, t.shape[0]) ])
    ax[i].plot(t, pt, 'k', label = "potential transpiration")  # potential transpiration
    y = np.maximum(y, 0)
    ax[i].plot(t, 10 * y / area, 'g', label = "actual transpiration")  # actual transpiration  according to soil model
    #if data2 is not None:
    #    t = data[i][0]
    #    c = data2[i]*-1
    #    ax[i].plot(t, c, 'r:', label = "carbon exudation [g/day]")  # actual transpiration  according to soil model
    #    dt = np.diff(t)
    #    cuc = np.cumsum(np.multiply(c[:-1], dt))
    #    print(str_[i], "cumulative carbon exudation", cuc[-1] , "g", c[0])

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
    #if data2 is not None:
    #    ax2.plot(t[1:],  cuc, 'r--', label = "cumulative carbon exudation [g]")  # cumulative transpiration (neumann)

    #ax2.legend(loc = 'center right')
    print("cumulative water uptake", cup[-1], "cm")

ax[i].set_xlabel("time [day]")

plt.tight_layout()
plt.show()