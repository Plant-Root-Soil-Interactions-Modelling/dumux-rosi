"""
total nitrate in soil over time
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

import van_genuchten as vg
import evapotranspiration as evap
import scenario_setup as scenario

from datetime import *

soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.maize(0)

""" pick... """

name = "maize"
str_ = "_sra2"
lai = evap.lai_maize
sim_time = 95

# name = "maize"
# str_ = "_cyl0"
# Kc = Kc_maize
# lai = evap.lai_maize
# ylim_ = None # -10

fname1 = "soil_" + name + str_
fname2 = "soilc_" + name + str_

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
water = np.load(path + fname1 + ".npy")
nitrate = np.load(path + fname2 + ".npy")
t_ = np.linspace(1., sim_time, water.shape[0])
z_ = np.linspace(-200, 0, 200)

s = vg.Parameters(soil_)
c_ = np.zeros((water.shape[0],))
for i in range(0, water.shape[0]):
    theta = vg.water_content(water[i,:], s)  ################################################################# -z_ has no effect
    c_[i] = np.sum(np.multiply(theta, nitrate[i,:]) * 1. * area)

plt.plot(t_, c_)
# plt.plot(water[0,:], z_)
plt.xlabel("time  [days]")
plt.ylabel("nitrate [g]")
plt.show()

