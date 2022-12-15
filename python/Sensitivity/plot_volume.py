"""
volume, surface, depth plot 
"""
import sys; sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

import numpy as np
import matplotlib.pyplot as plt

from xylem_flux import sinusoidal2
import evapotranspiration as evap

Kc_maize = 1.2
Kc_soybean = 1.15

name = "soybean"
str_ = ["_sra0"]

fnames_t = np.array(["transpiration_" + name + s for s in str_ ])
fnames = np.array(["vol_" + name + s for s in str_ ])
fnames2 = np.array(["surf_" + name + s for s in str_ ])
fnames3 = np.array(["krs_" + name + s for s in str_ ])
fnames4 = np.array(["depth_" + name + s for s in str_ ])
path = "results/"

SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title

""" volume surface depth plot """
dt_ = 360 / (24 * 3600)

# load data
times = [np.load(path + n_ + ".npy")  for n_ in fnames_t]
t = times[0][0][::10]

data = [np.load(path + n_ + ".npy")  for n_ in fnames]
data2 = [np.load(path + n_ + ".npy")  for n_ in fnames2]
data3 = [np.load(path + n_ + ".npy")  for n_ in fnames3]
data4 = [np.load(path + n_ + ".npy")  for n_ in fnames4]

fig, ax = plt.subplots(3, 1, figsize = (8, 18))

for i in range(0, data[0].shape[0]):  # Volume
    y = data[0][i,:]
    ax[0].plot(t, y, label = "subType " + str(i))
    ax[0].set_ylabel("Volume [cm3]")
yt = np.sum(data[0], axis = 0)
ax[0].plot(t, yt, label = "total")
print("total volume", yt[-1], "cm3")
ax[0].legend()

for i in range(0, data[0].shape[0]):  # Surface
    y = data2[0][i,:]
    ax[1].plot(t, y, label = "subType " + str(i))
    ax[1].set_ylabel("Surface [cm2]")
yt = np.sum(data2[0], axis = 0)
ax[1].plot(t, yt, label = "total")
print("total surface", yt[-1], "cm2")
ax[1].legend()

d_ = data4[0]
krs_ = data3[0]
ax[2].plot(t, d_)
ax[2].set_ylabel("Depth [cm]")
twin_axes = ax[2].twinx()
twin_axes.plot(t, krs_, "r")
twin_axes.set_ylabel("krs [cm2/day]")
print("max depth", np.max(y), "cm")
print("final Krs", np.max(krs_), "cm2/day")

ax[2].set_xlabel("Time [d]")
plt.tight_layout()
plt.show()
