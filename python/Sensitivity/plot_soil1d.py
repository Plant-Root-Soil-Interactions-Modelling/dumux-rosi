"""
Sink plot (noon and midnight), of a 1d soil 
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

# name = "soybean"
# str_ = "sra0"
# days = 0.25 * 87.5

name = "maize"
str_ = "sra0"
days = 0.25 * 95

# name = "maize"
# str_ = "cyl0"
# days = 0.25 * 95

fname = "soil_" + name + "_" + str_

l = 200  # cm soil depth
dx = 1  # cm resolution
ylim = 99.5

cell_volume = 4  # cm3
plot_times = [1., 5, 10, 15, 20]
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
data = np.load(path + fname + ".npy")
data = np.transpose(data)
data = data[:100:-1,:]

""" sink plot """
fig, ax = plt.subplots(1, 1, figsize = (18, 10))
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size = '5%', pad = 0.05)

cmap_reversed = matplotlib.cm.get_cmap('jet_r')
im = ax.imshow(data, cmap = cmap_reversed, aspect = 'auto', extent = [0, days, -200, 0])  #  interpolation = 'bicubic', interpolation = 'nearest',

cb = fig.colorbar(im, cax = cax, orientation = 'vertical')
cb.ax.get_yaxis().labelpad = 30
cb.set_label('Soil matric potential [cm]', rotation = 270)

ax.set_ylabel("depth [cm]")
ax.set_xlabel("time [days]")

print("range", np.min(data), np.max(data), "cm")

plt.tight_layout()
plt.show()

