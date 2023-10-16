"""
Sink plot (noon and midnight), of a 1d soil 
"""
import matplotlib.pyplot as plt
import numpy as np

name = "singleroot"
str_ = ["sra", "ups"]
# str_ = ["cyl", "sra", "agg"]
# str_ = ["cyl", "sra", "agg", "ups"]

# name = "rootsystem1d"
# str_ = ["sra", "ups"]

fnames = np.array(["sink_" + name + "_" + s for s in str_ ])

cmap = plt.get_cmap('Set1')
col = cmap([1, 0, 4, 3, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])  # adjust colors to jans plot

days = 21
l = 150  # cm soil depth
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
n = len(fnames)
data = [np.load(path + n_ + ".npy")  for n_ in fnames]

""" sink plot """
fig, ax = plt.subplots(1, 2, figsize = (18, 10))
ax[0].set_ylabel("depth [cm]")
ax[0].set_xlabel("sink term at noon [1/day]")
ax[1].set_xlabel("sink term at night [1/day]")
ax[0].plot([0, 0], [-l, 0.], "k:")
ax[1].plot([0, 0], [-l, 0.], "k:")
ls = ["-", "--", "-.", ":"]

""" noon """
for i in range(0, n):

    sink_ = data[i]
    soil_z_ = np.linspace(-l + dx / 2., -dx / 2., sink_.shape[1])  # segment mids

    peak_id = np.round(sink_.shape[0] / days * np.array([0.5 + i for i in plot_times]))
    peak_id = peak_id.astype(int)

    for ind, j in enumerate(peak_id):
        lstr = "{:g}d ({:s})".format(plot_times[ind], str_[i])
        ax[0].plot(sink_[j,:] / cell_volume, soil_z_, label = lstr, color = col[ind], linestyle = ls[i])

    ax[0].set_ylim([-ylim, 0.])
    ax[0].legend()

""" midnight """
for i in range(0, n):

    sink_ = data[i]
    soil_z_ = np.linspace(-l + dx / 2., -dx / 2., sink_.shape[1])  # segment mids

    redistribution_id = np.round(sink_.shape[0] / days * np.array([i for i in plot_times]))
    redistribution_id = redistribution_id.astype(int)

    for ind, j in enumerate(redistribution_id):
        lstr = "{:g}d ({:s})".format(plot_times[ind], str_[i])
        ax[1].plot(sink_[j,:] / cell_volume, soil_z_, label = lstr, color = col[ind], linestyle = ls[i])

    ax[1].set_ylim([-ylim, 0.])
    ax[1].legend()

plt.tight_layout()
plt.show()

