"""
single root plots - compares two sink in one axis

from xls result files (in results/) 
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

add_str = "_dry"

fnames = ["results/sink_" + "small_soil" + add_str,
          "results/sink_" + "small_sra" + add_str,
          "results/sink_" + "small_agg" + add_str]
titles = ["Cylindric", "Steady Rate Approach", "Parallel"]
days = 7.1

plot_times = range(0, 7)
L = 110  # soil depth

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

fig, ax = plt.subplots(1, 1, figsize = (14, 14))
ax = [ax]
ax[0].set_xlabel("sink [cm$^3$/day]")
ax[0].set_ylabel("depth [cm]")
ax[0].plot([0, 0], [-L, 0.], "k:")

ls = ['-', '-.']
for i in range(0, 2):
    df2 = pd.read_excel(fnames[i] + ".xls", header = None)  # open file
    sink_ = -df2.to_numpy()

    z_ = np.linspace(-L + 0.25, -0.25, sink_.shape[1])  # single root 100 segments, 0 - (-50) cm, segment mids

    peak_id = np.round(sink_.shape[0] / days * np.array([0.5 + i for i in plot_times]))
    peak_id = peak_id.astype(int)
    redistribution_id = np.round(sink_.shape[0] / days * np.array([i for i in plot_times]))
    redistribution_id = redistribution_id.astype(int)
    color_intensity = np.ones((sink_.shape[0]),) * 0.2 + np.linspace(1., 0., sink_.shape[0]) * 0.8

    for j in range(0, sink_.shape[0]):
        if j == peak_id[0]:
            ax[0].plot(sink_[j,:], z_, color = [color_intensity[j], 0., 0.], label = "peak " + titles[i], linestyle = ls[i])
        if j in peak_id[1:]:
            ax[0].plot(sink_[j,:], z_, color = [color_intensity[j], 0., 0.], linestyle = ls[i])
        if j == redistribution_id[0]:
            ax[0].plot(sink_[j,:], z_, 'b:', label = "initial " + titles[i])
        if j == redistribution_id[1]:
            ax[0].plot(sink_[j,:], z_, color = [0., color_intensity[j], 0.], label = "redistribution " + titles[i], linestyle = ls[i])
        if j in redistribution_id[2:]:
            ax[0].plot(sink_[j,:], z_, color = [0., color_intensity[j], 0.], linestyle = ls[i])

ax[0].legend()
plt.tight_layout()
plt.show()

