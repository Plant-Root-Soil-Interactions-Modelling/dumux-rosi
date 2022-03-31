"""
single root plots - soil potentials of 1-3 soils, with fancy background

from xls result files (in results/) 
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

add_str = "_dry"
fnames = ["results/soil_" + "small_sra" + add_str + ".xls",
          "results/soil_" + "small_agg" + add_str + ".xls"]
days = 7.1
titles = ["steady rate", "aggregated"]  # "steady rate", "aggregated", "rhizosphere",

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

fig, ax = plt.subplots(1, len(fnames), figsize = (14, 14))
if len(fnames) == 1:
    ax = [ax]
for i in range(0, len(fnames)):
    ax[i].set_xlabel("soil potential [cm]")
    ax[i].set_title(titles[i] + " (" + add_str[1:] + ")")
ax[0].set_ylabel("depth [cm]")

for i in range(0, len(fnames)):
    df2 = pd.read_excel(fnames[i], header = None)  # open file
    sink_ = df2.to_numpy()

    z_ = np.linspace(-L + 0.25, -0.25, sink_.shape[1])  # single root 100 segments, 0 - (-50) cm, segment mids

    peak_id = np.round(sink_.shape[0] / days * np.array([0.5 + i for i in range(0, 7)]))
    peak_id = peak_id.astype(int)
    redistribution_id = np.round(sink_.shape[0] / days * np.array([i for i in range(0, 7)]))
    redistribution_id = redistribution_id.astype(int)

    for j in range(0, sink_.shape[0]):
        ax[i].plot(sink_[j,:], z_, 'k', alpha = 0.01)

    for j in range(0, sink_.shape[0]):
        if j == peak_id[0]:
            ax[i].plot(sink_[j,:], z_, 'r:', label = "peak")
        if j in peak_id[1:]:
            ax[i].plot(sink_[j,:], z_, 'r:')
        if j == redistribution_id[0]:
            ax[i].plot(sink_[j,:], z_, 'b:', label = "initial")
        if j == redistribution_id[1]:
            ax[i].plot(sink_[j,:], z_, 'g', label = "redistribution")
        if j in redistribution_id[2:]:
            ax[i].plot(sink_[j,:], z_, 'g')

for i in range(0, len(fnames)):
    ax[i].legend()

plt.tight_layout()
plt.show()

