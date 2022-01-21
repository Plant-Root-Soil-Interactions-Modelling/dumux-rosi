"""
single root plots from xls result files (in results/) 
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

z_ = np.linspace(-49.75, -0.25, 100)  # single root 100 segments, 0 - (-50) cm, segment mids

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

path = "results/"

# """ DRY cylindric models """
# sink = ["soil_singleroot_cyl_dynamic_constkrkx_dry.xls", "soil_singleroot_sra_dynamic_constkrkx_dry.xls"]
sink = ["soil_singleroot_cyl_dynamic_constkrkx_dry.xls", "soil_singleroot_sra_dynamicA_constkrkx_dry.xls"]

# """ WET cylindric models """
# sink = ["soil_singleroot_cyl_dynamic_constkrkx_wet.xls", "soil_singleroot_sra_dynamic_constkrkx_wet.xls"]
# sink = ["soil_singleroot_cyl_dynamic_constkrkx_wet.xls", "soil_singleroot_sra_dynamicA_constkrkx_wet.xls"]

# labels = ["cyl peak", "cyl redistribution", "sra peak", "sra redistribution"]
labels = ["cyl peak", "cyl redistribution", "detached peak", "detached redistribution"]

fig, ax = plt.subplots(1, 2, figsize = (15, 10))
ax[0].set_ylabel("depth [cm]")
ax[0].set_xlabel("water content [1]")
ax[1].set_xlabel("water content [1]")

days = 7

reds = [[1. / (i / 3 + 1), 0., 0.] for i in range(0, days)]
greens = [[0., 1. / (i / 3 + 1), 0.] for i in range(0, days)]
ls = ['-', '-.']

for i in range(0, 2):

    df2 = pd.read_excel(path + sink[i], header = None)
    sink_ = df2.to_numpy()

    if i == 0:
        ax[1].plot(sink_[:, 0], z_, "k:", label = "initial")  # initial
    else:
        ax[1].plot(sink_[:, 0], z_, "k:")  # initial

    for j in range(0, days):
        if j == 0:  # label or not
            ax[0].plot(sink_[:, 2 + j * 4], z_, ls[i], label = labels[2 * i + 0], color = reds[j])
            ax[1].plot(sink_[:, 4 + j * 4], z_, ls[i], label = labels[2 * i + 1], color = greens[j])
        else:
            ax[0].plot(sink_[:, 2 + j * 4], z_, ls[i], color = reds[j])
            ax[1].plot(sink_[:, 4 + j * 4], z_, ls[i], color = greens[j])

ax[0].legend()
ax[1].legend()
plt.margins()
plt.show()

