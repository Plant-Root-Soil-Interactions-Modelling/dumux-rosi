"""
sink term, method comparison, two columns: peak tranpiration, and redistribution.

"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

L = 50  # cm length

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
# sink = ["sink_singleroot_cyl_dynamic_constkrkx_dry.xls", "sink_singleroot_sra_dynamic_constkrkx_dry.xls"]
sink = ["sink_singleroot_sra_dynamic_constkrkx_dry.xls", "sink_singleroot_cyl_dynamic_constkrkx_dry.xls"]

# """ WET cylindric models """
# sink = ["sink_singleroot_cyl_dynamic_constkrkx_wet.xls", "sink_singleroot_sra_dynamic_constkrkx_wet.xls"]
sink = ["sink_singleroot_sra_dynamic_constkrkx_wet.xls", "sink_singleroot_cyl_dynamic_constkrkx_wet.xls"]

# labels = ["cyl peak", "cyl redistribution", "sra peak", "sra redistribution"]
labels = ["cyl peak", "cyl redistribution", "detached peak", "detached redistribution"]

fig, ax = plt.subplots(1, 2, figsize = (15, 10))
ax[0].set_title("Peak transpiration")
ax[1].set_title("Redistribution")
ax[0].set_ylabel("depth [cm]")
ax[0].set_xlabel("root uptake [cm$^3$/day]")
ax[1].set_xlabel("root uptake [cm$^3$/day]")
# ax[1].plot(np.linspace(-300, -200, 100), z_, "k:", label = "initial") # initial wet
# ax[1].plot(np.linspace(-5000, -200, 100), z_, "k:", label = "initial") # initial dry

days = 7

reds = plt.get_cmap("Reds")
greens = plt.get_cmap("Greens")
ls = ['-', '-.']

for i in range(0, 2):

    df2 = pd.read_excel(path + sink[i], header = None)
    sink_ = df2.to_numpy()
    z_ = np.linspace(-0.25, -L + 0.25, sink_.shape[0])  # single root 100 segments, 0 - (-50) cm, segment mids

    ax[1].plot(sink_[:, 0], z_, "k:")  # initial

    for j in range(0, days * 4):
        ax[0].plot(sink_[:, j], z_, linestyle = ls[i], color = reds(0.1 * j / (days * 4)))
        ax[1].plot(sink_[:, j], z_, linestyle = ls[i], color = greens(0.1 * j / (days * 4)))

        if (j + 2) % 4 == 0:
            ax[0].plot(sink_[:, j], z_, linestyle = ls[i], color = "r")
        if j % 4 == 0:
            ax[1].plot(sink_[:, j], z_, linestyle = ls[i], color = "g")

ax[0].legend()
ax[1].legend()
# plt.colorbar()
plt.margins()
plt.show()

