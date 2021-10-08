"""
single root plots from xls result files (in results/) 
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

z_ = np.linspace(-0.25, -49.75, 100)  # single root 100 segments, 0 - (-50) cm, segment mids

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
sink = ["sink_singleroot_sra_dynamic_constkrkx_dry.xls", "sink_singleroot_cyl_dynamic_constkrkx_dry.xls"]
psi_x = ["psix_singleroot_sra_dynamic_constkrkx_dry.xls", "psix_singleroot_cyl_dynamic_constkrkx_dry.xls"]
psi_interface = ["psiinterface_singleroot_sra_dynamic_constkrkx_dry.xls", "psiinterface_singleroot_cyl_dynamic_constkrkx_dry.xls"]
labels = ["python sra", "python cylindrical"]

# """ WET cylindric models """
# sink = ["sink_singleroot_sra_dynamic_constkrkx_wet.xls", "sink_singleroot_cyl_dynamic_constkrkx_wet.xls"]
# psi_x = ["psix_singleroot_sra_dynamic_constkrkx_wet.xls", "psix_singleroot_cyl_dynamic_constkrkx_wet.xls"]
# psi_interface = ["psiinterface_singleroot_sra_dynamic_constkrkx_wet.xls", "psiinterface_singleroot_cyl_dynamic_constkrkx_wet.xls"]
# labels = ["python sra", "python cylindrical"]
#
# sink = ["sink_singleroot_sra_dynamic_constkrkx_wet2.xls", "sink_singleroot_cyl_dynamic_constkrkx_wet2.xls"]
# psi_x = ["psix_singleroot_sra_dynamic_constkrkx_wet2.xls", "psix_singleroot_cyl_dynamic_constkrkx_wet2.xls"]
# psi_interface = ["psiinterface_singleroot_sra_dynamic_constkrkx_wet2.xls", "psiinterface_singleroot_cyl_dynamic_constkrkx_wet2.xls"]
# labels = ["python sra", "python cylindrical"]

fig, ax = plt.subplots(1, 3, figsize = (15, 10))
ax[0].set_title("$\psi_x$")
ax[1].set_title("$\psi_{interface}$")
ax[2].set_title("Sink")
ax[0].set_ylabel("depth [cm]")
ax[0].set_xlabel("matric potential [cm]")
ax[1].set_xlabel("matric potential [cm]")
ax[2].set_xlabel("uptake per root segment [cm$^3$/day]")

# ax[1].plot(np.linspace(-300, -200, 100), z_, "k:", label = "classic sink") # initial wet
# ax[1].plot(np.linspace(-5000, -200, 100), z_, "k:", label = "classic sink") # initial dry

days = 7

labels = ["sra peak", "sra redistribution", "cyl peak", "cyl redistribution"]
reds = [[1. / (i / 3 + 1), 0., 0.] for i in range(0, days)]
greens = [[0., 1. / (i / 3 + 1), 0.] for i in range(0, days)]
ls = ['-', '-.']

for i in range(0, 2):

    print(i)
    print(colors[i])

    df0 = pd.read_excel(path + psi_x[i], header = None)
    psi_x_ = df0.to_numpy()
    for j in range(0, days):
        if j == 0:
            ax[0].plot(psi_x_[:, 2 + j * 4], z_, ls[i], label = labels[2 * i + 0], color = reds[j])
            ax[0].plot(psi_x_[:, 4 + j * 4], z_, ls[i], label = labels[2 * i + 1], color = greens[j])
        else:
            ax[0].plot(psi_x_[:, 2 + j * 4], z_, ls[i], color = reds[j])
            ax[0].plot(psi_x_[:, j * 4], z_, ls[i], color = greens[j])

    df1 = pd.read_excel(path + psi_interface[i], header = None)
    psi_interface_ = df1.to_numpy()
    for j in range(0, days):
        if j == 0:
            ax[1].plot(psi_interface_[:, 2 + j * 4], z_, ls[i], label = labels[2 * i + 0], color = reds[j])
            ax[1].plot(psi_interface_[:, 4 + j * 4], z_, ls[i], label = labels[2 * i + 1], color = greens[j])
        else:
            ax[1].plot(psi_interface_[:, 2 + j * 4], z_, ls[i], color = reds[j])
            ax[1].plot(psi_interface_[:, j * 4], z_, ls[i], color = greens[j])

    df2 = pd.read_excel(path + sink[i], header = None)
    sink_ = df2.to_numpy()
    for j in range(0, days):
        if j == 0:
            ax[2].plot(sink_[:, 2 + j * 4], z_, ls[i], label = labels[2 * i + 0], color = reds[j])
            ax[2].plot(sink_[:, 4 + j * 4], z_, ls[i], label = labels[2 * i + 1], color = greens[j])
        else:
            ax[2].plot(sink_[:, 2 + j * 4], z_, ls[i], color = reds[j])
            ax[2].plot(sink_[:, j * 4], z_, ls[i], color = greens[j])

ax[2].legend()
plt.margins()
plt.show()

