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

# # """ DRY """
# sink = ["sink_singleroot_constkrkx_dry.xls", "sink_singleroot_sra_constkrkx_dry.xls"]
# psi_x = ["psix_singleroot_constkrkx_dry.xls", "psix_singleroot_sra_constkrkx_dry.xls"]
# psi_interface = ["psiinterface_singleroot_constkrkx_dry.xls", "psiinterface_singleroot_sra_constkrkx_dry.xls"]
# labels = ["matlab", "python", "python table"]

# # dry 2 - same top, but in hydraulic equilibrium
# sink = ["sink_singleroot_sra_constkrkx_dry2.xls"]
# psi_x = ["psix_singleroot_sra_constkrkx_dry2.xls"]
# psi_interface = ["psiinterface_singleroot_sra_constkrkx_dry2.xls"]
# labels = ["python"]

# # """ DRY cylindric models """
sink = ["sink_singleroot_sra_constkrkx_dry.xls", "sink_singleroot_cyl_constkrkx_dry.xls"]
psi_x = ["psix_singleroot_sra_constkrkx_dry.xls", "psix_singleroot_cyl_constkrkx_dry.xls"]
psi_interface = ["psiinterface_singleroot_sra_constkrkx_dry.xls", "psiinterface_singleroot_cyl_constkrkx_dry.xls"]
labels = ["python sra", "python cylindrical"]

# sink = ["sink_singleroot_sra_constkrkx_dry.xls", "sink_singleroot_cyl_constkrkx_dry.xls", "sink_singleroot_cyl_constkrkx_dryHR.xls"]
# psi_x = ["psix_singleroot_sra_constkrkx_dry.xls", "psix_singleroot_cyl_constkrkx_dry.xls", "psix_singleroot_cyl_constkrkx_dryHR.xls"]
# psi_interface = ["psiinterface_singleroot_sra_constkrkx_dry.xls", "psiinterface_singleroot_cyl_constkrkx_dry.xls", "psiinterface_singleroot_cyl_constkrkx_dryHR.xls"]
# labels = ["python sra", "python cylindrical", "python hr"]

# """ WET """
# sink = ["sink_singleroot_constkrkx_wet.xls", "sink_singleroot_sra_constkrkx_wet.xls"]
# psi_x = ["psix_singleroot_constkrkx_wet.xls", "psix_singleroot_sra_constkrkx_wet.xls"]
# psi_interface = ["psiinterface_singleroot_constkrkx_wet.xls", "psiinterface_singleroot_sra_constkrkx_wet.xls"]
# labels = ["matlab", "python", "python table"]

# """ wet 2 - same top, but in hydraulic equilibrium """
# sink = ["sink_singleroot_sra_constkrkx_wet2.xls"]
# psi_x = ["psix_singleroot_sra_constkrkx_wet2.xls"]
# psi_interface = ["psiinterface_singleroot_sra_constkrkx_wet2.xls"]
# labels = ["python"]

# # """ WET cylindric models """
# sink = ["sink_singleroot_sra_constkrkx_wet.xls", "sink_singleroot_cyl_constkrkx_wet.xls"]
# psi_x = ["psix_singleroot_sra_constkrkx_wet.xls", "psix_singleroot_cyl_constkrkx_wet.xls"]
# psi_interface = ["psiinterface_singleroot_sra_constkrkx_wet.xls", "psiinterface_singleroot_cyl_constkrkx_wet.xls"]
# labels = ["python sra", "python cylindrical"]

fig, ax = plt.subplots(1, 3, figsize = (15, 10))
ax[0].set_title("$\psi_x$")
ax[1].set_title("$\psi_{interface}$")
ax[2].set_title("Sink")
ax[0].set_ylabel("depth [cm]")
ax[0].set_xlabel("matric potential [cm]")
ax[1].set_xlabel("matric potential [cm]")
ax[2].set_xlabel("uptake along root [cm$^3$/day]")

# ax[1].plot(np.linspace(-300, -200, 100), z_, "k:", label = "classic sink") # initial wet
# ax[1].plot(np.linspace(-5000, -200, 100), z_, "k:", label = "classic sink") # initial dry

for i in range(0, len(sink)):
    print(i)
    print(colors[i])
    df = pd.read_excel(path + psi_x[i], header = None)
    psi_x_ = df.to_numpy()
    ax[0].plot(psi_x_[:, 0], z_, '-*', color = colors[i])
    for j in range(1, min(psi_x_.shape[1], 7)):
        ax[0].plot(psi_x_[:, j], z_, ':', color = colors[i])

    df1 = pd.read_excel(path + psi_interface[i], header = None)
    psi_interface_ = df1.to_numpy()
    ax[1].plot(psi_interface_[:, 0], z_, '-*', color = colors[i])
    for j in range(1, min(psi_interface_.shape[1], 7)):
        ax[1].plot(psi_interface_[:, j], z_, ':', color = colors[i])

    df2 = pd.read_excel(path + sink[i], header = None)
    sink_ = df2.to_numpy()
    ax[2].plot(sink_[:, 0], z_, '-*', color = colors[i], label = labels[i])
    for j in range(1, min(sink_.shape[1], 7)):
        ax[2].plot(sink_[:, j], z_, ':', color = colors[i])

ax[2].legend()
plt.margins()
plt.show()

