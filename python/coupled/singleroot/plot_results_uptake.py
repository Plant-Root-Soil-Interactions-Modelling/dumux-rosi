"""
single root plots - compares 2 root uptake profiles in one plot (root xylem potential, soil-root interface potential, resulting sink)

from xls result files (in results/) 
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

add_str = "_comp"  # "_wet", "_dry"
# fnames = ["singleroot_cyl_constkrkx" + add_str + ".xls",
#           "singleroot_agg_constkrkx" + add_str + ".xls"]  #
# days = 0.51

fnames = ["singleroot_sra_dynamic_constkrkx" + add_str + ".xls",
          "singleroot_agg_dynamic_constkrkx" + add_str + ".xls"]  #
# fnames = ["singleroot_agg_dynamic_constkrkx" + add_str + ".xls"]
days = 7.1
titles = ["Steady rate", "Aggregated"]  # "Steady rate", , "Aggregated steady rate", "Rhizosphere"

plot_times = range(0, 7)
type_names = ["psix_", "psiinterface_", "sink_"]
L = 50  # cm root length
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

path = "results/"

# psi_interface = ["psiinterface_singleroot_sra_constkrkx_wet.xls", "psiinterface_singleroot_cyl_constkrkx_wet.xls"]
# labels = ["python sra", "python cylindrical"]

fig, ax = plt.subplots(1, 3, figsize = (15, 10))
ax[0].set_title("$\psi_x$ [cm]")
ax[1].set_title("$\psi_{interface}$ [cm]")
ax[2].set_title("sink [cm$^3$]")
ax[0].set_ylabel("depth [cm]")
ax[0].set_xlabel("xylem potential [cm]")
ax[1].set_xlabel("interface potential [cm]")
ax[2].set_xlabel("root uptake [cm$^3$/day]")

ls = ['-', '-.']
for i in range(0, len(fnames)):

    df = pd.read_excel(path + type_names[0] + fnames[i], header = None)
    psi_x_ = df.to_numpy()
    df1 = pd.read_excel(path + type_names[1] + fnames[i], header = None)
    psi_interface_ = df1.to_numpy()
    df2 = pd.read_excel(path + type_names[2] + fnames[i], header = None)
    sink_ = df2.to_numpy()
    z_ = np.linspace(-L + 0.25, -0.25, sink_.shape[1])  # single root 100 segments, 0 - (-50) cm, segment mids

    peak_id = np.round(sink_.shape[0] / days * np.array([0.5 + i for i in plot_times]))
    peak_id = peak_id.astype(int)
    redistribution_id = np.round(sink_.shape[0] / days * np.array([i for i in plot_times]))
    redistribution_id = redistribution_id.astype(int)
    color_intensity = np.ones((sink_.shape[0]),) * 0.2 + np.linspace(1., 0., sink_.shape[0]) * 0.8

    for j in range(0, sink_.shape[0]):
        if j == peak_id[0]:
            ax[0].plot(psi_x_[j,:], z_, color = [color_intensity[j], 0., 0.], linestyle = ls[i])
            ax[1].plot(psi_interface_[j,:], z_, color = [color_intensity[j], 0., 0.], linestyle = ls[i])
            ax[2].plot(sink_[j,:], z_, color = [color_intensity[j], 0., 0.], label = "peak " + titles[i], linestyle = ls[i])
        if j in peak_id[1:]:
            ax[0].plot(psi_x_[j,:], z_, color = [color_intensity[j], 0., 0.], linestyle = ls[i])
            ax[1].plot(psi_interface_[j,:], z_, color = [color_intensity[j], 0., 0.], linestyle = ls[i])
            ax[2].plot(sink_[j,:], z_, color = [color_intensity[j], 0., 0.], linestyle = ls[i])
        if j == redistribution_id[0]:
            ax[0].plot(psi_x_[j,:], z_, 'b', linestyle = ls[i])
            ax[1].plot(psi_interface_[j,:], z_, 'b:', linestyle = ls[i])
            ax[2].plot(sink_[j,:], z_, 'b', label = "initial " + titles[i])
        if j == redistribution_id[1]:
            ax[0].plot(psi_x_[j,:], z_, color = [0., color_intensity[j], 0.], linestyle = ls[i])
            ax[1].plot(psi_interface_[j,:], z_, color = [0., color_intensity[j], 0.], linestyle = ls[i])
            ax[2].plot(sink_[j,:], z_, color = [0., color_intensity[j], 0.], label = "redistribution " + titles[i], linestyle = ls[i])
        if j in redistribution_id[2:]:
            ax[0].plot(psi_x_[j,:], z_, color = [0., color_intensity[j], 0.], linestyle = ls[i])
            ax[1].plot(psi_interface_[j,:], z_, color = [0., color_intensity[j], 0.], linestyle = ls[i])
            ax[2].plot(sink_[j,:], z_, color = [0., color_intensity[j], 0.], linestyle = ls[i])

if add_str == "_dry":
    ax[2].set_xlim(-0.0025, 0.005)
ax[2].legend()
plt.tight_layout()
plt.show()

