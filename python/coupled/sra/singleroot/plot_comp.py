"""
single root plots - compares two sink in one axis

from xls result files (in results/) 
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

cmap = plt.get_cmap('Set1')
col = cmap([1, 0, 4, 3, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])  # adjust colors to jans plot

add_str = "_comp"  # "_wet", "_dry"

fnames = ["sink_singleroot_sra_dynamic_constkrkx" + add_str,
          "sink_singleroot_agg_dynamic_constkrkx" + add_str]  # for figure 1 & 2

fnames2 = ["psiinterface_singleroot_sra_dynamic_constkrkx" + add_str,
           "psiinterface_singleroot_agg_dynamic_constkrkx" + add_str]  # for figure 3 & 4

fnames3 = ["soil_singleroot_sra_dynamic_constkrkx" + add_str,
           "soil_singleroot_agg_dynamic_constkrkx" + add_str]  # for figure 3 & 4

fnames4 = ["psix_singleroot_sra_dynamic_constkrkx" + add_str,
           "psix_singleroot_sra_dynamic_constkrkx" + add_str]  # for figure 4

days = 21
titles = ["steady rate", "parallel"]  # "steady rate", , "aggregated" "rhizosphere"

plot_times = [1., 5, 10, 15, 20]
L = 100  # cm root length
dx = 1.  # cm root segment size
path = "results/"
n = len(fnames)

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

""" Sink at noon """
fig, ax = plt.subplots(1, n, figsize = (18, 8))

for i in range(0, 2):
    sink_ = np.load(path + fnames[i] + ".npy")
    z_ = np.linspace(-L + dx / 2., -dx / 2., sink_.shape[1])  # segment mids
    z_ = z_[::-1]

    peak_id = np.round(sink_.shape[0] / days * np.array([0.5 + i for i in plot_times]))
    peak_id = peak_id.astype(int)

    for j in range(0, sink_.shape[0]):
        if j in peak_id:
            ind = np.argwhere(j == peak_id)[0][0]
            ax[i].plot(z_, sink_[j,:], label = "{:g}d".format(plot_times[ind]), color = col[ind])
    ax[i].set_title(titles[i])
    ax[i].set_ylabel("sink term at noon [cm$^3$/day]")
    ax[i].set_xlabel("depth [cm]")
    ax[i].legend()
    ax[i].set_ylim(0., 0.2)

plt.tight_layout()
plt.show()

# """ Sink at midnight"""
# fig, ax = plt.subplots(1, n, figsize = (18, 8))
#
# for i in range(0, 2):
#     sink_ = np.load(path + fnames[i] + ".npy")
#     z_ = np.linspace(-L + dx / 2., -dx / 2., sink_.shape[1])  # segment mids
#     z_ = z_[::-1]
#
#     redistribution_id = np.round(sink_.shape[0] / days * np.array([i for i in plot_times]))
#     redistribution_id = redistribution_id.astype(int)
#
#     for j in range(0, sink_.shape[0]):
#         if j in redistribution_id:
#             ind = np.argwhere(j == redistribution_id)[0][0]
#             ax[i].plot(z_, sink_[j,:], label = "{:g}d".format(plot_times[ind]), color = col[ind])
#     ax[i].set_title(titles[i])
#     ax[i].set_ylabel("sink term at midnight [cm$^3$/day]")
#     ax[i].set_xlabel("depth [cm]")
#     ax[i].legend()
#     ax[i].set_ylim(-0.015, 0.05)
#
# plt.tight_layout()
# plt.show()

""" pressure head at interface and bulk """
fig, ax = plt.subplots(1, n, figsize = (18, 8))

for i in range(0, 2):
    psi_x = np.load(path + fnames2[i] + ".npy")
    psi_s = np.load(path + fnames3[i] + ".npy")

    z_ = np.linspace(-L + dx / 2., -dx / 2., psi_x.shape[1])  # segment mids
    z_ = z_[::-1]
    zs_ = np.linspace(-150 + 1. / 2., -1. / 2., psi_s.shape[1])  # segment mids
    # zs_ = z_[::-1]

    peak_id = np.round(psi_x.shape[0] / days * np.array([0.5 + i for i in plot_times]))
    peak_id = peak_id.astype(int)

    for j in range(0, psi_x.shape[0]):
        if j in peak_id:
            ind = np.argwhere(j == peak_id)[0][0]
            ax[i].plot(z_, psi_x[j,:], label = "{:g}d".format(plot_times[ind]), color = col[ind])
            ax[i].plot(zs_, psi_s[j,:], label = "bulk {:g}d".format(plot_times[ind]), color = col[ind], linestyle = "--")
    ax[i].set_title(titles[i])
    ax[i].set_ylabel("pressure head at interface [cm]")
    ax[i].set_xlabel("depth [cm]")
    ax[i].legend()
    ax[i].set_ylim(-15000., 0.)

plt.tight_layout()
plt.show()

""" hydraulic heads at root collar """
fig, ax = plt.subplots(1, n, figsize = (18, 8))

for i in range(0, 2):
    psi_i = np.load(path + fnames2[i] + ".npy")
    psi_x = np.load(path + fnames4[i] + ".npy")
    psi_s = np.load(path + fnames3[i] + ".npy")
    t_ = np.linspace(0, days, psi_x.shape[0])

    ax[i].plot(t_, psi_x[:, 0], label = "collar", color = col[0])
    ax[i].plot(t_, psi_i[:, 0], label = "soil-root interface", color = col[1])
    ax[i].plot(t_, psi_s[:, -1], label = "bulk soil", color = col[2])

    ax[i].set_title(titles[i])
    ax[i].set_ylabel("hydraulic head [cm]")
    ax[i].set_xlabel("depth [cm]")
    ax[i].legend()
    ax[i].set_ylim(-15000., 0.)

plt.tight_layout()
plt.show()
