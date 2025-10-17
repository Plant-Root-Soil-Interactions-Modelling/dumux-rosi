"""
    Dynamic: 

    Actual, potential and cumulative transpiration plot 
    manually set area (to calculate [mL/plant] to [mm])
    
    Daniel Leitner, 2025  
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import numpy as np
import matplotlib.pyplot as plt

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


def plot_transpiration(path, name, area = 76 * 3):
    "plots actual and potential transpiration"

    alldata = np.load(path + name + ".npz")
    times = alldata["times"]
    pot_trans = alldata["pot_trans"]
    act_trans = alldata["act_trans"]

    n = 1
    fig, ax = plt.subplots(n, 1, figsize = (18, 8))
    if n == 1:
        ax = [ax]

    for i in range(0, n):
        t = times
        ax[i].plot(t, -10 * pot_trans / area, 'k', label = "potential transpiration")  # potential transpiration, pot_trans [cm3/plant] * plant/cm2 * 10 == [mm]
        ax[i].plot(t, -10 * act_trans / area, 'g', label = "actual transpiration")  # actual transpiration  according to soil model
        ax[i].set_ylabel("transpiration [mm/day]")
        ax[i].legend(loc = 'upper left')
        ax2 = ax[i].twinx()
        dt = np.diff(t)
        cup = -10 * np.cumsum(np.multiply(act_trans[:-1], dt)) / area
        cum_pot = -10 * np.cumsum(np.multiply(pot_trans[:-1], dt)) / area
        ax2.plot(t[1:], cup, 'c--', label = "cumulative transpiration")  # cumulative transpiration (neumann)
        ax2.set_ylabel("cumulative [mm]")
        ax2.legend(loc = 'center right')
        print("Cumulative potential transpiration", cum_pot[-1])
        print("Cumulative plant water uptake", cup[-1], "cm", -10 * np.sum(np.multiply(act_trans[:-1], dt)) / area)

    ax[i].set_xlabel("time [day]")

    plt.tight_layout()
    if not os.path.exists("figures"):
        os.mkdir("figures")
    plt.savefig('figures/' + name + '_trans.png', dpi = 150, bbox_inches = 'tight')
    plt.show()


if __name__ == "__main__":

    # pick12 = [
    # "soybean_all14_e5288952625e69534d1987146abf0263c78dea24a4f8ebbcfd54b3ec02b33752",
    # "soybean_all14_8583e300050dcde7b7745e83d952833faac5258877ab5607822a1f0e4c639b85",
    # "soybean_all14_2d85b66aacdc462dee8c271e13c52d281aa78b54c32225c310526a7f9e0ec247",
    # "soybean_all14_fb6876837510555b1f91d490edd2707e5024bccd8b0cc89b8156e299648513f6",
    # "soybean_all14_6dcaf82cca09f5c136ca0b9da5c6f539197a85f733211ece6665d9afb69ec4e2",
    # "soybean_all14_6f0184db286042bf501cf953d9f82cbf2c6801c1ada76b8a0ae2911a57cfc189",
    # "soybean_all14_b95151025ffa954589fb8d863d5be5a9110ecd4538f23b5cc77750cbee997ee9",
    # "soybean_all14_a656e46e3de66c57292b1abc890a27cc0712984db218111d14d3ec513319ea70",
    # "soybean_all14_6d244e8c11a6a20ad6d4985c337944e9a074b7ce204cbf27c423cf9abed7973b",
    # "soybean_all14_a817d2cd48002b337996bd68fff0e653b5c9e677fae41fca4bab8f3d84f1205b",
    # "soybean_all14_b4d63e89a73c45e0ef4f338b53b5b1ea82434d1bd9b7de5b71e774f9ef9d5fd2",
    # "soybean_all14_3a7d79c45a73e419323d343f07ee3bec6c631bac462622ac21a73c3823c740d0"
    # ]
    # envirotypes = ["0", "1", "5", "36", "59"]

    # cluster_10_opt = [
    #     "soybean_all14_3564a0f636e9f600fe68bf96ffca4124135106ae4787b9a3332334b04abcdf1a",  # 213
    #     "soybean_all14_f7319dc5d83c72932fd39e4afbf6e50822c2f7bf13b27fc5749c1128642a95d2"  # 138
    #     ]

    path = "local_soybean_0/"
    name = "local_soybean_0_001"
    plot_transpiration(path, name)

    path = "local_soybean_free_0/"
    name = "local_soybean_free_0_001"
    plot_transpiration(path, name)

    path = "local_soybean_1/"
    name = "local_soybean_1_001"
    plot_transpiration(path, name)

    path = "local_soybean_free_1/"
    name = "local_soybean_free_1_001"
    plot_transpiration(path, name)

