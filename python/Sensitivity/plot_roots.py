"""
Dynamic:

    checks if segment analysers are stored

Daniel Leitner, 2025
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import run_SA
import scenario_setup as scenario
from plantbox.functional.PlantHydraulicModel import HydraulicModel_Doussan
from plantbox.functional.PlantHydraulicParameters import PlantHydraulicParameters


def plot_sa(path, name, param="kr", output_index=2):
    """
    param              "rx" (root matric potential), "rsx" (potential at root soil interface), "age", "kr", "kx", "axial_flux", "radial_flux"
    """

    ana = scenario.open_sa_numpy(path + "sa_" + name + f"_{output_index}" + ".npz")

    all = np.load(path + name + ".npz")  # corresponding dynamic file
    times_sa = all["times_sa"]
    print("Root system at time", times_sa[output_index])

    print("Bounds", ana.getMinBounds())

    vp.plot_plant(ana, param)


def plot_suf(path, name):

    all = np.load(path + name + ".npz")  # corresponding dynamic file
    times_sa = all["times_sa"]

    for oi, t in enumerate(times_sa):
        if oi == 0:  # skip first
            continue
        ana = scenario.open_sa_numpy(path + "sa_" + name + f"_{oi}" + ".npz")
        print(oi, t, "bounds", ana.getMinBounds())
        ms = pb.MappedSegments(ana.nodes, ana.segments, ana.data["radius"])
        params = PlantHydraulicParameters()
        params.setKrValues(ana.data["kr"])
        params.setKxValues(ana.data["kx"])
        hm = HydraulicModel_Doussan(ms, params)
        suf = hm.get_suf(t)
        ana.addData("suf", suf)
        n_layers = 100
        z_max = 0
        z_min = -100
        suf_dist = ana.distribution("suf", z_max, z_min, n_layers, True)
        plt.plot(suf_dist, np.linspace(z_max, z_min, n_layers), label="{:g} days".format(t))

    plt.legend()
    plt.show()


def print_sa_indics(path):
    """helps to find the right index in a local sensitivity analysis"""
    names, ranges = run_SA.read_ranges("local_soybean_0/")
    ind = 1
    for i, n in enumerate(names):
        print("Parameter {:s} from {:g} to {:g}".format(n, ind + 1, ind + len(ranges[i])))
        ind = ind + len(ranges[i])


if __name__ == "__main__":

    # path = "local_soybean_1/"
    # name = "local_soybean_1_001"

    plant1 = "ex_name_list1_0_120/soybean_all14_b6270e6280ec99e41746bd015dd7f6e20d4cd33ae33b9d346f5698bdb950/soybean_all14_b6270e6280ec99e41746bd015dd7f6e20d4cd33ae33b9d346f5698bdb95075ba.npz"
    plant2 = "exp_name_list_0_free/soybean_all14_b6270e6280ec99e41746bd015dd7f6e20d4cd33ae33b9d346f5698bdb950/soybean_all14_b6270e6280ec99e41746bd015dd7f6e20d4cd33ae33b9d346f5698bdb95075ba.npz"
    plant3 = "exp_name_list2_0_120/soybean_all14_b96c78843547dbff412faa758815378dde5a923d876fedcbf67d06985e14/soybean_all14_b96c78843547dbff412faa758815378dde5a923d876fedcbf67d06985e143dd7.npz"

    plant1 = "exp_name_list1_0_120/soybean_all14_b6270e6280ec99e41746bd015dd7f6e20d4cd33ae33b9d346f5698bdb950/soybean_all14_b6270e6280ec99e41746bd015dd7f6e20d4cd33ae33b9d346f5698bdb95075ba.npz"
    plant2 = "exp_name_list1_0_free/soybean_all14_0a1816e0417e4805333e8bf7648ee5ae65627703501ce043ed996bd0fe556b61/soybean_all14_0a1816e0417e4805333e8bf7648ee5ae65627703501ce043ed996bd0fe556b61.npz"
    plant3 = "exp_name_list2_0_120/soybean_all14_b96c78843547dbff412faa758815378dde5a923d876fedcbf67d06985e14/soybean_all14_b96c78843547dbff412faa758815378dde5a923d876fedcbf67d06985e143dd7.npz"
    plant4 = "exp_name_list2_0_free/soybean_all14_b292e892b188228d5a342b2ef503e21d2e5d72383bb7a32a7ef70e899b08d8bd/soybean_all14_b292e892b188228d5a342b2ef503e21d2e5d72383bb7a32a7ef70e899b08d8bd.npz"

    p = Path(plant3)
    path = str(p.parent) + "/"
    name = p.name[:-4]  # remove .npz extension

    print_sa_indics(path)  # to find a specific index

    plot_suf(path, name)
    plot_sa(path, name, param="subType")  # radial_flux, axial_flux
