"""
    Dynamic:
        
        checks if segment analysers are stored
        
    Daniel Leitner, 2025          
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb
import scenario_setup as scenario
import visualisation.vtk_plot as vp
from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.PlantHydraulicModel import HydraulicModel_Doussan

import run_SA


def plot_sa(path, name, param = "kr", output_index = 2):
    """
    param              "rx" (root matric potential), "rsx" (potential at root soil interfacee), "age", "kr", "kx", "axial_flux", "radial_flux"
    """

    ana = scenario.open_sa_numpy(path + "sa_" + name + f"_{output_index}" + ".npz")

    all = np.load(path + name + ".npz")  # corresponding dynamic file
    times_sa = all["times_sa"]
    print("Root system at time", times_sa[output_index])

    vp.plot_plant(ana, param)


def plot_suf(path, name):

    all = np.load(path + name + ".npz")  # corresponding dynamic file
    times_sa = all["times_sa"]

    for oi, t in enumerate(times_sa):
        print(oi, t)
        ana = scenario.open_sa_numpy(path + "sa_" + name + f"_{oi}" + ".npz")
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
        plt.plot(suf_dist, np.linspace(z_max, z_min, n_layers), label = "{:g} days".format(t))

    plt.legend()
    plt.show()


def print_sa_indics(path):
    """ helps to find the right index in a local sensitivity analysis """
    names, ranges = run_SA.read_ranges("local_soybean_0/")
    ind = 1
    for i, n in enumerate(names):
        print("Parameter {:s} from {:g} to {:g}".format(n, ind + 1, ind + len(ranges[i])))
        ind = ind + len(ranges[i])


if __name__ == "__main__":

    path = "local_soybean_1/"
    name = "local_soybean_1_001"

    print_sa_indics(path)  # to find a specific index

    plot_suf(path, name)
    plot_sa(path, name, param = "subType")  # radial_flux, axial_flux

