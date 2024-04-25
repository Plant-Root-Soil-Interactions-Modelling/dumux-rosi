""" 
    plots SUF of both root systems (Krs on the console)

"""
import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import plantbox as pb
from functional.Perirhizal import *
import visualisation.vtk_plot as vp
import scenario_setup as scenario

import vtk
import numpy as np
import matplotlib.pyplot as plt

SMALL_SIZE = 20
MEDIUM_SIZE = 20
BIGGER_SIZE = 20
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title


def suf_per_layer(plant, soil):
    outer_method = "surface"  # not used
    initial = -200  # not used
    dim = "1D"
    if plant == "springbarley":
        min_b, max_b, cell_number = scenario.springbarley_(dim)
    elif plant == "maize":
        min_b, max_b, cell_number = scenario.maize_(dim)
    else:
        raise "wrong plant name"

    r_, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = scenario.set_scenario(plant, dim, initial, soil, outer_method)
    r = r_.rs

    r_.update(rs_age)

    ana = pb.SegmentAnalyser(r.mappedSegments())
    krs, _ = r_.get_krs(rs_age)
    suf = r_.get_suf()
    ana.addData("SUF", suf)
    n = int(np.ceil(-min_b[2]))
    z_ = np.linspace(-0.5, -n + 0.5, n)
    suf_ = ana.distribution("SUF", 0., float(-n), int(n), False)

    rld_ = ana.distribution("surface", 0., float(-n), int(n), False)

    r_.rs.write(plant + ".rsml")

    layer_vol = np.prod(max_b - min_b) / n

    return np.array(suf_), np.array(rld_) / layer_vol, z_, krs


if __name__ == "__main__":

    for plant in ["springbarley", "maize"]:

        fig, ax = plt.subplots(1, 1, figsize = (8, 8))
        ax = [ax]

        for soil in ["hydrus_sandyloam"]:  # , "hydrus_clay", "hydrus_sandyloam"]: # result is soil independent

            suf, rld, z_, krs = suf_per_layer(plant, soil)

            ax_ = ax[0].twiny()

            ax[0].plot(suf, z_, "-*")
            ax_.plot(rld, z_, ":*", color = "g")

            ax[0].set_ylabel("depth [cm]")
            ax[0].set_xlabel("standard uptake fraction [1]")
            ax_.set_xlabel("root length density [cm cm$^{-3}]$")

            fig.legend(['SUF', 'RLD'], bbox_to_anchor = (0.95, 0.3))

            print()
            print(plant, soil, "SUF from", np.min(suf), "to", np.max(suf), "mean", np.mean(suf), "sum", np.sum(suf), suf.shape)
            print("krs", krs)
            print()

            plt.tight_layout()
            plt.show()
