""" 
Show how tap root length and branching density shapes the root system 
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import plantbox as pb

import visualisation.vtk_plot as vp

import scenario_setup as scenario

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
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(0)

fig, ax = plt.subplots(3, 3, figsize = (18, 14))


def plot_(i, j, simtime, mods):  # i row, j column

    # Open plant and root parameter from a file
    path = "../../../CPlantBox/modelparameter/structural/rootsystem/"
    name = "Glycine_max_Moraes2020_opt2"
    s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1)  # , times = x_, net_inf = y_
    xml_name = "data/" + name + "_modified" + ".xml"  # root growth model parameter file
    r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name, stochastic = False, mods = mods)  # pass parameter file for dynamic growth
    kr = 1.e-4
    kx = 1.e-3
    scenario.init_conductivities_const(r.params, kr, kx)
    rs = r.ms

    rs.setGeometry(pb.SDF_PlantBox(1.e6, 1.e6, 200.))
    rs.initialize()

    rs.simulate(simtime, True)

    ana = pb.SegmentAnalyser(rs)
    ana.addAge(simtime)

    orders = np.array(rs.getParameter("order"))
    print("\nnumber of roots", len(rs.getOrgans()))
    print("types", np.sum(orders == 1), np.sum(orders == 2), np.sum(orders == 3), np.sum(orders == 4), np.sum(orders == 5))
    print("number of nodes", len(ana.nodes))
    print("number of segments", len(ana.segments))
    print("volume", np.sum(ana.getParameter("volume")), "cm3")
    print("surface", np.sum(ana.getParameter("surface")), "cm2")
    print("Krs", r.get_krs(simtime)[0], "cm2/day")

    print("\nunconfined")
    print(ana.getMinBounds(), "-", ana.getMaxBounds())

    w = np.array(max_b) - np.array(min_b)
    ana.mapPeriodic(w[0], w[1])
    print("periodic")
    print(ana.getMinBounds(), "-", ana.getMaxBounds())

    # vp.plot_roots(ana, "age")

    dz = 0.5
    exact = False
    domain_size = np.array(max_b) - np.array(min_b)
    slice_volume = domain_size[0] * domain_size[1] * 1  # cm3
    z_ = np.linspace(max_b[2] - dz, min_b[2] + dz, cell_number[2])
    rld = np.array(ana.distribution("length", 0., min_b[2], cell_number[2], exact)) / slice_volume
    rsd = np.array(ana.distribution("surface", 0., min_b[2], cell_number[2], exact)) / slice_volume

    xlim = [0.3, 0.9, 1.3]
    ax[i, j].plot(rld, z_)
    ax[i, j].set_xlabel("root length density [cm / cm3]")
    ax[i, j].set_ylabel("depth [cm]")
    ax[i, j].set_xlim(0., xlim[j])
    # ax[1].plot(rsd, z_)
    # ax[1].set_xlabel("root surface density [cm2 / cm3]")
    # ax[1].set_ylabel("depth [cm]")

    ana.write("results/soybean{:g}_{:g}.vtp".format(simtime, i))
    w = np.array(max_b) - np.array(min_b)
    rs.setGeometry(pb.SDF_PlantBox(w[0], w[1], w[2]))
    rs.write("results/soybean_box.py")

# mods = {"lmax145":1., "lmax2":1., "lmax3":1., "theta45":1.5708, "r145":1., "r2":1., "a":1.}
    return r.get_krs(simtime)[0]


krs_ = np.zeros((3, 3))
simtimes = [20, 50, 87.5]

for j in range(0, 3):
    mods = {}
    krs_[0, j] = plot_(0, j, simtimes[j], mods)

for j in range(0, 3):
    mods = {"lmax": 0.666, "r": 0.666}  # "lmax145": 0.5,
    krs_[1, j] = plot_(1, j, simtimes[j], mods)

for j in range(0, 3):
    mods = {"ln": 1.5}
    krs_[2, j] = plot_(2, j, simtimes[j], mods)

print(krs_)

plt.tight_layout()
plt.show()

