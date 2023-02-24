
""" 
Analyse and choose (modify) maize root system parameters 
"""
import sys; sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb
import vtk_plot as vp
from xylem_flux import *
import scenario_setup as scenario

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

# typical domain for maize
soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.maize(0)
simtime = 30  # 1.* 95  # between 90-100 days

rs = pb.MappedRootSystem()  # RootSystem
rs.setRectangularGrid(pb.Vector3d(100 * min_b[0], 100 * min_b[1], min_b[2]), pb.Vector3d(100 * max_b[0], 100 * max_b[1], max_b[2]),
                      pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), cut = False)

# Open plant and root parameter from a file
path = "../../../CPlantBox/modelparameter/rootsystem/"
name = "Zeamays_synMRI"
# Zea_mays_1_Leitner_2010
# Zea_mays_2_Pagès_2014
# Zea_mays_3_Postma_2011,
# Zea_mays_4_Leitner_2014
# Zea_mays_5_Leitner_2014
# Zea_mays_6_Leitner_2014
# Zea_mays_Streuber_2020 (groß, viele shoot borne überirdisch)
# Zeamays_synMRI
rs.readParameters(path + name + ".xml")

srp = rs.getOrganRandomParameter(pb.OrganTypes.seed)
print(srp[0])
# print("original seminals")  # between two and six seminal root
# print("first", srp[0].firstB)
# print("delay", srp[0].delayB)
# print("max", srp[0].maxB)
# print("original shootborne")  # The maize rootstock develops ∼70 shoot-borne roots in the course of development (Hoppe et al., 1986)
# print("first", srp[0].firstSB)
# print("delay", srp[0].delaySB)
# print("delayRC", srp[0].delayRC)
# print("number of crowns", srp[0].nC)
# print("shift in crowns", srp[0].nz)
# print()

srp[0].seedPos.z = -3.

""" Seminal roots """
srp[0].firstB = 0.5
srp[0].delayB = 0.1
srp[0].maxB = 6  # between two and six seminal root, Andrea 0-10

""" Shoot borne roots (brace roots, nodal roots) """
srp[0].firstSB = 1
srp[0].delaySB = 1.5
srp[0].delayRC = 15
srp[0].nC = 11  # number of roots per root crown
srp[0].nz = 0.1

""" root parameters """
# add shoot borne parameter set (copy subType 4)
rrp = rs.getOrganRandomParameter(pb.OrganTypes.root)
rp5 = rrp[4].copy(rs)
rp5.subType = 5
rp5.name = "shootborne"
rs.setOrganRandomParameter(rp5)
print(rp5)

print(rrp[2].lmax, rrp[2].lmaxs)
rrp[2].lmaxs = rrp[2].lmax
rrp[2].lb = 1.4
# rrp[4].lmax = 75  # wild guess
rrp[4].theta = 85. / 180 * np.pi  # wild guess
# rrp[4].ln = 0.2
# rrp[4].r = 0.5
#
# rrp[2].ldelay = 6
#
# rrp[2].r = 0.5
# rrp[3].r = 0.5

for i, p in enumerate(rrp[1:]):
    p.la = p.ln  # Never use la for delay based

scenario.set_all_sd(rs, 0)
rrp = rs.getOrganRandomParameter(pb.OrganTypes.root)

print(rrp[1].tropismN, rrp[1].tropismS)
print(rrp[4].tropismN, rrp[1].tropismS)
print(rrp[5].tropismN, rrp[1].tropismS)
rrp[1].tropismN = 0.1
rrp[4].tropismN = 0.1
rrp[5].tropismN = 0.1

for i, p in enumerate(rrp[1:]):
    if i == 0:
        p.theta = 0
        p.thetas = 0
    print("\nSubType", p.subType, p.name)
    print("Radius", p.a)
    print("lmax", p.lmax)
    print("ln", p.ln)
    print("lb", p.lb)
    print("la", p.la)
    # print("radius", p.a)
    print("initial growth rate", p.r)
    # print("theta", p.theta / np.pi * 180)
    print("delay", p.ldelay)
    print("lnk", p.lnk)
    p.dxMin = 0.005
    # print(p.dxMin)
    # print("successor", p.successor)
    p.dx = 0.5  # probably enough
    # p.dxMin = 0.05

""" initialize """
rs.setGeometry(pb.SDF_PlantBox(1.e6, 1.e6, np.abs(min_b[2])))
rs.initializeDB()
rs.writeParameters("data/" + name + "_modified" + ".xml")  # remember the modifications
print("'Resulting number of shoot segments", len(rs.getShootSegments()))  # due to crown root nodes
r = XylemFluxPython(rs)
scenario.init_maize_conductivities(r)

""" simulation """
rs.simulate(simtime, True)

# """ animation """
# param_name = "subType"  # kx, kr, cell_id
# anim = vp.AnimateRoots(rs)
# anim.min = np.array([-20., -20, -50])
# anim.max = np.array([20, 20, 0.])
# anim.res = np.array([1, 1, 1])
# anim.avi_name = "results/maize_anim"
# anim.start()
# dt = 0.1  # days
# N = round(simtime / dt)
# for i in range(0, N):
#     rs.simulate(dt, False)
#     ana = pb.SegmentAnalyser(rs.mappedSegments())
#     ana.addConductivities(r, i * dt + 1)
#     ana.addCellIds(rs)
#     anim.rootsystem = ana
#     anim.root_name = "cell_id"
#     anim.update()

""" Analyse """
ana = pb.SegmentAnalyser(rs.mappedSegments())
ana.addConductivities(r, simtime)
ana.addAge(simtime)
ana.addCellIds(rs.mappedSegments())

vp.plot_roots(ana, "subType", name)
# vp.plot_roots(ana, "cell_id", name)
# vp.plot_roots(ana, "kr", name + "_kr")
# vp.plot_roots(ana, "kx", name + "_kx")
# s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1)
# data = ana.data["cell_id"]
# print(np.min(data), np.max(data))
# # vp.plot_roots_and_soil(rs.mappedSegments(), "cell_id", data, s, True, min_b, max_b, cell_number, "dummy")
# vp.plot_roots(ana, "cell_id")

orders = np.array(rs.getParameter("subType"))
print("\nnumber of roots", len(rs.getRoots()))
print("types", np.sum(orders == 1), np.sum(orders == 2), np.sum(orders == 3), np.sum(orders == 4), np.sum(orders == 5))
print("number of nodes", len(ana.nodes))
print("number of segments", len(ana.segments))
print("volume", np.sum(ana.getParameter("volume")), "cm3")
print("surface", np.sum(ana.getParameter("surface")), "cm2")
print("\nunconfined", ana.getMinBounds(), "-", ana.getMaxBounds())

# w = np.array(max_b) - np.array(min_b)
# ana.mapPeriodic(w[0], w[1])
# print("periodic")
# print(ana.getMinBounds(), "-", ana.getMaxBounds())

""" RLD """
dz = 0.5
exact = False
domain_size = np.array(max_b) - np.array(min_b)
slice_volume = domain_size[0] * domain_size[1] * 1  # cm3
z_ = np.linspace(max_b[2] - dz, min_b[2] + dz, cell_number[2])
rld = np.array(ana.distribution("length", 0., min_b[2], cell_number[2], exact)) / slice_volume
rsd = np.array(ana.distribution("surface", 0., min_b[2], cell_number[2], exact)) / slice_volume
fig, ax = plt.subplots(1, 1, figsize = (10, 10))
ax = [ax]
ax[0].plot(rld, z_)
ax[0].set_xlabel("root length density [cm / cm3]")
ax[0].set_ylabel("depth [cm]")
# ax[1].plot(rsd, z_)
# ax[1].set_xlabel("root surface density [cm2 / cm3]")
# ax[1].set_ylabel("depth [cm]")
plt.tight_layout()
plt.savefig("results/" + name + "_rld" + ".png")
plt.show()  # e.g. compare to Zhuang et al. 2001

""" export """
ana.write("results/maize.vtp")  # rs.write writes polylines, ana.write writes segments
rs.write("results/maize.rsml")
w = np.array(max_b) - np.array(min_b)
print(w)
rs.setGeometry(pb.SDF_PlantBox(w[0], w[1], w[2]))
rs.write("results/maize_box.py")
ana.mapPeriodic(w[0], w[1])
ana.write("results/maize_periodic.vtp")

#
# print("fin")
