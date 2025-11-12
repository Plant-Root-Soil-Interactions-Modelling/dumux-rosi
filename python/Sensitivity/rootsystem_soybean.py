""" 
    Pre-study:

    Modify initial soybean system parameters (to fix final initial parametrisation) 
    plots SUF & RLD vs. depth, and other macroscopic parameters after 87.5 days
    
    Daniel Leitner, 2025      
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import plantbox as pb
import visualisation.vtk_plot as vp

import scenario_setup as scenario
import soil_model
import hydraulic_model

import numpy as np
import matplotlib.pyplot as plt

SMALL_SIZE = 24
MEDIUM_SIZE = 24
BIGGER_SIZE = 24
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

""" 1. Modify parameter xml file  ******************************************************** """
# typical domain for soybean
soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(0)
simtime = 87.5  # between 75-100 days

# Open plant and root parameter from a file
rs = pb.Plant()
path = "../../../CPlantBox/modelparameter/structural/rootsystem/"
name = "Glycine_max_Moraes2020_opt2"
rs.readParameters(path + name + ".xml")

srp = rs.getOrganRandomParameter(pb.OrganTypes.seed)
src = srp[0].maxB  # 14
print("maxB", srp[0].maxB)
print("firstB", srp[0].firstB)
print("delayB", srp[0].delayB)

rrp = rs.getOrganRandomParameter(pb.OrganTypes.root)
for p in rrp:
    print("\nSubType", p.subType)
    print("Radius", p.a)
    print("dx", p.dx, "cm")
    print("theta", p.theta, "cm")
    print("lmax", p.lmax, "cm")
    print("changed to 0.5 cm to be faster...")
    p.dx = 0.5  # probably enough
    # print(p.hairsElongation)
    # print(p.hairsZone)
    # print(p.hairsLength)
    # p.hairsZone = 1.7
    # p.hairsLength = 0.1
    # p.hairsElongation = 0.3

rrp[1].theta = 0.8 * rrp[1].theta  # otherwise the initial peak in RLD is a bit too high
rrp[1].thetas = 0.1 * rrp[1].theta  # 10% std
rs.writeParameters("data/" + name + "_modified3" + ".xml")  # remember the modifications ################################### not touching the original _modified

""" 2. Analyse ******************************************************** """
p = np.array([1.* 2 ** x for x in np.linspace(-2., 2., 9)])

s = soil_model.create_richards(soil_, min_b, max_b, cell_number)  # , times = x_, net_inf = y_

xml_name = "data/" + name + "_modified3" + ".xml"  # root growth model parameter file
mods = {"lmax145":1., "lmax2":1., "lmax3":1., "theta45":1.5708, "src":src}
mods = {}

r, params = hydraulic_model.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name, stochastic = False, mods = mods)  # pass parameter file for dynamic growth
# scenario.set_conductivities(params, mods, cdata)

rs = r.ms
# kr = 1.e-4
# kx = 1.e-3
# scenario.init_conductivities_const(r.params, kr, kx)
scenario.init_lupine_conductivities(params, 0.05, 1.)

# Simulate
rs.simulate(simtime, True)
# rs.calcExchangeZoneCoefs()

# Analyse
# vp.plot_roots(rs, "age")
ana = pb.SegmentAnalyser(rs)
ana.addAge(simtime)

orders = np.array(rs.getParameter("radius"))
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

n = len(rs.radii)
radii = np.array([rs.getEffectiveRadius(i) for i in range(0, n)])
print("radii max", n, np.max(radii))

rs.radii = radii
ana = pb.SegmentAnalyser(rs.mappedSegments())
ana.addData("distanceTip", rs.distanceTip)
# vp.plot_roots(ana, "distanceTip")

suf = r.get_suf(simtime)
ana.addData("suf", suf)

dz = 0.5
exact = False
domain_size = np.array(max_b) - np.array(min_b)
slice_volume = domain_size[0] * domain_size[1] * 1  # cm3
z_ = np.linspace(max_b[2] - dz, min_b[2] + dz, cell_number[2])
rld = np.array(ana.distribution("length", 0., min_b[2], cell_number[2], exact)) / slice_volume
rsd = np.array(ana.distribution("surface", 0., min_b[2], cell_number[2], exact)) / slice_volume
suf = np.array(ana.distribution("suf", 0., min_b[2], cell_number[2], exact))
fig, ax = plt.subplots(1, 1, figsize = (8, 10))
ax = [ax]
ax[0].plot(rld, z_, "r", label = "RLD")
ax2 = ax[0].twiny()
ax2.plot(suf, z_, "g", label = "SUF")
ax2.set_xlabel("Standard uptake fraction SUF [1]", color = "green")
ax[0].set_xlabel("Root length density RLD [cm/cm3]", color = "red")
ax[0].set_ylabel("Depth [cm]")
# ax[1].plot(rsd, z_)
# ax[1].set_xlabel("root surface density [cm2 / cm3]")
# ax[1].set_ylabel("depth [cm]")
plt.tight_layout()
plt.show()

