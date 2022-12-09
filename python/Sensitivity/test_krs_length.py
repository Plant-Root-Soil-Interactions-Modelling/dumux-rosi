""" 
Check Krs for soy and maize under different hydraulic conductivities 
"""
import sys; sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb
import vtk_plot as vp
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

""" Soybean """
soil_, table_name, p_top, min_b, max_b, cell_number, area, Kc = scenario.soybean(0)
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = p_top, p_bot = (p_top + 200), type = 1)
xml_name = "Glycine_max_Moraes2020_opt2_modified.xml"
p1 = np.array([1.* 2 ** x for x in np.linspace(-1., 1., 9)])
sim_time = 40

# Case: Meunier et al. (2018)
vol, surf, len_, krs1 = [], [], [], []
for p in p1:
    mods = {"lmax145": p, "r145":p}
    r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name, mods = mods)  # pass parameter file for dynamic growth
    scenario.init_lupine_conductivities(r)
    r.rs.simulate(sim_time)
    ana = pb.SegmentAnalyser(r.rs.mappedSegments())
    vol.append(ana.getSummed("volume"))
    surf.append(ana.getSummed("surface"))
    len_.append(ana.getSummed("length"))
    krs, _ = r.get_krs(1 + sim_time)
    krs1.append(krs)

plt.plot(p1, krs1)
plt.xlabel("length and growth rate altered")
plt.ylabel("Krs [cm2/day]")
plt.title("Soybean")
plt.legend()
print(krs1)
print(vol)
print(surf)
print(len_)
plt.show()

# ana = pb.SegmentAnalyser(r.rs.mappedSegments())
# ana.addData("suf", r.get_suf(1 + sim_time))
# vp.plot_roots(ana, "suf")

# """ Maize """
# soil_, table_name, p_top, min_b, max_b, cell_number, area, Kc = scenario.maize(0)
# s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = p_top, p_bot = (p_top + 200), type = 1)
# xml_name = "Zeamays_synMRI_modified.xml"
# p1 = np.array([1.* 2 ** x for x in np.linspace(-1., 1., 9)])
# sim_time = 40
#
# # Case: Meunier et al. (2018)
# vol, surf, len_, krs1 = [], [], [], []
# for p in p1:
#     mods = {"lmax145": p, "r145":p}
#     r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name, mods = mods)  # pass parameter file for dynamic growth
#     scenario.init_maize_conductivities(r)
#     r.rs.simulate(sim_time)
#     ana = pb.SegmentAnalyser(r.rs.mappedSegments())
#     vol.append(ana.getSummed("volume"))
#     surf.append(ana.getSummed("surface"))
#     len_.append(ana.getSummed("length"))
#     krs, _ = r.get_krs(1 + sim_time)
#     krs1.append(krs)
#
# plt.plot(p1, krs1)
# plt.xlabel("days")
# plt.ylabel("Krs [cm2/day]")
# plt.title("Maize")
# plt.legend()
# print(krs1)
# print(vol)
# print(surf)
# print(len_)
# plt.show()

# ana = pb.SegmentAnalyser(r.rs.mappedSegments())
# ana.addData("suf", r.get_suf(1 + sim_time))
# vp.plot_roots(ana, "suf")
