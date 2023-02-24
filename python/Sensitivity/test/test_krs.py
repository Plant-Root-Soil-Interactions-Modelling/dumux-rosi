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

sim_time = 95

# """ Soybean """
# soil_, table_name, p_top, min_b, max_b, cell_number, area, Kc = scenario.soybean(0)
# s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = p_top, p_bot = (p_top + 200), type = 1)
# xml_name = "Glycine_max_Moraes2020_opt2_modified.xml"
#  # r.plot_conductivities(monocot = False, axes_ind = [1, 4], lateral_ind = [2, 3])  # for soy
# times_ = np.linspace(0, sim_time, 100)
# dt_ = np.diff(times_)
#
# # Case: constant
# age = 1.
# r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
# scenario.init_conductivities_const(r, kr_const = 1.e-4, kx_const = 1.e-3)  # kx_const =0.1
# krs0 = [0.]
# for dt in dt_:
#     print("*", end = '')
#     age += dt
#     r.rs.simulate(dt)
#     krs, _ = r.get_krs(age)
#     krs0.append(krs)
# final_suf0 = r.get_suf(age)
#
# # Case: Zarebanadkouki et al. (2016)
# age = 1.
# r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
# scenario.init_lupine_conductivities(r)
# krs1 = [0.]
# for dt in dt_:
#     print("*", end = '')
#     age += dt
#     r.rs.simulate(dt)
#     krs, _ = r.get_krs(age)
#     krs1.append(krs)
# final_suf1 = r.get_suf(age)
#
# # ana = pb.SegmentAnalyser(r.rs.mappedSegments())
# # ana.addData("suf", r.get_suf(age))
# # ana.addConductivities(r, age)
# # d = ana.data
# # print(d["kx"][0:10])
# # d["kx"][0] = 0.
# # ana.data = d
# # ana.addAge(age)
# # vp.plot_roots(ana, "age")
#
# # Case: Meunier et al. (2018)
# age = 1.
# r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
# scenario.init_lupine_conductivities2(r)
# krs2 = [0.]
# for dt in dt_:
#     print("*", end = '')
#     age += dt
#     r.rs.simulate(dt)
#     krs, _ = r.get_krs(age)
#     krs2.append(krs)
# final_suf2 = r.get_suf(age)
#
# plt.plot(times_, krs0, label = "Constant (kr=1e-4, kx=1e-3)")
# plt.plot(times_, krs1, label = "Zarebanadkouki et al. (2016)")
# plt.plot(times_, krs2, label = "Meunier et al. (2018)")
# plt.xlabel("days")
# plt.ylabel("Krs [cm2/day]")
# print("\n\nkrs", krs0[-1], krs1[-1], krs2[-1])
# plt.title("Soybean")
# print("SUF")
# print(np.min(final_suf0), np.max(final_suf0), np.sum(final_suf0), np.sum(final_suf0 > 0), np.sum(final_suf0 < 0))
# print(np.min(final_suf1), np.max(final_suf1), np.sum(final_suf1), np.sum(final_suf1 > 0), np.sum(final_suf1 < 0))
# print(np.min(final_suf2), np.max(final_suf2), np.sum(final_suf2), np.sum(final_suf2 > 0), np.sum(final_suf2 < 0))
# plt.legend()
# plt.show()
# # krs 0.004279242936161055 0.2644226136609699 0.00614224506711642

""" Maize """
soil_, table_name, p_top, min_b, max_b, cell_number, area, Kc = scenario.maize(0)
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = p_top, p_bot = (p_top + 200), type = 1)
xml_name = "Zeamays_synMRI_modified.xml"
# r.plot_conductivities(monocot = False, axes_ind = [1, 4], lateral_ind = [2, 3])  # for soy
times_ = np.linspace(0, sim_time, 100)
dt_ = np.diff(times_)

# Case: constant:  kr = 1.8e-4, kx = 0.1
age = 1.
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
rs = r.rs
segs = rs.segments
mapping = rs.getSegmentMapper()  # because seg2cell is a dict
for i in range(0, len(segs)):
    if segs[i].x == 0:
        collar_ind = i  # segment index of root collar
        break
scenario.init_conductivities_const(r, kr_const = 1.e-4, kx_const = 1.e-3)  # kx_const =0.1
krs0 = [0.]
for dt in dt_:
    print("*", end = '')
    age += dt
    r.rs.simulate(dt)
    krs, _ = r.get_krs(age, [collar_ind])
    krs0.append(krs)
final_suf0 = r.get_suf(age)

# Case: Doussan et al. (1998)
age = 1.
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
scenario.init_maize_conductivities(r)
krs1 = [0.]
for dt in dt_:
    print("*", end = '')
    age += dt
    r.rs.simulate(dt)
    krs, _ = r.get_krs(age, [collar_ind])
    krs1.append(krs)
final_suf1 = r.get_suf(age)

# Case: Meunier et al. (2018)
age = 1.
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
scenario.init_maize_conductivities2(r)
krs2 = [0.]
for dt in dt_:
    print("*", end = '')
    age += dt
    r.rs.simulate(dt)
    krs, _ = r.get_krs(age, [collar_ind])
    krs2.append(krs)
final_suf2 = r.get_suf(age)

plt.plot(times_, krs0, label = "Constant (kr=1e-4, kx=1e-3)")
plt.plot(times_, krs1, label = "Doussan et al. (1998)")
plt.plot(times_, krs2, label = "Meunier et al. (2018)")
plt.xlabel("days")
plt.ylabel("Krs [cm2/day]")
print("\n\nkrs", krs0[-1], krs1[-1], krs2[-1])
print("SUF")
print(np.min(final_suf0), np.max(final_suf0), np.sum(final_suf0), np.sum(final_suf0 > 0), np.sum(final_suf0 < 0))
print(np.min(final_suf1), np.max(final_suf1), np.sum(final_suf1), np.sum(final_suf1 > 0), np.sum(final_suf1 < 0))
print(np.min(final_suf2), np.max(final_suf2), np.sum(final_suf2), np.sum(final_suf2 > 0), np.sum(final_suf2 < 0))
plt.title("Maize")
plt.legend()
plt.show()
