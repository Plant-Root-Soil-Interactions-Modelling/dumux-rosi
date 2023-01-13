import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver

import plantbox as pb
import vtk_plot as vp
import scenario_setup as scenario
import sra

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def sigmoid(x, L , x0, k, b):
    y = L / (1 + np.exp(-k * (x - x0))) + b
    return y


def exp2(x, L , x0, k, b):
    y = L * np.exp(-k * (x - x0) * (x - x0)) + b
    return y


typename = "length"

""" maize """
soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.maize(0)
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1)  # , times = x_, net_inf = y_
sra_table_lookup = sra.open_sra_lookup("data/" + table_name)
xml_name = "data/Zeamays_synMRI_modified.xml"  # root growth model parameter file
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
maize_vol_ = [0.]
for i in range(0, 60):
    r.rs.simulate(1.)
    maize_vol_.append(np.sum(r.rs.getParameter(typename)))
maize_vol_ = np.array(maize_vol_)
maize_vol_ = maize_vol_ / maize_vol_[-1]

""" soybean """
soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(0)
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1)  # , times = x_, net_inf = y_
sra_table_lookup = sra.open_sra_lookup("data/" + table_name)
xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
soy_vol_ = [0.]
for i in range(0, 80):
    r.rs.simulate(1.)
    soy_vol_.append(np.sum(r.rs.getParameter(typename)))
soy_vol_ = np.array(soy_vol_)
soy_vol_ = soy_vol_ / soy_vol_[-1]

""" data estimated from literature """
times_soybean = [0., 20., 27, 38, 51, 45, 60, 65, 73, 80, 88, 95, 105]
lai_soybean = [0., 0.2, 0.5, 1., 2.9, 2.2, 4.1, 7.2, 7.2, 9.3, 10., 5.9, 5.9 ]  # from Priscila et al. (2013)

times_maize = [0., 10, 20, 26, 32, 39, 50, 60, 65, 75, 81., 93, 102]
lai_maize = [0.01, 0.2, 0.5, 1.3, 2.2, 3.1, 3.9, 3.9, 3.8, 3.7, 3.6, 2.8, 2.2 ]  # from Boedhram et al. (2001)

print(1. / np.std(times_maize))

p0_soy = [max(lai_soybean), 90, 0.01, min(lai_soybean)]  # this is an mandatory initial guess
p0_maize = [max(lai_maize), np.median(times_maize), 0.02, min(lai_maize)]

popt1, pcov = curve_fit(sigmoid, times_soybean, lai_soybean, p0_soy, method = 'dogbox')
popt1b, pcov = curve_fit(exp2, times_soybean, lai_soybean, p0_soy, method = 'dogbox')

popt2, pcov = curve_fit(sigmoid, times_maize, lai_maize, p0_maize, method = 'dogbox')
popt3, pcov = curve_fit(exp2, times_maize, lai_maize, p0_maize, method = 'dogbox')

# popt2 = p0_maize

x1_ = np.linspace(np.min(times_soybean), np.max(times_soybean), 100)
y1_ = np.array([sigmoid(x, popt1[0], popt1[1], popt1[2], popt1[3]) for x in x1_])
y1b_ = np.array([exp2(x, popt1b[0], popt1b[1], popt1b[2], popt1b[3]) for x in x1_])
x12_ = np.linspace(0, 80, 61)
y12_ = maize_vol_ * np.max(lai_soybean)

x2_ = np.linspace(np.min(times_maize), np.max(times_maize), 100)
y2_ = np.array([sigmoid(x, popt2[0], popt2[1], popt2[2], popt2[3]) for x in x2_])
y3_ = np.array([exp2(x, popt3[0], popt3[1], popt3[2], popt3[3]) for x in x2_])
x4_ = np.linspace(0, 60, 61)
y4_ = maize_vol_ * np.max(lai_maize)

fig, ax = plt.subplots(1, 2, figsize = (18, 6))

ax[0].plot(times_soybean, lai_soybean, "*")
# ax[0].plot(x1_, y1_, label = "sigmoid")
ax[0].plot(x1_, y1b_, label = "normal distribution")
# ax[0].plot(x12_, y12_, 'r:', label = "proportional")
ax[0].set_title("Soybean")
ax[0].set_xlabel("time [day]")
ax[0].legend()
ax[0].set_ylabel("leaf area index [1]")

ax[1].plot(times_maize, lai_maize, "*")
ax[1].set_title("Maize")
# ax[1].plot(x2_, y2_, label = "sigmoid")
ax[1].plot(x2_, np.maximum(y3_, 0), label = "normal distribution")  # cut off negative values
# ax[1].plot(x4_, y4_, 'r:', label = "proportional")
ax[1].set_xlabel("time [day]")
ax[1].legend()
plt.show()

print("\nSoybean")
print("sigmoidal")
print(popt1)
print("neg exp^2")
print(popt1b)

print("\nMaize")
print("sigmoidal")
print(popt2)
print("neg exp^2")
print(popt3)
