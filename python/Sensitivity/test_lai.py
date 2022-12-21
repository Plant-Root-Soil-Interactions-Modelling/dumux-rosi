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
    return (y)


def exp2(x, L , x0, k, b):
    y = L * np.exp(-k * (x - x0) * (x - x0)) + b
    return (y)


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

times_soybean = [0., 20., 30, 50, 60, 70, 80, 90 ]
lai_soybean = [0., 0.2, 0.3, 3., 4.1, 6.5, 8., 7. ]  # Priscila et al. (2013)

times_maize = [0., 20, 40, 60, 80, 100]
lai_maize = [0.01, 0.5, 3.1, 4., 3.5, 2.5 ]  # Boedhram et al. (2001)

p0_soy = [max(lai_soybean), np.median(times_soybean), 1, min(lai_soybean)]  # this is an mandatory initial guess
p0_maize = [max(lai_maize), np.mean(times_maize), 1. / np.std(times_maize), min(lai_maize)]

popt1, pcov = curve_fit(sigmoid, times_soybean, lai_soybean, p0_soy, method = 'dogbox')
popt2, pcov = curve_fit(sigmoid, times_maize, lai_maize, p0_maize, method = 'dogbox')
popt3, pcov = curve_fit(exp2, times_maize, lai_maize, p0_maize, method = 'dogbox')

# popt2 = p0_maize

x1_ = np.linspace(np.min(times_soybean), np.max(times_soybean), 100)
y1_ = np.array([sigmoid(x, popt1[0], popt1[1], popt1[2], popt1[3]) for x in x1_])
x12_ = np.linspace(0, 80, 61)
y12_ = maize_vol_ * np.max(lai_soybean)

x2_ = np.linspace(np.min(times_maize), np.max(times_maize), 100)
y2_ = np.array([sigmoid(x, popt2[0], popt2[1], popt2[2], popt2[3]) for x in x2_])
y3_ = np.array([exp2(x, popt3[0], popt3[1], popt3[2], popt3[3]) for x in x2_])
x4_ = np.linspace(0, 60, 61)
y4_ = maize_vol_ * np.max(lai_maize)

fig, ax = plt.subplots(2, 1, figsize = (8, 16))

ax[0].plot(times_soybean, lai_soybean, "*")
ax[0].plot(x1_, y1_, label = "sigmoid")
ax[0].plot(x12_, y12_, 'r:', label = "proportional")
ax[0].set_title("Soybean")
# ax[0].xlabel("Time")
ax[0].legend()
ax[0].set_ylabel("LAI")

ax[1].plot(times_maize, lai_maize, "*")
ax[1].set_title("Maize")
ax[1].plot(x2_, y2_, label = "sigmoid")
ax[1].plot(x2_, y3_, label = "normal distribution")
ax[1].plot(x4_, y4_, 'r:', label = "proportional")
ax[1].set_xlabel("Time")
ax[1].set_ylabel("LAI")
ax[1].legend()
plt.show()

print("soybean")
print(popt1)

print("maize")
print("sigmoidal")
print(popt2)
print("neg exp^2")
print(popt3)
