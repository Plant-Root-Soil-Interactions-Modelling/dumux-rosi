""" 
Check conductivities for soy and maize 
"""
import sys; sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb
import vtk_plot as vp
import scenario_setup as scenario

soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.maize(0)
fname = "data/Zeamays_synMRI_modified.xml"

# soil_, table_name, p_top, min_b, max_b, cell_number, area, Kc = scenario.soybean(0)
# fname = "data/Glycine_max_Moraes2020_opt2_modified.xml"

s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1)
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, fname)  # pass parameter file for dynamic growth

# # kr = 1.e-4
# # kx = 1.e-3

# scenario.init_lupine_conductivities(r)
#  scenario.init_dynamic_simple_growth(r, kr0, kr1, kx0, kx1, dt0 = 30., dt1 = 17., kr_f = 0.25, kx_f = 5.)

# # scenario.init_lupine_conductivities2(r)
# r.plot_conductivities(monocot = False, axes_ind = [1, 4], lateral_ind = [2, 3])  # for soy

scenario.init_maize_conductivities(r)
r.plot_conductivities(monocot = True)  # for maize

