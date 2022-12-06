""" 
Check Krs for soy and maize under different hydraulic conductivities 
"""
import sys; sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb
import vtk_plot as vp
import scenario_setup as scenario

# soil_, table_name, p_top, min_b, max_b, cell_number, area, Kc = scenario.maize(0)
# name = "Zeamays_synMRI_modified.xml"

soil_, table_name, p_top, min_b, max_b, cell_number, area, Kc = scenario.soybean(0)
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = p_top, p_bot = (p_top + 200), type = 1)
name = "Glycine_max_Moraes2020_opt2"
xml_name = name + "_modified" + ".xml"  # root growth model parameter file
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth

scenario.init_lupine_conductivities2(r)
r.plot_conductivities(monocot = False, axes_ind = [1, 4], lateral_ind = [2, 3])  # for soy
