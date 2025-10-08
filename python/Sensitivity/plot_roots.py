"""
    Dynamic:
        
        checks if segment analysers are stored
        
    Daniel Leitner, 2025          
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb
import scenario_setup as scenario
import visualisation.vtk_plot as vp

path = "results/"
name = "soybean_test_0"

all = np.load(path + name + ".npz")

times_sa = all["times_sa"]
print(times_sa)

ana = scenario.open_sa_numpy(path + "sa_"+name +"_1"+ ".npz")

vp.plot_plant(ana, "kr")