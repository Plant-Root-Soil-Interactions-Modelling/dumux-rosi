import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

"""
   TODO recreates the root system using the xml files and visualizes simulation results on that
"""
import matplotlib.pyplot as plt
import numpy as np

from rosi_richardsnc import RichardsNCSP  # C++ part (Dumux binding), macroscopic soil model
from rosi_richards import RichardsSP
# from rosi_richards import RichardsSPnum as  RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
import plantbox as pb  # CPlantBox

sim_time = 40

initial_age = 1

name = "local_singleroot_conductivities64_27"
name = "soybean_test_0"

alldata = np.load("results/" + name + ".npz")
times = alldata["times"]
# pot_trans = alldata["pot_trans"]
# act_trans = alldata["act_trans"]

# mimic hydraulic_model
rs = pb.MappedPlant()
rs.setSeed(1)
rs.readParameters("results/" + name + ".xml")
rs.setGeometry(pb.SDF_PlantBox(1.e6, 1.e6, -200))
rs.initializeLB(4, 5)
rs.simulate(initial_age, False)

# mimic simulation loop
dt = 360 / (24 * 3600)
N = int(np.ceil(sim_time / dt))  # number of iterations
for i in range(0, N):
    rs.simulate(dt, False)

print("fin")
