""" 
calculates the equivalent soil water potential 

1. SUF from a CPlantBox simulated root system
2. soil matric potential from a .vtu (from a prior dumux-rosi simulation)
3. calculate the ESWP 

TODO cut segments along grid
"""
import sys; sys.path.append("../../../python/modules/"); sys.path.append("../../../../CPlantBox/")
sys.path.append("../../../build-cmake/cpp/python_binding/")

from xylem_flux import XylemFluxPython  # Python hybrid solver

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import plantbox as pb
import vtk_plot as vp
import vtk_tools as vt

import numpy as np
import matplotlib.pyplot as plt

import time

""" 1. SUF """

""" Parameters """
kz = 4.32e-2  # axial conductivity [cm^3/day]
kr = 1.728e-4  # radial conductivity [1/day]
simtime = 154  # [day] for task b

""" root system """
rs = pb.MappedRootSystem()
p_s = np.linspace(-200, -500, 3001)  # 3 meter down, from -200 to -500, resolution in mm
soil_index = lambda x, y, z : int(-10 * z)  # maps to p_s (hydrostatic equilibirum)
rs.setSoilGrid(soil_index)

soilcore = pb.SDF_PlantBox(1e6, 1e6, 149.9)
rs.setGeometry(soilcore)

path = "../../../../CPlantBox//modelparameter/rootsystem/"
name = "Glycine_max"  # Zea_mays_1_Leitner_2010
rs.setSeed(1)
rs.readParameters(path + name + ".xml")
rs.initialize()
rs.simulate(simtime, False)

""" set up xylem parameters """
r = XylemFluxPython(rs)
r.setKr([kr])  # or use setKrTables, see XylemFlux.h
r.setKx([kz])

""" numerical solution of transpiration -1 cm3/day"""
rx = r.solve_neumann(simtime, -1, p_s, True)  # True: matric potential given per cell (not per segment)
print("solved")

fluxes = r.segFluxes(simtime, rx, p_s, False, True)  # cm3/day (double simTime,  rx,  sx,  approx, cells
print("Transpiration", r.collar_flux(simtime, rx, p_s), np.sum(fluxes), "cm3/day")
suf = np.array(fluxes) / -1.  # [1]

# """ Additional vtk plot """
ana = pb.SegmentAnalyser(r.rs)
ana.addData("SUF", np.minimum(suf, 1.e-2))  # cut off for vizualisation
# vp.plot_roots(ana, "SUF", "Soil uptake fraction (cm3 day)")  # "fluxes"

""" 2. SOIL MATRIC POTENTIAL """

name = "soybean_Honly-00001"

# Open .vtu
pd = vp.read_vtu(name + ".vtu")
print(pd.GetBounds())  # xmin, xmax, ymin, ymax, zmin, zmax
print("Number of cells", vt.np_cells(pd).shape[0])

data, _ = vt.np_data(pd, 9, True)  # grid, data_index, cell data
print("Data range from {:g} to {:g}".format(np.min(data), np.max(data)))

min_ = np.array([-18.5, -3, -150])
max_ = np.array([18.5, 3, 0.])
res_ = np.array([18, 3, 75])
periodic = True

s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_, max_, res_, periodic)  # [cm]
loam = [0.08, 0.43, 0.04, 1.6, 50]  # we do not plan to calculate a thing, but we need parameters for initialisation
s.setVGParameters([loam])
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.initializeProblem()
s.setInitialCondition(data)  # put data to the grid

""" Coupling (map indices) """
picker = lambda x, y, z : s.pick([x, y, z])
r.rs.setSoilGrid(picker)  # maps segments

# """ 3. EQUIVALENT SOIL WATER POTENTIAL """
t = time.time()
eswp = 0.
n = len(ana.segments)
seg2cell = r.rs.seg2cell
for i in range(0, n):
    d = 0.
    try:
        d = suf[i] * data[seg2cell[i]]
    except:  # out of domain!
        d = 0.
    eswp += d

print(eswp)
elapsed = time.time() - t
print("Time elapsed: ", elapsed)
