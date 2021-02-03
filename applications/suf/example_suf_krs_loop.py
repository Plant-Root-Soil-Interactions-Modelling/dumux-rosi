""" 
soil uptake fraction of a root system (soil is in hydrostatic equilibrium) 
"""
import sys; sys.path.append("../../python/modules/"); sys.path.append("../../../CPlantBox/"); 
from astropy.units import dd
sys.path.append("../../../CPlantBox/src/python_modules")

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt

""" Parameters """
# kr0 = np.array([[0, 2.1e-7], [100, 2.1e-7]])
kr0 = np.array([[0, 2.55e-6], [12.5, 2.55e-6], [20.9, 8.9e-7], [44.6, 8.9e-7], [62.7, 2.1e-7], [100, 2.1e-7]])
kr1 = np.array([[0, 2.55e-6], [12.5, 2.55e-6], [20.9, 8.9e-7], [44.6, 8.9e-7], [62.7, 2.1e-7], [100, 2.1e-7]])
# kr1 = np.array([[0, 2.e-4], [10, 2.e-4], [15, 3.e-5], [20, 3.e-5]])
# kr2 = np.array([[0, 2.55e-6], [12.5, 2.55e-6], [20.9, 8.9e-7], [44.6, 8.9e-7], [62.7, 2.1e-7], [100, 2.1e-7]])
kr2 = np.array([[0, 2.e-4], [10, 2.e-4], [15, 3.e-5], [20, 3.e-5]])
kr3 = np.array([[0, 2.e-4], [10, 2.e-4], [15, 3.e-5], [20, 3.e-5]])
kr4 = np.array([[0, 2.55e-6], [12.5, 2.55e-6], [20.9, 8.9e-7], [44.6, 8.9e-7], [62.7, 2.1e-7], [100, 2.1e-7]])
kr5 = np.array([[0, 2.55e-6], [12.5, 2.55e-6], [20.9, 8.9e-7], [44.6, 8.9e-7], [62.7, 2.1e-7], [100, 2.1e-7]])
kr0 = np.array([[-154, 0.], [0, 1.e-12], [1e20, 1.e-12]])

# kz0 = np.array([[0, 5.752278e-2], [100, 5.752278e-2]])
kz0 = np.array([[0, 2.3148e-4], [18.3, 2.3148e-4], [27.8, 4.0509e-3], [36.4, 4.0509e-3], [51.1, 5.752278e-2], [100, 5.752278e-2]])
kz1 = np.array([[0, 2.3148e-4], [18.3, 2.3148e-4], [27.8, 4.0509e-3], [36.4, 4.0509e-3], [51.1, 5.752278e-2], [100, 5.752278e-2]])
# kz1 = np.array([[0, 1.e-6], [9, 2.e-4], [13, 6.e-4], [20, 6.e-4]])
# kz2 = np.array([[0, 2.3148e-4], [18.3, 2.3148e-4], [27.8, 4.0509e-3], [36.4, 4.0509e-3], [51.1, 5.752278e-2], [100, 5.752278e-2]])
kz2 = np.array([[0, 1.e-6], [9, 2.e-4], [13, 6.e-4], [20, 6.e-4]])
kz3 = np.array([[0, 1.e-6], [9, 2.e-4], [13, 6.e-4], [20, 6.e-4]])
kz4 = np.array([[0, 2.3148e-4], [18.3, 2.3148e-4], [27.8, 4.0509e-3], [36.4, 4.0509e-3], [51.1, 5.752278e-2], [100, 5.752278e-2]])
kz5 = np.array([[0, 2.3148e-4], [18.3, 2.3148e-4], [27.8, 4.0509e-3], [36.4, 4.0509e-3], [51.1, 5.752278e-2], [100, 5.752278e-2]])
kz0 = np.array([[-154, 0.], [0, 1.], [1e20, 1.]])

# kr0 = np.array([[-154, 0.], [0, 1.728e-4], [1e20, 1.728e-4]])
# kr1 = np.array([[-154, 0.], [0, 1.728e-4], [1e20, 1.728e-4]])
# kr2 = np.array([[-154, 0.], [0, 1.728e-4], [1e20, 1.728e-4]])
# kr3 = np.array([[-154, 0.], [0, 1.728e-4], [1e20, 1.728e-4]])
# 
# kz0 = np.array([[-154, 0.], [0, 4.32e-2], [1e20, 4.32e-2]])
# kz1 = np.array([[-154, 0.], [0, 4.32e-2], [1e20, 4.32e-2]])
# kz2 = np.array([[-154, 0.], [0, 4.32e-2], [1e20, 4.32e-2]])
# kz3 = np.array([[-154, 0.], [0, 4.32e-2], [1e20, 4.32e-2]])

simtime = 60  # [day] for task b

""" root system """
rs = pb.MappedRootSystem()
rs.setGeometry(pb.SDF_PlantBox(1.e6, 1.e6, 1.e6))  # not allowed to grow upwards out of soil
p_s = np.linspace(-500, -200, 3001)  #  -200.*np.ones((2001, 1))   # 3 meter down, from -200 to -500, resolution in mm
soil_index = lambda x, y, z: int(-10 * z)  # maps to p_s (hydrostatic equilibirum)
rs.setSoilGrid(soil_index)

path = "../../../CPlantBox//modelparameter/rootsystem/"
name = "Zea_mays_1_Leitner_2010"  # "Glycine_max"  # "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
rs.setSeed(2)
rs.readParameters(path + name + ".xml")
rs.getRootSystemParameter().seedPos.z = -3

# for p in rs.getRootRandomParameter():
#     p.dx = 0.01
rs.initialize()

rs.simulate(simtime, False)
# rs.simulate(9, False)

""" set up xylem parameters """
r = XylemFluxPython(rs)
r.setKrTables([kr0[:, 1], kr1[:, 1], kr2[:, 1], kr3[:, 1], kr4[:, 1], kr5[:, 1]], [kr0[:, 0], kr1[:, 0], kr2[:, 0], kr3[:, 0], kr4[:, 0], kr5[:, 0]])
r.setKxTables([kz0[:, 1], kz1[:, 1], kz2[:, 1], kz3[:, 1], kz4[:, 1], kz5[:, 1]], [kz0[:, 0], kz1[:, 0], kz2[:, 0], kz3[:, 0], kz4[:, 0], kz5[:, 0]])

# alter shoot subType
shoot_segs = rs.getShootSegments()
print("Shoot segments", [str(s) for s in shoot_segs])
print("Shoot type", rs.subTypes[0])

r.test()

Krs_ = []
l_ = []
eswp_ = []
""" numerical solution of transpiration -1 cm3/day"""
for j in range(10, simtime):
    
    # rs.simulate(1., False)
    
    rx = r.solve_neumann(j, -1.e5, p_s, True)  # True: matric potential given per cell (not per segment) high number to recuce spurious fluxes
    print("solved")

    fluxes = r.segFluxes(j, rx, p_s, False, True)  # cm3/day, simTime,  rx,  sx,  approx, cells
    print("Transpiration", r.collar_flux(j, rx, p_s), np.sum(fluxes), "cm3/day")
    suf = np.array(fluxes) / -1.e5  # [1]
    print("Sum of SUF", np.sum(suf), "from", np.min(suf), "to", np.max(suf), "summed positive", np.sum(suf[suf >= 0]))

    eswp = 0.
    n = len(r.rs.segments)
    seg2cell = r.rs.seg2cell
    for i in range(0, n):
        eswp += suf[i] * p_s[seg2cell[i]]

    Krs = 1.e5 / (eswp - rx[1])  # 
    print("Krs: ", Krs, "rx[1]:", rx[1])

    Krs_.append(Krs)
    eswp_.append(eswp) 
    l_.append(rs.getSummed("length"))

time = np.linspace(10, simtime, simtime - 10)
plt.plot(time, Krs_)
plt.xlabel("Number of days")
plt.ylabel("Global root system conductance $[cm^3 hPa^{-1} d^{-1}]$")
# plt.ylim(0., 0.0016)
plt.savefig("Krs_const_soybean.pdf")
plt.show()

