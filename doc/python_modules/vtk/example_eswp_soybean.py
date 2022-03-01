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
import glob
import time

""" 1. SUF """

""" Parameters """
kr0 = np.array([[-154, 0.0006736268], [0, 0.0006736268], [2, 0.00109], [4, 0.00103], [6, 0.000983], [8, 0.000935], [10, 0.00089], [12, 0.000847], [14, 0.000806], [16, 0.000767], [18, 0.00073], [20, 0.000695], [22, 0.000662], [24, 0.00063], [26, 0.000599], [28, 0.00057], [30, 0.000543], [32, 0.000517], [60, 0.0008], [1e20, 0.0008]])
kr1 = np.array([[-154, 0.0006736268], [0, 0.0006736268], [2, 0.00109], [4, 0.00103], [6, 0.000983], [8, 0.000935], [10, 0.00089], [12, 0.000847], [14, 0.000806], [16, 0.000767], [18, 0.00073], [20, 0.000695], [22, 0.000662], [24, 0.00063], [26, 0.000599], [28, 0.00057], [30, 0.000543], [32, 0.000517], [60, 0.0008], [1e20, 0.0008]])
kr2 = np.array([[-154, 0.0006736268], [0, 0.0006736268], [2, 0.00389], [4, 0.00367], [6, 0.00347], [8, 0.00328], [10, 0.0031], [12, 0.00293], [14, 0.00277], [16, 0.00262], [18, 0.00248], [20, 0.00234], [22, 0.00221], [24, 0.00209], [26, 0.00198], [28, 0.00187], [30, 0.00177], [32, 0.00167], [60, 0.0018], [1e20, 0.0018]])
kr3 = np.array([[-154, 0.0006736268], [0, 0.0006736268], [2, 0.00389], [4, 0.00367], [6, 0.00347], [8, 0.00328], [10, 0.0031], [12, 0.00293], [14, 0.00277], [16, 0.00262], [18, 0.00248], [20, 0.00234], [22, 0.00221], [24, 0.00209], [26, 0.00198], [28, 0.00187], [30, 0.00177], [32, 0.00167], [60, 0.0018], [1e20, 0.0018]])

kz0 = np.array([[-154, 0.0006736268], [0, 0.0006736268], [2, 0.074759246], [4, 0.08296797], [6, 0.09207803], [8, 0.102188394], [10, 0.113408897], [12, 0.125861436], [14, 0.13968129], [16, 0.155018593], [18, 0.172039965], [20, 0.190930319], [22, 0.211894875], [24, 0.235161384], [26, 0.260982605], [28, 0.289639051], [30, 0.321442035], [32, 0.356737056], [60, 4.3], [1e20, 4.3]])
kz1 = np.array([[-154, 0.0006736268], [0, 0.0006736268], [2, 0.074759246], [4, 0.08296797], [6, 0.09207803], [8, 0.102188394], [10, 0.113408897], [12, 0.125861436], [14, 0.13968129], [16, 0.155018593], [18, 0.172039965], [20, 0.190930319], [22, 0.211894875], [24, 0.235161384], [26, 0.260982605], [28, 0.289639051], [30, 0.321442035], [32, 0.356737056], [60, 4.3], [1e20, 4.3]])
kz2 = np.array([[-154, 0.000407], [0, 0.000407], [1, 0.0005], [2, 0.000615], [4, 0.000756], [6, 0.00093], [8, 0.00114], [10, 0.00141], [12, 0.00173], [14, 0.00212], [16, 0.00261], [18, 0.00321], [20, 0.00395], [22, 0.00486], [24, 0.00597], [26, 0.00734], [28, 0.00903], [30, 0.0111], [32, 0.0136], [60, 0.43], [1e20, 0.43]])
kz3 = np.array([[-154, 0.000407], [0, 0.000407], [1, 0.0005], [2, 0.000615], [4, 0.000756], [6, 0.00093], [8, 0.00114], [10, 0.00141], [12, 0.00173], [14, 0.00212], [16, 0.00261], [18, 0.00321], [20, 0.00395], [22, 0.00486], [24, 0.00597], [26, 0.00734], [28, 0.00903], [30, 0.0111], [32, 0.0136], [60, 0.43], [1e20, 0.43]])
simtime = 154  # [day] for task b

""" root system """
rs = pb.MappedRootSystem()
p_s = np.linspace(-500, -200, 3001)  # 3 meter down, from -200 to -500, resolution in mm
soil_index = lambda x, y, z : int(-10 * z)  # maps to p_s (hydrostatic equilibirum)
rs.setSoilGrid(soil_index)

soilcore = pb.SDF_PlantBox(1e6, 1e6, 149.9)
rs.setGeometry(soilcore)

path = "../../../../CPlantBox/modelparameter/rootsystem/"
name = "Glycine_max"  # Zea_mays_1_Leitner_2010
rs.setSeed(2)
rs.readParameters(path + name + ".xml")
rs.initialize()
rs.simulate(simtime, False)

""" set up xylem parameters """
r = XylemFluxPython(rs)
r.setKrTables([kr0[:, 1], kr1[:, 1], kr2[:, 1], kr3[:, 1]], [kr0[:, 0], kr1[:, 0], kr2[:, 0], kr3[:, 0]])
r.setKxTables([kz0[:, 1], kz1[:, 1], kz2[:, 1], kz3[:, 1]], [kz0[:, 0], kz1[:, 0], kz2[:, 0], kz3[:, 0]])

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
name = "soybean_Honly_2003"
Hseq_t = []
# Open .vtu
filelist = glob.iglob(r'../../../build-cmake/cpp/coupled_1p_richards/results_' + name + '/*.vtu')
for filepath in sorted(filelist):
    pd = vp.read_vtu(filepath)
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
    """ 3. EQUIVALENT SOIL WATER POTENTIAL """
    
    t = time.time()
    eswp = 0.
    n = len(ana.segments)
    seg2cell = r.rs.seg2cell
    
    for i in range(0, n):
    	eswp += suf[i] * data[seg2cell[i]]
    print("\nEquivalent soil water potential", eswp)
    elapsed = time.time() - t
    print("\nTime elapsed: ", elapsed)
    
    Hseq_t.append(eswp)
    print(filepath)
    
#    Figure
fig, ax1 = plt.subplots()
time = np.linspace(0,154,3697)						# (initial day, final day, number of vtps)
ax1.plot(time, Hseq_t, 'b-')    						# time vs equivalent soil pressure using SUF as a weighing factor 
ax1.set_xlabel("Time [days]")
ax1.set_ylabel("Equivalent soil water potential [cm]")
ax1.legend(loc = 'upper left')
fig.savefig("../../../build-cmake/cpp/coupled_1p_richards/results_" + name + "/" + name + "_equivalent_p.pdf", bbox_inches='tight')
plt.show()
