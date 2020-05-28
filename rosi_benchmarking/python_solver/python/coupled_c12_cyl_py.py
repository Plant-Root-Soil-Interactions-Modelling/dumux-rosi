import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

from solver.xylem_flux import XylemFluxPython  # Python hybrid solver
import solver.plantbox as pb
import solver.rsml_reader as rsml
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding)
from solver.richards import RichardsWrapper  # Python part

import van_genuchten as vg

from math import *
import numpy as np
import matplotlib.pyplot as plt
import timeit


def sinusoidal(t):
    return np.sin(2. * pi * np.array(t) - 0.5 * pi) + 1.

""" 
Benchmark M1.2 static root system in soil, coupled to cylindrical richards

also works parallel with mpiexec (only slightly faster, due to overhead)
"""

""" Parameters """
min_b = [-5., -5., -20.]
max_b = [5., 5., 0.]
cell_number = [10, 10, 20]
loam = [0.08, 0.43, 0.04, 1.6, 50]
initial = -100

kx = 4.32e-2 
kr = 1.728e-4
trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

NC = 10  # dof of the cylindrical problem
logbase = 1.5

periodic = False  
 
sim_time = 1  # [day]
NT = 1000  # iteration

sim_time = 1  # [day] for task b

""" Initialize macroscopic soil model """
cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()
s.createGrid(min_b, max_b, cell_number)  # [cm]
s.setHomogeneousIC(initial, False)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([loam])
s.initializeProblem()
s.ddt = 1.e-5  # [day] initial Dumux time step 

""" Initialize xylem model (a)"""
r = XylemFluxPython("../grids/RootSystem.rsml")
r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]))  
r.setKr([kr])  # [cm^3/day]
r.setKx([kx])  # [1/day]

picker = lambda x, y, z : s.pick([x, y, z])
r.rs.setSoilGrid(picker)  # maps segments
cci = picker(nodes[0, 0], nodes[0, 1], nodes[0, 2])  # collar cell index

inner_radii = np.array(r.rs.radii)
outer_radii = r.segOuterRadii()
seg_length = r.segLength()
nodes = r.get_nodes()

""" Cylindrical models (around each root segment) """
cyls = []
ns = len(seg_length)  # number of segments 
print(ns); input()
for i in range(0, ns): 
    a_in = inner_radii[i]
    a_out = outer_radii[i]    
    if a_in < a_out:        
        cpp_base = RichardsCylFoam()
        cyl = RichardsWrapper(cpp_base)    
        cyl.initialize()    
        # plt.plot(points, [0.] * len(points), "b*"); plt.show()
        points = np.logspace(np.log(a_in) / np.log(logbase), np.log(a_out) / np.log(logbase), NC, base=logbase)
        cyl.createGrid1d(points)                 
        cyl.setHomogeneousIC(initial)  # cm pressure head
        cyl.setVGParameters([loam])
        cyl.setOuterBC("fluxCyl", 0.)  # [cm/day]
        cyl.setInnerBC("fluxCyl", 0.)  # [cm/day]         
        cyl.initializeProblem()
        cyl.setCriticalPressure(wilting_point_cyl)  # cm pressure head                                             
        cyls.append(cyl)  
    else:
        cyls.append([])
        print("Segment", i, "[", a_in, a_out, "]")  # this happens if elements are not within the domain
        input()

""" Simulation """
rsx = np.zeros((ns,))  # xylem pressure at the root soil interface
dt = sim_time / NT 

water_domain = []

cell_volumes = s.getCellVolumes()
net_flux = np.zeros(cell_volumes.shape)

for i in range(0, NT):

    t = i * dt

    # solutions of previous time step
    sx = s.getSolutionHead()  # [cm]         
    for j, rc in enumerate(cyls):  # for each segment
        rsx[j] = rc.getInnerHead()  # [cm]            
                       
     # solves root model                      
    rx = r.solve(t, -trans * sinusoidal(t), sx[cci], rsx, False, wilting_point)  # [cm]   
#     min_rx.append(np.min(np.array(rx)))                
#     print("Minimal root xylem pressure", min_rx[-1])  # working
                  
    # fluxes per segment according to root system model and cylindrical models                       
    seg_fluxes = r.segFluxes(0., rx, rsx, approx=False)  # [cm3/day] 
                                        
    # water 
    soil_water = np.multiply(np.array(s.getWaterContent()), cell_volumes)
    water_domain.append(np.sum(soil_water))
            
    seg_outer_fluxes = r.splitSoilFluxes(net_flux / dt)
            
    # run cylindrical models           
    for j, rc in enumerate(cyls):  
        l = seg_length[j]
        rc.setInnerFluxCyl(seg_fluxes[j] / (2 * np.pi * inner_radii[j] * l))  # [cm3/day] -> [cm /day]                                                                    
        rc.setOuterFluxCyl(seg_outer_fluxes[j] / (2 * np.pi * outer_radii[j] * l))  # [cm3/day] -> [cm /day]                                                                    
        rc.ddt = 1.e-5  # [day] initial time step  
        try:
            rc.solve(dt)
        except:
            x = rc.getDofCoordinates()
            y = rc.getSolutionHead()
            plt.plot(x, y)
            plt.xlabel("x (cm)")
            plt.ylabel("pressure (cm)")
            plt.show()
                
        seg_fluxes[j] = -rc.getInnerFlux() * (2 * np.pi * inner_radii[j] * l) / inner_radii  # [cm/day] -> [cm3/day], ('/inner_radii' comes from cylindrical implementation)   
    
    print("one iteration done")
    input()

    soil_fluxes = r.sumSoilFluxes(seg_fluxes)  # [cm3/day] 
    print("here we go"); input()
        
    # run macroscopic soil model
    s.setSource(soil_fluxes.copy())  # [cm3/day], richards.py
    print("here we go"); input()
    
    s.solve(dt)    
    print("here we go"); input()

    # calculate net fluxes
    net_flux = (np.multiply(np.array(s.getWaterContent()), cell_volumes) - soil_water) 
    for k, v in soil_fluxes.items():
        net_flux[k] -= v * dt;    
    print("sum", np.sum(net_flux))

