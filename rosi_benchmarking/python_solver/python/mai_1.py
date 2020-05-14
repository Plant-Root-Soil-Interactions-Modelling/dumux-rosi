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

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

""" 
Mai et al (2019) scenario 1 water movement  
"""
N = 3  # number of cells in each dimension 
loam = [0.08, 0.43, 0.04, 1.6, 50]
initial = -100.

r_root = 0.02  # [cm]
kr = 2.e-13 * 1000 * 9.81  # [m / (Pa s)] -> [ 1 / s ]
kx = 5.e-17 * 1000 * 9.81  # [m^4 / (Pa s)] -> [m3 / s] 
kr = kr * 24 * 3600  # [ 1 / s ] -> [1/day]
kx = kx * 1.e6 * 24 * 3600  # [ m3 / s ] -> [cm3/day]
# print(kr, kx)

NC = 10  # dof of the cylindrical problem
logbase = 1.5

q_r = 1.e-5 * 24 * 3600  # [cm / s] -> [cm / day] 
sim_time = 20  # 20  # days
NT = 100  # iteration

""" Soil problem """
cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()    
s.setHomogeneousIC(initial)  # cm pressure head
s.setTopBC("noflux")
s.setBotBC("noflux")
s.createGrid([-1.5, -1.5, -3.], [1.5, 1.5, 0.], [N, N, N])  # [cm] 3x3x3
s.setVGParameters([loam])
s.initializeProblem()
s.setCriticalPressure(-15000)
s.ddt = 1.e-5  # initial Dumux time step [days]
 
""" Root problem"""
n, segs = [], []
for i in range(0, 4):  # nodes
    n.append(pb.Vector3d(0, 0, -float(i)))
for i in range(0, len(n) - 1):  # segments
    segs.append(pb.Vector2i(i, i + 1))
rs = pb.MappedSegments(n, segs, [r_root] * len(segs))  # a single root
rs.setRectangularGrid(pb.Vector3d(-1.5, -1.5, -3.), pb.Vector3d(1.5, 1.5, 0.), pb.Vector3d(N, N, N))
r = XylemFluxPython(rs)
r.setKr([kr])
r.setKx([kx])

r_outer = r.segOuterRadii()
seg_length = r.segLength()
# print(r_outer); print(seg_length)

""" Coupling (map indices) """
picker = lambda x, y, z : s.pick([x, y, z])
r.rs.setSoilGrid(picker)
cci = picker(0, 0, 0)  # collar cell index

""" Cylindrical problems """
cyls = []
ns = len(seg_length)  # number of segments 
for i in range(0, ns):     
    cpp_base = RichardsCylFoam()
    cyl = RichardsWrapper(cpp_base)    
    cyl.initialize()    
    points = np.logspace(np.log(r_root) / np.log(logbase), np.log(r_outer[i]) / np.log(logbase), NC, base=logbase)
    # plt.plot(points, [0.] * len(points), "b*"); plt.show()
    cyl.createGrid1d(points)                 
    cyl.setHomogeneousIC(initial)  # cm pressure head
    cyl.setVGParameters([loam])
    cyl.setOuterBC("fluxCyl", 0.)  # [cm/day]
    cyl.setInnerBC("fluxCyl", 0.)  # [cm/day]         
    cyl.initializeProblem()
    cyl.setCriticalPressure(-15000.)  # cm pressure head                
    cyl.setRegularisation(1.e-4, 1.e-4)
    cyl.ddt = 1.e-5  # [day] initial time step                                
    cyls.append(cyl)    

""" Simulation """
rsx = np.zeros((ns,))  # xylem pressure at the root soil interface
dt = sim_time / NT 
uptake = []

for _ in range(0, NT):

    # solutions of previous time step
    sx = s.getSolutionHead()         
    for j, rc in enumerate(cyls):  # for each segment
        rsx[j] = rc.getInnerHead()            
                       
     # solves root model                  
    print(rsx)     
    rx = r.solve(0., -q_r, sx[cci], rsx, False, -15000)                    
    print("Minimal root xylem pressure", np.min(np.array(rx)))  # working
                  
    # fluxes between 1d - 3d                         
    seg_fluxes = r.segFluxes(0., rx, rsx, approx=False)  
    soil_fluxes = r.sumSoilFluxes(seg_fluxes)                                   
             
    sum_flux = 0.
    for f in soil_fluxes.values():
        sum_flux += f
    print("Summed fluxes {:g}, predescribed {:g}".format(sum_flux, -q_r))
    uptake.append(sum_flux)
            
    # run cylindrical model            
    for j, rc in enumerate(cyls):  # set cylindrical model fluxes
        l = seg_length[j]
        print("Set inner flux to", seg_fluxes[j] / (2 * np.pi * r_root * l), "[cm day-1]")  
        rc.setInnerFluxCyl(seg_fluxes[j] / (2 * np.pi * r_root * l))    
                                            
    for j, rc in enumerate(cyls):  # simualte time step
        print("Solving segment", j)
        rc.ddt = 1.e-5  
        rc.solve(dt)
    
    # run macroscopic soil model
    s.setSource(soil_fluxes)  # g day-1, richards.py
    s.solve(dt)    
                
    x0 = s.getSolutionHeadAt(cci)
    print("Matric potential in collar cell" , x0, " cm")

plt.title("Water uptake")
plt.plot(np.linspace(0, sim_time, NT), -np.array(uptake))
plt.xlabel("Time (days)")
plt.ylabel("Uptake (cm/day)")
plt.show()                 

print("fin")
