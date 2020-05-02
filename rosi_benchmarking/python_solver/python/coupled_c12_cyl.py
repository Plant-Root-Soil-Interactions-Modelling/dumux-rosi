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


def sinusoidal(t):
    return np.sin(2. * pi * np.array(t) - 0.5 * pi) + 1.

""" 
Benchmark M1.2 static root system in soil, coupled to cylindrical richards

also works parallel with mpiexec (only slightly faster, due to overhead)
"""

""" Parameters """
sim_time = 7  # [day] for task b
trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -10000  # cm
loam = [0.08, 0.43, 0.04, 1.6, 50]
periodic = False  

""" Root problem (a) or (b)"""
r = XylemFluxPython("../grids/RootSystem.rsml")
r.setKr([ 1.728e-4])  # [cm^3/day]
r.setKx([4.32e-2 ])  # [1/day]
nodes = r.get_nodes()

""" Soil problem """
cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()
min_b = [-4., -4., -20.]
max_b = [4., 4., 0.]
cell_number = [8, 8, 20]
s.createGrid(min_b, max_b, cell_number)  # [cm]
cell_volume = np.prod(np.array(max_b) - np.array(min_b)) / (np.prod(cell_number))
print("cell volume:", cell_volume, "cm3")
r.rs.setRectangularGrid(pb.Vector3d(-4., -4., -20.), pb.Vector3d(4., 4., 0.), pb.Vector3d(8, 8, 20))  # cut root segments to grid (segments are not mapped after)

s.setHomogeneousIC(-100, False)  # cm pressure head, equilibrium
# s.setHomogeneousIC(-669.8 - 10, True)  # cm pressure head, equilibrium

s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([loam])
s.initializeProblem()

""" Coupling (map indices) """
picker = lambda x, y, z : s.pick([x, y, z])
r.rs.setSoilGrid(picker)  # maps segments
cci = picker(nodes[0, 0], nodes[0, 1], nodes[0, 2])  # collar cell index

""" Set up cylindrical problems """
segments = r.rs.segments
cell2seg = r.rs.cell2seg  # cell to segments mapper
nodes = r.get_nodes()

rich_cyls = []
N = 20 

seg_radii = np.array(r.rs.radii)
seg_length = np.zeros(seg_radii.shape)
seg_prop = np.zeros(seg_radii.shape)  # proportionality factor
seg_outer_radii = np.zeros(seg_radii.shape)

# interate over the cells
for cellId in cell2seg.keys():  # create seg volumes per cell
    segs = cell2seg[cellId]    
    # print(cellId, ":", segs)
    v = 0.  # volume of all segments
    for i in segs:
        n1 = nodes[int(segments[i].x), :]
        n2 = nodes[int(segments[i].y), :]
        l = np.linalg.norm(n1 - n2)
        seg_length[i] = l
        a = seg_radii[i]
        v += 2 * np.pi * a * l
    for i in segs:        
        l = seg_length[i] 
        a = seg_radii[i]
        seg_prop[i] = (2 * np.pi * a * l) / v
        tv = seg_prop[i] * cell_volume  # target volume                
        a_out = (tv + (2 * np.pi * a * l)) / (2 * np.pi * l) 
        # print("target:", tv, "inner", a, "outer", a_out)
        seg_outer_radii[i] = a_out

for i, _ in enumerate(segments):    
    
    a_in = seg_radii[i]
    a_out = seg_outer_radii[i]
    
    cpp_base = RichardsCylFoam()
    rcyl = RichardsWrapper(cpp_base)
    
    rcyl.initialize([""], False)  # [""], False
    try:
        rcyl.createGrid([a_in], [a_out], [N])  # [cm]
        rcyl.setHomogeneousIC(-100.)  # cm pressure head
        rcyl.setOuterBC("fluxCyl", 0.)  #  [cm/day]
        rcyl.setInnerBC("fluxCyl", 0.)  # [cm/day] 
        rcyl.setVGParameters([loam])
        rcyl.initializeProblem()
        rcyl.setCriticalPressure(-15000)  # cm pressure head                
        rich_cyls.append(rcyl)    
    except:
        rich_cyls.append([])  # in case a_out = 0, e.g. if segment does not lie within the domain
        print("\nError", a_in, a_out, "\n")        

""" Numerical solution (a) """
start_time = timeit.default_timer()
x_, y_, w_, cpx, cps = [], [], [], [], []
sx = s.getSolutionHead()  # inital condition, solverbase.py
rsx = np.zeros(seg_radii.shape)  # root soil interface

dt = 120. / (24 * 3600)  # [days] Time step must be very small
N = sim_time * round(1. / dt)
t = 0.

for i in range(0, N):

    if rank == 0:  # Root part is not parallel
        
        for j, rc in enumerate(rich_cyls):
            if not isinstance(rc, list):                
                rsx[j] = rc.getInnerHead()
        
        # user rsx as soil pressures around the individual segments (cells = False)
        rx = r.solve(t, -trans * sinusoidal(t), sx[cci], rsx, False, wilting_point)  # xylem_flux.py        
        
        # For the soil model 
        seg_fluxes = r.segFluxes(t, rx, rsx, approx=False)  # class XylemFlux is defined in CPlantBox XylemFlux.h
        soil_fluxes = r.sumSoilFluxes(seg_fluxes)  # class XylemFlux is defined in CPlantBox XylemFlux.h

        # set seg_fluxes
        for j, rc in enumerate(rich_cyls):
            if not isinstance(rc, list):                
                rsx[j] = rc.setInnerBC("flux", seg_fluxes[j])  # TODO units

        sum_flux = 0.
        for f in soil_fluxes.values():
            sum_flux += f
        print("Fuxes ", sum_flux, "= prescribed", -trans * sinusoidal(t) , "= collar flux", r.collar_flux(0., rx, sx))

    else:
        fluxes = None

    fluxes = comm.bcast(fluxes, root=0)  # Soil part coult runs parallel

    s.setSource(soil_fluxes)  # richards.py
    s.solve(dt)

    sx = s.getSolutionHead()  # richards.py
    water = s.getWaterVolume()

    if rank == 0:
        n = round(float(i) / float(N) * 100.)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "]")
        f = float(r.collar_flux(t, rx, sx))  # exact root collar flux
        x_.append(t)
        y_.append(f)
        w_.append(water)
        cpx.append(rx[0])
        cps.append(float(sx[cci]))

        # print("Time:", t, ", collar flux", f, "cm^3/day at", rx[0], "cm xylem ", float(sx_old[cci]), "cm soil", "; domain water", s.getWaterVolume(), "cm3")

    t += dt

s.writeDumuxVTK("c12_final")

""" Plot """
if rank == 0:
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")
    fig, ax1 = plt.subplots()
    ax1.plot(x_, trans * sinusoidal(x_), 'k')  # potential transpiration
    ax1.plot(x_, -np.array(y_), 'g')  # actual transpiration (neumann)
    ax2 = ax1.twinx()
    ax2.plot(x_, np.cumsum(-np.array(y_) * dt), 'c--')  # cumulative transpiration (neumann)
    ax1.set_xlabel("Time [d]")
    ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
    ax1.legend(['Potential', 'Actual', 'Cumulative'], loc='upper left')
    plt.show()

