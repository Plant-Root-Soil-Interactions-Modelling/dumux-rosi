import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src/python_modules/") 
from cmath import isnan
sys.path.append("../../build-cmake/cpp/python_binding/")

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml_reader as rsml
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import vtk_plot as vp
import van_genuchten as vg
from root_conductivities import *

from math import *
import numpy as np
import matplotlib.pyplot as plt
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()


def sinusoidal(t):
    return np.sin(2. * pi * np.array(t) - 0.5 * pi) + 1.

""" 
Benchmark M1.2 static root system in soil

The classic sink is decoupled from movement, 
i.e. water movement is calculated first, in a second step water is taken up by the roots (classical sink)
python implementation for developing RichardsWrapper.applySink 

MPI todo
"""

""" Parameters """
min_b = [-4., -4., -15.]
max_b = [4., 4., 0.]
cell_number = [7, 7, 15]  #  [8, 8, 15]  # [16, 16, 30]  # [32, 32, 60]  # [8, 8, 15]
# [7, 7, 15], [14,14,30] works ???
periodic = False

name = "DuMux_1cm"
sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil = loam

initial = -659.8 + 7.5  # -659.8

trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

sim_time = 3  # [day] for task b
age_dependent = False  # conductivities
dt = 120 / (24 * 3600)  # [days] Time step must be very small
skip = 1  

""" Initialize macroscopic soil model """
sp = vg.Parameters(soil)  
s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([soil])
s.setParameter("Newton.EnableChop", "True")
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.initializeProblem()
s.setCriticalPressure(wilting_point)

""" Initialize xylem model (a) or (b)"""
# fname = "../../../CPlantBox/tutorial/examples/python/results/example_3c.rsml"   
fname = "../../grids/RootSystem8.rsml"
r = XylemFluxPython(fname)
r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)  # cutting
init_conductivities(r, age_dependent)
# r.plot_conductivities()
r.test()  # sanity checks
rs_age = np.max(r.get_ages())
# print("press any key"); input()

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])    
cci = picker(r.rs.nodes[0].x, r.rs.nodes[0].y, r.rs.nodes[0].z)  # collar cell index
r.rs.setSoilGrid(picker)  # maps segment

""" Numerical solution (a) """
start_time = timeit.default_timer()
x_, y_, y2_, cpx, cps, cf = [], [], [], [], [], []
sx = s.getSolutionHead()  # inital condition, solverbase.py

N = round(sim_time / dt)

t = 0.
vols = s.getCellVolumes()  # cell volumes [cm3] 

for i in range(0, N):

    # old time step
    old_water_volume = s.getWaterVolume()
    sx_k = s.getSolutionHead()
    w_k = np.multiply(vols, s.getWaterContent())

    sx = sx_k.copy()  # initial value
    w = w_k.copy()
    
    iter = 0
    while True:

        # f(sx,trans)
        if rank == 0: 
            rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., sx, True, wilting_point, [])  # xylem_flux.py, cells = True
            fluxes = r.soilFluxes(rs_age + t, rx, sx, False)  # class XylemFlux is defined in MappedOrganism.h, approx = True
        else:
            fluxes = None     
        fluxes = comm.bcast(fluxes, root=0)  # Soil part runs parallel

        # sink(sx, trans)
        sink = np.zeros(w.shape)  # convert fluxes to sink    
        for key, value in fluxes.items():
            sink[key] += value
        
        w_old = w          
        w = w_k + dt * sink         
        eps = np.linalg.norm(w_old - w)  # other norm?!
        # eps = np.max(np.abs(w_old - w))
                          
        wc = np.divide(w, vols)
        sx = vg.pressure_head(wc, sp)           
        sx = np.maximum(sx, wilting_point)  # only to avoid nan 
            
        iter += 1
        print(eps)
        if eps < 1.e-5 or iter > 10: 
            if isnan(eps):
                input() 
            
            break
        
    sx = np.maximum(sx, wilting_point)    
    s.setInitialCondition(sx)  # need to revise method
    s.solve(dt)    

    uptake = float((s.getWaterVolume() - old_water_volume) / dt)  # the actual uptake in this time step 

    if rank == 0 and i % skip == 0:
        x_.append(t)
        min_sx = np.min(sx)
        max_sx = np.max(sx)
        min_rx = np.min(rx)        
        max_rx = np.max(rx)   
        sum_flux = 0.
        for f in fluxes.values():
            sum_flux += f                             
        print("Uptake {:g} summed fluxes {:g} = collar flux {:g} = prescribed {:g}"
              .format(uptake, sum_flux, float(r.collar_flux(rs_age + t, rx, sx)), -trans * sinusoidal(t)))
        y_.append(uptake)  # cm3/day
        y2_.append(sum_flux)  # cm3/day
        cf.append(float(r.collar_flux(rs_age + t, rx, sx)))  # cm3/day
        cpx.append(rx[0])  # cm
        cps.append(float(sx[cci]))  # cm
        n = round(float(i) / float(N) * 100.)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days {:g}\n"
              .format(min_sx, max_sx, min_rx, max_rx, s.simTime, rx[0]))

    t += dt

s.writeDumuxVTK(name)

""" Plot """
if rank == 0:
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")
    # vp.plot_roots_and_soil(r.rs, "pressure head", rx, s, periodic, min_b, max_b, cell_number, name)  # VTK vizualisation
    fig, ax1 = plt.subplots()
    ax1.plot(x_, trans * sinusoidal(x_), 'k')  # potential transpiration
    ax1.plot(x_, -np.array(y_), 'g')  # actual transpiration (neumann)
    ax1.plot(x_, -np.array(y2_), 'b:')  # actual transpiration (neumann)
    ax2 = ax1.twinx()
    ax2.plot(x_, np.cumsum(-np.array(y_) * dt), 'c--')  # cumulative transpiration (neumann)
    ax1.set_xlabel("Time [d]")
    ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
    ax1.legend(['Potential', 'Actual', 'Root fluxes', 'Cumulative'], loc='upper left')
    np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter=';')
    plt.show()

