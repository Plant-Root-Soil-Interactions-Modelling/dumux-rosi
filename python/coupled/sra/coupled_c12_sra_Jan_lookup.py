import sys; sys.path.append("../../modules/"); sys.path.append("../../../../CPlantBox/");  sys.path.append("../../../../CPlantBox/src/python_modules")
sys.path.append("../../../build-cmake/cpp/python_binding/"); sys.path.append("../../modules/fv/");

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml_reader as rsml
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import vtk_plot as vp
import van_genuchten as vg
from root_conductivities import *
from rhizo_models import plot_transpiration

from math import *
import numpy as np
import matplotlib.pyplot as plt
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
from scipy.optimize import fsolve

from sra_table_lookup import *  # <------- new


def sinusoidal(t):
    return np.sin(2. * pi * np.array(t) - 0.5 * pi) + 1.


def soil_root_interface_table(rx, sx, r, sp, outer_r, f):
    assert rx.shape == sx.shape          
    inner_ = np.zeros(rx.shape)
    outer_ = np.zeros(rx.shape)
    for i in range(0, len(rx)):
        inner_[i] = max(min(r.rs.radii[i] , 0.2), 0.01)
        outer_[i] = max(min(outer_r[i] , 20), 0.1)      
    rsx = f((rx, sx, inner_ , outer_))                                                      
    return rsx

""" 
Benchmark M1.2 static root system in soil (with Jans new SRA sink)
"""

""" Parameters """
min_b = [-4., -4., -15.]
max_b = [4., 4., 0.]
cell_number = [7, 7, 15]  #  [8, 8, 15]  # [16, 16, 30]  # [32, 32, 60]  # [8, 8, 15]
periodic = False
fname = "../../../grids/RootSystem8.rsml"

name = "jan_lookup2"
loam = [0.08, 0.43, 0.04, 1.6, 50]
soil = loam

initial = -659.8 + 7.5  # -659.8

trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

sim_time = 7  # [day] for task b
dt = 60. / (24 * 3600)  # [days] Time step must be very small

skip = 1  # output

""" Initialize macroscopic soil model """
sp = vg.Parameters(soil)  
vg.create_mfp_lookup(sp, -1.e5, 1000)
s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([soil])
s.setParameter("Newton.EnableChop", "True")
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
# s.setParameter("Soil.SourceSlope", "1000")  # turns regularisation of the source term on
s.initializeProblem()
s.setCriticalPressure(wilting_point)

""" Initialize xylem model (a) or (b)"""
r = XylemFluxPython(fname)
r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)  # cutting
init_conductivities(r, False)

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])    
cci = picker(r.rs.nodes[0].x, r.rs.nodes[0].y, r.rs.nodes[0].z)  # collar cell index
r.rs.setSoilGrid(picker)  # maps segment
outer_r = r.rs.segOuterRadii() 

""" sanity checks """
# r.plot_conductivities()
r.test()  # sanity checks
rs_age = np.max(r.get_ages())
seg_length = r.segLength()
ns = len(seg_length)
# print("press any key"); input()
print("outer radii", np.min(outer_r) , np.max(outer_r))
radii = r.rs.radii
print("inner radii", np.min(radii) , np.max(radii))

sra_table_lookup = open_sra_lookup("table")

# quick check
# rsx2 = soil_root_inerface(np.array([-15000]), np.array([-700]), r, sp, outer_r)   
# print(r.rs.radii[0])
# print(outer_r[0])
# rsx3 = sra_table_lookup((-15000, -700, 0.1679960056208074, 0.6952332821448589))   
# print(rsx2, rsx3)        

""" Numerical solution (a) """
start_time = timeit.default_timer()
x_, y_, w_, cpx, cps, cf = [], [], [], [], [], []
sx = s.getSolutionHead()  # inital condition, solverbase.py
hsb = np.array([sx[r.rs.seg2cell[j]][0] for j in range(0, ns)])  # soil bulk matric potential per segment
rsx = hsb.copy()  # initial values for fix point iteration 

N = round(sim_time / dt)
t = 0.

for i in range(0, N):

    if rank == 0:  # root part is not parallel
 
        if i == 0:  # only first time
            rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, [])  # xylem_flux.py, cells = False        
 
        rsx = soil_root_interface_table(rx[1:], hsb, r, sp, outer_r, sra_table_lookup)
        rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, [])  # xylem_flux.py, cells = False        
        rx_old = rx.copy()
        
        err = 1.e6
        c = 1
        while err > 1. and c < 100:
            rsx = soil_root_interface_table(rx[1:], hsb, r, sp, outer_r, sra_table_lookup)
            rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, [])  # xylem_flux.py, cells = False
            err = np.linalg.norm(rx - rx_old)
            # print(err)
            rx_old = rx.copy()
            c += 1
        
        print(c, "iterations")
        
        # rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, [])  # xylem_flux.py, cells = False
        fluxes = r.segFluxes(rs_age + t, rx, rsx, False)  # class XylemFlux is defined in MappedOrganism.h, approx = True
        min_rsx = np.min(rsx)  # for console output
        max_rsx = np.max(rsx)        
 
#         sum_flux = 0.
#         for f in fluxes.values():
#             sum_flux += f        
#         print(sum_flux, r.collar_flux(rs_age + t, rx, rsx))    
 
    else:
        fluxes = None
     
    fluxes = comm.bcast(r.sumSoilFluxes(fluxes), root=0)  # Soil part runs parallel
    s.setSource(fluxes.copy())  # richards.py 
    s.solve(dt)
    sx = s.getSolutionHead()  # richards.py
    hsb = np.array([sx[r.rs.seg2cell[j]][0] for j in range(0, ns)])  # soil bulk matric potential per segment    
    water = s.getWaterVolume()
 
    if rank == 0 and i % skip == 0:
        min_sx = np.min(sx)
        max_sx = np.max(sx)
        min_rx = np.min(rx)
        max_rx = np.max(rx)                
        x_.append(t)
        sum_flux = 0.
        for f in fluxes.values():
            sum_flux += f
        cf_ = r.collar_flux(rs_age + t, rx, rsx, k_soil=[], cells=False)
        print("Summed fluxes ", sum_flux, "= collar flux", cf_, "= prescribed", -trans * sinusoidal(t))
        y_.append(sum_flux)  # cm4/day
        w_.append(water)  # cm3
        cf.append(cf_)  # cm3/day
        cpx.append(rx[0])  # cm
        cps.append(float(sx[cci]))  # cm
        n = round(float(i) / float(N) * 100.)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil sx [{:g}, {:g}], interface [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days\n"
              .format(min_sx, max_sx, min_rsx, max_rsx, min_rx, max_rx, s.simTime))

    t += dt

s.writeDumuxVTK(name)

""" Plot """
if rank == 0:
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")
    
    vp.plot_roots_and_soil(r.rs, "pressure head", rx, s, periodic, min_b, max_b, cell_number, name)  # VTK vizualisation
    
    ana = pb.SegmentAnalyser(r.rs)
    ana.addData("pressure", rx)
    vp.plot_roots(ana, "pressure")
    
    plot_transpiration(x_, y_, cf, lambda t: trans * sinusoidal(t))
    np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter=';')

