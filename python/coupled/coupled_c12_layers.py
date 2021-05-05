import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src/python_modules/") 
sys.path.append("../../build-cmake/cpp/python_binding/")

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


def sinusoidal(t):
    return np.sin(2. * pi * np.array(t) - 0.5 * pi) + 1.

""" 
Benchmark M1.2 static root system in soil (with the classic sink)

also works parallel with mpiexec (slower, due to overhead?)
"""

""" Parameters """
min_b = [-4., -4., -15.]
max_b = [4., 4., 0.]
cell_number = [7, 7, 15]  #  [8, 8, 15]  # [16, 16, 30]  # [32, 32, 60]  # [8, 8, 15]
periodic = False
fname = "../../grids/RootSystem8.rsml"

# min_b = [-6.25, -1.5, -180.]  # cm
# max_b = [6.25, 1.5, 1]  # cm
# cell_number = [13, 3, 180]
# periodic = True
# fname = "spring_barley_CF12_107d.rsml"

name = "DuMux_1cm"
sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil = loam

initial = -200  # -659.8

trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

sim_time = 3  # [day] for task b
age_dependent = False  # conductivities
dt = 30. / (24 * 3600)  # [days] Time step must be very small
skip = 1

""" Initialize macroscopic soil model """
sp = vg.Parameters(soil)  # for debugging
cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setLayersZ([2, 2, 1, 1], [-15, -7., -7., 0.])  # sample points ([1], [cm])
s.setVGParameters([loam,sand])
s.setParameter("Newton.EnableChop", "True")
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Soil.SourceSlope", "2000")  # turns regularisation of the source term on
s.initializeProblem()
s.setCriticalPressure(wilting_point)

""" Initialize xylem model (a) or (b)"""

r = XylemFluxPython(fname)
r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)  # cutting
init_conductivities(r, age_dependent)

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])    
cci = picker(r.rs.nodes[0].x, r.rs.nodes[0].y, r.rs.nodes[0].z)  # collar cell index
r.rs.setSoilGrid(picker)  # maps segment

""" sanity checks """
# r.plot_conductivities()
r.test()  # sanity checks
rs_age = np.max(r.get_ages())
# print("press any key"); input()

""" Numerical solution (a) """
start_time = timeit.default_timer()
x_, y_, w_, cpx, cps, cf = [], [], [], [], [], []
sx = s.getSolutionHead()  # inital condition, solverbase.py

N = round(sim_time / dt)
t = 0.

for i in range(0, N):

    if rank == 0:  # Root part is not parallel
 
        rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., sx, True, wilting_point, [])  # xylem_flux.py, cells = True
        fluxes = r.soilFluxes(rs_age + t, rx, sx, False)  # class XylemFlux is defined in MappedOrganism.h, approx = True
 
    else:
        fluxes = None
 
    fluxes = comm.bcast(fluxes, root=0)  # Soil part runs parallel
    s.setSource(fluxes.copy())  # richards.py 
    s.solve(dt)
    sx = s.getSolutionHead()  # richards.py
    water = s.getWaterVolume()
 
    if rank == 0 and i % skip == 0:
        min_sx = np.min(sx)
        min_rx = np.min(rx)
        max_sx = np.max(sx)
        max_rx = np.max(rx)                
        x_.append(t)
        sum_flux = 0.
        for f in fluxes.values():
            sum_flux += f
        print("Summed fluxes ", sum_flux, "= collar flux", r.collar_flux(rs_age + t, rx, sx), "= prescribed", -trans * sinusoidal(t))
        y_.append(sum_flux)  # cm4/day
        w_.append(water)  # cm3
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
    vp.plot_roots_and_soil(r.rs, "pressure head", rx, s, periodic, min_b, max_b, cell_number, name)  # VTK vizualisation
    plot_transpiration(x_, y_, cf, lambda t: trans * sinusoidal(t))
    np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter=';')

