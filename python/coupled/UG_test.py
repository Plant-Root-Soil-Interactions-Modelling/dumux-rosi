import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

from functional.xylem_flux import *  # Python hybrid solver
import plantbox as pb
import rsml.rsml_reader as rsml
from rosi_richards import RichardsUG  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import visualisation.vtk_plot as vp
import functional.van_genuchten as vg
from functional.root_conductivities import *

import numpy as np
import matplotlib.pyplot as plt
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

""" 
Benchmark M1.2 static root system in soil (with the classic sink)

using an unstructured grid
"""

""" Parameters """
# not needed for UG, only for vg-plots
# min_b = [-4., -4., -15]
# max_b = [4., 4., 0.]
# cell_number = [8, 8, 15]  # [8, 8, 15]  # [8, 8, 15]  # [16, 16, 30]  # [32, 32, 60]  # [8, 8, 15]
periodic = False

name = "UG_test"
sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil = loam
initial = -659.8 + 7.5  # -659.8

trans = 6.4  # cm3 /day (sinusoidal) TESTING, DEFAULT IS 6.4!
wilting_point = -15000  # cm

sim_time = 0.1  # [day] for task b
age_dependent = False  # conductivities
dt = 120. / (24 * 3600)  # [days] Time step must be very small

""" Initialize macroscopic soil model """
sp = vg.Parameters(soil)  # for debugging
cpp_base = RichardsUG()
s = RichardsWrapper(cpp_base)
s.initialize()

# s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
# s.readGrid("../../grids/cylinder_tobi_fast_r2.81.msh")  # [cm] MRI_default_file (old dumux dimenions, dimensionsin vtu % vtp off by factor 100, root system is very large)
s.readGrid("../../grids/cylinder_tobi_fast_r2.81(x100).msh")  # [cm] MRI_default_file scaled x100, dimensions of vtu & vtp in paraview match, but initial pressure distribution (with equilibrium) is totally off

s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium Homogenious IC does not work with cylinder_tobi_fast_r2.81(x100).msh
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([soil])
s.setParameter("Newton.EnableChop", "True")
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.initializeProblem()
s.setCriticalPressure(wilting_point)

""" Initialize xylem model (a) or (b)"""
r = XylemFluxPython("../../grids/RootSystem8.rsml")
# r = XylemFluxPython("../../grids/III_Soil_1W-00000.rsml")") # vtp2rsml-script rsml does not work yet (too many orders, emergence time way off, etc)

# r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),  	# only for vp plot
#                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)  # cutting
kz = 4.32e-2
kr = 1.728e-4
r.setKr([kr])
r.setKx([kz])
# init_conductivities(r, age_dependent)

r.rs.sort()  # ensures segment is located at index s.y-1
r.test()  # sanity checks
nodes = r.get_nodes()
rs_age = np.max(r.get_ages())

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])
r.rs.setSoilGrid(picker)  # maps segments
cci = picker(nodes[0, 0], nodes[0, 1], nodes[0, 2])  # collar cell index

""" Numerical solution (a) """
start_time = timeit.default_timer()
x_, y_, w_, cpx, cps = [], [], [], [], []
sx = s.getSolutionHead()  # inital condition, solverbase.py
print(sx)
N = round(sim_time / dt)
t = 0.
test = 0
for i in range(0, N):

    if rank == 0:  # Root part is not parallel

        rx = r.solve(rs_age + t, -trans * sinusoidal(t), sx[cci], sx, True, wilting_point, [])  # xylem_flux.py, cells = True

        fluxes = r.soilFluxes(rs_age + t, rx, sx, False)  # class XylemFlux is defined in MappedOrganism.h, approx = True

        sum_flux = 0.
        for f in fluxes.values():
            sum_flux += f
        print("Summed fluxes ", sum_flux, "= collar flux", r.collar_flux(rs_age + t, rx, sx), "= prescribed", -trans * sinusoidal(t))

    else:
        fluxes = None

    fluxes = comm.bcast(fluxes, root = 0)  # Soil part runs parallel
    s.setSource(fluxes)  # richards.py

   #  s.ddt = dt / 10
    s.solve(dt)

    sx = s.getSolutionHead()  # richards.py
    water = s.getWaterVolume()

    if rank == 0:
        n = round(float(i) / float(N) * 100.)
        min_sx = np.min(sx)
        min_rx = np.min(rx)
        max_sx = np.max(sx)
        max_rx = np.max(rx)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g} days {:g}"
              .format(min_sx, max_sx, min_rx, max_rx, s.simTime, rx[0]))
        f = float(r.collar_flux(rs_age + t, rx, sx))  # exact root collar flux
        x_.append(t)
        y_.append(sum_flux)
        w_.append(water)
        cpx.append(rx[0])
        cps.append(float(sx[cci]))
        # print("Time:", t, ", collar flux", f, "cm^3/day at", rx[0], "cm xylem ", float(sx_old[cci]), "cm soil", "; domain water", s.getWaterVolume(), "cm3")

    test += 1
    t += dt

    s.writeDumuxVTK(name + str(test))
# s.writeDumuxVTK(name)

# vp-plot not functional or useful with UG simulation
"""
if rank == 0:
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")
    vp.plot_roots_and_soil(r.rs, "pressure head", rx, s, periodic, min_b, max_b, cell_number, name)  # VTK vizualisation
    fig, ax1 = plt.subplots()
    ax1.plot(x_, trans * sinusoidal(x_), 'k')  # potential transpiration
    ax1.plot(x_, -np.array(y_), 'g')  # actual transpiration (neumann)
    ax2 = ax1.twinx()
    ax2.plot(x_, np.cumsum(-np.array(y_) * dt), 'c--')  # cumulative transpiration (neumann)
    ax1.set_xlabel("Time [d]")
    ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
    ax1.legend(['Potential', 'Actual', 'Cumulative'], loc = 'upper left')
    np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter = ';')
    plt.show()
"""
