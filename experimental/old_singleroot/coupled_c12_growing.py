import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/"); 

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
similar to Benchmark M1.2 but with simulatied growing root system in soil (with the classic sink)

NOT working in parallel, needs revisions

also works parallel with mpiexec (only slightly faster, due to overhead)
"""

""" Parameters """
min_b = [-4., -4., -25.]
max_b = [4., 4., 0.]
cell_number = [7, 7, 53]  # [8, 8, 15]  # [16, 16, 30]  # [32, 32, 60]  # [8, 8, 15]
periodic = False

path = "../../../CPlantBox//modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
loam = [0.08, 0.43, 0.04, 1.6, 50]
initial = -659.8 + 7.5  # -659.8

trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

sim_time = 2  # [day] for task b
rs_age = 7
age_dependent = False  # conductivities
dt = 120. / (24 * 3600)  # [days] Time step must be very small

""" Initialize macroscopic soil model """
s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([loam])
s.initializeProblem()
s.setCriticalPressure(wilting_point)

""" Initialize root system """
rs = pb.MappedRootSystem()
rs.readParameters(path + name + ".xml")
if not periodic:
    sdf = pb.SDF_PlantBox(0.99 * (max_b[0] - min_b[0]), 0.99 * (max_b[1] - min_b[1]), max_b[2] - min_b[2])
    rs.setGeometry(sdf)

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])
rs.setSoilGrid(picker)  # maps segments
rs.setRectangularGrid(pb.Vector3d(min_b), pb.Vector3d(max_b), pb.Vector3d(cell_number), True)

""" initialize xylem model """
rs.initialize()
rs.simulate(rs_age, False)
r = XylemFluxPython(rs)
init_conductivities(r, age_dependent)
nodes = r.get_nodes()
cci = picker(nodes[0, 0], nodes[0, 1], nodes[0, 2])  # collar cell index

""" Numerical solution (a) """
start_time = timeit.default_timer()
x_, y_, w_, cpx, cps = [], [], [], [], []
sx = s.getSolutionHead()  # inital condition, solverbase.py
N = round(sim_time / dt)
t = 0.

# anim = vp.AnimateRoots(rs)

for i in range(0, N):

    if rank == 0:  # Root part is not parallel
        rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., sx, True, wilting_point, [])  # xylem_flux.py
        fluxes = r.soilFluxes(rs_age + t, rx, sx, approx=False)  # class XylemFlux is defined in MappedOrganism.h
#         x_.append(t)
#         y_.append(float(r.collar_flux(rs_age + t, rx, sx)))  # exact root collar flux
    else:
        fluxes = None

    fluxes = comm.bcast(fluxes, root=0)  # Soil part runs parallel
    s.setSource(fluxes)  # richards.py
    s.solve(dt)

    sx = s.getSolutionHead()  # richards.py
    water = s.getWaterVolume()

    if rank == 0:
        rs.simulate(dt, False)

#     if i % 10 == 0:
#         anim.update()

    if rank == 0:  # just output
        n = round(float(i) / float(N) * 100.)
        min_sx, min_rx, max_sx, max_rx = np.min(sx), np.min(rx), np.max(sx), np.max(rx)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g} days {:g}"
              .format(min_sx, max_sx, min_rx, max_rx, s.simTime, rx[0]))
        f = float(r.collar_flux(rs_age + t, rx, sx))  # exact root collar flux
        x_.append(t)
        y_.append(f)
        w_.append(water)
        cpx.append(rx[0])
        cps.append(float(sx[cci]))

        # print("Time:", t, ", collar flux", f, "cm^3/day at", rx[0], "cm xylem ", float(sx_old[cci]), "cm soil", "; domain water", s.getWaterVolume(), "cm3")

    t += dt

s.writeDumuxVTK(name)

""" Plot """
if rank == 0:
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    # VTK vizualisation
    vp.plot_roots_and_soil(r.rs, "pressure head", rx, s, periodic, min_b, max_b, cell_number, name)

    fig, ax1 = plt.subplots()
    ax1.plot(x_, trans * sinusoidal(x_), 'k')  # potential transpiration
    ax1.plot(x_, -np.array(y_), 'g')  # actual transpiration (neumann)
    ax2 = ax1.twinx()
    ax2.plot(x_, np.cumsum(-np.array(y_) * dt), 'c--')  # cumulative transpiration (neumann)
    ax1.set_xlabel("Time [d]")
    ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
    ax1.legend(['Potential', 'Actual', 'Cumulative'], loc='upper left')
    np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter=';')
    plt.show()

