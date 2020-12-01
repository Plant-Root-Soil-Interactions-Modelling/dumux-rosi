import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src/python_modules/") 
sys.path.append("../../build-cmake/cpp/python_binding/")

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml_reader as rsml
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import vtk_plot as vp
import vtk_tools as vt
import van_genuchten as vg
from root_conductivities import *

from math import *
import numpy as np
import matplotlib.pyplot as plt
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()


""" 
Benchmark M1.2 static root system in soil (with the classic sink)

also works parallel with mpiexec (only slightly faster, due to overhead)
"""

"""Functions"""
def sinusoidal(t):
    return np.sin(2. * pi * np.array(t) - 0.5 * pi) + 1.

def ESWP():
    """SUF"""
    r = XylemFluxPython("../../grids/RootSystem8.rsml")
    p_s = np.linspace(-500, -200, 3001)  # 3 meter down, from -200 to -500, resolution in mm
    init_conductivities(r, age_dependent)
    r.rs.sort()  # ensures segment is located at index s.y-1
    nodes = r.get_nodes()
    rs_age = np.max(r.get_ages())
    soil_index = lambda x, y, z : int(-10 * z)  # maps to p_s (hydrostatic equilibirum)
    r.rs.setSoilGrid(soil_index)
    """ numerical solution of transpiration -10000 cm3/day"""
    rx = r.solve_neumann(sim_time, -10000, p_s, True)  # True: matric potential given per cell (not per segment)
    fluxes = r.segFluxes(sim_time, rx, p_s, False, True)  # cm3/day (double simTime,  rx,  sx,  approx, cells
    print("Transpiration", r.collar_flux(sim_time, rx, p_s), np.sum(fluxes), "cm3/day")
    suf = np.array(fluxes) / -10000.  # [1]
    print("Sum of suf", np.sum(suf), "from", np.min(suf), "to", np.max(suf))


    """Equivalent soil water potential"""
    min_b = [-4., -4., -15.]
    max_b = [4., 4., 0.]
    cell_number = [8, 8, 15]  
    periodic=False;

    s = RichardsWrapper(RichardsSP())
    s.initialize()
    s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
    loam = [0.08, 0.43, 0.04, 1.6, 50]  # we do not plan to calculate a thing, but we need parameters for initialisation
    s.setVGParameters([loam])
    s.setTopBC("noFlux")
    s.setBotBC("noFlux")
    s.initializeProblem()
    # Open .vtu
    pd = vp.read_vtu("DuMux_1cm-00000.vtu")
    print(pd.GetBounds())  # xmin, xmax, ymin, ymax, zmin, zmax
    print("Number of cells", vt.np_cells(pd).shape[0])

    data, _ = vt.np_data(pd,0, True)  # grid, data_index, cell data
    data = s.to_head(data)
    print("Data range from {:g} to {:g}".format(np.min(data), np.max(data)))
    s.setInitialCondition(data)  # put data to the grid


    """ Coupling (map indices) """
    picker = lambda x, y, z : s.pick([x, y, z])
    r.rs.setSoilGrid(picker)  # maps segments

    """ 3. EQUIVALENT SOIL WATER POTENTIAL """
    eswp = 0.
    ana = pb.SegmentAnalyser(r.rs)
    n = len(ana.segments)
    seg2cell_ = r.rs.seg2cell
    for i in range(0, n):
        eswp += suf[i] * data[seg2cell_[i]]
    print("\nEquivalent soil water potential", eswp)
    return eswp


""" Parameters """
min_b = [-4., -4., -15.]
max_b = [4., 4., 0.]
cell_number = [8, 8, 15]  
periodic = False

name = "DuMux_1cm"
sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil = clay

initial = -1659.8*2 + 7.5  # -659.8

trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

sim_time = 1  # [day] for task b
age_dependent = False  # conductivities
dt = 120. / (24 * 3600)  # [days] Time step must be very small

""" Initialize macroscopic soil model """
sp = vg.Parameters(soil)  # for debugging
cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
s.setTopBC("atmospheric", 0.5, [[0., 0.5,1.,1.e10], [0.1,0.1,1., 1.]])
s.setBotBC("noFlux")
s.setVGParameters([soil])
s.setParameter("Newton.EnableChop", "True")
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.initializeProblem()
s.setCriticalPressure(wilting_point)

""" Initialize xylem model (a) or (b)"""
r = XylemFluxPython("../../grids/RootSystem8.rsml")
r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)  # cutting
init_conductivities(r, age_dependent)
r.rs.sort()  # ensures segment is located at index s.y-1
r.test()  # sanity checks
nodes = r.get_nodes()
rs_age = np.max(r.get_ages())

""" Coupling (map indices) """
picker = lambda x, y, z : s.pick([x, y, z])
r.rs.setSoilGrid(picker)  # maps segments
cci = picker(nodes[0, 0], nodes[0, 1], nodes[0, 2])  # collar cell index

""" Numerical solution (a) """
start_time = timeit.default_timer()
x_, y_, w_, cpx, cps = [], [], [], [], []
sx = s.getSolutionHead()  # inital condition, solverbase.py

N = round(sim_time / dt)
t = 0.

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

    t += dt

s.writeDumuxVTK(name)
print("Root collar potential: ",rx[0]," eswp: ",ESWP())


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


