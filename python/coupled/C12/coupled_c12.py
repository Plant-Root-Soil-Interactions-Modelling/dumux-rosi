"""
Benchmark M1.2 static root system in soil (root hydraulics with Doussan or Meunier using the classic sink)

Fro age dependent conductivities:
*) are very unstable, work better for larger spatial resolutions only (e.g. [4, 4, 15]), small temporal dt = 60 s
*) cache = False for Meunier, otherwise wrong

run times (for 7.1 days)
Meunier               195s
Meunier cached        200s (slower than direct...)
Doussan               213
Doussan cached        186

"""

import timeit

import matplotlib.pyplot as plt
import numpy as np
import plantbox as pb
import plantbox.functional.van_genuchten as vg
import plantbox.rsml.rsml_reader as rsml
import plantbox.visualisation.vtk_plot as vp
from mpi4py import MPI
from plantbox.functional.PlantHydraulicModel import (
    HydraulicModel_Doussan,
    HydraulicModel_Meunier,
)
from plantbox.functional.PlantHydraulicParameters import PlantHydraulicParameters
from plantbox.functional.root_conductivities import *
from rosi.richards import RichardsWrapper  # Python part
from rosi.rosi_richards import RichardsSPSSORC as RichardsSP  # C++ part (Dumux binding)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def plot_transpiration(t, soil_uptake, root_uptake, potential_trans, title=""):
    """plots potential and actual transpiration over time

    depending on discretisation soil_uptake and root_uptake might differ

    @param t                  times [day]
    @param soil_uptake        actual transpiration [cm3/day] of soil
    @param root_uptake        actual transpiration [cm3/day] according to root model
    @param potential_trans    function in t stating the potential transpiration [cm3/day]
    """
    fig, ax1 = plt.subplots()
    ax1.plot(t, [potential_trans(t_) for t_ in t], "k", label="potential transpiration")  # potential transpiration
    ax1.plot(t, -np.array(soil_uptake), "g", label="soil uptake")  # actual transpiration  according to soil model
    ax1.plot(t, -np.array(root_uptake), "r:", label="root system uptake")  # actual transpiration according root model
    ax1.set_xlabel("Time [d]")
    ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
    ax2 = ax1.twinx()
    dt = np.diff(t)
    so = np.array(soil_uptake)
    cum_transpiration = np.cumsum(-np.multiply(so[:-1], dt))
    ax2.plot(t[1:], cum_transpiration, "c--")  # cumulative transpiration (neumann)
    ax2.set_ylabel("Cumulative soil uptake $[cm^3]$")
    print("Cumulative soil uptake", cum_transpiration[-1], "[cm^3]")
    fig.legend()
    plt.title(title)
    if title:
        plt.savefig(title + ".png")
    plt.show()


""" Parameters """
min_b = [-4.0, -4.0, -15.0]
max_b = [4.0, 4.0, 0.0]
cell_number = [4, 4, 15]  # [8, 8, 15]
periodic = False
path = ""
fname = "../../../grids/RootSystem8.rsml"

name = "c12"
sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil = loam

initial = -659.8 + 7.5  # -659.8

trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

sim_time = 3.1  # 7.1  # [day] for task b
age_dependent = True  # conductivities
dt = 60.0 / (24 * 3600)  # [days] Time step must be very small
skip = 10

""" Initialize macroscopic soil model """
sp = vg.Parameters(soil)  # for debugging
s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([soil])
s.setParameter("Newton.EnableChop", "True")
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Soil.SourceSlope", "1000")  # turns regularisation of the source term on, will change the shape of actual transpiration...
s.initializeProblem()
s.setCriticalPressure(wilting_point)

""" Initialize root hydraulic model (a) or (b)"""
sinusoidal = HydraulicModel_Meunier.sinusoidal  # rename
params = PlantHydraulicParameters()
init_conductivities(params, age_dependent)

r = HydraulicModel_Doussan(fname, params, cached=True)  # or HydraulicModel_Doussan, HydraulicModel_Meunier

r.ms.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]), pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), cut=False)  # cutting

""" sanity checks """
# r.params.plot_conductivities()
r.test()  # sanity checks
rs_age = np.max(r.get_ages())
# print("press any key"); input()

""" Numerical solution (a) """
start_time = timeit.default_timer()

x_, y_, z_ = [], [], []
sink1d = []
sx = s.getSolutionHead()  # inital condition, solverbase.py
N = round(sim_time / dt)
t = 0.0


for i in range(0, N):

    if rank == 0:  # Root part is not parallel
        if i == 0:
            rx = r.solve(rs_age, -trans * sinusoidal(t), sx, cells=True)
        rx = r.solve_again(rs_age + t, -trans * sinusoidal(t), sx, cells=True)
        fluxes = r.soil_fluxes(rs_age + t, rx, sx)

    else:
        fluxes = None

    fluxes = comm.bcast(fluxes, root=0)  # Soil part runs parallel

    water = s.getWaterVolume()

    s.setSource(fluxes.copy())
    s.solve(dt)

    if rank == 0:
        old_sx = sx.copy()

    sx = s.getSolutionHead()

    soil_water = (s.getWaterVolume() - water) / dt  # since no-flux bc everywhere, change in water volume should equal transpirational flux

    if rank == 0 and i % skip == 0:
        x_.append(t)
        sum_flux = 0.0
        for f in fluxes.values():
            sum_flux += f
        act_trans = r.get_transpiration(rs_age + t, rx, old_sx, cells=True)
        print("Summed fluxes ", sum_flux, "= collar flux", act_trans, "= prescribed", -trans * sinusoidal(t))
        # y_.append(soil_water)  # cm3/day (soil uptake)
        y_.append(sum_flux)  # cm3/day (root system uptake)
        z_.append(act_trans)  # cm3/day
        n = round(float(i) / float(N) * 100.0)
        print("[" + "".join(["*"]) * n + "".join([" "]) * (100 - n) + "], soil [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days {:g}\n".format(np.min(sx), np.max(sx), np.min(rx), np.max(rx), s.simTime, rx[0]))

        # """ Additional sink plot """
        # if i % 60 == 0:  # every 6h
        #     ana = pb.SegmentAnalyser(r.ms)
        #     fluxes = r.radial_fluxes(rs_age + t, rx, old_sx, cells = True)
        #     ana.addData("fluxes", fluxes)  # cut off for vizualisation
        #     flux1d = ana.distribution("fluxes", max_b[2], min_b[2], 15, False)
        #     sink1d.append(np.array(flux1d))
        #     print("\nSink integral = ", np.sum(np.array(flux1d)), "\n")

    t += dt

s.writeDumuxVTK(name)

""" Plot """
if rank == 0:
    print("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    # vp.plot_roots_and_soil(r.ms, "pressure head", rx, s, periodic, min_b, max_b, cell_number, name)  # VTK vizualisation

    plot_transpiration(x_, y_, z_, lambda t: trans * sinusoidal(t))
    # TODO getWaterVolume() does not work with MPI

    # np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter = ';')
    # sink1d = np.array(sink1d)
    # np.save("sink1d", sink1d)
    # print(sink1d.shape)
    # print(sink1d[-1])
