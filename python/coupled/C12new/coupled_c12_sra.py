""" 
Benchmark M1.2 static root system in soil, root hydrualics with Meunier with rhizosphere (using steady rate approach)
"""

import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import plantbox as pb

from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.PlantHydraulicModel import HydraulicModel_Doussan
from functional.PlantHydraulicModel import HydraulicModel_Meunier
from functional.Perirhizal import PerirhizalPython  # Steady rate helper

import functional.van_genuchten as vg
from functional.root_conductivities import *
import rsml.rsml_reader as rsml
import visualisation.vtk_plot as vp

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
from rhizo_models import plot_transpiration

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import timeit

""" Parameters """
min_b = [-4., -4., -15.]
max_b = [4., 4., 0.]
cell_number = [8, 8, 15]  #  [8, 8, 15]  # [16, 16, 30]  # [32, 32, 60]  # [8, 8, 15]
periodic = False
fname = "../../../grids/RootSystem8.rsml"
age_dependent = False

name = "c12_sra"
loam = [0.08, 0.43, 0.04, 1.6, 50]
soil = loam

initial = -659.8 + 7.5  # -659.8

trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

sim_time = 7.1  # [day] for task b
dt = 360. / (24 * 3600)  # [days] Time step must be very small
skip = 1  # output

""" Initialize macroscopic soil model """
sp = vg.Parameters(soil)
# vg.create_mfp_lookup(sp, -1.e5, 1000)
s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([soil])
s.setParameter("Newton.EnableChop", "True")
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Soil.SourceSlope", "1000")  # turns regularisation of the source term on
s.initializeProblem()
s.setCriticalPressure(wilting_point)

""" Initialize root hydraulic model (a) or (b)"""
sinusoidal = HydraulicModel_Doussan.sinusoidal  # rename
params = PlantHydraulicParameters()
init_conductivities(params, age_dependent)

r = HydraulicModel_Meunier(fname, params, cached = True)  # or HydraulicModel_Doussan, HydraulicModel_Meunier

r.ms.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)  # cutting

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])
r.ms.setSoilGrid(picker)  # maps segment
outer_r = r.ms.segOuterRadii()
inner_r = r.ms.radii

""" sanity checks """
# r.plot_conductivities()
r.test()  # sanity checks
rs_age = np.max(r.get_ages())
print("Root system age ", rs_age)
seg_length = r.ms.segLength()
print("outer radii", np.min(outer_r) , np.max(outer_r))
print("inner radii", np.min(inner_r) , np.max(inner_r))
print()

peri = PerirhizalPython()
peri.open_lookup("../table_loam")

""" Numerical solution (a) """
start_time = timeit.default_timer()
x_, y_, z_ = [], [], []

mapping = r.ms.getSegmentMapper()
sx = s.getSolutionHead_()  # richards.py
hsb = np.array([sx[j] for j in mapping])
rsx = hsb.copy()  # initial values for fix point iteration
rho_ = np.divide(outer_r, np.array(inner_r))

N = round(sim_time / dt)
t = 0.

for i in range(0, N):

    rx = r.solve(rs_age + t, -trans * sinusoidal(t), rsx, cells = False)
    rx_old = rx.copy()

    kr_ = r.params.getKr(rs_age + t)
    inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up; here const

    err = 1.e6
    c = 0
    while err > 1. and c < 100:

        """ interpolation """
        rsx = peri.soil_root_interface_potentials(rx[1:], hsb, inner_kr_, rho_, None)

        """ xylem matric potential """
        rx = r.solve_again(rs_age + t, -trans * sinusoidal(t), rsx, cells = False)
        err = np.linalg.norm(rx - rx_old)

        rx_old = rx.copy()
        c += 1

    fluxes = r.radial_fluxes(rs_age + t, rx, rsx)
    collar_flux = r.get_transpiration(rs_age + t, rx.copy(), rsx.copy())
    err = np.linalg.norm(np.sum(fluxes) - collar_flux)
    if err > 1.e-6:
        print("error: summed root surface fluxes and root collar flux differ" , err, r.neumann_ind, collar_flux, np.sum(fluxes))
    err2 = np.linalg.norm(-trans * sinusoidal(t) - collar_flux)
    if r.last == "neumann":
        if err2 > 1.e-6:
            print("error: potential transpiration differs root collar flux in Neumann case" , err2)

    water = s.getWaterVolume()
    soil_fluxes = r.sumSegFluxes(fluxes)
    s.setSource(soil_fluxes.copy())  # richards.py
    s.solve(dt)
    soil_water = (s.getWaterVolume() - water) / dt

    sx = s.getSolutionHead_()  # richards.py
    hsb = np.array([sx[j] for j in mapping])  # soil bulk matric potential per segment

    if i % skip == 0:

        x_.append(t)
        sum_flux = 0.
        for f in soil_fluxes.values():
            sum_flux += f
        cf_ = r.get_transpiration(rs_age + t, rx, rsx)
        print("Summed fluxes ", sum_flux, "= collar flux", cf_, "= prescribed", -trans * sinusoidal(t))
        y_.append(sum_flux)  # cm4/day
        z_.append(soil_water)  # cm3
        n = round(float(i) / float(N) * 100.)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil sx [{:g}, {:g}], interface [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days\n"
              .format(np.min(sx), np.max(sx), np.min(rsx), np.max(rsx), np.min(rx), np.max(rx), s.simTime))

    t += dt

s.writeDumuxVTK(name)

""" Plot """
print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

# vp.plot_roots_and_soil(r.ms, "pressure head", rx, s, periodic, min_b, max_b, cell_number, name)  # VTK vizualisation
# ana = pb.SegmentAnalyser(r.rs)
# ana.addData("pressure", rx)
# vp.plot_roots(ana, "pressure")

plot_transpiration(x_, y_, z_, lambda t: trans * HydraulicModel_Doussan.sinusoidal(t))
np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter = ';')

