""" 
couple_c12 without mpi 
"""
import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import plantbox as pb
from functional.PlantHydraulicModel import HydraulicModel_Meunier
from functional.PlantHydraulicParameters import PlantHydraulicParameters
import functional.van_genuchten as vg
from functional.root_conductivities import *
import rsml.rsml_reader as rsml
import visualisation.vtk_plot as vp

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
from rhizo_models import plot_transpiration

import numpy as np
import matplotlib.pyplot as plt
import timeit

""" Parameters """
min_b = [-4., -4., -15.]
max_b = [4., 4., 0.]
cell_number = [8, 8, 15]
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

sim_time = 7.1  # 7.1  # [day] for task b
age_dependent = False  # conductivities
dt = 120. / (24 * 3600)  # [days] Time step must be very small
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
s.setParameter("Soil.SourceSlope", "100")  # turns regularisation of the source term on, will change the shape of actual transpiration...
s.initializeProblem()
s.setCriticalPressure(wilting_point)

""" Initialize root hydraulic model (a) or (b)"""
sinusoidal = HydraulicModel_Meunier.sinusoidal  # rename
params = PlantHydraulicParameters()
init_conductivities(params, age_dependent)

r = HydraulicModel_Meunier(fname, params, cached = True)  # or HydraulicModel_Doussan, HydraulicModel_Meunier

r.ms.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), cut = False)  

""" Numerical solution (a) """
start_time = timeit.default_timer()

x_, y_, z_ = [], [], []
N = round(sim_time / dt)
t = 0.
rs_age = np.max(r.get_ages())

sx = s.getSolutionHead()  
rx = r.solve(rs_age, -trans * sinusoidal(t), sx, cells = True)

for i in range(0, N):

    rx = r.solve_again(rs_age + t, -trans * sinusoidal(t), sx, cells = True)
    
    fluxes = r.soil_fluxes(rs_age + t, rx, sx)
    s.setSource(fluxes.copy())
    s.solve(dt)
    
    old_sx = sx.copy()
    sx = s.getSolutionHead()

    if i % skip == 0:
        x_.append(t)
        sum_flux = 0.
        for f in fluxes.values():
            sum_flux += f
        act_trans = r.get_transpiration(rs_age + t, rx, old_sx, cells = True)
        print("Summed fluxes ", sum_flux, "= collar flux", act_trans, "= prescribed", -trans * sinusoidal(t))
        y_.append(sum_flux)  # cm3/day (root system uptake)
        z_.append(act_trans)  # cm3/day
        n = round(float(i) / float(N) * 100.)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days {:g}\n"
              .format(np.min(sx), np.max(sx), np.min(rx), np.max(rx), s.simTime, rx[0]))

    t += dt

""" Plot """
print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")
plot_transpiration(x_, y_, z_, lambda t: trans * sinusoidal(t))


