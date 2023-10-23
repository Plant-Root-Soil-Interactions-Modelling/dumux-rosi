import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src/python_modules/")
sys.path.append("../../build-cmake/cpp/python_binding/")

import numpy as np
import matplotlib.pyplot as plt

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml_reader as rsml
import vtk_plot as vp
import van_genuchten as vg
from root_conductivities import *

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
from rhizo_models import plot_transpiration


def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.

""" 
Benchmark M1.2 static root system in soil (with the classic sink)
"""

""" Parameters """
min_b = [-4., -4., -15.]
max_b = [4., 4., 0.]
cell_number = [7, 7, 15]  #  [8, 8, 15]  # [1"../../../grids/RootSystem8.rsml"6, 16, 30]  # [32, 32, 60]  # [8, 8, 15]
periodic = False
path = ""
fname = "../../grids/RootSystem8.rsml"

# min_b = [-6.25, -1.5, -180.]  # cm
# max_b = [6.25, 1.5, 1]  # cm
# cell_number = [13, 3, 180]
# periodic = True
# path = "../../../CPlantBox/modelparameter/rootsystem/"
# fname = "Lupinus_albus_Leitner_2014"  # "spring_barley_CF12_107d.rsml"
# rs = pb.RootSystem()
# rs.readParameters(path + fname + ".xml")
# rs.setGeometry(pb.SDF_PlantBox(1e6, 1e6, -min_b[2]))  # to not let roots grow out of soil
# rs.initialize()
# rs.simulate(8, True)
# rs.write(fname + ".rsml"); fname = fname + ".rsml"

""" Soil """
name = "classical360"
sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil = loam

initial = -300 + 7.5  # -659.8

trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

sim_time = 7  # [day] for task b
age_dependent = False  # conductivities
dt = 360. / (24 * 3600)  # [days] Time step must be very small
skip = 1

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

""" Initialize xylem model """
r = XylemFluxPython(fname)
r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)  # cutting
init_conductivities(r, age_dependent)
rs_age = np.max(r.get_ages())  # initial root system age

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])
cci = picker(r.rs.nodes[0].x, r.rs.nodes[0].y, r.rs.nodes[0].z)  # collar cell index
r.rs.setSoilGrid(picker)  # maps segment

""" sanity checks """
# r.plot_conductivities()
r.test()  # sanity checks
# print("press any key"); input()

""" Numerical solution """
times_, soil_flux_, collar_flux_, sink1d_ = [], [], [], []
sx = s.getSolutionHead()  # inital condition, solverbase.py

N = round(sim_time / dt)
t = 0.

for i in range(0, N):

    rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., sx, True, wilting_point, [])  # xylem_flux.py, cells = True
    fluxes = r.soilFluxes(rs_age + t, rx, sx, False)  # class XylemFlux is defined in MappedOrganism.h, approx = True

    s.setSource(fluxes.copy())  # richards.py
    s.solve(dt)
    sx = s.getSolutionHead()  # richards.py

    if i % skip == 0:
        times_.append(t)
        sum_flux = 0.
        for f in fluxes.values():
            sum_flux += f
        soil_flux_.append(sum_flux)  # cm3/day
        collar_flux_.append(float(r.collar_flux(rs_age + t, rx, sx)))  # cm3/day
        print("Summed fluxes ", sum_flux, "= collar flux", collar_flux_[-1], "= prescribed", -trans * sinusoidal(t))
        min_sx, max_sx = np.min(sx), np.max(sx)
        min_rx, max_rx = np.min(rx), np.max(rx)
        n = round(float(i) / float(N) * 100.)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days {:g}\n"
              .format(min_sx, max_sx, min_rx, max_rx, s.simTime, rx[0]))

        """ Additional sink plot """
        if i % 60 == 0:  # every 6h
            ana = pb.SegmentAnalyser(r.rs)
            fluxes = r.segFluxes(rs_age + t, rx, sx, False, cells = True)  # cm3/day
            ana.addData("fluxes", fluxes)  # cut off for vizualisation
            flux1d = ana.distribution("fluxes", max_b[2], min_b[2], 15, True)
            sink1d_.append(np.array(flux1d))

    t += dt

""" Save results and plot """
s.writeDumuxVTK(name)
vp.plot_roots_and_soil(r.rs, "pressure head", rx, s, periodic, min_b, max_b, cell_number, name)  # VTK vizualisation
plot_transpiration(times_, soil_flux_, collar_flux_, lambda t: trans * sinusoidal(t))
np.savetxt(name, np.vstack((times_, -np.array(soil_flux_))), delimiter = ';')
np.save("sink1d", np.array(sink1d_))

