import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src/python_modules/")
sys.path.append("../../build-cmake/cpp/python_binding/")

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml_reader as rsml
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import vtk_plot as vp
import vtk_tools as vt
import van_genuchten as vg

from math import *
import numpy as np
import matplotlib.pyplot as plt
import timeit


def sinusoidal(t):
    return np.sin(2. * pi * np.array(t) - 0.5 * pi) + 1.

""" 
Benchmark M1.2 static root system in soil, coupled to cylindrical richards (Dumux solver)

also works parallel with mpiexec (only slightly faster, due to overhead)
"""

""" Parameters """
min_b = [-4., -4., -15.]
max_b = [4., 4., 0.]
cell_number = [8, 8, 15]  # [8, 8, 15]  # [16, 16, 30]  # [32, 32, 60]  # [8, 8, 15]
loam = [0.08, 0.43, 0.04, 1.6, 50]
soil = vg.Parameters(loam)
initial = -659.8 + 7.5  # -659.8

kx = 4.32e-2
kr = 1.73e-4
trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

NC = 10  # dof of the cylindrical problem
logbase = 1.5

sim_time = 0.25  # [day]

split_type = 0  # type 0 == volume, type 1 == surface, type 2 == length

NT = round(10 * sim_time * 24 * 3600 / 1200)
domain_volume = np.prod(np.array(max_b) - np.array(min_b))

name = "dumux_c12_025cm"

""" Initialize macroscopic soil model """
cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()
s.createGrid(min_b, max_b, cell_number)  # [cm]
s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([loam])
s.initializeProblem()
s.setCriticalPressure(wilting_point)  # new source term regularisation
s.setRegularisation(1.e-4, 1.e-4)
s.ddt = 1.e-5  # [day] initial Dumux time step

""" Initialize xylem model (a)"""
r = XylemFluxPython("../../grids/RootSystem8.rsml")
r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)
r.setKr([kr])  # [cm^3/day] # todo check order of age (segs are sorted)
r.setKx([kx])  # [1/day]
picker = lambda x, y, z : s.pick([x, y, z])
r.rs.setSoilGrid(picker)  # maps segments
r.rs.sort()  # <- ensures segment is located at index s.y-1

nodes = r.get_nodes()
segs = r.get_segments()
cci = picker(nodes[0, 0], nodes[0, 1], nodes[0, 2])  # collar cell index

inner_radii = np.array(r.rs.radii)
outer_radii = r.segOuterRadii(split_type)
seg_length = r.segLength()

r.test()  # sanity checks
print("Initial pressure head", s.getSolutionHeadAt(cci), s.getSolutionHeadAt(picker(0., 0., min_b[2])))

""" Initialize local soil models (around each root segment) """
grids = []
cyls = []
ns = len(seg_length)  # number of segments
for i in range(0, ns):
    a_in = inner_radii[i]
    a_out = outer_radii[i]
    if a_in < a_out:
        cpp_base = RichardsCylFoam()
        cyl = RichardsWrapper(cpp_base)
        cyl.initialize()
        points = np.logspace(np.log(a_in) / np.log(logbase), np.log(a_out) / np.log(logbase), NC, base = logbase)
        grids.append(points)  # to remember
        cyl.createGrid1d(points)
        cyl.setHomogeneousIC(initial)  # cm pressure head
        cyl.setVGParameters([loam])
        cyl.setOuterBC("fluxCyl", 0.)  # [cm/day]
        cyl.setInnerBC("fluxCyl", 0.)  # [cm/day]
        cyl.initializeProblem()
        cyl.setCriticalPressure(wilting_point)  # cm pressure head
        cyls.append(cyl)
    else:
        cyls.append([])
        print("Segment not in domain", i, "[", a_in, a_out, "]")  # this happens if elements are not within the domain
        input()

""" Simulation """
print("Starting simulation")
start_time = timeit.default_timer()

rsx = np.zeros((ns,))  # xylem pressure at the root soil interface
dt = sim_time / NT

min_rx, min_rsx, collar_sx, collar_flux = [], [], [], []  # cm
water_uptake, water_collar_cell, water_cyl, water_domain = [], [], [], []  # cm3

rsx = np.zeros((ns,))  # matric potential at the root soil interface [cm]
cell_volumes = s.getCellVolumes()
inital_soil_water = np.sum(np.multiply(np.array(s.getWaterContent()), cell_volumes))

net_flux = np.zeros(cell_volumes.shape)
realized_inner_fluxes = np.zeros((len(cyls),))

for i in range(0, NT):

    t = i * dt

    """ 
    Xylem model 
    """
    csx = s.getSolutionHeadAt(cci)  # [cm]
    for j, rc in enumerate(cyls):  # for each segment
        rsx[j] = rc.getInnerHead()  # [cm]
#         if rsx[j] > s.getSolutionHeadAt(r.rs.seg2cell[j]) + 1:
#             print("strange segment", j, "in cell", r.rs.seg2cell[j], "root soil interface", rsx[j], "vs macro soil", s.getSolutionHeadAt(r.rs.seg2cell[j]))

    soil_k = np.divide(vg.hydraulic_conductivity(rsx, soil), inner_radii)  # only valid for homogenous soil
    print("Conductivities", np.min(soil_k), kr)
    rx = r.solve(t, -trans * sinusoidal(t) , csx, rsx, False, wilting_point, soil_k)  # [cm]   * sinusoidal(t)
    collar_flux.append(r.collar_flux(t, rx, rsx, soil_k, False))

    min_rsx.append(np.min(np.array(rsx)))
    collar_sx.append(csx)
    min_rx.append(np.min(np.array(rx)))
    print("Minimum of cylindrical model {:g} cm, soil cell {:g} cm, root xylem {:g} cm".format(min_rsx[-1], np.min(s.getSolutionHead()), min_rx[-1]))

    """
    Local soil model
    """
    proposed_inner_fluxes = r.segFluxes(0., rx, rsx, approx = False)  # [cm3/day]
    proposed_outer_fluxes = r.splitSoilFluxes(net_flux / dt, split_type)
    local_models_time = timeit.default_timer()
    for j, cyl in enumerate(cyls):  # run cylindrical models
        l = seg_length[j]
        cyl.setInnerFluxCyl(proposed_inner_fluxes[j] / (2 * np.pi * inner_radii[j] * l))  # [cm3/day] -> [cm /day]
        cyl.setOuterFluxCyl(proposed_outer_fluxes[j] / (2 * np.pi * outer_radii[j] * l))  # [cm3/day] -> [cm /day]
        cyl.ddt = 1.e-5  # [day] initial time step
        try:
            cyl.solve(dt)
        except:
            x = cyl.getDofCoordinates()
            y = cyl.getSolutionHead()
            plt.plot(x, y)
            plt.xlabel("x (cm)")
            plt.ylabel("pressure (cm)")
            plt.show()
        realized_inner_fluxes[j] = -float(cyl.getInnerFlux()) * (2 * np.pi * inner_radii[j] * l) / inner_radii[j]  # [cm/day] -> [cm3/day], ('/inner_radii' comes from cylindrical implementation)

    print ("Local models solved in ", timeit.default_timer() - local_models_time, " s")

    """
    Macroscopic soil model
    """
    soil_water = np.multiply(np.array(s.getWaterContent()), cell_volumes)  # water per cell [cm3]
    soil_fluxes = r.sumSoilFluxes(realized_inner_fluxes)  # [cm3/day]
    s.setSource(soil_fluxes.copy())  # [cm3/day], richards.py

    summed_soil_fluxes = 0.
    for k, v in soil_fluxes.items():
        summed_soil_fluxes += v

    print("Fluxes: realized per segment", summed_soil_fluxes, np.sum(realized_inner_fluxes), "predescribed: ", collar_flux[-1], -trans * sinusoidal(t))

    s.solve(dt)

    new_soil_water = np.multiply(np.array(s.getWaterContent()), cell_volumes)
    net_flux = new_soil_water - soil_water  # change in water per cell [cm3]
    for k, root_flux in soil_fluxes.items():
        net_flux[k] -= root_flux * dt
    print("Summed net flux {:g}, min movement {:g}, max {:g} cm3".format(np.sum(net_flux), np.min(net_flux), np.max(net_flux)))  # summed fluxes should equal zero
    soil_water = new_soil_water

    """ 
    Water (for output)
    """
    water_domain.append(np.min(soil_water))  # from previous time step
    water_collar_cell.append(soil_water[cci])
    sum_flux = 0.
    for k, f in soil_fluxes.items():
        sum_flux += f
    water_uptake.append(sum_flux)  # cm3/day

    cyl_water = 0.
    for k in r.rs.cell2seg[cci]:
        cyl_water_content = cyls[k].getWaterContent()  # segment 0
        for j, wc in enumerate(cyl_water_content):
            r1 = grids[k][j]
            r2 = grids[k][j + 1]
            cyl_water += np.pi * (r2 * r2 - r1 * r1) * seg_length[k] * wc

    print("Water volume cylindric", cyl_water, "soil", soil_water[cci], cyl_water / soil_water[cci], cci)
    water_cyl.append(cyl_water)

    n = round(float(i) / float(NT) * 100.)
    print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], {:g} days".format(s.simTime))

print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

# VTK vizualisation
ana = pb.SegmentAnalyser(r.rs)
ana.addData("soil pressure", rsx)
ana.addData("xylem pressure", rx)
pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", "soil pressure", "xylem pressure"])
rootActor, rootCBar = vp.plot_roots(pd, "soil pressure" , "", False)

soil_grid = vp.uniform_grid(np.array(min_b), np.array(max_b), np.array(cell_number))
soil_water_content = vt.vtk_data(np.array(s.getWaterContent()))
soil_water_content.SetName("water content")
soil_grid.GetCellData().AddArray(soil_water_content)
soil_pressure = vt.vtk_data(np.array(s.getSolutionHead()))
soil_pressure.SetName("pressure head")  # in macroscopic soil
soil_grid.GetCellData().AddArray(soil_pressure)
meshActors, meshCBar = vp.plot_mesh_cuts(soil_grid, "pressure head", 5, "", False)

lut = meshActors[-1].GetMapper().GetLookupTable()  # same same
rootActor.GetMapper().SetLookupTable(lut)
meshActors.extend([rootActor])
vp.render_window(meshActors, name + " at {:g} days".format(sim_time) , meshCBar, soil_grid.GetBounds()).Start()

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

ax1.set_title("Water amount")
ax1.plot(np.linspace(0, sim_time, NT), np.array(water_collar_cell), label = "water cell")
ax1.plot(np.linspace(0, sim_time, NT), np.array(water_cyl), label = "water cylindric")
ax1.legend()
ax1.set_xlabel("Time (days)")
ax1.set_ylabel("(cm3)")

ax2.set_title("Pressure")
ax2.plot(np.linspace(0, sim_time, NT), np.array(collar_sx), label = "soil at root collar")
ax2.plot(np.linspace(0, sim_time, NT), np.array(min_rx), label = "root collar")
ax2.plot(np.linspace(0, sim_time, NT), np.array(min_rsx), label = "1d model at root surface")
ax2.legend()
ax2.set_xlabel("Time (days)")
ax2.set_ylabel("Matric potential (cm)")
# plt.ylim(-15000, 0)

ax3.set_title("Water uptake")
ax3.plot(np.linspace(0, sim_time, NT), -np.array(water_uptake))
ax3.set_xlabel("Time (days)")
ax3.set_ylabel("Uptake (cm/day)")

ax4.set_title("Water in domain")
ax4.plot(np.linspace(0, sim_time, NT), np.array(water_domain))
ax4.set_xlabel("Time (days)")
ax4.set_ylabel("cm3")
plt.show()

fig, ax1 = plt.subplots()
x_ = np.linspace(0, sim_time, NT)
ax1.plot(x_, trans * sinusoidal(x_), 'k', label = "potential")  # potential transpiration * sinusoidal(x_)
ax1.plot(x_, -np.array(water_uptake), 'g', label = "actual")  # actual transpiration (neumann)
ax1.plot(x_, -np.array(collar_flux), 'r:', label = "collar flux")  # actual transpiration (neumann)
ax2 = ax1.twinx()
ax2.plot(np.linspace(0, sim_time, NT), np.array(min_rx), label = "root collar")
ax1.set_xlabel("Time [d]")
ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
ax2.set_ylabel("Min collar pressure $[cm]$")
fig.legend()

np.savetxt(name, np.vstack((x_, -np.array(collar_flux), -np.array(water_uptake))), delimiter = ';')

plt.show()
