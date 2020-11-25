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
from root_conductivities import *

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
periodic = False

name = "dumux_c12_1cm_cyl_dumux"
sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil_ = loam
soil = vg.Parameters(soil_)
initial = -659.8 + (max_b[2]-min_b[2])/2 

kx = 4.32e-2
kr = 1.73e-4
trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

NC = 10  # dof of the cylindrical problem
logbase = 1.5

sim_time = 7  #  [day]
age_dependent = False  # conductivities
predefined_growth = False  # growth by setting radial conductivities
rs_age = 8 * (not predefined_growth) + 1 * predefined_growth  # rs_age = 0 in case of growth, else 8 days

split_type = 1  # type 0 == volume, type 1 == surface, type 2 == length

dt = 60/(24*3600)   # time step
NT = int(np.ceil(sim_time / dt))
skip = 1  # for output and results, skip iteration

domain_volume = np.prod(np.array(max_b) - np.array(min_b))

""" Initialize macroscopic soil model """
cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()
s.createGrid(min_b, max_b, cell_number)  # [cm]
s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([soil_])
#s.setParameter("Newton.EnableChop", "True")
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.initializeProblem()
s.setCriticalPressure(wilting_point)  # new source term regularisation
#s.setRegularisation(1.e-4, 1.e-4)
s.ddt = 1.e-5  # [day] initial Dumux time step

""" Initialize xylem model (a)"""
r = XylemFluxPython("../../grids/RootSystem8.rsml")
r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)
#r.setKr([kr])  # [cm^3/day] # todo check order of age (segs are sorted)
#r.setKx([kx])  # [1/day]
init_conductivities(r, age_dependent)   #age_dependent is a boolean, root conductivies are given in the file /root_conductivities.py
picker = lambda x, y, z : s.pick([x, y, z])
r.rs.setSoilGrid(picker)  # maps segments
r.rs.sort()  # <- ensures segment is located at index s.y-1

nodes = r.get_nodes()
cci = picker(nodes[0, 0], nodes[0, 1], nodes[0, 2])  # collar cell index

segs = r.get_segments()
seg_ages = r.get_ages(rs_age)
seg_types = r.rs.types
seg_length = r.segLength()
inner_radii = np.array(r.rs.radii)
outer_radii = r.segOuterRadii(split_type)

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
        cyl.setVGParameters([soil_])
        cyl.setOuterBC("fluxCyl", 0.)  # [cm/day]
        cyl.setInnerBC("fluxCyl", 0.)  # [cm/day]
        #cyl.setParameter("Newton.EnableChop", "True")
        cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
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

min_rx, min_rsx, collar_sx, collar_flux, out_times = [], [], [], [], []  # cm
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
    # print("Conductivities", np.min(soil_k), kr)
    rx = r.solve(rs_age+t, -trans * sinusoidal(t) , csx, rsx, False, wilting_point, soil_k)  # [cm]   * sinusoidal(t)

    if i % skip == 0:
        out_times.append(t)
        collar_flux.append(r.collar_flux(rs_age + t, rx, rsx, soil_k, cells = False))
        min_rsx.append(np.min(np.array(rsx)))
        collar_sx.append(csx)
        min_rx.append(np.min(np.array(rx)))
        # print("Conductivities", np.min(soil_k), kr_f(0., 0))
        print("Minimum of cylindrical model {:g} cm, soil cell {:g} cm, root xylem {:g} cm".format(min_rsx[-1], np.min(s.getSolutionHead()), min_rx[-1]))

    """
    Local soil model
    """
    proposed_inner_fluxes = r.segFluxes(rs_age + t, rx, rsx, approx = False, cells = False, soil_k = soil_k)  # [cm3/day]
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
            print("in",proposed_inner_fluxes[j] / (2 * np.pi * inner_radii[j] * l),"out", proposed_outer_fluxes[j] / (2 * np.pi * outer_radii[j] * l) )
            print("inner", inner_radii[j], "outer", outer_radii[j], "l", l )
            plt.show()
            raise
            #plt.pause(3)
            #plt.close()
        realized_inner_fluxes[j] = -float(cyl.getInnerFlux()) * (2 * np.pi * inner_radii[j] * l) / inner_radii[j]  # [cm/day] -> [cm3/day], ('/inner_radii' comes from cylindrical implementation)

    # print ("Local models solved in ", timeit.default_timer() - local_models_time, " s")

    """
    Macroscopic soil model
    """
    soil_water = np.multiply(np.array(s.getWaterContent()), cell_volumes)  # water per cell [cm3]
    soil_fluxes = r.sumSoilFluxes(realized_inner_fluxes)  # [cm3/day]  per soil cell

    s.setSource(soil_fluxes.copy())  # [cm3/day], in richards.py
    s.solve(dt)  # in solverbase.py

    new_soil_water = np.multiply(np.array(s.getWaterContent()), cell_volumes)
    net_flux = new_soil_water - soil_water  # change in water per cell [cm3]
    for k, root_flux in soil_fluxes.items():
        net_flux[k] -= root_flux * dt
    soil_water = new_soil_water

    if i % skip == 0:
        min_soil_fluxes, max_soil_fluxes, summed_soil_fluxes = 1.e9, -1.e9, 0.
        for k, v in soil_fluxes.items():
            summed_soil_fluxes += v
            if max_soil_fluxes < v:
                max_soil_fluxes = v
            if min_soil_fluxes > v:
                min_soil_fluxes = v
        print("Fluxes: realized", summed_soil_fluxes, "proposed", np.sum(proposed_inner_fluxes), 
              "predescribed", -trans * sinusoidal(t), "collar flux", collar_flux[-1])    
         # print("      : min {:g}, max {:g}".format(min_soil_fluxes, max_soil_fluxes))

    """ 
    Water (for output)
    """
    if i % skip == 0:
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

        #print("Water volume cylindric", cyl_water, "soil", soil_water[cci], cyl_water / soil_water[cci], cci)
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
ax1.plot(out_times, np.array(water_collar_cell), label = "water cell")
ax1.plot(out_times, np.array(water_cyl), label = "water cylindric")
ax1.legend()
ax1.set_xlabel("Time (days)")
ax1.set_ylabel("(cm3)")

ax2.set_title("Pressure")
ax2.plot(out_times, np.array(collar_sx), label = "soil at root collar")
ax2.plot(out_times, np.array(min_rx), label = "root collar")
ax2.plot(out_times, np.array(min_rsx), label = "1d model at root surface")
ax2.legend()
ax2.set_xlabel("Time (days)")
ax2.set_ylabel("Matric potential (cm)")
# plt.ylim(-15000, 0)

ax3.set_title("Water uptake")
ax3.plot(out_times, -np.array(water_uptake))
ax3.set_xlabel("Time (days)")
ax3.set_ylabel("Uptake (cm/day)")

ax4.set_title("Water in domain")
ax4.plot(out_times, np.array(water_domain))
ax4.set_xlabel("Time (days)")
ax4.set_ylabel("cm3")
plt.show()

fig, ax1 = plt.subplots()
ax1.plot(out_times, trans * sinusoidal(out_times), 'k', label = "potential")  # potential transpiration * sinusoidal(x_)
ax1.plot(out_times, -np.array(water_uptake), 'g', label = "actual")  # actual transpiration (neumann)
ax1.plot(out_times, -np.array(collar_flux), 'r:', label = "collar flux")  # actual transpiration (neumann)
ax2 = ax1.twinx()
ax2.plot(out_times, np.array(min_rx), label = "root collar")
ax1.set_xlabel("Time [d]")
ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
ax2.set_ylabel("Min collar pressure $[cm]$")
fig.legend()

np.savetxt("results/" + name, np.vstack((out_times, -np.array(collar_flux), -np.array(water_uptake))), delimiter = ';')

plt.show()
