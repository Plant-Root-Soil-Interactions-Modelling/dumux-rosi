import sys; sys.path.append("../modules/"); sys.path.append("../modules/fv/"); sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src/python_modules/")
sys.path.append("../../build-cmake/cpp/python_binding/")

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb  # CPlantBox
import rsml_reader as rsml
from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part
from fv_grid import *
import fv_richards as rich  # local cylindrical models
import vtk_plot as vp
import vtk_tools as vt
import van_genuchten as vg
from root_conductivities import *

from math import *
import numpy as np
import matplotlib.pyplot as plt
import timeit
from multiprocessing import Pool


def sinusoidal(t):
    """ sinusoidal function (used for transpiration) """
    return np.sin(2. * pi * np.array(t) - 0.5 * pi) + 1.   # TODO: for growing root systems: adapt to either LAI or root system size

""" 
Benchmark M1.2 static root system in soil, coupled to cylindrical richards (Python solver)
"""

""" Parameters """
min_b = [-4., -4., -15.]
max_b = [4., 4., 0.]
cell_number = [8, 8, 15]  # [8, 8, 15]  # [16, 16, 30]  # [32, 32, 60]  # [8, 8, 15]
periodic = False

name = "dumux_c12_2cm"
sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil_ = sand
soil = vg.Parameters(soil_)
initial = -659.8 + (max_b[2]-min_b[2])/2   # -659.8 + 7.5 because -659.8 is the value at the top, but we need the average value in the domain

trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

sim_time = 7  # 0.65  # 0.25  # [day]
age_dependent = False  # conductivities
predefined_growth = False  # growth by setting radial conductivities
rs_age = 8 * (not predefined_growth) + 1 * predefined_growth  # rs_age = 0 in case of growth, else 8 days

NC = 10  # dof+1
logbase = 1.5  # according to Mai et al. (2019)
split_type = 0  # type 0 == volume, type 1 == surface, type 2 == length

dt = 120/(24*3600)   # time step
NT = int(np.ceil(sim_time / dt))
skip = 10  # for output and results, skip iteration
domain_volume = np.prod(np.array(max_b) - np.array(min_b))

""" Initialize macroscopic soil model """
cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
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

""" Initialize xylem model (a) or (b)"""
r = XylemFluxPython("../../grids/RootSystem8.rsml")
print("number of segments", len(r.get_segments()))
r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)  # True: root segments are cut according to the soil grid so that each segment is completely within one soil control element(works only for rectangular grids so far)
init_conductivities(r, age_dependent)   #age_dependent is a boolean, root conductivies are given in the file /root_conductivities.py
picker = lambda x, y, z : s.pick([x, y, z])   #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
r.rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions 
r.rs.sort()  # <- ensures segment is located at index s.y-1

nodes = r.get_nodes()
cci = picker(nodes[0, 0], nodes[0, 1], nodes[0, 2])  # collar cell index
segs = r.get_segments()
seg_ages = r.get_ages(rs_age)
seg_types = r.rs.subTypes
seg_length = r.segLength()
inner_radii = np.array(r.rs.radii)
outer_radii = r.segOuterRadii(split_type)

# # For debugging
# ana2 = pb.SegmentAnalyser(r.rs.nodes, r.rs.segments, r.rs.nodeCTs[1:], r.rs.radii)
# types = np.array(r.rs.types, dtype = np.float64)
# ana2.addData("subType", types)
# ana2.addData("age", r.get_ages())
# pd = vp.segs_to_polydata(ana2, 1., ["radius", "subType", "creationTime", "age"])
# vp.plot_roots(pd, "creationTime")

r.test()  # sanity checks
print("Initial root system age ", rs_age)
print("Initial pressure head", s.getSolutionHeadAt(cci), s.getSolutionHeadAt(picker(0., 0., min_b[2])))
# input()

""" Initialize local soil models (around each root segment) """
ns = len(seg_length)  # number of segments
cyls = [None] * ns
ndof = NC - 1    # rhizosphere grid (dofs live on the cells)


def initialize_cyl(i):
    """ initialization of local cylindrical model """
    a_in = inner_radii[i]
    a_out = outer_radii[i]
    if a_in < a_out:
        z = 0.5 * (nodes[segs[i, 0], 2] + nodes[segs[i, 1], 2])
        points = np.logspace(np.log(a_in) / np.log(logbase), np.log(a_out) / np.log(logbase), NC, base = logbase)
        # points = np.linspace(a_in, a_out, NC)
        grid = FVGrid1Dcyl(points)
        cyl = rich.FVRichards1D(grid, soil_)
        cyl.x0 = np.ones((ndof,)) * (initial - (max_b[2]-min_b[2])/2 - z)
        return cyl
    else:
        print("Segment", i, "[", a_in, a_out, "]")  # this happens if elements are not within the domain
        return []


def simulate_cyl(cyl):
    try:
        # cyl.solve_single(sim_time / NT, False)
        cyl.solve([sim_time / NT], sim_time / NT, False)
    except:
        x = cyl.grid.centers()
        y = cyl.x0
        plt.plot(x, y)
        plt.xlabel("x (cm)")
        plt.ylabel("pressure (cm)")
        plt.show()
        print(cyl.bc[(0, 0)])
        input()
    return cyl


pool = Pool()  # defaults to number of available CPU's
start_time = timeit.default_timer()
cyls = pool.map(initialize_cyl, range(ns))
print ("Initialized in", timeit.default_timer() - start_time, " s")

""" Simulation """
print("Starting simulation")
start_time = timeit.default_timer()

min_rx, min_rsx, collar_sx, collar_flux, out_times = [], [], [], [], []  # cm
water_uptake, water_collar_cell, water_cyl, water_domain = [], [], [], []  # cm3

rsx = np.zeros((ns,))  # matric potential at the root soil interface [cm]
cell_volumes = s.getCellVolumes()
inital_soil_water = np.sum(np.multiply(np.array(s.getWaterContent()), cell_volumes))

net_flux = np.zeros(cell_volumes.shape)
realized_inner_fluxes = np.zeros((len(cyls),))

for i in range(0, NT):

    wall_iteration = timeit.default_timer()
    t = i * dt

    """ 
    Xylem model 
    """
    csx = s.getSolutionHeadAt(cci)  # [cm]
    for j, rc in enumerate(cyls):  # for each segment
        rsx[j] = rc.getInnerHead()  # [cm]
#         if rsx[j] > s.getSolutionHeadAt(r.rs.seg2cell[j]) + 1:
#             print("strange segment", j, "in cell", r.rs.seg2cell[j], "root soil interface", rsx[j], "vs macro soil", s.getSolutionHeadAt(r.rs.seg2cell[j]))

    wall_root_model = timeit.default_timer()
    soil_k = np.divide(vg.hydraulic_conductivity(rsx, soil), inner_radii)  # only valid for homogenous soil
    rx = r.solve(rs_age + t, -trans * sinusoidal(t), csx, rsx, False, wilting_point, soil_k)  # [cm]   False means that rsx is given per root segment not per soil cell
    wall_root_model = timeit.default_timer() - wall_root_model

    if i % skip == 0:
        out_times.append(t)
        collar_flux.append(r.collar_flux(rs_age + t, rx, rsx, soil_k, False))
        min_rsx.append(np.min(np.array(rsx)))
        collar_sx.append(csx)
        min_rx.append(np.min(np.array(rx)))
        # print("Conductivities", np.min(soil_k), kr_f(0., 0))
        print("Minimum of cylindrical model {:g} cm, soil cell {:g} cm, root xylem {:g} cm".format(min_rsx[-1], np.min(s.getSolutionHead()), min_rx[-1]))

    """
    Local soil model
    """
    # proposed_inner_fluxes = r.segFluxes(0., rx, rsx, approx=False)  # [cm3/day]
    proposed_outer_fluxes = r.splitSoilFluxes(net_flux / dt, split_type)

    for j, cyl in enumerate(cyls):  # boundary condtions
        l = seg_length[j]
        age = seg_ages[j] + t  # age  = (maxCT -nodeCT[j]) + t = (maxCT+t) - nodeCT[j] = rs_sim_time - nodeCT[j]
        type_ = seg_types[j]
#         rx_approx = 0.5 * (rx[segs[j][0]] + rx[segs[j][1]])
#         cyl.bc[(0, 0)] = ("rootsystem", [rx_approx, kr_f(age, type_)])
        # print(r.kr_f(age, type_), r.kx_f(age, type_))
        cyl.bc[(0, 0)] = ("rootsystem_exact", [rx[segs[j][0]], rx[segs[j][1]], r.kr_f(age, type_), r.kx_f(age, type_), inner_radii[j], l ])
        dx_outer = cyl.grid.nodes[ndof] - cyl.grid.center(ndof - 1)
        q_outer = proposed_outer_fluxes[j] / (2 * np.pi * outer_radii[j] * l)
        cyl.bc[(ndof - 1, 1)] = ("flux_in_out", [q_outer , wilting_point, dx_outer])

    wall_rhizo_models = timeit.default_timer()

    cyls = pool.map(simulate_cyl, cyls)  # simulate
#     for cyl in cyls:
#         cyl.solve([sim_time / NT], sim_time / NT, False)

    wall_rhizo_models = timeit.default_timer() - wall_rhizo_models

    for j, cyl in enumerate(cyls):  # res
        realized_inner_fluxes[j] = cyl.getInnerFlux() * (2 * np.pi * inner_radii[j] * seg_length[j]) / dt # divide by dt is correct here! getInnerFlux only gives the source in cm3/cm2

    """
    Macroscopic soil model
    """
    soil_water = np.multiply(np.array(s.getWaterContent()), cell_volumes)  # water per cell [cm3]
    soil_fluxes = r.sumSoilFluxes(realized_inner_fluxes)  # [cm3/day]  per soil cell

    wall_soil_model = timeit.default_timer()
    s.setSource(soil_fluxes.copy())  # [cm3/day], in richards.py
    s.solve(dt)  # in solverbase.py
    wall_soil_model = timeit.default_timer() - wall_soil_model

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
        print("Fluxes: summed", summed_soil_fluxes, "predescribed", -trans * sinusoidal(t), 
              "collar flux", collar_flux[-1])
        # print("      : min {:g}, max {:g}".format(min_soil_fluxes, max_soil_fluxes))
        # print("Summed net flux {:g}, min movement {:g}, max {:g} cm3".format(np.sum(net_flux), np.min(net_flux), np.max(net_flux)))  # summed fluxes should equal zero

    """ 
    Water (for output only)
    """
    wall_iteration = timeit.default_timer() - wall_iteration
    if i % skip == 0:
        water_domain.append(np.min(soil_water))  # from previous time step
        water_collar_cell.append(soil_water[cci])
        water_uptake.append(summed_soil_fluxes)  # cm3/day
        cyl_water = 0.
        for k in r.rs.cell2seg[cci]:
            cyl_water_content = cyls[k].getWaterContent()  # segment 0
            for j, wc in enumerate(cyl_water_content):
                r1 = cyls[k].grid.nodes[j]
                r2 = cyls[k].grid.nodes[j + 1]
                cyl_water += np.pi * (r2 * r2 - r1 * r1) * seg_length[k] * wc
        # print("Water volume cylindric", cyl_water, "soil", soil_water[cci], cyl_water / soil_water[cci], cci)
        water_cyl.append(cyl_water)
        print("Iteration {:g} s, rhizo {:g} s, {:g}% root, {:g}% rhizo {:g}% soil".
              format(wall_iteration, wall_rhizo_models, wall_root_model / wall_iteration, wall_rhizo_models / wall_iteration, wall_soil_model / wall_iteration))
        n = round(float(i) / float(NT) * 100.)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], {:g} days".format(s.simTime))
        print()

print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

vp.plot_roots_and_soil(r.rs, "pressure head", rsx, s, periodic, min_b, max_b, cell_number, name)  # VTK vizualisation

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
x_ = out_times

ax1.set_title("Water amount")
ax1.plot(x_, np.array(water_collar_cell), label = "water cell")
ax1.plot(x_, np.array(water_cyl), label = "water cylindric")
ax1.legend()
ax1.set_xlabel("Time (days)")
ax1.set_ylabel("(cm3)")

ax2.set_title("Pressure")
ax2.plot(x_, np.array(collar_sx), label = "soil at root collar")
ax2.plot(x_, np.array(min_rx), label = "root collar")
ax2.plot(x_, np.array(min_rsx), label = "1d model at root surface")
ax2.legend()
ax2.set_xlabel("Time (days)")
ax2.set_ylabel("Matric potential (cm)")
# plt.ylim(-15000, 0)

ax3.set_title("Water uptake")
ax3.plot(x_, -np.array(water_uptake))
ax3.set_xlabel("Time (days)")
ax3.set_ylabel("Uptake (cm/day)")

ax4.set_title("Water in domain")
ax4.plot(x_, np.array(water_domain))
ax4.set_xlabel("Time (days)")
ax4.set_ylabel("cm3")
plt.show()

fig, ax1 = plt.subplots()
ax1.plot(x_, trans * sinusoidal(x_), 'k', label = "potential")  # potential transpiration * sinusoidal(x_)
ax1.plot(x_, -np.array(water_uptake), 'g', label = "actual")  # actual transpiration (neumann)
ax1.plot(x_, -np.array(collar_flux), 'r:', label = "collar flux")  # actual transpiration (neumann)
ax2 = ax1.twinx()
ax2.plot(x_, np.array(min_rx), label = "root collar")
ax1.set_xlabel("Time [d]")
ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
ax2.set_ylabel("Min collar pressure $[cm]$")
fig.legend()

np.savetxt("results/" + name, np.vstack((x_, -np.array(collar_flux), -np.array(water_uptake))), delimiter = ';')

plt.show()

# old_rs = XylemFluxPython("../grids/RootSystem_big.rsml")
# ana = pb.SegmentAnalyser(old_rs.rs)
# ana.filter("creationTime", 0., 8)
# ana.crop(pb.SDF_PlantBox(7.76, 7.76, 14.76))
# ana.pack()
# segCT = ana.data["creationTime"]  # per segment
# nodeCT = np.zeros((len(segCT) + 1))  # convert segCT to nodeCT
# for i, seg in enumerate(ana.segments):
#     nodeCT[seg.y] = segCT[i]
# subType = np.array(ana.data["subType"], dtype = np.int64)  # convert to int
# rs = pb.MappedSegments(ana.nodes, nodeCT, ana.segments, ana.data["radius"], subType)
# r = XylemFluxPython(rs)
