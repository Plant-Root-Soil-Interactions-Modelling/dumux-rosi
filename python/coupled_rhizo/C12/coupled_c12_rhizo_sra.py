import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src");

import plantbox as pb  # CPlantBox
from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
from functional.xylem_flux import *  # root system Python hybrid solver
from rhizo_models import *  # Helper class for cylindrical rhizosphere models

import visualisation.vtk_plot as vp
import functional.van_genuchten as vg
from functional.root_conductivities import *

import numpy as np
import timeit
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

""" 
Benchmark M1.2 static root system in soil

Coupled to cylindrical rhizosphere models using 1d richards equation (DUMUX solver)

Method description: 
We use a Dirichlet-Neumann boundary condition, and prescribe the pressure head at the soil-root interface, which we compute via iterations from the xylem pressure head and the bulk soil pressure head. Contrary to the first method, this method does not require a very small time step. The  result is the same as with the neumann-neumann boundary condition, but tiem steps can be chosen larger and it is more stable 
Before running this approach, a lookup table for the pressure heads must be created by running create_sri_lookup_tables.py

complicated MPI support (a non-mpi version of richards_cyl is needed, see script dumux3_nompi.sh)
"""


def open_sri_lookup(filename):
    """ opens the look-up table from a file, to quickly find soil root interface potential """
    sri_table = np.load(filename + ".npy")
    x = np.load(filename + "_.npy", allow_pickle = True)
    rx_ = x[0]
    sx_ = x[1]
    inner_ = x[2]
    outer_ = x[3]
    return RegularGridInterpolator((rx_, sx_, inner_, outer_), sri_table)  # method = "nearest" fill_value = None , bounds_error=False


def soil_root_interface_table(rx, sx, inner_kr_, rho_, f):
    """
    finds potential at the soil root interface
        
    rx             xylem matric potential [cm]
    sx             bulk soil matric potential [cm]
    inner_kr       root radius times hydraulic conductivity [cm/day] 
    rho            geometry factor [1]
    f              function to look up the potentials
    """
    try:
        rsx = f((rx, sx, inner_kr_ , rho_))
    except:
        print("rx", np.min(rx), np.max(rx))  # 0, -16000
        print("sx", np.min(sx), np.max(sx))  # 0, -16000
        print("inner_kr", np.min(inner_kr_), np.max(inner_kr_))  # 1.e-7 - 1.e-4
        print("rho", np.min(rho_), np.max(rho_))  # 1. - 200.
    return rsx

""" 
Parameters  
"""
name = "c12_rhizo_1cm_sra"  # scenario name, to save results

""" soil """
min_b = [-4., -4., -15.]  # cm
max_b = [4., 4., 0.]  # cm
cell_number = [7, 7, 15]  # [8, 8, 15]  # [16, 16, 30]  # [32, 32, 60]  # [8, 8, 15] # [1]
periodic = False
fname = "../../../grids/RootSystem8.rsml"

domain_volume = np.prod(np.array(max_b) - np.array(min_b))

sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil_ = loam
soil = vg.Parameters(soil_)
initial = -659.8 + (max_b[2] - min_b[2]) / 2  # -659.8 + 7.5 because -659.8 is the value at the top, but we need the average value in the domain

sri_table_lookup = open_sri_lookup('../../coupled/table_loam')
root_interface = soil_root_interface_table
""" root system """
trans = 6.4  # average per day [cm3 /day] (sinusoidal)
wilting_point = -15000  # [cm]
age_dependent = False  # conductivities
predefined_growth = False  # root growth by setting radial conductivities
rs_age = 8 * (not predefined_growth) + 1 * predefined_growth  # rs_age = 0 in case of growth, else 8 days

""" rhizosphere models """
mode = "dumux_dirichlet"  # or "dumux_exact"
NC = 10  # dof+1
logbase = 1.5  # according to Mai et al. (2019)
split_type = 0  # type 0 == volume, type 1 == surface, type 2 == length

""" simulation time """
sim_time = 3  # 7  # 0.65  # 0.25  # [day]
dt = 30 / (24 * 3600)  # time step [day], 120 schwankt stark
NT = int(np.ceil(sim_time / dt))  # number of iterations
skip = 1  # for output and results, skip iteration
max_iter = 1000  # maximum for fix point iteration

""" 
Initialize macroscopic soil model (Dumux binding)
"""
s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([soil_])
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
# s.setParameter("Soil.SourceSlope", "1000")  # turns regularisation of the source term on
s.initializeProblem()
s.setCriticalPressure(wilting_point)  # new source term regularisation
s.ddt = 1.e-5  # [day] initial Dumux time step
print()

""" 
Initialize xylem model 
"""
rs = RhizoMappedSegments(fname, wilting_point, NC, logbase, mode)
rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)
r = XylemFluxPython(rs)  # wrap the xylem model around the MappedSegments
init_conductivities_growth(r, age_dependent, 0.05)  # age_dependent is a boolean, root conductivies are given in the file src/python_modules/root_conductivities.py

picker = lambda x, y, z: s.pick([x, y, z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
rs.set_xylem_flux(r)

# For debugging
# r.plot_conductivities()
r.test()  # sanity checks (todo need improvements...)

""" 
Initialize local soil models (around each root segment) 
"""
start_time = timeit.default_timer()
x = s.getSolutionHead()[:]  # initial condition of soil [cm]
x = comm.bcast(x, root = 0)  # Soil part runs parallel
nodes = rs.nodes
segs = rs.segments
ns = len(segs)
seg2cell = rs.seg2cell
cell2seg = rs.cell2seg
mapping = rs.getSegmentMapper()

for i in range(0, len(segs)):
    if segs[i].x == 0:
        collar_ind = i  # segment index of root collar
        break

hsb = np.array([x[j] for j in mapping])  # soil bulk matric potential per segment
cell_centers = s.getCellCenters_()
cell_centers_z = np.array([cell_centers[j][2] for j in mapping])
seg_centers_z = rs.getSegmentZ()
ns = len(rs.segments)
dcyl = int(np.floor(ns / max_rank))
if rank + 1 == max_rank:
    rs.initialize(soil_, x, np.array(range(rank * dcyl, ns)))
    print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, rank * dcyl, ns, timeit.default_timer() - start_time))
else:
    rs.initialize(soil_, x, np.array(range(rank * dcyl, (rank + 1) * dcyl)))
    print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, rank * dcyl, (rank + 1) * dcyl, timeit.default_timer() - start_time))
# print("press any key"); input()

outer_r = rs.segOuterRadii()
inner_r = rs.radii
types = rs.subTypes
rho_ = np.divide(outer_r, np.array(inner_r))
rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200)

""" 
Simulation 

loop
1. xylem model
2. local soil models
3. macroscopic soil model 
"""
print("Starting simulation")
start_time = timeit.default_timer()

# for post processing
min_sx, min_rx, min_rsx, collar_sx, collar_flux = [], [], [], [], []  # cm
water_uptake, water_collar_cell, water_cyl, water_domain = [], [], [], []  # cm3
sink1d = []
sink1d2 = []
out_times = []  # days
cci = picker(rs.nodes[0].x, rs.nodes[0].y, rs.nodes[0].z)  # collar cell index
cell_volumes = s.getCellVolumes()  # cm3
cell_volumes = comm.bcast(cell_volumes, root = 0)
net_flux = np.zeros(cell_volumes.shape)

for i in range(0, NT):

    wall_iteration = timeit.default_timer()

    t = i * dt  # current simulation time
    wall_root_model = timeit.default_timer()

    sx = s.getSolutionHead_()  # richards.py
    hsb = np.array([sx[j] for j in mapping])  # soil bulk matric potential per segment
    rsx = hsb.copy()  # initial values for fix point iteration

    kr_ = r.getKr(rs_age + t)
    inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up; here const
    inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)  ############################################ (too keep within table)
    inner_kr_ = np.minimum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-4)  ############################################ (too keep within table)

    wall_iteration = timeit.default_timer()
    wall_fixpoint = timeit.default_timer()

    err = 1.e6  # cm
    c = 0

    # r.init_solve_static(rs_age + t, rsx, False, wilting_point, soil_k = [])  # LU factorisation for speed up
    rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, soil_k = [])
    rx_old = rx.copy()

    hsb_ = hsb - cell_centers_z  # from total matric potential to matric potential
    hsb_ = np.maximum(hsb_, np.ones(hsb_.shape) * -15000.)  ############################################ (too keep within table)
    hsb_ = np.minimum(hsb_, np.zeros(hsb_.shape))  ############################################ (too keep within table)
    while err > 1 and c < max_iter:

        """ interpolation """
        wall_interpolation = timeit.default_timer()
        rx_ = rx[1:] - seg_centers_z  # from total matric potential to matric potential
        rx_ = np.maximum(rx_, np.ones(rx_.shape) * -15000.)  ############################################ (too keep within table)
        rsx = root_interface(rx_ , hsb_, inner_kr_, rho_, sri_table_lookup)
        rsx = rsx + seg_centers_z  # from matric potential to total matric potential

        """ xylem matric potential """
        wall_xylem = timeit.default_timer()
        # print("Segment size from Python ", len(r.rs.segments), ns)
        rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, soil_k = [])  # xylem_flux.py, cells = False
        err = np.linalg.norm(rx - rx_old)
        wall_xylem = timeit.default_timer() - wall_xylem

        rx_old = rx.copy()
        c += 1
        print('number of iterations', c, err)

    proposed_inner_fluxes = r.segFluxes(rs_age + t, rx, rsx, approx = False, cells = False, soil_k = [])  # [cm3/day]
    seg_fluxes = np.array(proposed_inner_fluxes)
    soil_k = np.divide(vg.hydraulic_conductivity(rsx, soil), rs.radii)  # only valid for homogenous soil

    seg_fluxes = comm.bcast(seg_fluxes, root = 0)
    wall_root_model = timeit.default_timer() - wall_root_model

    """ 2. local soil models """
    wall_rhizo_models = timeit.default_timer()
    if rank == 0:
        proposed_outer_fluxes = r.splitSoilFluxes(net_flux / dt, split_type)
    else:
        proposed_outer_fluxes = None
    proposed_outer_fluxes = comm.bcast(proposed_outer_fluxes, root = 0)

    rs.solve(dt, rsx, proposed_outer_fluxes)
    # rs.solve(dt, seg_rx, proposed_outer_fluxes)
    # realized_inner_fluxes = rs.get_inner_fluxes()
    # realized_inner_fluxes = comm.bcast(realized_inner_fluxes, root = 0)
    wall_rhizo_models = timeit.default_timer() - wall_rhizo_models

    """ 3a. macroscopic soil model """
    wall_soil_model = timeit.default_timer()
    water_content = np.array(s.getWaterContent())
    water_content = comm.bcast(water_content, root = 0)
    soil_water = np.multiply(water_content, cell_volumes)  # water per cell [cm3]
    soil_fluxes = r.sumSegFluxes(seg_fluxes)  # [cm3/day]  per soil cell
    s.setSource(soil_fluxes.copy())  # [cm3/day], in richards.py
    s.solve(dt)  # in solverbase.py

    """ 3b. calculate net fluxes """
    water_content = np.array(s.getWaterContent())
    water_content = comm.bcast(water_content, root = 0)
    new_soil_water = np.multiply(water_content, cell_volumes)  # calculate net flux
    net_flux = new_soil_water - soil_water  # change in water per cell [cm3]
    for k, root_flux in soil_fluxes.items():
        net_flux[k] -= root_flux * dt
    soil_water = new_soil_water
    wall_soil_model = timeit.default_timer() - wall_soil_model

    wall_iteration = timeit.default_timer() - wall_iteration

    if i % skip == 0:
        collar_sx.append(s.getSolutionHeadAt(cci))
        min_sx.append(np.min(s.getSolutionHead()))
        water_cyl.append(np.sum(rs.get_water_volume()))  # cm3

        if rank == 0:
            out_times.append(t)
            collar_flux.append(r.collar_flux(rs_age + t, rx, rsx, [], False))
            min_rsx.append(np.min(np.array(rsx)))
            min_rx.append(np.min(np.array(rx)))
            print("Cylindrical model: minimum root soil interface {:g} cm, soil {:g} cm, root xylem {:g} cm".format(min_rsx[-1], min_sx[-1], min_rx[-1]))
            min_soil_fluxes, max_soil_fluxes, summed_soil_fluxes = 1.e9, -1.e9, 0.
            for k, v in soil_fluxes.items():
                summed_soil_fluxes += v
                if max_soil_fluxes < v:
                    max_soil_fluxes = v
                if min_soil_fluxes > v:
                    min_soil_fluxes = v
            print("Fluxes: summed local fluxes {:g}, collar flux {:g}, predescribed {:g}".format(summed_soil_fluxes, collar_flux[-1], -trans * sinusoidal(t)))
            water_domain.append(np.min(soil_water))  # cm3
            water_collar_cell.append(soil_water[cci])  # cm3
            water_uptake.append(summed_soil_fluxes)  # cm3/day
            n = round(float(i) / float(NT) * 100.)
            print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], {:g} days".format(s.simTime))
            print("Iteration {:g} took {:g} seconds [{:g}% root, {:g}% rhizo {:g}% soil ]\n".
                  format(i, wall_iteration, wall_root_model / wall_iteration, wall_rhizo_models / wall_iteration, wall_soil_model / wall_iteration))
#             if min_rsx[-1] < -16000:
#                 print("breaksim time ", sim_time)
#                 break

            """ Additional sink plot """
            if i % (60 * 12) == 0:  # every 6h
                ana = pb.SegmentAnalyser(r.rs)
                fluxes = r.segFluxes(rs_age + t, rx, rsx, approx = False, cells = False, soil_k = soil_k)  # cm3/day
                ana.addData("fluxes", fluxes)  # cut off for vizualisation
                ana.addData("fluxes2", seg_fluxes)  # cut off for vizualisation
                flux1d = ana.distribution("fluxes", max_b[2], min_b[2], 15, False)
                flux1d2 = ana.distribution("fluxes2", max_b[2], min_b[2], 15, False)
                sink1d.append(np.array(flux1d))
                sink1d2.append(np.array(flux1d2))
                # realized_inner_fluxes!!!!!!!!!!!!

""" plots and output """
if rank == 0:
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    vp.plot_roots_and_soil(r.rs, "pressure head", rsx, s, periodic, min_b, max_b, cell_number, name)  # VTK vizualisation

#     rsx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models [cm]
#     soil_k = np.divide(vg.hydraulic_conductivity(rsx, soil), rs.radii)  # only valid for homogenous soil
#     rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, soil_k)
#     fluxes = r.segFluxes(rs_age + t, rx, rsx, approx=False, cells=False, soil_k=soil_k)
#     ana = pb.SegmentAnalyser(r.rs)
#     ana.addData("rsx", rsx)
#     ana.addData("rx", rx)
#     ana.addData("fluxes", fluxes)
#     vp.plot_roots(ana, "rsx")  # VTK vizualisation
#     vp.plot_roots(ana, "fluxes")  # VTK vizualisation

    crit_min_i, crit_max_i, crit_min_o, crit_max_o = rs.plot_cylinders()
    # print(crit_min_i)
    rs.plot_cylinder(crit_min_i, name)
    print(rs.radii[crit_min_i])

    plot_transpiration(out_times, water_uptake, collar_flux, lambda t: trans * sinusoidal(t), name)  # in rhizo_models.py
    # plot_info(out_times, water_collar_cell, water_cyl, collar_sx, min_sx, min_rx, min_rsx, water_uptake, water_domain)  # in rhizo_models.py

    np.savetxt('results/' + name, np.vstack((out_times, -np.array(collar_flux), -np.array(water_uptake))), delimiter = ';')

    sink1d = np.array(sink1d)
    np.save('results/' + name + "sink1d_rhizo", sink1d)
    np.save('results/' + name + "sink1d2_rhizo", sink1d2)

    print(sink1d.shape)
