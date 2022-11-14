""" 
Single root scenario: Soil depletion due to constant Dirichlet collar potential

using steady rate approach and fix point iteration
"""

import sys; sys.path.append("../../../modules/"); sys.path.append("../../../../../CPlantBox/");  sys.path.append("../../../../../CPlantBox/src/python_modules")
sys.path.append("../../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../modules/fv/");
sys.path.append("../"); sys.path.append("../scenarios/");

import plantbox as pb  # CPlantBox
from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
from xylem_flux import *  # root system Python hybrid solver
from rhizo_models import *  # Helper class for cylindrical rhizosphere models

import vtk_plot as vp
import van_genuchten as vg
from sra_table_lookup import *
import aggregated_rs as agg

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve


def soil_root_interface(rx, sx, inner_kr, rho, sp):
    """
    rx             xylem matric potential [cm]
    sx             bulk soil matric potential [cm]
    inner_kr       root radius times hydraulic conductivity [cm/day] 
    rho            geometry factor [1]
    sp             soil van Genuchten parameters (type vg.Parameters)
    """

    k_soilfun = lambda hsoil, hint: (vg.fast_mfp[sp](hsoil) - vg.fast_mfp[sp](hint)) / (hsoil - hint)

    # rho = outer_r / inner_r  # Eqn [5]
    rho2 = rho * rho  # rho squared
    # b = 2 * (rho2 - 1) / (1 + 2 * rho2 * (np.log(rho) - 0.5))  # Eqn [4]
    b = 2 * (rho2 - 1) / (1 - 0.53 * 0.53 * rho2 + 2 * rho2 * (np.log(rho) + np.log(0.53)))  # Eqn [7]

    fun = lambda x: (inner_kr * rx + b * sx * k_soilfun(sx, x)) / (b * k_soilfun(sx, x) + inner_kr) - x
    rsx = fsolve(fun, rx)

    return rsx


def soil_root_interface_table2(rx, sx, inner_kr_, rho_, f):
    assert rx.shape == sx.shape
    rsx = f((rx, sx, inner_kr_ , rho_))
    return rsx

""" 
Parameters  
"""

""" soil """
p_top = -300  # -5000 (_dry), -300 (_wet)
p_bot = -200  #
sstr = "_wet"  # <---------------------------------------------------------- (dry or wet)

min_b = [-0.5, -0.5, -50.]
max_b = [0.5, 0.5, 0.]
cell_number = [1, 1, 100]
periodic = False

alpha = 0.018  # (cm-1) # spoil
n = 1.8
Ks = 28.46  # (cm d-1)
loam = [0.08, 0.43, alpha, n, Ks]
soil_ = loam
soil = vg.Parameters(soil_)
vg.create_mfp_lookup(soil, -1.e5, 1000)  # creates the matrix flux potential look up table (in case for exact)
sra_table_lookup = open_sra_lookup("../table_jan2")  # opens the precomputed soil root interface potentials

""" root system """
collar = -8000  # dirichlet
radius = 0.05  # cm
wilting_point = -15000

""" simulation time """
sim_time = 0.51  # 0.65  # 0.25  # [day]
dt = 60 / (24 * 3600)  # time step [day], 120 schwankt stark
NT = int(np.ceil(sim_time / dt))  # number of iterations
skip = 1  # for output and results, skip iteration

""" 
Initialize macroscopic soil model (Dumux binding)
"""
s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setLinearIC(p_top, p_bot)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([soil_])
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.initializeProblem()
s.setCriticalPressure(wilting_point)  # new source term regularisation
s.ddt = 1.e-5  # [day] initial Dumux time step

water0 = s.getWaterVolume()

""" 
Initialize xylem model 
"""
ns = 100
r = agg.create_singleroot(ns = ns, l = 50, a = radius)
r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), cut = False)
picker = lambda x, y, z: s.pick([x, y, z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
r.rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
agg.init_comp_conductivities_const(r)
nodes = r.rs.nodes
segs = r.rs.segments

""" sanity checks """
r.test()  # sanity checks
mapping = np.array([r.rs.seg2cell[j] for j in range(0, ns)])
outer_r = r.rs.segOuterRadii()
inner_r = r.rs.radii
types = r.rs.subTypes
rho_ = np.divide(outer_r, np.array(inner_r))
print("Krs", r.get_krs(0.))

""" Numerical solution (a) """
start_time = timeit.default_timer()

# for post processing
psi_x_ = []
psi_s_ = []
psi_s2_ = []
sink_ = []
x_, y_ = [], []

sx = s.getSolutionHead()  # inital condition, solverbase.py
cell_centers = s.getCellCenters()
cell_centers_z = np.array([cell_centers[mapping[j]][2] for j in range(0, ns)])
seg_centers_z = np.array([0.5 * (nodes[s.x].z + nodes[s.y].z) for s in segs])
hsb = np.array([sx[mapping[j]][0] for j in range(0, ns)])  # soil bulk matric potential per segment
rsx = hsb.copy()  # initial values for fix point iteration

t = 0.
rs_age = 0.

for i in range(0, NT):

    t = i * dt  # current simulation time

    wall_iteration = timeit.default_timer()
    wall_fixpoint = timeit.default_timer()

    if i == 0:  # only first time
        rx = r.solve_dirichlet(rs_age + t, [collar], 0., rsx, cells = False, soil_k = [])
        rx_old = rx.copy()

    kr_ = np.array([r.kr_f(rs_age + t, types[j]) for j in range(0, len(outer_r))])
    inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up

    err = 1.e6
    c = 0
    while err > 1 and c < 100:

        """ interpolation """
        wall_interpolation = timeit.default_timer()
        rx_ = rx[1:] - np.array([nodes[ii].z for ii in range(1, ns + 1)])  # from total matric potential to matric potential
        hsb_ = hsb - cell_centers_z  # from total matric potential to matric potential
        rsx = soil_root_interface_table2(rx_ , hsb_, inner_kr_, rho_, sra_table_lookup)  # rsx = soil_root_interface(rx[1:] , hsb, inner_kr_, rho_, soil)
        rsx = rsx + np.array([nodes[ii].z for ii in range(1, ns + 1)])  # from matric potential to total matric potential
        wall_interpolation = timeit.default_timer() - wall_interpolation

        """ xylem matric potential """
        wall_xylem = timeit.default_timer()
        rx = r.solve(0., [collar], 0., rsx, False, wilting_point, soil_k = [])  # xylem_flux.py, cells = False
        err = np.linalg.norm(rx - rx_old)
        wall_xylem = timeit.default_timer() - wall_xylem

        rx_old = rx.copy()
        c += 1

    print(i, c, "iterations", wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem))

    # wall_fixpoint = timeit.default_timer() - wall_fixpoint

    fluxes = r.segFluxes(rs_age + t, rx, rsx, False)
    collar_flux = r.collar_flux(rs_age + t, rx.copy(), rsx.copy(), k_soil = [], cells = False)
    err = np.linalg.norm(np.sum(fluxes) - collar_flux)
    if err > 1.e-10:
        print("error: summed root surface fluxes and root collar flux differ" , err)
        raise

    wall_soil = timeit.default_timer()

    soil_fluxes = r.sumSegFluxes(fluxes)
    s.setSource(soil_fluxes.copy())  # richards.py
    s.solve(dt)
    sx = s.getSolutionHead()  # richards.py
    hsb = np.array([sx[mapping[j]][0] for j in range(0, ns)])  # soil bulk matric potential per segment
    wall_soil = timeit.default_timer() - wall_soil

    wall_iteration = timeit.default_timer() - wall_iteration

    """ remember results ... """
    if i % skip == 0:
        x_.append(t)
        sum_flux = 0.
        for f in soil_fluxes.values():
            sum_flux += f
        y_.append(sum_flux)  # cm3/day
        psi_x_.append(rx.copy())
        psi_s_.append(rsx.copy())
        dd = np.array(sx)
        psi_s2_.append(dd[:, 0])
        sink_.append(fluxes.copy())

print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

water_end = s.getWaterVolume()
print("\ntotal uptake", water0 - water_end, "cm3")

""" file output """
file1 = 'results/psix_singleroot_sra_constkrkx' + sstr
np.save(file1, np.array(psi_x_))  # xylem pressure head per segment [cm]

file2 = 'results/psiinterface_singleroot_sra_constkrkx' + sstr
np.save(file2, np.array(psi_s_))  # pressure head at interface per segment [cm]

file3 = 'results/sink_singleroot_sra_constkrkx' + sstr
np.save(file3, -np.array(sink_))  # sink per segment [cm3/day]

file4 = 'results/transpiration_singleroot_sra_constkrkx' + sstr
np.save(file4, np.vstack((x_, -np.array(y_))))  # time [day], transpiration [cm3/day]

file5 = 'results/soil_singleroot_sra_constkrkx' + sstr
np.save(file5, np.array(psi_s2_))  # soil potential per cell [cm]

print("fin")

