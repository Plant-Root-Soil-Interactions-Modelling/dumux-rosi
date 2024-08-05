""" 
Benchmark M1.2 static root system in soil, root hydrualics with Meunier with rhizosphere (using steady rate approach)

also works parallel with mpiexec (only slightly faster, due to overhead)
"""

import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import plantbox as pb
from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
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


def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


def open_sra_lookup(filename):
    """ opens the look-up table from a file, to quickly find soil root interface potential """
    sra_table = np.load(filename + ".npy")
    x = np.load(filename + "_.npy", allow_pickle = True)
    rx_ = x[0]
    sx_ = x[1]
    inner_ = x[2]
    outer_ = x[3]
    return RegularGridInterpolator((rx_, sx_, inner_, outer_), sra_table)  # method = "nearest" fill_value = None , bounds_error=False


def soil_root_interface(rx, sx, inner_kr, rho, sp):
    """
    finds potential at the soil root interface
    
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
        raise

    return rsx


""" Parameters """
min_b = [-4., -4., -15.]
max_b = [4., 4., 0.]
cell_number = [8, 8, 15]  #  [8, 8, 15]  # [16, 16, 30]  # [32, 32, 60]  # [8, 8, 15]
periodic = False
fname = "../../../grids/RootSystem8.rsml"

name = "c12_sra"
loam = [0.08, 0.43, 0.04, 1.6, 50]  # do not change, or adapt look up table L144
soil = loam

initial = -659.8 + 7.5  # -659.8

trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

sim_time = 7.1  # [day] for task b
dt = 360. / (24 * 3600)  # [days] Time step must be very small
skip = 1  # output

""" Initialize macroscopic soil model """
sp = vg.Parameters(soil)
vg.create_mfp_lookup(sp, -1.e5, 1000)
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

""" Initialize xylem model (a) or (b)"""
r = XylemFluxPython(fname)
r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)  # cutting
init_conductivities(r, False)

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])
r.rs.setSoilGrid(picker)  # maps segment
outer_r = r.rs.segOuterRadii()
inner_r = r.rs.radii

""" sanity checks """
# r.plot_conductivities()
r.test()  # sanity checks
rs_age = np.max(r.get_ages())
print("Root system age ", rs_age)
seg_length = r.rs.segLength()
ns = len(seg_length)
# print("press any key"); input()
print("outer radii", np.min(outer_r) , np.max(outer_r))
print("inner radii", np.min(inner_r) , np.max(inner_r))
print()

sra_table_lookup = open_sra_lookup("table_loam")

""" Numerical solution (a) """
start_time = timeit.default_timer()
x_, y_, z_ = [], [], []

mapping = r.rs.getSegmentMapper()
sx = s.getSolutionHead_()  # richards.py
hsb = np.array([sx[j] for j in mapping])
rsx = hsb.copy()  # initial values for fix point iteration
rho_ = np.divide(outer_r, np.array(inner_r))

N = round(sim_time / dt)
t = 0.

for i in range(0, N):

    rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, [])  # xylem_flux.py, cells = False
    rx_old = rx.copy()

    kr_ = r.getKr(rs_age + t)
    inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up; here const

    err = 1.e6
    c = 0
    while err > 1. and c < 100:

        """ interpolation """
        rsx = soil_root_interface_table(rx[1:], hsb, inner_kr_, rho_, sra_table_lookup)

        """ xylem matric potential """
        # print("solve", rs_age + t, -trans * sinusoidal(t), rsx.shape, wilting_point)
        rx = r.solve(rs_age + t, -trans * sinusoidal(t), 0., rsx, False, wilting_point, [])  # xylem_flux.py, cells = False
        err = np.linalg.norm(rx - rx_old)

        rx_old = rx.copy()
        c += 1

    fluxes = r.segFluxes(rs_age + t, rx, rsx, approx = False, cells = False)
    collar_flux = r.collar_flux(rs_age + t, rx.copy(), rsx.copy(), k_soil = [], cells = False)  # validity checks
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
        cf_ = r.collar_flux(rs_age + t, rx, rsx, k_soil = [], cells = False)
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

# vp.plot_roots_and_soil(r.rs, "pressure head", rx, s, periodic, min_b, max_b, cell_number, name)  # VTK vizualisation
# ana = pb.SegmentAnalyser(r.rs)
# ana.addData("pressure", rx)
# vp.plot_roots(ana, "pressure")

plot_transpiration(x_, y_, z_, lambda t: trans * sinusoidal(t))
np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter = ';')

