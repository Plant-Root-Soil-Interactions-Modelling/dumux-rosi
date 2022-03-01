import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src/python_modules/") 
sys.path.append("../../build-cmake/cpp/python_binding/")
import os

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml_reader as rsml
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import vtk_plot as vp
import van_genuchten as vg
from root_conductivities import *
from shutil import copyfile

from math import *
import numpy as np
import matplotlib.pyplot as plt
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()


def sinusoidal(t):
    return np.sin(2. * pi * np.array(t) - 0.5 * pi) + 1.

""" 
Benchmark M1.2 static root system in soil (with the classic sink)

also works parallel with mpiexec (only slightly faster, due to overhead)
"""

""" Parameters """
min_b = [-4., -4., -15.]
max_b = [4., 4., 0.]

periodic = False

approach = "27-05-2021_CStudy_classic_sink_"

# KOMMA!!!!!!!!!!!!!!!!!!!!!!! "CStudy_loam_dry_4,0cm"
name = "clay_dry_0,8cm"


#name = "CStudy_loam_wet_0,4cm"

dirname= approach +name
cell_number = [10, 10, 19]  # [8, 8, 15]  # [8, 8, 15]  # [16, 16, 30]  # [32, 32, 60]  # [8, 8, 15]

#[2, 2, 4] 	# 4   cm 
#[3, 3, 5] 	# 3   cm 
#[4, 4, 8]	# 2   cm 
#[5, 5, 10]	# 1.5 cm 
#[8, 8, 15]	# 1   cm 
#[10, 10, 19]	# 0.8 cm 
#[13, 13, 25]	# 0.6 cm 
#[20, 20, 38]	# 0.4 cm 



#0.6 cm MPI DONE AND RESULTS OK, for loam wet time step = 30 required. 0.4 grid gives wrong results with dt = 30;

sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil = clay

kz = 4.32e-2  # axial conductivity [cm^3/day]
kr = 1.728e-4  # radial conductivity [1/day]

# dry scenario for clay & loam
initial =  -659.8 + 7.5  # -659.8

#wet scenario for clay & loam
#initial =  -300 + 7.5  # -659.8

# sand scenarios
#initial =  -100 + (max_b[2]-min_b[2])/2 

trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

sim_time = 14 # [day] for task b
age_dependent = False  # conductivities
dt = 120 / (24 * 3600)  # [days] Time step must be very small

""" Initialize macroscopic soil model """
sp = vg.Parameters(soil)  # for debugging
cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
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
#s.setRegularisation(1.e-3, 1e-3)

""" Initialize xylem model (a) or (b)"""
fname = "../../grids/RootSystem8.rsml"
r = XylemFluxPython(fname)
r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)  # cutting
init_conductivities(r, age_dependent)
r.rs.sort()  # ensures segment is located at index s.y-1
r.test()  # sanity checks
nodes = r.get_nodes()
rs_age = np.max(r.get_ages())
r.setKr([kr])  # or use setKrTables, see XylemFlux.h
r.setKx([kz])

""" DOCUMENTATION """
if not os.path.exists("convergence_study_27-05"):
    os.mkdir("convergence_study_27-05")

if not os.path.exists("convergence_study_27-05/"+dirname):
    os.mkdir("convergence_study_27-05/"+dirname)
os.chdir("convergence_study_27-05/"+dirname)

#copyfile("/home/tobi/Schreibtisch/DUMUX_up2date/dumux-rosi/python/coupled/convergence_c12.py", "/home/tobi/Schreibtisch/DUMUX_up2date/dumux-rosi/python/coupled/convergence_study_27-05/"+dirname+"/documentation_"+dirname+".py")

""" Coupling (map indices) """
picker = lambda x, y, z : s.pick([x, y, z])
r.rs.setSoilGrid(picker)  # maps segments
cci = picker(nodes[0, 0], nodes[0, 1], nodes[0, 2])  # collar cell index

""" Numerical solution (a) """
start_time = timeit.default_timer()
x_, y_, w_, cpx, cps, mean_bulk_soil_pressure, min_bulk_soil_pressure, max_bulk_soil_pressure, ESWP, mean_xylem_pressure , ESWP_check = [], [], [], [], [], [] , [], [] , [] , [] , []
sx = s.getSolutionHead()  # inital condition, solverbase.py

N = round(sim_time / dt)
t = 0.


#Tobi
suf = r.get_suf(0.)
print("Sum of SUF", np.sum(suf), "from", np.min(suf), "to", np.max(suf), "summed positive", np.sum(suf[suf >= 0]))

#breakpoint()

for i in range(0, N):

    if rank == 0:  # Root part is not parallel

        rx = r.solve(rs_age + t, -trans * sinusoidal(t), sx[cci], sx, True, wilting_point, [])  # xylem_flux.py, cells = True

        fluxes = r.soilFluxes(rs_age + t, rx, sx, False)  # class XylemFlux is defined in MappedOrganism.h, approx = True

        sum_flux = 0.
        for f in fluxes.values():
            sum_flux += f
        print("Summed fluxes ", sum_flux, "= collar flux", r.collar_flux(rs_age + t, rx, sx), "= prescribed", -trans * sinusoidal(t))

    else:
        fluxes = None

    fluxes = comm.bcast(fluxes, root = 0)  # Soil part runs parallel
    s.setSource(fluxes)  # richards.py

    s.solve(dt)

    sx = s.getSolutionHead()  # richards.py
    water = s.getWaterVolume()


    if rank == 0:
        """ Equivalent soil water potential """
        eswp = r.get_eswp(t, sx)
        print("root collar pressure             " +  str(rx[0]))
        print("Equivalent soil water potential", eswp)

        n = round(float(i) / float(N) * 100.)
        min_sx = np.min(sx)
        min_rx = np.min(rx)
        max_sx = np.max(sx)
        max_rx = np.max(rx)
        mean_bulk_soil_pressure.append(np.mean(sx))
        mean_xylem_pressure.append(np.mean(rx))
        
        print("mean bulk soil pressure          " + str(np.mean(sx)))
        min_bulk_soil_pressure.append(min_sx)
        max_bulk_soil_pressure.append(max_sx)
        ESWP.append(float(eswp))
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "\n" + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g} days {:g}"
              .format(min_sx, max_sx, min_rx, max_rx, s.simTime, rx[0]))
        f = float(r.collar_flux(rs_age + t, rx, sx))  # exact root collar flux
        x_.append(t)
        y_.append(sum_flux)
        w_.append(water)
        cpx.append(rx[0])
        cps.append(float(sx[cci]))


        eswp_check = eswp - rx[0]
        if (eswp_check < 0):
             ESWP_check.append(eswp_check)
    t += dt

s.writeDumuxVTK(dirname)



""" Plot """
if rank == 0:
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")
    #vp.plot_roots_and_soil(r.rs, "pressure head", rx, s, periodic, min_b, max_b, cell_number, name)  # VTK vizualisation

    print(" ESWP < collar pressure?  ", ESWP_check)
    WCT= round(timeit.default_timer() - start_time) # grab wall-clock-time for plot

    # VTK vizualisation
    #os.chdir("convergence_study/"+dirname)
    fig, ax1 = plt.subplots()
    ax1.plot(x_, trans * sinusoidal(x_), 'k')  # potential transpiration
    ax1.plot(x_, -np.array(y_), 'g')  # actual transpiration (neumann)
    ax2 = ax1.twinx()
    ax2.plot(x_, np.cumsum(-np.array(y_) * dt), 'c--')  # cumulative transpiration (neumann)
    ax1.set_xlabel("Time [d]")
    ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
    ax1.legend(['Potential', 'Actual', 'Cumulative'], loc = 'upper left')
    # default np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter = ';')
    # Simtime, Tact, Tpot, TCum
    fig.suptitle(approach + name + "  WCT= " + str(WCT)+"s", )#y=1)
    np.savetxt(dirname, np.vstack((x_, -np.array(y_), trans * sinusoidal(x_), np.cumsum(-np.array(y_) * dt) )), delimiter = ';')
    fig.savefig(str(dirname))
    plt.show()

 
    fig, ax1 = plt.subplots()
    ax1.plot(x_, mean_bulk_soil_pressure, "chocolate")
    ax1.plot(x_, ESWP, "b--")
    #ax1.plot(x_, mean_xylem_pressure, "g")
    ax1.plot(x_, cpx, "g")
    #ax1.plot(x_, min_bulk_soil_pressure)
    ax1.set_xlabel("Time [d]")
    ax1.set_ylabel("soil pressure [cm]")
    ax1.legend(["mean bulk soil pressure", "ESWP", "root collar pressure"], loc = "upper right") 
    fig.suptitle(approach + name + "  WCT= " + str(WCT)+"s", )#y=1)
    fig.set_tight_layout(True)
    fig.savefig(str(dirname+ "_pressures"))
    np.savetxt(dirname + "_pressures", np.vstack((x_, mean_bulk_soil_pressure, ESWP, cpx, max_bulk_soil_pressure, min_bulk_soil_pressure, mean_xylem_pressure)), delimiter = ';')
    plt.show()


    fig, ax1 = plt.subplots()
    ax1.plot(x_, mean_bulk_soil_pressure, "chocolate")
    ax1.plot(x_, ESWP, "b--")
    ax1.plot(x_, mean_xylem_pressure, "g")
    ax1.plot(x_, cpx)
    ax1.plot(x_, min_bulk_soil_pressure)
    ax1.set_xlabel("Time [d]")
    ax1.set_ylabel("soil pressure [cm]")
    ax1.legend(["mean bulk soil pressure", "ESWP", "mean xylem pressure", "cpx", "min bulk soil pressure"], loc = "upper right") 
    fig.suptitle(approach + name + "  WCT= " + str(WCT)+"s", )#y=1)
    fig.set_tight_layout(True)
    fig.savefig(str(dirname+ "_pressures_detailed"))
    plt.show()
