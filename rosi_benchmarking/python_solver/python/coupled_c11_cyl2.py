import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

from solver.xylem_flux import XylemFluxPython  # Python hybrid solver
import solver.plantbox as pb
import solver.rsml_reader as rsml
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding)
from solver.richards import RichardsWrapper  # Python part

import van_genuchten as vg

from math import *
import numpy as np
import matplotlib.pyplot as plt
import timeit

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

""" 
Benchmark M1.1 single root in in soil

Prestudy for cylindrical models around segments:
soil sink is based on the cylindrical model, no information from mascroscopic soil to cylindrical model (yet)

also works parallel with mpiexec (only slightly faster, due to overhead)
"""

NC = 20 
r_root = 0.02  # cm
r_out = 0.6


def solve(soil, simtimes, q_r, N):
    """ 
    soil            soil type
    simtimes        simulation output times
    q_r             root flux [cm/day]   
    N               spatial resolution
    """
    
    q_r = q_r * 1  # * 1 g/cm3 = g/cm2/day
    q_r = q_r * 2 * r_root * np.pi * 1.  # g/day
    print("Qr as sink", q_r, "g day-1")
    trans = q_r
    l = np.sqrt((r_out * r_out - r_root * r_root) * np.pi) / 2  # same area as cylindrical

    """ Root problem"""
    n1 = pb.Vector3d(0, 0, 0)
    n2 = pb.Vector3d(0, 0, -0.5)
    n3 = pb.Vector3d(0, 0, -1)
    seg1 = pb.Vector2i(0, 1)
    seg2 = pb.Vector2i(1, 2)
    # rs = pb.MappedSegments([n1, n3], [seg1], [r_root])  # a single root
    rs = pb.MappedSegments([n1, n2, n3], [seg1, seg2], [r_root, r_root])  # a single root
    rs.setRectangularGrid(pb.Vector3d(-l, -l, -1.), pb.Vector3d(l, l, 0.), pb.Vector3d(N, N, 1))
    r = XylemFluxPython(rs)
    r.setKr([1.e-7])
    r.setKx([1.e-7])
    r_outer = r.segOuterRadii()
    seg_length = r.segLength()

    """ Soil problem """
    cpp_base = RichardsSP()
    s = RichardsWrapper(cpp_base)
    s.initialize()    
    s.setHomogeneousIC(-100)  # cm pressure head
    s.setTopBC("noflux")
    s.setBotBC("noflux")
    s.createGrid([-l, -l, -1.], [l, l, 0.], [N, N, 1])  # [cm]
    s.setVGParameters([soil[0:5]])
    s.initializeProblem()
    s.setCriticalPressure(-15000) 
    # s.setRegularisation(1.e-6, 1.e-6) 

    """ Cylindrical problems """
    rich_cyls = []
    ns = len(r.rs.segments)
    rsx = np.zeros((ns,))  # root soil interface
    for i in range(0, ns):     
        cpp_base = RichardsCylFoam()
        rcyl = RichardsWrapper(cpp_base)    
        rcyl.initialize([""], False)  # [""], False
        rcyl.createGrid([r_root], [r_outer[i]], [NC])  # [cm]        
        rcyl.setHomogeneousIC(-100.)  # cm pressure head
        rcyl.setVGParameters([soil[0:5]])
        rcyl.setOuterBC("fluxCyl", 0.)  #  [cm/day]
        # rcyl.setOuterBC("pressure", -100)  # [cm/day]  
        rcyl.setInnerBC("fluxCyl", 0.)  # [cm/day]         
        rcyl.initializeProblem()
        rcyl.setCriticalPressure(-15000)  # cm pressure head                
        rcyl.setRegularisation(1.e-4, 1.e-4)
        rcyl.ddt = 1.e-5  # days                                
        rich_cyls.append(rcyl)    

    """ Coupling (map indices) """
    picker = lambda x, y, z : s.pick([x, y, z])
    r.rs.setSoilGrid(picker)
    cci = picker(0, 0, 0)  # collar cell index

    """ Numerical solution """
    start_time = timeit.default_timer()
    nsp = N  # number of sample points for output
    y_ = np.linspace(0, l, nsp)
    y_ = np.expand_dims(y_, 1)
    x_ = np.hstack((np.zeros((nsp, 1)), y_, np.zeros((nsp, 1))))

    dt_ = np.diff(simtimes)
    s.ddt = 1.e-5  # initial Dumux time step [days]

    for dt in dt_:

        sx = s.getSolutionHead()

        if rank == 0:  # Root part is not parallel
            
            for j, rc in enumerate(rich_cyls):
                rsx[j] = rc.getInnerHead()            
                       
#             for j, rc in enumerate(rich_cyls):  # set outer BC 
#                 print("set to", sx[r.rs.seg2cell[j]], "cm")
#                 rc.setOuterPressure(sx[r.rs.seg2cell[j]])
            
            # rx = r.solve_neumann(0., -trans, rsx, False)  # xylem_flux.py, True: soil matric potential per cell, False: soil matric potential per segment
            
            rx = r.solve(0., -trans, sx[cci], rsx, False, -15000)  # xylem_flux.py                    
            # print("Minimal root xylem pressure", np.min(np.array(rx))) # working
            
            # For the soil model 
            seg_fluxes = r.segFluxes(0., rx, rsx, approx=False)  # class XylemFlux is defined in CPlantBox XylemFlux.h
            soil_fluxes = r.sumSoilFluxes(seg_fluxes)  # class XylemFlux is defined in CPlantBox XylemFlux.h                        
            
            fluxes_exact = r.soilFluxes(0., rx, sx, False)  # class XylemFlux is defined in MappedOrganism.h            
            collar_flux = r.collar_flux(0., rx, sx)
            print("fluxes at ", cci, "exact", fluxes_exact[cci], "collar flux", collar_flux, "[g day-1]")  # == cm3 / day
            print("fluxes at ", cci, "exact", soil_fluxes[cci])
                        
            fluxes = soil_fluxes
            sum_flux = 0.
            for f in fluxes.values():
                sum_flux += f
            print("Summed fluxes {:g}".format(sum_flux + trans))
            
            # run cylindrical model            
            for j, rc in enumerate(rich_cyls):  # set sources
                l = seg_length[j]
                print("Set inner flux to", seg_fluxes[j] / (2 * np.pi * r_root * l), "[cm day-1]")  
                rc.setInnerFluxCyl(seg_fluxes[j] / (2 * np.pi * r_root * l))  # /  
                print("set to", sx[r.rs.seg2cell[j]], "cm")
                rc.setOuterBC("pressure", rc.to_pa(sx[r.rs.seg2cell[j]]))
            
            for j, rc in enumerate(rich_cyls):  # simualte time step
                print("Solving segment", j)
                try:
                    rc.solve(dt)
                except:
                    x = rc.getDofCoordinates()
                    y = rc.getSolutionHead()
                    plt.plot(x, y)
                    plt.show()                      
        else:
            fluxes = None

        fluxes = comm.bcast(fluxes, root=0)  # Soil part runs parallel
        s.setSource(fluxes)  # g day-1, richards.py

        s.solve(dt)

        x0 = s.getSolutionHeadAt(cci)
        if x0 < -14000:
            if rank == 0:
                print("Simulation time at -15000 cm > {:.3f} cm after {:.3f} days".format(float(x0), s.simTime))
            y = s.to_head(s.interpolateNN(x_))
            yc = []
            for rc in rich_cyls:
                yc.append(rc.getSolutionHead())
            return y, yc, x_[:, 1], s.simTime

    y = s.to_head(s.interpolateNN(x_))
    yc = []
    for rc in rich_cyls:
        yc.append(rc.getSolutionHead())
    return y, yc, x_[:, 1], s.simTime


if __name__ == "__main__":

    N = 20

    sand = [0.045, 0.43, 0.15, 3, 1000, "Sand"]
    loam = [0.08, 0.43, 0.04, 1.6, 50, "Loam"]
    clay = [0.1, 0.4, 0.01, 1.1, 10, "Clay"]

    sim_times = np.linspace(0, 25, 250)  # temporal resolution of 0.1 d

    if rank == 0:
        fig, ax = plt.subplots(2, 3, figsize=(14, 14))
        t0 = timeit.default_timer()

    jobs = ([sand, 0.1, 0, 0], [loam, 0.1, 0, 1], [clay, 0.1, 0, 2], [sand, 0.05, 1, 0], [loam, 0.05, 1, 1], [clay, 0.05, 1, 2])

    for soil, qj, i, j in jobs:
        print("**********************************")
        print(soil[5])
        print("**********************************")
        y, yc, x, t = solve(soil, sim_times, qj, N)
        if rank == 0:
            ax[i, j].plot(x, y, "r*", label="sink based")
            for yc_ in yc:
                ax[i, j].plot(np.linspace(r_root, r_out, NC), yc_, "b", label="cylindrical")
            ax[i, j].set_xlabel("r (cm)")
            ax[i, j].set_ylabel("water potential (cm)")
            ax[i, j].title.set_text(soil[5] + ", q = {:g} cm/d, final: {:g} d".format(qj, t))
            # ax[i, j].set_ylim(-15000, 0.)
            ax[i, j].legend()                        

    if rank == 0:
        print("Elapsed time: ", timeit.default_timer() - t0, "s")
        plt.show()
