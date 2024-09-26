import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

from functional.xylem_flux import *  # Python hybrid solver
import plantbox as pb
import rsml.rsml_reader as rsml
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding)
from richards import RichardsWrapper  # for soil part wrappring RichardsSP
from richards_no_mpi import RichardsNoMPIWrapper  # for cylindrical part, RichardsCylFoam must be compiled disabeling MPI

from math import *
import numpy as np
import matplotlib.pyplot as plt
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

""" 
THINK
"""
r_root = 0.01  # cm root radius
r_out = 1 # cm, perirhizal radius


def solve(soil, simtimes, q_r, N):
    """ 
    soil            soil type
    simtimes        simulation output times
    q_r             root flux [cm3/day]   
    N               spatial resolution
    """
    rc = RichardsNoMPIWrapper(RichardsCylFoam())
    rc.initialize([""], False) # )  # [""], False
    rc.createGrid([r_root], [r_out], [N])  # [cm]
    rc.setHomogeneousIC(-10.)  # cm pressure head
    rc.setVGParameters([soil[0:5]])
    rc.setOuterBC("fluxCyl", 0.)  #  [cm/day]
    rc.setInnerBC("fluxCyl", -q_r)  # [cm/day]
    rc.initializeProblem()
    rc.setCriticalPressure(-15000)  # cm pressure head

    """ Numerical solution """
    start_time = timeit.default_timer()
    dt_ = np.diff(simtimes)


    # initial_water = rc.getWaterVolume() # cm3
    # print("initial_water", initial_water)
    t_ = []
    t =0 
    cum_qr = 0
    for i, dt in enumerate(dt_):
        
        t += dt 
        t_.append(t)
        rc.setInnerBC("fluxCyl", -q_r*sinusoidal2(t,dt))  # [cm/day]       
        rc.solve(dt)
        qr_act_dict = rc.getAllNeumann()  
        qr_act = qr_act_dict[0] * 1.e4  
        rsx = rc.getInnerHead()
        cum_qr += qr_act*dt

        # summary
        # print("time", simtimes[i], "rsx", rsx, "qr_pot", -q_r, "qr_act", qr_act)
        
        # Questions:  *) does water volume equals intial water volume - cumulative uptake
        #             *) qr_act < q_pot only if rsx<-15000 (that is not the case, some regularisation?)
        
    # final_water = rc.getWaterVolume() # cm3
    # print("final_water", initial_water)    
    
    y = rc.getSolutionHead()  
    plt.clf()
    plt.plot(np.linspace(r_root, r_out, N), y, "b", label = "cylindrical")
    plt.title("qr {:g} [cm/day]".format(qr))
    plt.ylabel("matric potential [cm]")
    plt.xlabel("radius [cm]")
    # plt.show()    
    plt.savefig(soil[-1]+"_dynamic_peak_{:g}.png".format(qr))
    
    
    # plt.plot(t_, [q_r*sinusoidal2(t,dt) for t in t_])
    # plt.show()

    return cum_qr


if __name__ == "__main__":


    N = 101

    # replace by hydrus soils?
    sand = [0.045, 0.43, 0.15, 3, 1000, "Sand"]
    loam = [0.08, 0.43, 0.04, 1.6, 50, "Loam"]
    clay = [0.1, 0.4, 0.01, 1.1, 10, "Clay"] 

    sim_times = np.linspace(0, 6.5, 24*6*7)  # temporal resolution of 0.1 h
    QRN = 4

    qr_ = 100*np.linspace(0.01, 0.1, QRN) / 7 # [cm/day]    
    soil = sand
    cum_qr = np.zeros(qr_.shape)
    for i, qr in enumerate(qr_):
        print(i, "qr", qr)
        cum_qr[i] = solve(soil, sim_times, qr, N)
    print("Sand (dynamic)")
    print("qr_", np.min(qr_), np.max(qr_))#
    print("cum_qr", np.min(cum_qr), np.max(cum_qr))
    plt.title("Sand (dynamic)")
    plt.plot(qr_,cum_qr)
    plt.xlabel("qr [cm/day]")
    plt.ylabel("cum_qr [cm]")
    plt.show()
    
    qr_ = 100*np.linspace(0.05, 0.5, QRN) / 7 # [cm/day]
    soil = loam
    cum_qr = np.zeros(qr_.shape)
    for i, qr in enumerate(qr_):
        cum_qr[i] = solve(soil, sim_times, qr, N)
    print("Loam (dynamic)")
    print("qr_", np.min(qr_), np.max(qr_))#
    print("cum_qr", np.min(cum_qr), np.max(cum_qr))
    plt.title("Loam (dynamic)")
    plt.plot(qr_,cum_qr) 
    plt.xlabel("qr [cm/day]")
    plt.ylabel("cum_qr [cm]")
    plt.show()
    
    
    qr_ = 100*np.linspace(0.01, 0.1, QRN)/ 7 # [cm/day]
    soil = clay
    cum_qr = np.zeros(qr_.shape)
    for i, qr in enumerate(qr_):
        cum_qr[i] = solve(soil, sim_times, qr, N)
    print("Clay (dynamic)")
    print("qr_", np.min(qr_), np.max(qr_))#
    print("cum_qr", np.min(cum_qr), np.max(cum_qr))
    plt.title("Clay (dynamic)")
    plt.plot(qr_,cum_qr)
    plt.xlabel("qr [cm/day]")
    plt.ylabel("cum_qr [cm]")
    plt.show()

