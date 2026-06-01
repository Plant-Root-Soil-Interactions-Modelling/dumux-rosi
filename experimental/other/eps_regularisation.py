import sys; sys.path.append("../../python/modules");  sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi.richards import RichardsWrapper  # Python part |\label{l61:paths_a}|
from rosi.rosi_richards import RichardsSP  # C++ part (Dumux binding) |\label{paths_e}|
from rosi.rosi_richardsnc_cyl import RichardsNCCylFoam

import plantbox.functional.van_genuchten as vg
import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import numpy as np

""" 
select the eps_reg for the sand param set
"""


def solve():
    x= [-0.00087584, -0.00166197, -0.00229157, -0.00284242, -0.00333934, -0.00379047,
                -0.00419244, -0.00452672, -0.00474684]
    soil = [0.03, 0.414, 0.038, 2, 1864]
    s = RichardsWrapper(RichardsNCCylFoam())
    logbase = 0.5
    NC = 10  # spatial resolution (1D model)
    a_in = 0.02  # cm
    a_out = 0.6  # cm
    points = np.logspace(np.log(a_in) / np.log(logbase), np.log(a_out) / np.log(logbase), NC, base=logbase)
    s.initialize()
    s.createGrid1d(points)
    s.setHomogeneousIC(x, False)  # cm pressure head
    s.setBotBC("constantFlux", 0)
    s.setTopBC("constantFlux", 0)  #  [cm/day]
    s.setParameter("Component.MolarMass", str(12))
    s.setParameter("Soil.MolarMass", str(12))
    s.setParameter("Soil.solidDensity", str(1000))
    s.setParameter("Flux.UpwindWeight", "1")
    s.setVGParameters([soil])
    s.initializeProblem()

    print('x at',x, end=', ')
    cm_init = np.full(NC - 1, -x)
    Pa_init = np.array(s.to_pa(cm_init))
    s.base.setSolution(Pa_init, 0)
    Pa_divide = Pa_init
    
    Pa_divide[Pa_divide == 0] = 1.
    
    thetainit = s.getWaterContent()
    pHead1 = s.getSolution()
    pHead = s.getSolutionHead()
    thetainitth = np.array([vg.water_content( p_mean_, s.vg_soils[0]) for
                            p_mean_ in pHead])
    theta_divide = np.where(thetainitth!=0,thetainitth,1)
    assert len(pHead) == (NC - 1)
    #print('pHead',pHead,pHead1,'Pa_init',Pa_init)
    assert (np.logical_or( (abs((pHead1 - Pa_init)/Pa_divide)*100 < 1e-5) , 
                           (abs(pHead1- Pa_init) < 1e-9) )).all()
    #print('thetainit',thetainit,'thetainitth',thetainitth)
    assert (np.logical_or( (abs((thetainit - thetainitth)/theta_divide)*100 < 1e-5) ,
                           (abs(thetainit- thetainitth) < 1e-9) )).all()
    print(thetainitth)
        
if __name__ == "__main__":

    solve()

