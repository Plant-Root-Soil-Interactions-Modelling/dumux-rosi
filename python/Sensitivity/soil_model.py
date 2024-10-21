""" 
Functions to simplify setup of the scenarios for the INARI project
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

# from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); size = comm.Get_size()
import numpy as np
# import timeit
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
# from datetime import *

from rosi_richardsnc import RichardsNCSP  # C++ part (Dumux binding), macroscopic soil model
from rosi_richards import RichardsSPnum as  RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model

# import plantbox as pb  # CPlantBox
import functional.van_genuchten as vg
# from functional.PlantHydraulicParameters import PlantHydraulicParameters
# from functional.PlantHydraulicModel import HydraulicModel_Doussan
# from functional.PlantHydraulicModel import HydraulicModel_Meunier
# import evapotranspiration as evap


def create_richards(soil_, min_b , max_b , cell_number, times = None, net_inf = None, bot_bc = "noFlux", bot_value = 0.):
    """
        Creates a soil domain from @param min_b to @param max_b with resolution @param cell_number
        
        soil type is fixed and homogeneous !!! 
        
        domain is periodic (if 2d or 3d)
        
        returns soil_model (RichardsWrapper(RichardsSP()))
    """
    wet = False  # wet scenario

    s = RichardsWrapper(RichardsSP())  # water only
    s.initialize()
    s.createGrid(min_b, max_b, cell_number, periodic = False)  # [cm]  # periodicity is not needed, since coupling via s.pick(0., 0., z.)

    # BC
    if times is not None:
        if wet:  # wet scenario increases precipitation and decreases evaporation
            net_inf[net_inf > 0] = net_inf[net_inf > 0] * 1.2  # increase precipitation for 20%
            net_inf[net_inf < 0] = net_inf[net_inf < 0] * 0.8  # decrease evaporation for 20%
        s.setTopBC("atmospheric", 0.5, [times, net_inf])  # 0.5 is dummy value
    else:
        s.setTopBC("noFlux")
    # s.setBotBC("freeDrainage")
    # s.setBotBC("noFlux")
    # s.setBotBC("potential", 80)
    s.setBotBC(bot_bc, bot_value)

    # Parameters
    s.setVGParameters([soil_])
    s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
    s.setParameter("Soil.SourceSlope", "1000")  # turns regularisation of the source term on
    if type == 2:
        s.setParameter("Component.MolarMass", "6.2e-2")  # nitrate 62,0049 g/mol
        s.setParameter("Component.LiquidDiffusionCoefficient", "1.7e-9")  # m2 s-1 # nitrate = 1700 um^2/sec

    s.setHomogeneousIC(-1000)

    s.initializeProblem()
    wilting_point = -15000
    s.setCriticalPressure(wilting_point)  # for boundary conditions constantFlow, constantFlowCyl, and atmospheric
    s.ddt = 1.e-5  # [day] initial Dumux time step

    # # IC
    # h = np.load("data/initial_potential.npy")
    # s.setInitialConditionHead(h)  # cm

    # h = s.getSolutionHead()
    # plt.plot(h, np.linspace(-200., 0., h.shape[0]))
    # plt.xlabel("soil matric potential [cm]")
    # plt.ylabel("depth (cm)")
    # plt.tight_layout()
    # plt.show()
    # plt.plot(c, np.linspace(-200, 0., c.shape[0]))
    # plt.xlabel("nitrate concentration [g/cm3]")
    # plt.ylabel("depth (cm)")
    # plt.tight_layout()
    # plt.show()

    return s

