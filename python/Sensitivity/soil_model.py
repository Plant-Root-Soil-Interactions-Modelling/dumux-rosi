""" 
    Sets up the soil model 
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
# from rosi_richards import RichardsSP
from richards import RichardsWrapper  # Python part, macroscopic soil model

import functional.van_genuchten as vg


def create_richards(soil_, min_b , max_b , cell_number, times = None, net_inf = None, bot_bc = "noFlux", bot_value = 0., initial_totalpotential = -100):
    """
        Creates a soil domain from @param min_b to @param max_b with resolution @param cell_number     
        
        domain is periodic (if 2d or 3d)
        
        returns soil_model (RichardsWrapper(RichardsSP()))
    """
    s = RichardsWrapper(RichardsSP())  # water only
    s.initialize()
    s.createGrid(min_b, max_b, cell_number, periodic = False)  # [cm]  # periodicity is not needed, since coupling via s.pick(0., 0., z.)

    # BC
    if times is not None:
        s.setTopBC("atmospheric", 0.5, [times, net_inf])  # 0.5 is dummy value
    else:
        s.setTopBC("noFlux")
    s.setBotBC(bot_bc, bot_value)

    # Parameters
    s.setVGParameters([soil_])
    s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
    s.setParameter("Soil.SourceSlope", "1000")  # turns regularisation of the source term on
    if type == 2:
        s.setParameter("Component.MolarMass", "6.2e-2")  # nitrate 62,0049 g/mol
        s.setParameter("Component.LiquidDiffusionCoefficient", "1.7e-9")  # m2 s-1 # nitrate = 1700 um^2/sec

    s.setHomogeneousIC(initial_totalpotential, equilibrium = True)  # cm

    s.initializeProblem()
    s.setBotBC(bot_bc, bot_value)

    wilting_point = -15000
    s.setCriticalPressure(wilting_point)  # for boundary conditions constantFlow, constantFlowCyl, and atmospheric
    s.ddt = 1.e-5  # [day] initial Dumux time step

    return s


def set_initial_conditions(s, filename):
    """ TODO set constant total potential or from file """
    if  isinstance(filename, str):
       h = np.load(filename)
       s.setInitialConditionHead(h)  # cm
    else:
        s.setHomogeneousIC(-1000)


def set_bottom_bc_watertable(s, times, values):
    """ TODO set variable water table """
    pass


def plot_solution1D(s, depth = -200):
    """ for debugging e.g inital conditions"""
    h = s.getSolutionHead()
    plt.plot(h, np.linspace(-200., 0., h.shape[0]))
    plt.xlabel("soil matric potential [cm]")
    plt.ylabel("depth (cm)")
    plt.tight_layout()
    plt.show()

