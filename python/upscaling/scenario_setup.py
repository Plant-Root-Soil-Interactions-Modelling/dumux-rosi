""" some functions to simplify setup """

import sys; sys.path.append("../../build-cmake/cpp/python_binding/"); sys.path.append("../modules/");
sys.path.append("../../../CPlantBox/src/python_modules"); sys.path.append("../../../CPlantBox/");

from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
import van_genuchten as vg
import aggregated_rs as agg


def create_soil_model(min_b , max_b , cell_number, p_top, p_bot,):
    """
        Creates a soil domain from @param min_b and max_b with resolution @param cell_number
        soil type is fixed  
        soil is periodic (if 2d or 3d)
        initial potentials are linear from @param p_top to @param p_bot
        
        returns soil_model and soil parameter
    """
    alpha = 0.0383  # (cm-1) soil
    n = 1.3774
    Ks = 60.  # (cm d-1)
    soil_ = [0.025, 0.403, alpha, n, Ks]  # corresponding soil root interface potentials table is  "../table_jan_comp"
    soil = vg.Parameters(soil_)
    vg.create_mfp_lookup(soil, -1.e5, 1000)
    if (cell_number[0] == 1 and cell_number[1] == 1):  # 1D
        periodic = False
    else:  # 2D or 3D
        periodic = True
    s = RichardsWrapper(RichardsSP())
    s.initialize()
    s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
    s.setLinearIC(p_top, p_bot)  # cm pressure head, equilibrium
    s.setTopBC("noFlux")
    s.setBotBC("noFlux")
    s.setVGParameters([soil_])
    s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
    s.initializeProblem()
    wilting_point = -15000
    s.setCriticalPressure(wilting_point)  # for boundary conditions constantFlow, constantFlowCyl, and atmospheric
    s.ddt = 1.e-5  # [day] initial Dumux time step

    return s, soil


def create_mapped_rootsystem(soil_model):

    return r


def create_mapped_singleroot(soil_model):

    return r


if __name__ == '__main__':

        s, soil = create_soil_model([1, 1, 55], -310, -200)

        print()
        print(s)
        print(soil)
