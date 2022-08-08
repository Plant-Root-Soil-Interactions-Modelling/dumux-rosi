""" functions to simplify setup of the scenarios """

import sys; sys.path.append("../../build-cmake/cpp/python_binding/"); sys.path.append("../modules/");
sys.path.append("../../../CPlantBox/src/python_modules"); sys.path.append("../../../CPlantBox/");

import numpy as np

import plantbox as pb  # CPlantBox
import van_genuchten as vg
from xylem_flux import XylemFluxPython

from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model



def create_soil_model(soil_, min_b , max_b , cell_number, p_top, p_bot):
    """
        Creates a soil domain from @param min_b to @param max_b with resolution @param cell_number
        soil type is fixed and homogeneous 
        domain is periodic (if 2d or 3d)
        initial potentials are linear from @param p_top to @param p_bot
        
        returns soil_model (RichardsWrapper(RichardsSP())) and soil parameter (vg.Parameters)
    """  
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
    s.setBotBC("freeDrainage")
    s.setVGParameters([soil_])
    s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
    s.initializeProblem()
    wilting_point = -15000
    s.setCriticalPressure(wilting_point)  # for boundary conditions constantFlow, constantFlowCyl, and atmospheric
    s.ddt = 1.e-5  # [day] initial Dumux time step

    return s, soil


def init_conductivities_const(r):
    """ Hydraulic conductivities - for Jans scenarios, but constant """
    # Scenario 1
    kr_const = 1.8e-4  # [1/day]
    kx_const = 0.1  # [cm3/day]
    r.setKr([kr_const])
    r.setKx([kx_const])


def create_mapped_singleroot(min_b , max_b , cell_number, soil_model, ns = 100, l = 50 , a = 0.05):
    """ creates a single root mapped to a soil with @param ns segments, length l, and radius a """
    r = create_singleroot(ns, l, a)
    r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), cut = False)
    picker = lambda x, y, z: soil_model.pick([x, y, z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
    r.rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
    init_conductivities_const(r)
    return r


def create_singleroot(ns = 100, l = 50 , a = 0.05):
    """ creates a single root with @param ns segments, length l, and radius a """
    radii = np.array([a] * ns)
    nodes = [pb.Vector3d(0, 0, 0)]
    segs = []
    dx = l / ns
    z_ = np.linspace(-dx, -l , ns)
    for i in range(0, ns):
        nodes.append(pb.Vector3d(0, 0, z_[i]))
        segs.append(pb.Vector2i(i, i + 1))
    rs = pb.MappedSegments(nodes, segs, radii)
    return XylemFluxPython(rs)


def create_mapped_rootsystem(soil_model):

    return r


def write_files(file_name, psi_x, psi_i, sink, times, trans, psi_s):
    """  saves numpy arrays ass npy files """
    np.save('results/psix_' + file_name, np.array(psi_x))  # xylem pressure head per segment [cm]
    np.save('results/psiinterface_' + file_name, np.array(psi_i))  # pressure head at interface per segment [cm]
    np.save('results/sink_' + file_name, -np.array(sink))  # sink per segment [cm3/day]
    np.save('results/transpiration_' + file_name, np.vstack((times, -np.array(trans))))  # time [day], transpiration [cm3/day]
    np.save('results/soil_' + file_name, np.array(psi_s))  # soil potential per cell [cm]


if __name__ == '__main__':

        s, soil = create_soil_model([-1, -1, -150.], [1, 1, 0.], [1, 1, 55], -310, -200)

        print()
        print(s)
        print(soil)
        
        """ TODO: tests would be nice, or a minimal example setup ... """
