import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src/python_modules/")
sys.path.append("../../build-cmake/cpp/python_binding/")

import plantbox as pb
import xylem_flux  # to use its rsml reader 

from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding) of cylindrcial model
from richards import RichardsWrapper  # Python part of cylindrcial model

import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool


class RhizoMappedSegments(pb.MappedSegments):
    """
        Adds 1-dimensional rhizosphere models to each root segment of a MappedSegments (or later MappedPlant)
    """
    
    # todo copy mapped segments constructors ...
        
    def __init__(self, rs, wilting_point, NC, logbase):
        """ @param rs is either a pb.MappedRootSystem, pb.MappedSegments, or a string containing a rsml filename"""
        ms = xylem_flux.XylemFluxPython.read_rsml(rs)
        super().__init__(ms.nodes, ms.nodeCTs, ms.segments, ms.radii, ms.subTypes) 
        self.mode = "none"        
        self.cyls = []        
        self.wilting_point = wilting_point
        self.NC = NC
        self.logbase = logbase
                
    def initialize(self, soil, x, mode="dumux"):
        """ calls the specific initializer according to @param mode (no parallelisation)  
        @param soil      van genuchten parameters as list        
        @param x         is the solution (or initial condition) of the soil model
        @param mode      currently "dumux", later additionally other solvers
        """
        self.mode = mode  # more precise RichardsCylFoam
        self.outer_radii = self.segOuterRadii()         
        self.seg_length = self.segLength()
        self.cyls = []
        if self.mode == "dumux":
            for i in range(0, len(self.segments)): 
                self.cyls.append(self.initialize_dumux_(i, soil, x))        

    def initialize_mp(self, soil, x, mode="dumux"):
        """ calls the specific initializer according to @param mode (using multiple threads)  
        @param soil      van genuchten parameters as list        
        @param x         is the solution (or initial condition) of the soil model
        @param mode      currently "dumux", later additionally other solvers
        """
        self.soil = soil
        self.x = x
        self.mode = mode  # more precise RichardsCylFoam
        self.outer_radii = self.segOuterRadii()  # single threaded
        self.seg_length = self.segLength()  # single threaded
        pool = Pool()  # defaults to number of available CPU's    
        cyls = pool.map(self.initialize_dumux2_, range(len(self.outer_radii)))

    def initialize_dumux2_(self, i):
        self.initialize(self.soil, x, mode) 

    def initialize_dumux_(self, i, soil, x): 
        """ see initialize_dumux_cyl(...) """
        a_in = self.radii[i]
        a_out = self.outer_radii[i]
        if a_in < a_out: 
            cpp_base = RichardsCylFoam()
            cyl = RichardsWrapper(cpp_base)
            cyl.initialize()                
            lb = self.logbase     
            points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), self.NC, base=lb)
            cyl.createGrid1d(points)                                          
            cyl.setHomogeneousIC(x[self.seg2cell[i]])  # cm pressure head                
            cyl.setVGParameters([soil])
            cyl.setOuterBC("fluxCyl", 0.)  # [cm/day]
            cyl.setInnerBC("fluxCyl", 0.)  # [cm/day]
            cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
            cyl.initializeProblem()
            cyl.setCriticalPressure(self.wilting_point)  # cm pressure head   
            return cyl             
        else:
            print("RhizoMappedSegments.initialize_dumux_cyl_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))  
            input()
            return []

    def get_inner_heads(self):
        """ """
        rsx = np.zeros((len(self.cyls),))
#         if self.mode == "python_fv":
#             for j, cyl in enumerate(self.cyls):  # for each segment
#                 rsx[j] = cyl.getInnerHead()  # [cm], needed for conductiv
#             return rsx                
        if self.mode == "dumux":
            for j, cyl in enumerate(self.cyls):  # for each segment
                rsx[j] = cyl.getInnerHead()  # [cm], needed for conductiv
            return rsx                
        print("RhizoMappedSegments.get_inner_heads: Warning, mode {:s} unknown".format(mode))

    def get_inner_fluxes(self):
        """ """
        fluxes = np.zeros((len(self.cyls),))
        if self.mode == "dumux":
            for j, cyl in enumerate(self.cyls):  # for each segment
                fluxes[j] = -float(cyl.getInnerFlux()) * (2 * np.pi * self.radii[j] * self.seg_length[j]) / self.radii[j]  # [cm/day] -> [cm3/day], ('/inner_radii' comes from cylindrical implementation)
            return fluxes                
        print("RhizoMappedSegments.get_inner_fluxes: Warning, mode {:s} unknown".format(mode))
       
    def solve_dumux(self, dt, proposed_inner_fluxes, proposed_outer_fluxes): 
        """ """
        for j, cyl in enumerate(self.cyls):  # run cylindrical models
            l = self.seg_length[j]
            cyl.setInnerFluxCyl(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l))  # [cm3/day] -> [cm /day]
            cyl.setOuterFluxCyl(proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l))  # [cm3/day] -> [cm /day]
            # cyl.ddt = 1.e-5  # [day] initial time step
            try:
                cyl.solve(dt)
            except:
                x = cyl.getDofCoordinates()
                y = cyl.getSolutionHead()
                plt.plot(x, y)
                plt.xlabel("x (cm)")
                plt.ylabel("pressure (cm)")
                print("in", proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l), "out", proposed_outer_fluxes[j] / (2 * np.pi * self.radii[j] * l))
                print("inner", self.radii[j], "outer", self.outer_radii[j], "l", l)
                plt.show()
                raise

    def get_water_volume(self):
        """ returns the water volume in the cylindrical models [cm3] """
        volumes = []
        for i, cyl in enumerate(self.cyls):  # run cylindrical models
            cyl_water = 0.
            cyl_water_content = cyl.getWaterContent()  # segment 0
            nodes = cyl.getPoints()
            for j, wc in enumerate(cyl_water_content):
                r1 = nodes[j]
                r2 = nodes[j + 1]
                cyl_water += np.pi * (r2 * r2 - r1 * r1) * self.seg_length[i] * wc
            volumes.append(cyl_water)
        return volumes        


class RhizoRootSystem(RhizoMappedSegments):
    """
        ??? Either RhizoMappedSegments, or MappedRootSystem, can we use both?
    
        Adds 1-dimensional rhizosphere models to each root segment of a MappedRootSystem (or later MappedPlant)
    """

    def __init__(self):
        
        # Todo which super to call?
        
        self.mode = "none"
        self.cyls = []
        
    def simulate_cyl_(self):
        pass
    
    def initialize_dumux_cyl(self, soil, x):
        """ x is the solution (or initial condition) of the soil model
        """
        self.mode = "dumux"  # more precise RichardsCylFoam     


def plot_transpiration(t, water_uptake, collar_flux, trans):
    """ potential and actual transpiration over time """
    fig, ax1 = plt.subplots()
    ax1.plot(t, trans(t), 'k', label="potential")  # potential transpiration
    ax1.plot(t, -np.array(water_uptake), 'g', label="actual")  # actual transpiration 
    ax1.plot(t, -np.array(collar_flux), 'r:', label="collar flux")  # actual transpiration
    ax1.set_xlabel("Time [d]")
    ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
    fig.legend()
    plt.show()


def plot_info(x_, water_collar_cell, water_cyl, collar_sx, min_sx, min_rx, min_rsx, water_uptake, water_domain): 
    """ 2x2 plot with additional information """    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)    
    ax1.set_title("Water amount")
    ax1.plot(x_, np.array(water_collar_cell), label="water cell")
    ax1.plot(x_, np.array(water_cyl), label="water cylindric")
    ax1.legend()
    ax1.set_xlabel("Time (days)")
    ax1.set_ylabel("(cm3)")    
    ax2.set_title("Pressure")
    ax2.plot(x_, np.array(collar_sx), label="soil at root collar")
    ax2.plot(x_, np.array(min_sx), label="min soil")
    ax2.plot(x_, np.array(min_rx), label="min xylem")
    ax2.plot(x_, np.array(min_rsx), label="min 1d at root surface")
    ax2.legend()
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Matric potential (cm)")    
    ax3.set_title("Water uptake")
    ax3.plot(x_, -np.array(water_uptake))
    ax3.set_xlabel("Time (days)")
    ax3.set_ylabel("Uptake (cm/day)")
    ax4.set_title("Water in domain")
    ax4.plot(x_, np.array(water_domain))
    ax4.set_xlabel("Time (days)")
    ax4.set_ylabel("cm3")
    plt.show()        
    
# print(len(rs.nodeCTs), len(rs.segments))
# ana2 = pb.SegmentAnalyser(r.rs.nodes, r.rs.segments, r.rs.nodeCTs[1:], r.rs.radii)
# types = np.array(r.rs.subTypes, dtype=np.float64)
# ana2.addData("subType", types)
# ana2.addData("age", r.get_ages())
# pd = vp.segs_to_polydata(ana2, 1., ["radius", "subType", "creationTime", "age"])
# vp.plot_roots(pd, "creationTime")

