import sys; sys.path.append("../modules/"); sys.path.append("../modules/fv/"); sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src/python_modules/")
sys.path.append("../../build-cmake/cpp/python_binding/")

import plantbox as pb
import xylem_flux  # to use its rsml reader 

from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding) of cylindrcial model
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial model (a single cylindrical model is not allowed to run in parallel)
from fv_grid import *
import fv_richards as rich  # local pure Python cylindrical models
import van_genuchten as vg

import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()


class RhizoMappedSegments(pb.MappedSegments):
    """
        Adds 1-dimensional rhizosphere models to each root segment of a MappedSegments (or later MappedPlant)
        
        do not use MPI with DUMUX rhizosphere models (i have no idea how to turn off mpi for 10 dof cylindrical models)
        
    """
    
    # todo copy mapped segments constructors ...
        
    def __init__(self, rs, wilting_point, NC, logbase, mode="dumux"):
        """ @param rs is either a pb.MappedRootSystem, pb.MappedSegments, or a string containing a rsml filename"""
        ms = xylem_flux.XylemFluxPython.read_rsml(rs)
        super().__init__(ms.nodes, ms.nodeCTs, ms.segments, ms.radii, ms.subTypes) 
        self.cyls = []        
        self.wilting_point = wilting_point
        self.NC = NC
        self.logbase = logbase
        self.mode = mode  # more precise RichardsCylFoam, mode="dumux"
        self.last_dt = 0.        
                
    def initialize(self, soil, x, eidx=None):
        """ calls the specific initializer according to @param mode (no parallelisation)  
        @param soil     van genuchten parameters as list        
        @param x        is the solution (or initial condition) of the soil model
        @param mode     currently "dumux", later additionally other solvers
        @þaram eidx     in case of mpi the cylinder indices per process
        """        
        if eidx is None: 
            eidx = np.array(range(0, len(self.radii)), np.int64)  # segment indices for process
        self.eidx = np.array(eidx, np.int64)
        self.outer_radii = np.array(self.segOuterRadii())  # not parallel yet        
        self.seg_length = self.segLength()
        self.soil = soil
        self.vg_soil = vg.Parameters(soil) 
        self.cyls = []
        if self.mode == "dumux":
            for i in eidx: 
                self.cyls.append(self.initialize_dumux_(i, x[self.seg2cell[i]]))
        elif self.mode == "python": 
            for i in eidx: 
                x0 = x[self.seg2cell[i]]
                self.cyls.append(self.initialize_python_(i, x0))
        else:
            raise "RhizoMappedSegments.initialize: unknown solver {}".format(mode)

    def initialize_dumux_(self, i, x): 
        """ Dumux RichardsCylFoam solver"""
        a_in = self.radii[i]
        a_out = self.outer_radii[i]
        if a_in < a_out: 
            cyl = RichardsNoMPIWrapper(RichardsCylFoam())  # only works for RichardsCylFoam compiled without MPI
            cyl.initialize()                
            lb = self.logbase     
            points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), self.NC, base=lb)
            cyl.createGrid1d(points)                                          
            cyl.setHomogeneousIC(x)  # cm pressure head                
            cyl.setVGParameters([self.soil])
            cyl.setOuterBC("fluxCyl", 0.)  # [cm/day]
            cyl.setInnerBC("fluxCyl", 0.)  # [cm/day]
            cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
            cyl.initializeProblem()
            cyl.setCriticalPressure(self.wilting_point)  # cm pressure head   
            return cyl             
        else:
            print("RhizoMappedSegments.initialize_dumux_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))  
            return []

    def initialize_python_(self, i, x): 
        """ Python home grown richards fv solver"""
        a_in = self.radii[i]
        a_out = self.outer_radii[i]
        if a_in < a_out:
            lb = self.logbase     
            points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), self.NC, base=lb)
            grid = FVGrid1Dcyl(points)
            cyl = rich.FVRichards1D(grid, self.soil)
            cyl.x0 = np.ones((self.NC - 1,)) * x            
            cyl.solver_initialize()  # solve() would call that, but we use solve_single_step                           
            return cyl
        else:
            print("RhizoMappedSegments.initialize_python_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))  
            return []

    def set_xylem_flux(self, rs):
        """ we need XylemFlux for the kr and kx call back functions (for python)"""
        self.rs = rs

    def get_inner_heads(self):
        """ matric potential at the root surface interface [cm]"""
        rsx = np.zeros((len(self.cyls),))           
        if self.mode == "dumux":
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                rsx[i] = cyl.getInnerHead()  # [cm], needed for conductiv                          
        elif self.mode == "python":
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                rsx[i] = cyl.getInnerHead()  # [cm], needed for conductiv
        else:
            print("RhizoMappedSegments.get_inner_heads: Warning, mode {:s} unknown".format(self.mode))
        return self._map(self._flat0(comm.gather(rsx, root=0)))  # gathers and maps correctly 

    def get_inner_fluxes(self):
        """ fluxes [cm3/day] at the root surface, i.e. inner boundary  """
        fluxes = np.zeros((len(self.cyls),))
        if self.mode == "dumux":
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                fluxes[i] = -float(cyl.getInnerFlux()) * (2 * np.pi * self.radii[j] * self.seg_length[j]) / self.radii[j]  # [cm/day] -> [cm3/day], ('/inner_radii' comes from cylindrical implementation)            
        elif self.mode == "python":
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                fluxes[i] = cyl.getInnerFlux() * (2 * np.pi * self.radii[j] * self.seg_length[j]) / self.last_dt  # divide by dt is correct here! getInnerFlux only gives the source in cm3/cm2
        else:
            print("RhizoMappedSegments.get_inner_fluxes: Warning, mode {:s} unknown".format(self.mode))        
        return self._map(self._flat0(comm.gather(fluxes, root=0)))  # gathers and maps correctly 
       
    def solve(self, dt, *argv): 
        """ set bc for cyl and solves it """
        self.last_dt = dt
        if self.mode == "dumux":
            proposed_inner_fluxes = argv[0] 
            proposed_outer_fluxes = argv[1]
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                l = self.seg_length[j]
                cyl.setInnerFluxCyl(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l))  # [cm3/day] -> [cm /day]
                cyl.setOuterFluxCyl(proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l))  # [cm3/day] -> [cm /day]
                try:
                    cyl.solve(dt)
                except:
                    str = "RhizoMappedSegments.solve: dumux exception with boundaries in flow {:g} cm3/day, out flow {:g} cm3/day, segment radii [{:g}-{:g}] cm"                
                    str = str.format(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l), proposed_outer_fluxes[j] / (2 * np.pi * self.radii[j] * l), self.radii[j], self.outer_radii[j])
                    raise str                
        elif self.mode == "python": 
            rx = argv[0]        
            proposed_outer_fluxes = argv[1]    
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                l = self.seg_length[j]                
                seg_ct = self.nodeCTs[self.segments[j].y]
                age = 1.  # self.seg_ages[j] + t  # age  = (maxCT -nodeCT[j]) + t = (maxCT+t) - nodeCT[j] = rs_sim_time - nodeCT[j]
                type_ = self.subTypes[j]
                x0 = rx[self.segments[j].x]
                x1 = rx[self.segments[j].y]
                #         rx_approx = 0.5 * (x0 + x1)
                #         cyl.bc[(0, 0)] = ("rootsystem", [rx_approx, self.kr_f(age, type_)])
                cyl.bc[(0, 0)] = ("rootsystem_exact", [x0, x1, self.rs.kr_f(age, type_), self.rs.kx_f(age, type_), self.radii[j], l])                
                ndof = self.NC - 1
                dx_outer = cyl.grid.nodes[ndof] - cyl.grid.center(ndof - 1)
                q_outer = proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l)
                cyl.bc[(ndof - 1, 1)] = ("flux_in_out", [q_outer , self.wilting_point, dx_outer])
                cyl.solve([dt], dt, False)
                                        
        else:
            raise "RhizoMappedSegments.initialize: unknown solver {}".format(self.mode)
        
    def get_water_volume(self):
        """ returns the water volume of the cylindrical models [cm3] """
        volumes = np.zeros((len(self.cyls),))
        if self.mode == "dumux":
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                cyl_water = 0.
                cyl_water_content = cyl.getWaterContent()  # segment 0
                nodes = cyl.getPoints()
                for k, wc in enumerate(cyl_water_content):
                    r1 = nodes[k]
                    r2 = nodes[k + 1]
                    cyl_water += np.pi * (r2 * r2 - r1 * r1) * self.seg_length[j] * wc
                volumes[i] = cyl_water
        elif self.mode == "python":
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                cyl_water = 0.
                cyl_water_content = cyl.getWaterContent()  # segment 0
                for k, wc in enumerate(cyl_water_content):
                    r1 = cyl.grid.nodes[k]
                    r2 = cyl.grid.nodes[k + 1]
                    cyl_water += np.pi * (r2 * r2 - r1 * r1) * self.seg_length[j] * wc            
                volumes[i] = cyl_water
        return self._map(self._flat0(comm.gather(volumes, root=0)))  # gathers and maps correctly 
    
    def _map(self, x):
        """Converts @param x to a numpy array and maps it to the right indices                 """
        indices = self._flat0(comm.gather(self.eidx, root=0))  # segment indices
        if indices:  # only for rank 0 it is not empty
            assert len(indices) == len(x), "RhizoMappedSegments._map: indices and values have different length"            
            p = np.zeros((len(x),), dtype=np.float64)
            for i in range(0, len(indices)):  #
                p[indices[i]] = x[i]
            return p
        else:
            return 0

    def _flat0(self, xx):
        """flattens the gathered list in rank 0, empty list for other ranks """
        if rank == 0:
            return [item for sublist in xx for item in sublist]
        else:
            return []
    

def plot_transpiration(t, soil_uptake, root_uptake, potential_trans):
    """ plots potential and actual transpiration over time 
    
    depending on discretisation soil_uptake and root_uptake might differ  
    
    @param t                  times [day]
    @param soil_uptake        actual transpiration [cm3/day] of soil 
    @param root_uptake        actual transpiration [cm3/day] according to root model 
    @param potential_trans    function in t stating the potential transpiration [cm3/day]
    """
    fig, ax1 = plt.subplots()
    ax1.plot(t, potential_trans(t), 'k', label="potential transpiration")  # potential transpiration
    ax1.plot(t, -np.array(soil_uptake), 'g', label="soil uptake")  # actual transpiration  according to soil model
    ax1.plot(t, -np.array(root_uptake), 'r:', label="root system uptake")  # actual transpiration according root model
    ax1.set_xlabel("Time [d]")
    ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
    ax2 = ax1.twinx()
    dt = np.diff(t)
    so = np.array(soil_uptake)
    ax2.plot(t[1:], np.cumsum(-np.multiply(so[:-1], dt)), 'c--')  # cumulative transpiration (neumann)
    ax2.set_ylabel("Cumulative soil uptake $[cm^3]$")    
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

