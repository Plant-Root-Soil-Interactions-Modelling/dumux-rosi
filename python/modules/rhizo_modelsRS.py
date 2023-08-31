# import sys; sys.path.append("../modules/"); sys.path.append("../modules/fv/"); sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src");
# sys.path.append("../../build-cmake/cpp/python_binding/")

import plantbox as pb
import functional.xylem_flux as xylem_flux
import sys
from functional.xylem_flux import XylemFluxPython
from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding) of cylindrcial model
from rosi_richards2c_cyl import Richards2CCylFoam  # C++ part (Dumux binding)
from rosi_richardsnc_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)
from rosi_richards10c_cyl import Richards10CCylFoam # C++ part (Dumux binding)
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial model (a single cylindrical model is not allowed to run in parallel)
from fv.fv_grid import *
import fv.fv_richards as rich  # local pure Python cylindrical models
import functional.van_genuchten as vg

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import norm
from scipy.interpolate import griddata as gd
from mpl_toolkits.mplot3d import Axes3D
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import multiprocessing
from multiprocessing import Process, active_children
import psutil
from threading import Thread 
from air_modelsPlant import AirSegment
from scipy.interpolate import PchipInterpolator,  CubicSpline

class RhizoMappedSegments(pb.MappedRootSystem):#XylemFluxPython):#
    """
        Adds 1-dimensional rhizospere models to each root segment of a MappedSegments (or later MappedPlant)    
        
        modes:        
        "dumux"              inner boundary is FluxCyl,                    <- currently best
        "dumux_exact"        inner boundary is rootSystemExact (same as "dumux", but slower, for experiments...) 
        "python"             inner boundary is flux_in,
        "python_exact"       inner boundary is rootsystem_exact (slower, for experiments...)
    """

    # TODO copy mapped segments constructors (!)...

    def __init__(self, wilting_point, NC, logbase, mode, soil, recreateComsol_):
        """ @param file_name is either a pb.MappedRootSystem, pb.MappedSegments, or a string containing a rsml filename"""
        super().__init__()
        self.wilting_point = wilting_point
        self.NC = NC #dof +1
        self.logbase = logbase
        self.mode = mode  # more precise RichardsCylFoam, mode="dumux"
        #print(self.mode)

        # changes with growing plant (to update)
        self.cyls = []
        self.outer_radii = None
        self.seg_length = np.array([]) 
        self.eidx = []
        self.dx2 = []
        self.soilModel = soil
        self.recreateComsol = recreateComsol_

        # additional variables
        self.last_dt = 0.
        
        
    
    
    def initializeRhizo(self, soil, x, eidx = None, cc = None):
        """ calls the specific initializer according to @param mode (no parallelisation)  
        @param soil     van genuchten parameters as list        
        @param x        is the solution (or initial condition) of the soil model
        @þaram eidx     in case of mpi the cylinder indices per process
        """
        self.soil = soil
        self.vg_soil = vg.Parameters(soil)
        vg.create_mfp_lookup(self.vg_soil, -1.e5, 1000)
        self.cyls = []
        self.dx2 = []
        self.eidx = np.array([], dtype=np.int64)
        
        
        self.solidDensity = 2700 # [kg/m^3 solid]
        self.solidMolarMass = 60.08e-3 # [kg/mol] 
            
        # [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
        self.solidMolDensity = self.solidDensity/self.solidMolarMass
        # [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
        self.bulkDensity_m3 = self.solidMolDensity*(1.- self.soil[1]) #porosity == theta_s
        self.molarMassWat = 18. # [g/mol]
        self.densityWat = 1. #[g/cm3]
        # [mol/cm3] = [g/cm3] /  [g/mol] 
        self.molarDensityWat =  self.densityWat /self.molarMassWat # [mol/cm3]    
        self.numFluidComp = 2
        self.numComp = 8
        self.update(x, eidx, cc)
        
    def update(self, x, newEidx = None, cc = None):
        """ calls the specific initializer according to @param mode (no parallelisation)  
        @param soil     van genuchten parameters as list        
        @param x        is the solution (or initial condition) of the soil model
        @þaram newEidx  in case of mpi the cylinder indices per process
        """
        if len(x.shape)!=1:
            raise Exception
        if newEidx is None:
            newEidx = np.array(range(len(self.eidx), len(self.radii)), np.int64)  # segment indices for process
        self.eidx = np.concatenate((self.eidx,np.array(newEidx, dtype = np.int64)), dtype = np.int64)
        
        if self.recreateComsol:
            self.outer_radii = np.full(len(self.radii), 0.6)#np.array(self.segOuterRadii())  # in the future, segment radius might vary with time. TODO: how to handle that?
        else:
            self.outer_radii = np.array(self.segOuterRadii())  # in the future, segment radius might vary with time. TODO: how to handle that?
            
        self.seg_length = self.segLength()#node might have shifte: new length for pre-existing segments
        for i, cyl in enumerate(self.cyls):
            self.updateOld(i, cyl, x, cc)
        for i in newEidx:#only initialize the new eidx
            self.initialize_(i,x,cc)
    
    def updateOld(self, i,  cyl, x, cc = None):
        """ update distribution of cylinder if it s volume has changed
        """
        oldPoints = np.array(cyl.getPoints()).flatten()
        lb = self.logbase            
        a_in = self.radii[i] 
        a_out = self.outer_radii[i] * 2.
        points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), self.NC, base = lb)
        if (len(points) != len(oldPoints)) or (max(abs(points - oldPoints)) > 1e-2):#new vertex distribution
            cellsOld = np.array(cyl.getCellCenters()).flatten()
            oldHead = np.array(cyl.getSolutionHead()).flatten()
            cAllOld = [np.array(cyl.getSolution_(i+1)) for i in range(self.numComp)]
            # add a point far away which is the mean value of the voxel
            
            # # because I don t want to use the linear interpolation of dumux:
            # # but nor necessary => outer_radii will only ever decrease
            
            # chip = PchipInterpolator(cellsOld, cssOld)#for interpolation
            # spl =  CubicSpline(cellsOld, cssOld, bc_type='not-a-knot')
            # def getValUpdate(newX):
                # if (newX >= min(cellsOld)) and (newX <= max(cellsOld)): #interpolation
                    # return chip(newX)
                # else:#extrapolation
                    # return spl(newX)
            # print("new points",points)
            # print("old points",oldPoints)
            # print(cssOld)
            # cssNew = np.array([getValUpdate(xx) for xx in points])
            
            print("cAllOld",cAllOld, cyl.getCellCenters())
            
            self.cyls[i] = self.initialize_dumux_nc_( i, 
                                                 x = oldHead, 
                                                 c1 = cAllOld[0], 
                                                 c2 = cAllOld[1], 
                                                 c3 = cAllOld[2], 
                                                 c4 = cAllOld[3], 
                                                 c5 = cAllOld[4], 
                                                 c6 = cAllOld[5], 
                                                 c7 = cAllOld[6], 
                                                 c8 = cAllOld[7], 
                                                 oldCells = cellsOld)
        
    def initialize_(self,i,x,cc ):
        if self.seg2cell[i] > 0:
            if self.mode == "dumux":
                self.cyls.append(self.initialize_dumux_(i, x[self.seg2cell[i]], False, False))
            elif self.mode == "dumux_exact":
                self.cyls.append(self.initialize_dumux_(i, x[self.seg2cell[i]], True, False))
            elif self.mode == "dumux_dirichlet":
                self.cyls.append(self.initialize_dumux_(i, x[self.seg2cell[i]], False, True))
            elif self.mode == "dumux_dirichlet_2c":
                self.cyls.append(self.initialize_dumux_nc_(i, x[self.seg2cell[i]], cc[self.seg2cell[i]]))
            elif (self.mode == "dumux_dirichlet_nc") or(self.mode == "dumux_dirichlet_10c") or (self.mode == "dumux_10c") :
                self.cyls.append(self.initialize_dumux_nc_(i, x[self.seg2cell[i]], cc[self.seg2cell[i]]))
            elif self.mode == "dumux_nc":
                self.cyls.append(self.initialize_dumux_(i, x[self.seg2cell[i]], False, True))
            elif self.mode == "python" or self.mode == "python_exact":
                x0 = x[self.seg2cell[i]]
                self.cyls.append(self.initialize_python_(i, x0))
            else:
                print("let s see", self.mode,"dumux_dirichlet_2c",self.mode == "dumux_dirichlet_2c" )
                raise Exception("RhizoMappedSegments.initialize: unknown solver {}".format(self.mode))
        else:
            a_in = self.radii[i]#get Perimeter instead? not used for now anyway
            self.cyls.append(AirSegment(a_in)) #psi_air computed by the photosynthesis module.
            #print("seg",i,"is out of the soil domain and has no rhizosphere")        
    

    def initialize_dumux_nc_(self, i, x, c1 = 0.1 / 18, 
                                         c2 = 10/ 18, 
                                         c3 =  0.011 * 1e6/ (2700/ 60.08e-3* (1. - 0.43)), 
                                         c4 = 0.05 * 1e6/(2700/ 60.08e-3* (1. - 0.43)), 
                                         c5 = 0.011 * 1e6/(2700/ 60.08e-3* (1. - 0.43)), 
                                         c6 =  0.05 * 1e6/(2700/ 60.08e-3* (1. - 0.43)), 
                                         c7 = 0./(2700/ 60.08e-3* (1. - 0.43)), 
                                         c8 = 0./(2700/ 60.08e-3* (1. - 0.43)), 
                                         oldCells = []):
        cAll = [c1, c2 , c3 , c4 , c5 , c6 , c7 , c8 ]
        a_in = self.radii[i]
        a_out = self.outer_radii[i]
        if a_in < a_out:
            if self.mode == "dumux_dirichlet_2c":
                cyl = RichardsNoMPIWrapper(Richards2CCylFoam())  # only works for RichardsCylFoam compiled without MPI
            if self.mode == "dumux_dirichlet_nc":
                cyl = RichardsNoMPIWrapper(RichardsNCCylFoam())  # only works for RichardsCylFoam compiled without MPI
            if self.mode == "dumux_dirichlet_10c":
                cyl = RichardsNoMPIWrapper(Richards10CCylFoam())  # only works for RichardsCylFoam compiled without MPI
            if self.mode == "dumux_10c":
                cyl = RichardsNoMPIWrapper(Richards10CCylFoam())  # only works for RichardsCylFoam compiled without MPI
            cyl.initialize(verbose = False)
            cyl.setVGParameters([self.soil])
            lb = self.logbase
            
            if self.recreateComsol:
                nCells = self.NC
                cyl.createGrid([0.02], [0.6], [nCells])
            else:
                points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), self.NC, base = lb)
                cyl.createGrid1d(points)
                if len(self.dx2) > i:
                    self.dx2[i] = 0.5 * (points[1] - points[0]) #when updating
                else:
                    self.dx2.append(0.5 * (points[1] - points[0]))#what is this?
            # print(a_in,a_out,lb)
            # print( a_in)
            if i==1:
                cyl.setParameter("Problem.verbose", "0")
            else:
                cyl.setParameter("Problem.verbose", "-1")
                
            cyl.setParameter( "Soil.Grid.Cells",str( self.NC)) # -1 to go from vertices to cell
            if self.recreateComsol:
                cyl.setHomogeneousIC(-100.)#x)  # cm pressure head
            else:
                cyl.setHomogeneousIC(x)  # cm pressure head
                cyl.setParameter("Soil.IC.P", cyl.dumux_str(x))
            # cyl.setICZ_solute(c)  # [kg/m2] 
            
            #default: no flux
            cyl.setInnerBC("fluxCyl", 0.)  # [cm/day] #Y 0 pressure?
            #cyl.setInnerBC_solute("fluxCyl", 0.)  # [kg/m2], that s fair
            cyl.setOuterBC("fluxCyl", 0.)
            #cyl.setOuterBC_solute("fluxCyl", 0.)
            
            Ds = 1e-8 # m^2/s
            Dl = 1e-9
            
            cyl.setParameter( "Soil.MolarMass", str(self.solidMolarMass))
            cyl.setParameter( "Soil.solidDensity", str(self.solidDensity))
            
            cyl.setParameter("Soil.betaC", str(0.001 ))
            cyl.setParameter("Soil.betaO", str(0.1 ))
            cyl.setParameter("Soil.C_S_W_thresC", str(0.1 )) #mol/cm3
            cyl.setParameter("Soil.C_S_W_thresO", str(0.05 )) #mol/cm3
            cyl.setParameter("Soil.k_decay", str(0.2 ))
            cyl.setParameter("Soil.k_decay2", str(0.6 ))
            #cyl.setParameter("Soil.k_decay3", str(1 ))
            cyl.setParameter("Soil.k_DC", str(1  )) # 1/d
            cyl.setParameter("Soil.k_DO", str(1  )) # 1/d
            #cyl.setParameter("Soil.k_growthBis", str(1 )) #bool
            cyl.setParameter("Soil.k_growthC", str(0.2 ))
            cyl.setParameter("Soil.k_growthO", str(0.2 ))
            cyl.setParameter("Soil.K_L", str(8.3))#[mol/cm3]
            cyl.setParameter("Soil.k_phi", str(0.1 ))
            cyl.setParameter("Soil.k_RC", str(0.1 ))
            cyl.setParameter("Soil.k_RO", str(0.1 ))

            k_SC = 1
            k_SO = 10
            cyl.setParameter("Soil.k_SC", str(k_SC )) #cm^3/mol/d
            #cyl.setParameter("Soil.k_SCBis", str(k_SC )) #cm^3/mol/d
            cyl.setParameter("Soil.k_SO", str(k_SO )) #cm^3/mol/d
            #cyl.setParameter("Soil.k_SOBis", str(k_SO )) #cm^3/mol/d

            m_maxC = 0.1 
            m_maxO = 0.02 
            cyl.setParameter("Soil.m_maxC", str(m_maxC  ))# 1/d
            cyl.setParameter("Soil.m_maxO", str(m_maxO  ))# 1/d
            cyl.setParameter("Soil.micro_maxC", str(2 ))# 1/d
            cyl.setParameter("Soil.micro_maxO", str(0.01 ))# 1/d
            cyl.setParameter("Soil.v_maxL", str(1.5 ))#[d-1]

            k_sorp = 0.2*100
            cyl.setParameter("Soil.k_sorp", str(k_sorp)) # mol / cm3
            f_sorp = 0.9
            cyl.setParameter("Soil.f_sorp", str(f_sorp)) #[-]
            CSSmax = 1e-4*10000
            cyl.setParameter("Soil.CSSmax", str(CSSmax)) #[mol/cm3 scv]
            cyl.setParameter("Soil.alpha", str(0.1)) #[1/d]


            cyl.setParameter( "Soil.BC.Bot.C1Type", str(3)) #put free flow later
            cyl.setParameter( "Soil.BC.Top.C1Type", str(3))
            cyl.setParameter( "Soil.BC.Bot.C1Value", str(0)) 
            cyl.setParameter( "Soil.BC.Top.C1Value", str(0 )) 

            cyl.setParameter("1.Component.LiquidDiffusionCoefficient", str(Ds)) #m^2/s

            cyl.setParameter( "Soil.BC.Bot.C2Type", str(3))
            cyl.setParameter( "Soil.BC.Top.C2Type", str(3))
            cyl.setParameter( "Soil.BC.Bot.C2Value", str(0)) 
            cyl.setParameter( "Soil.BC.Top.C2Value", str(0 )) 
            cyl.setParameter("2.Component.LiquidDiffusionCoefficient", str(Dl)) #m^2/s

            # cyl.decay = 0. #1.e-5
            
            for j in range(self.numFluidComp + 1, self.numComp+1):
                cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Type", str(3))
                cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Type", str(3))
                cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Value", str(0)) 
                cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Value", str(0 )) 
            print("cAll",cAll)                
            for j in range( 1, self.numComp+1):       
                # print("cAll[j-1]", cAll[j-1],cyl.dumux_str(cAll[j-1])  )
                cyl.setParameter("Soil.IC.C"+str(j), cyl.dumux_str(cAll[j-1]) ) 
            # C_S =cyl.dumux_str(c1) #mol/cm3 wat
            # cyl.setParameter("Soil.IC.C1", str(C_S) )  #mol/cm3 / mol/cm3 = mol/mol                 
            # C_L = cyl.dumux_str(c2) #mol/cm3 wat
            # cyl.setParameter("Soil.IC.C2", str(C_L) )  #mol/cm3 / mol/cm3 = mol/mol
            # COa = cyl.dumux_str(c3 ) #mol C / m3 space
            # cyl.setParameter("Soil.IC.C3",str(COa)) #mol C / mol Soil 
            # #cyl.setParameter("Soil.IC.C3", str(0.009127163)) #[mol/mol soil] == 233.8 [mol/m3 bulk soil]
            # COd =cyl.dumux_str(c4)#mol C / m3 space
            # cyl.setParameter("Soil.IC.C4", str(COd))
            # CCa = cyl.dumux_str(c5)#mol C / m3 space
            # cyl.setParameter("Soil.IC.C5", str(CCa)) 
            # CCd = cyl.dumux_str(c6)  #mol C / m3 space
            # cyl.setParameter("Soil.IC.C6", str(CCd))
            # CSS2_init =cyl.dumux_str( c7 )#CSSmax*1e6 * (C_S/(C_S+ k_sorp)) * (1 - f_sorp) #mol C/ m3 scv

            # cyl.setParameter("Soil.IC.C7", str(CSS2_init))#CSS2_init/bulkDensity_m3))#[1:(len(str(CSS2_init/bulkDensity_m3))-1)] ) #mol C / mol scv
            # cyl.setParameter("Soil.IC.C8", str(cyl.dumux_str(c8) ))
            
            
            if len(oldCells) > 0:#in case we update a cylinder, to allow dumux to do interpolation
                oldCellsStr = cyl.dumux_str(oldCells/100)
                cyl.setParameter("Soil.IC.Z",oldCellsStr)
                if len(oldCells)!= len(x):
                    print("oldCells, x",oldCells, x, len(oldCells), len(x))
                    raise Exception
                for j in range( 1, self.numComp+1):
                    cyl.setParameter("Soil.IC.C"+str(j)+"Z",oldCellsStr)
                    if len(oldCells)!= len( cAll[j-1]):
                        print("oldCells,  cAll[j-1]",oldCells,  cAll[j-1], 
                                len(oldCells), len(cAll[j-1]), j)
                        raise Exception
            
            cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
            cyl.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
                        
            cyl.setParameter("Problem.EnableGravity", "false")
            cyl.setParameter("Problem.verbose", "1")
            cyl.setParameter("Flux.UpwindWeight", "0.5")
            cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
            cyl.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
            cyl.setParameter("Newton.EnableChop", "true")
            cyl.setParameter("Newton.EnableResidualCriterion", "true")
            cyl.setParameter("Newton.EnableShiftCriterion", "true")
            cyl.setParameter("Newton.MaxAbsoluteResidual", "1e-10")

            cyl.setParameter("Newton.MaxRelativeShift", "1e-10")

            cyl.setParameter("Newton.MaxSteps", "30")
            cyl.setParameter("Newton.ResidualReduction", "1e-10")
            cyl.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
            cyl.setParameter("Newton.TargetSteps", "10")
            cyl.setParameter("Newton.UseLineSearch", "false")
            cyl.setParameter("Newton.EnablePartialReassembly", "true")
            cyl.setParameter("Grid.Overlap", "0")  #no effec5

            cyl.initializeProblem()
            cyl.setCriticalPressure(self.wilting_point)  # cm pressure head
            
            return cyl
        else:
            print("RhizoMappedSegments.initialize_dumux_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))
            return []

    def initialize_dumux_(self, i, x, exact, dirichlet = False):
        """ Dumux RichardsCylFoam solver"""
        raise Exception
        a_in = self.radii[i]
        a_out = self.outer_radii[i]
        if a_in < a_out:
            cyl = RichardsNoMPIWrapper(RichardsCylFoam())  # only works for RichardsCylFoam compiled without MPI
            cyl.initialize()
            lb = self.logbase
            points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), self.NC, base = lb)
            print(points)
            cyl.createGrid1d(points)
            self.dx2.append(0.5 * (points[1] - points[0]))
            cyl.setHomogeneousIC(x)  # cm pressure head
            cyl.setVGParameters([self.soil])
            cyl.setOuterBC("fluxCyl", 0.)  # [cm/day]
            if exact:
                cyl.setInnerBC("rootSystemExact")  # parameters are passed with cyl.setRootSystemBC (defined in richards_cyl.hh)
            else:
                if dirichlet:
                    cyl.setInnerBC("pressure", 0.)  # [cm/day]
                else:
                    cyl.setInnerBC("fluxCyl", 0.)  # [cm/day]
            cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
            cyl.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
            cyl.initializeProblem()
            cyl.setCriticalPressure(self.wilting_point)  # cm pressure head
            return cyl
        else:
            raise Exception("RhizoMappedSegments.initialize_dumux_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))
            return []

    def initialize_python_(self, i, x):
        """ Python home grown richards fv solver"""
        a_in = self.radii[i]
        a_out = self.outer_radii[i]
        if a_in < a_out:
            lb = self.logbase
            points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), self.NC, base = lb)
            grid = FVGrid1Dcyl(points)
            cyl = rich.FVRichards1D(grid, self.soil)
            cyl.x0 = np.ones((self.NC - 1,)) * x
            return cyl
        else:
            raise Exception("RhizoMappedSegments.initialize_python_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))
            return []

    def set_phloem_flux(self, rs):
        """ we need XylemFlux for the kr and kx call back functions (for python)"""
        self.rs = rs

    def get_inner_heads(self, shift = 0, weather:dict={}):#        
        """ matric potential at the root surface interface [cm]"""
        #print(self.mode, weather)
        rsx = np.zeros((len(self.cyls),))
        for i, cyl in enumerate(self.cyls):  # run cylindrical models
            if isinstance(cyl, AirSegment):
                rsx[i] = cyl.get_inner_head(val = self.rs.getPsiAir(weather["ea"]/weather["es"], weather["TairC"]))  # [cm]
            elif self.mode.startswith("dumux"):
                rsx[i] = cyl.getInnerHead(shift)  # [cm] (in richards.py, then richards_cyl.hh)
            elif self.mode.startswith("python"):
                rsx[i] = cyl.get_inner_head()  # [cm]
            else:
                raise Exception("RhizoMappedSegments.get_inner_heads: Warning, mode {:s} unknown".format(self.mode))
        return self._map(self._flat0(comm.gather( rsx, root = 0))) # gathers and maps correctly

    def get_inner_solutes(self, shift = 0):
        """ matric potential at the root surface interface [cm]"""
        rsx = np.zeros((len(self.cyls),))
        for i, cyl in enumerate(self.cyls):  # run cylindrical models
            rsx[i] = cyl.getInnerSolutes(shift)  # [cm]
            #print('rsx', rsx[i])
        return self._map(self._flat0(comm.gather(rsx, root = 0)))  # gathers and maps correctly

    def get_soil_k(self, rx):
        """ TODO """
        soil_k = np.zeros((len(self.cyls),))
        if self.mode.startswith("dumux"):
            for i, idx in enumerate(self.eidx):  # run cylindrical models
                rsx = self.cyls[i].getInnerHead()
                nidx = idx + 1  # segment index+1 = node index
                try:
                    soil_k[i] = ((vg.fast_mfp[self.vg_soil](rsx) - vg.fast_mfp[self.vg_soil](rx[nidx])) / (rsx - rx[nidx])) / self.dx2[i]
                except:
                    print(rsx, rx[nidx])
        else:
            raise Exception("RhizoMappedSegments.get_soil_k: Warning, mode {:s} not implemented or unknown".format(self.mode))
        return self._map(self._flat0(comm.gather(soil_k, root = 0)))

    def get_dx2(self):
        """ TODO doc me AND only for mode="dumux" yet (set in initialize)"""
        return self._map(self._flat0(comm.gather(self.dx2, root = 0)))

    def get_inner_fluxes(self):
        """ fluxes [cm3/day] at the root surface, i.e. inner boundary  """
        fluxes = np.zeros((len(self.cyls),))#
        for i, cyl in enumerate(self.cyls):  # run cylindrical models
            if isinstance(cyl, AirSegment):   
                j = self.eidx[i]  # for one process j == i
                fluxes[i] = float(cyl.getInnerFlux()) 
            elif self.mode.startswith("dumux"):            
                j = self.eidx[i]  # for one process j == i
                fluxes[i] = -float(cyl.getInnerFlux()) * (2 * np.pi * self.radii[j] * self.seg_length[j]) / self.radii[j]  # [cm/day] -> [cm3/day], ('/inner_radii' comes from cylindrical implementation)
            elif self.mode.startswith("python"):
                j = self.eidx[i]  # for one process j == i
                fluxes[i] = cyl.get_inner_flux() * (2 * np.pi * self.radii[j] * self.seg_length[j]) / self.last_dt  # divide by dt is correct here! getInnerFlux only gives the source in cm3/cm2
            else:
                raise Exception("RhizoMappedSegments.get_inner_fluxes: Warning, mode {:s} unknown".format(self.mode))
        return self._map(self._flat0(comm.gather(fluxes, root = 0)))  # gathers and maps correctly

    def get_inner_concentrations(self):  # TODO
        raise Exception("do not use get_inner_concentrations yet")
        """ solute concentration at the root surface interface [g / cm3]"""
        rsx = np.zeros((len(self.cyls),))
        if self.mode.startswith("dumux"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                 rsx[i] = cyl.getInnerHead()  # [cm] ?!
            #pass  # currently zero flux !!!!!!
        else:
            raise Exception("RhizoMappedSegments.get_inner_concentrations: Warning, mode {:s} unknown".format(self.mode))
        return self._map(self._flat0(comm.gather(rsx, root = 0)))  # gathers and maps correctly

    def get_inner_mass_fluxes(self):  # TODO check
        """ fluxes [cm3/day] at the root surface, i.e. inner boundary  """
        fluxes = np.zeros((len(self.cyls),))
        if self.mode.startswith("dumux"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                fluxes[i] = -float(cyl.getInnerFlux(1)) * (2 * np.pi * self.radii[j] * self.seg_length[j]) / self.radii[j]  # [cm/day] -> [cm3/day], ('/inner_radii' comes from cylindrical implementation)
        else:
            raise Exception("RhizoMappedSegments.get_inner_fluxes: Warning, mode {:s} unknown".format(self.mode))
        return self._map(self._flat0(comm.gather(fluxes, root = 0)))  # gathers and maps correctly

    def get_water_volumes(self):
        """ fluxes [cm3/day] at the root surface, i.e. inner boundary  """
        volumes = np.zeros((len(self.cyls),))
        if self.mode.startswith("dumux"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                volumes[i] = cyl.getWaterVolume()
        else:
            raise Exception("RhizoMappedSegments.get_water_volumes: Warning, mode {:s} unknown".format(self.mode))
        return self._map(self._flat0(comm.gather(volumes, root = 0)))  # gathers and maps correctly

    def solve(self, dt, *argv):#dt, seg_rx, proposed_outer_fluxes,kex,proposed_outer_sol_fluxes
        """ set bc for cyl and solves it """
        self.last_dt = dt
        
        if self.mode == "dumux":
            proposed_inner_fluxes = argv[0]
            proposed_outer_fluxes = argv[1]
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                l = self.seg_length[j]
                if isinstance(cyl, AirSegment):  
                    cyl.setInnerFluxCyl(proposed_inner_fluxes[j])
                else:                
                    cyl.setInnerFluxCyl(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l))  # [cm3/day] -> [cm /day]
                    cyl.setOuterFluxCyl(proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l))  # [cm3/day] -> [cm /day]
                    try:
                        cyl.solve(dt)
                    except:
                        myStr = "RhizoMappedSegments.solve: dumux exception with boundaries in flow {:g} cm3/day, out flow {:g} cm3/day, segment radii [{:g}-{:g}] cm"
                        myStr = myStr.format(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l), proposed_outer_fluxes[j] / (2 * np.pi * self.radii[j] * l), self.radii[j], self.outer_radii[j])
                        print("node ", self.nodes[self.segments[j].y])
                        self.plot_cylinder(j)
                        self.plot_cylinders()
                        raise Exception(myStr)
        elif ((self.mode == "dumux_dirichlet_nc") or (self.mode == "dumux_dirichlet_2c")
                or (self.mode == "dumux_dirichlet_10c")):            
            rx = argv[0]
            proposed_outer_fluxes = argv[1]
            proposed_inner_fluxes = argv[2]
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                    
                j = self.eidx[i]  # for one process j == i
                if isinstance(cyl, AirSegment):  
                    cyl.setInnerFluxCyl(proposed_inner_fluxes[j])
                    #cyl.setParameter( "Soil.BC.Bot.C1Value", str(exud)) 
                else:
                    l = self.seg_length[j]
                    cyl.setInnerMatricPotential(rx[i])
                    #inner = bot = plant
                    #outer = top = soil
                    
                    cyl.setParameter( "Soil.BC.Top.C1Value", proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l))  # [cm3/day] -> [cm /day]
                    
                    if not isinstance(argv[2],float): 
                        cyl.setInnerBC_solute("constantFluxCyl", argv[2][i])
                        #setBotBC_solute
                        cyl.setParameter(self.param_group + "BC.Bot.CValue", str(value_bot))
                    else:
                        cyl.setInnerBC_solute("constantFluxCyl", argv[2])
                    cyl.setOuterBC_solute("constantFluxCyl", argv[3])
                    try:
                        cyl.solve(dt)
                    except:
                        # str = "RhizoMappedSegments.solve: dumux exception with boundaries in flow {:g} cm3/day, out flow {:g} cm3/day, segment radii [{:g}-{:g}] cm"
                        # str = str.format(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l), proposed_outer_fluxes[j] / (2 * np.pi * self.radii[j] * l), self.radii[j], self.outer_radii[j])
                        # print("node ", self.nodes[self.segments[j].y])
                        self.plot_cylinder(j)
                        self.plot_cylinders()
                        self.plot_cylinders_solute()
                        raise Exception
        elif ((self.mode == "dumux_nc") or (self.mode == "dumux_10c")):
            #inner = bot = plant
            #outer = top = soil
            proposed_inner_fluxes = argv[0]
            proposed_outer_fluxes = argv[1]
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                if isinstance(cyl, AirSegment):  
                    cyl.setInnerFluxCyl(proposed_inner_fluxes[j])
                else:
                    l = self.seg_length[j]
                    if self.recreateComsol:
                        cyl.setInnerFluxCyl(-0.26)#proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l))  # [cm3/day] -> [cm /day]
                        cyl.setOuterFluxCyl(0)#proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l))  # [cm3/day] -> [cm /day]
                    else:   
                        cyl.setInnerFluxCyl(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l))  # [cm3/day] -> [cm /day]
                        cyl.setOuterFluxCyl(proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l))  # [cm3/day] -> [cm /day] 
                         
                    if not isinstance(argv[2],float): 
                        # cyl.setInnerBC_solute("constantFluxCyl", argv[2][i])
                        botVal = argv[2][i]
                        topVal = argv[3][i]
                    else:
                        # cyl.setInnerBC_solute("constantFluxCyl", argv[2])
                        botVal = argv[2]
                        topVal = argv[3]
                    
                    # cyl.setInnerFluxCyl(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l))  # [cm3/day] -> [cm /day]
                    #cyl.setParameter( "Soil.BC.Bot.C1Value",str(1))# str(botVal / (2 * np.pi * self.radii[j] * l)))  # [cm3/day] -> [cm /day]
                    #cyl.setParameter( "Soil.BC.Top.C1Value",str(1))#str(topVal / (2 * np.pi * self.outer_radii[j] * l)))  # [cm3/day] -> [cm /day]
                    typeBC = np.full(self.numComp,3)
                    
                    valueBC = np.full(self.numComp,0.)
                    if not self.recreateComsol:
                        valueBC[0] = topVal / (2 * np.pi * self.outer_radii[j] * l) # [cm3/day] -> [cm /day]
                        valueBC[1] = topVal / (2 * np.pi * self.outer_radii[j] * l) # [cm3/day] -> [cm /day]
                    cyl.setSoluteTopBC(typeBC, valueBC)
                    if self.recreateComsol:
                        valueBC[0] = 1.;valueBC[1] = 1.
                    else :
                        valueBC[0] = botVal / (2 * np.pi * self.radii[j] * l) # [cm3/day] -> [cm /day]
                        valueBC[1] = botVal / (2 * np.pi * self.radii[j] * l) # [cm3/day] -> [cm /day]
                    cyl.setSoluteBotBC(typeBC, valueBC)
                    try:
                        cyl.solve(dt, maxDt = 2500/(24*3600))
                    except:
                        myStr = "RhizoMappedSegments.solve: dumux exception with boundaries in flow {:g} cm3/day, out flow {:g} cm3/day, segment radii [{:g}-{:g}] cm"
                        myStr = myStr.format(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l), proposed_outer_fluxes[j] / (2 * np.pi * self.radii[j] * l), self.radii[j], self.outer_radii[j])
                        raise Exception(myStr)
        else:
            print(self.mode)
            raise Exception("RhizoMappedSegments.initialize: unknown solver {}".format(self.mode))

    def get_water_volume(self):
        """ returns the water volume of the cylindrical models [cm3] """
        volumes = np.zeros((len(self.cyls),))
        if self.mode.startswith("dumux"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                cyl_water = 0.
                cyl_water_content = cyl.getWaterContent()  # getWaterContent() in cpp/pyhton_binding/richards.hh
                nodes = cyl.getPoints()
                for k, wc in enumerate(cyl_water_content):
                    r1 = nodes[k]
                    r2 = nodes[k + 1]
                    cyl_water += np.pi * (r2 * r2 - r1 * r1) * self.seg_length[j] * wc
                volumes[i] = cyl_water
        elif self.mode.startswith("python"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                cyl_water = 0.
                cyl_water_content = cyl.get_water_content()  # segment 0
                for k, wc in enumerate(cyl_water_content):
                    r1 = cyl.grid.nodes[k]
                    r2 = cyl.grid.nodes[k + 1]
                    cyl_water += np.pi * (r2 * r2 - r1 * r1) * self.seg_length[j] * wc
                volumes[i] = cyl_water
        else:
            raise Exception("RhizoMappedSegments.get_water_volume: unknown solver {}".format(self.mode))
        return self._map(self._flat0(comm.gather(volumes, root = 0)))  # gathers and maps correctly

    def _map(self, x):
        """Converts @param x to a numpy array and maps it to the right indices                 """
        indices = self._flat0(comm.gather(self.eidx, root = 0))  # gather segment indices from all threads
        if indices:  # only for rank 0 it is not empty
            assert len(indices) == len(x), "RhizoMappedSegments._map: indices and values have different length"
            p = np.zeros((len(x),), dtype = np.float64)
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

    def plot_cylinder(self, i):
        """ plots a specific cylinder (DUMUX only, TODO) """
        cyl = self.cyls[i]
        x_ = cyl.getDofCoordinates()
        y_ = cyl.getSolutionHead()

        SMALL_SIZE = 22
        MEDIUM_SIZE = 22
        BIGGER_SIZE = 22
        plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
        plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
        plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
        plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
        plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
        plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
        plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title
        # plt.xlim([0., 0.2 ])
        plt.plot(x_, y_)
        plt.xlabel("distance [cm]")
        plt.ylabel("matric potential [cm]")
        plt.show()

    def plot_cylinders(self):
        """ plots a specific cylinder (DUMUX only, TODO) """
        inner, outer = [], []
        zz = -self.minBound.z
        for i, cyl in enumerate(self.cyls):
            if not isinstance(cyl, AirSegment):  
                x_ = cyl.getDofCoordinates()
                y_ = cyl.getSolutionHead()
                inner.append(y_[0])
                outer.append(y_[-1])
                j = self.segments[i].y
                z = self.nodes[j].z
                col_i = int(-z / zz * 255.)
                c_ = '#%02x%02x%02x' % (col_i, col_i, 64)
                plt.plot(x_, y_, alpha = 0.1, c = c_)
        plt.xlabel("distance [cm], deeper roots are yellow")
        plt.ylabel("matric potential [cm]")
        # plt.xlim([0.05, 0.6])
        # plt.ylim([-8500, 0. ])
        plt.show()
        return  np.argmin(inner), np.argmax(inner), np.argmin(outer), np.argmax(inner)

    def plot_cylinders_solute(self):
        """ plots a specific cylinder (DUMUX only, TODO) """
        inner, outer = [], []
        zz = -self.minBound.z
        for i, cyl in enumerate(self.cyls):
            if not isinstance(cyl, AirSegment):  
                x_ = cyl.getDofCoordinates()
                y_ = cyl.getSolution_(1)
                #y_ = cyl.getWaterContent() works
                inner.append(y_[0])
                outer.append(y_[-1])
                j = self.segments[i].y
                z = self.nodes[j].z
                col_i = int(-z / zz * 255.)
                c_ = '#%02x%02x%02x' % (col_i, col_i, 64)
                plt.plot(x_-x_[0], y_, alpha = 0.1, c = c_)
        plt.xlabel("distance [cm], deeper roots are yellow")
        plt.ylabel('solute concentration (g/cm3)')
        # plt.xlim([0.05, 0.6])
        # plt.ylim([-8500, 0. ])
        plt.show()
        return  np.argmin(inner), np.argmax(inner), np.argmin(outer), np.argmax(inner)


    def map_cylinders_solute(self, XX,YY,ZZ,name = ""):
        """ maps cylinders to soil grid """

        shape = np.shape(XX)
        conc = np.zeros((shape))
        for i, cyl in enumerate(self.cyls):
            if not isinstance(cyl, AirSegment):  
                R_ = cyl.getDofCoordinates()
                vv_ = cyl.getSolution_(1)
                
                p0 = np.array(self.nodes[self.segments[i].x])
                p1 = np.array(self.nodes[self.segments[i].y])
                
                v = p1 - p0
                mag = norm(v)
                v = v / mag
                not_v = np.array([1, 0, 0])
                if (v == not_v).all():
                    not_v = np.array([0, 1, 0])
                n1 = np.cross(v, not_v)
                n1 /= norm(n1)
                n2 = np.cross(v, n1)
                t = np.linspace(0, mag, 20)
                theta = np.linspace(0, 2 * np.pi, 20)
                t, theta = np.meshgrid(t, theta)
                
                x = []
                y = []
                z = []
                vv = []
                
                for j in range(0,len(R_)):
                    R = R_[j]
                    X, Y, Z = [p0[i] + v[i] * t + R * np.sin(theta) * n1[i] + R * np.cos(theta) * n2[i] for i in [0, 1, 2]]
                    x.extend(X.flatten())
                    y.extend(Y.flatten())
                    z.extend(Z.flatten())
                    vv.extend(np.ones((len(X.flatten())))*vv_[j])

                # interpolate "data.v" on new grid "inter_mesh"
                V = gd((x,y,z), vv, (XX.flatten(),YY.flatten(),ZZ.flatten()), method='linear')
                V[np.isnan(V)] = 0
                V = np.array(V.reshape(shape))
                conc = conc+V
                print('cyl '+str(i)+' of ' + str(len(self.cyls))+ ' is finished!')

        return conc
        
    def calc_model(number,XX,YY,ZZ,q):
        q.put([(self.calculate_cyls(number,XX,YY,ZZ,))])
        
    def calculate_cyls(i,XX,YY,ZZ): 
        cyl = self.cyls[i]
        R_ = cyl.getDofCoordinates()
        vv_ = cyl.getSolution_(1)
        
        p0 = np.array(self.nodes[self.segments[i].x])
        p1 = np.array(self.nodes[self.segments[i].y])
        
        v = p1 - p0
        mag = norm(v)
        v = v / mag
        not_v = np.array([1, 0, 0])
        if (v == not_v).all():
            not_v = np.array([0, 1, 0])
        n1 = np.cross(v, not_v)
        n1 /= norm(n1)
        n2 = np.cross(v, n1)
        t = np.linspace(0, mag, 20)
        theta = np.linspace(0, 2 * np.pi, 20)
        t, theta = np.meshgrid(t, theta)
        
        x = []
        y = []
        z = []
        vv = []
        
        for j in range(0,len(R_)):
            R = R_[j]
            X, Y, Z = [p0[i] + v[i] * t + R * np.sin(theta) * n1[i] + R * np.cos(theta) * n2[i] for i in [0, 1, 2]]
            x.extend(X.flatten())
            y.extend(Y.flatten())
            z.extend(Z.flatten())
            vv.extend(np.ones((len(X.flatten())))*vv_[j])

        # interpolate "data.v" on new grid "inter_mesh"
        V = gd((x,y,z), vv, (XX.flatten(),YY.flatten(),ZZ.flatten()), method='linear')
        V[np.isnan(V)] = 0
        V = np.array(V.reshape(shape))
        return V

    def map_cylinders_solute_parallel(self, XX,YY,ZZ,name = ""):
        """ maps cylinders to soil grid """
        
        shape = np.shape(XX)
        conc = np.zeros((shape))
        q = multiprocessing.Queue()
        for i in range(0,len(self.cyls)):
            if not isinstance(cyl, AirSegment):              
                p = multiprocessing.Process(target=self.calc_model, args=(i,XX,YY,ZZ,q,))
                p.start()
            
        for i in range(0,len(self.cyls)):
            conc = np.add(conc,np.array(q.get()))
        return conc
                  
    def collect_cylinder_solute_data(self):
        """ collects solute data from cylinders  """
        
        l_all = self.segLength()
        a_all = self.radii
        dist = []
        conc = []
        l = []
        a = []
        for i, cyl in enumerate(self.cyls):
            if not isinstance(cyl, AirSegment):  
                x = cyl.getDofCoordinates()
                x_ = x-x[0]
                dist.append(x_)
                conc.append(cyl.getSolution_(1))
                l.append(l_all[i]*np.ones((len(cyl.getDofCoordinates()))))
                #print('shape dist', np.shape(dist))

        return  dist, conc, l
    
def plot_transpiration(t, soil_uptake, realised_trans):
    """ plots potential and actual transpiration over time 
    
    depending on discretisation soil_uptake and root_uptake might differ  
    
    @param t                  times [day]
    @param soil_uptake        actual transpiration [cm3/day] of soil 
    @param realised_trans    [cm3/day]
    """
    fig, ax1 = plt.subplots()
    ax1.plot(t, realised_trans, 'k', label = "realised transpiration")  # realised transpiration
    ax1.plot(t, np.array(soil_uptake), 'g', label = "soil uptake")  # actual transpiration  according to soil model
    ax1.set_xlabel("Time [d]")
    ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
    ax2 = ax1.twinx()
    dt = np.diff(t)
    so = np.array(soil_uptake)
    cum_transpiration = np.cumsum(np.multiply(so[:-1], dt))
    ax2.plot(t[1:], cum_transpiration, 'c--', label = "Cumulative soil uptake")  # cumulative transpiration (neumann)
    ax2.set_ylabel("Cumulative soil uptake $[cm^3]$")
    print("Cumulative soil uptake", cum_transpiration[-1], "[cm^3]")
    fig.legend()
    plt.show()


def plot_info(x_, water_collar_cell, water_cyl, collar_sx, min_sx, min_rx, min_rsx, water_uptake, water_domain):
    """ 2x2 plot with additional information """
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    ax1.set_title("Water amount")
    ax1.plot(x_, np.array(water_collar_cell), label = "water cell")
    ax1.plot(x_, np.array(water_cyl), label = "water cylindric")
    ax1.legend()
    ax1.set_xlabel("Time (days)")
    ax1.set_ylabel("(cm3)")
    ax2.set_title("Pressure")
    ax2.plot(x_, np.array(collar_sx), label = "soil at root collar")
    ax2.plot(x_, np.array(min_sx), label = "min soil")
    ax2.plot(x_, np.array(min_rx), label = "min xylem")
    ax2.plot(x_, np.array(min_rsx), label = "min 1d at root surface")
    ax2.set_ylim([-15000, 0])
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

