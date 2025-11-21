""" 
Functions to simplify setup of the scenarios for the INARI project
"""

import sys;
sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../modules/");
sys.path.append("/data");
sys.path.append("../../../../CPlantBox/modelparameter/functional/plant_hydraulics");

import numpy as np
import pandas as pd
import timeit
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
                                                                                            
from rosi_richards4c_cyl import Richards4CCylFoam # C++ part (Dumux binding)
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial model (a single cylindrical model is not allowed to run in parallel)
from rosi_richards4c import Richards4CSPILU as Richards4CSP  # C++ part (Dumux binding), macroscopic soil model

# from rosi_richards2c import Richards2CSP  # C++ part (Dumux binding), macroscopic soil model
#from rosi_richards10c_cyl import RichardsNCCylFoam # C++ part (Dumux binding)
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial model (a single cylindrical model is not allowed to run in parallel)
#from rosi_richards10c import RichardsNCSPILU as RichardsNCSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
from functional.phloem_flux import PhloemFluxPython  # root system Python hybrid solver



import plantbox as pb  # CPlantBox
import functional.van_genuchten as vg
from functional.xylem_flux import *
from datetime import *
import wheat1997_conductivities 
from helpfull import *
import weatherFunctions 
from PhloemPhotosynthesis import *




def getBiochemParam(s,paramIdx):    
    """ define TraiRhizo biochemical parameters 
        @param: the dumux soil object
        @ param: index of the TraiRhizo parameter set to use
        
    """
    s.doSimpleReaction = 1 #only diffusion, decay, sorption and not Mona's complete model 
    
    s.Ds =7e-10  # for P from 10.1093/aob/mcaa120, Fig S7, S8 and eq 40 of article m^2/s, and from https://www.aqion.de/site/diffusion-coefficients
    s.Dl = 1.9e-9 #m^2/s https://www.aqion.de/site/diffusion-coefficients
    #s.Vmax_decay = 7.32e-5 #mol C / m^3 scv / s #decay rate from Nideggen et al. 
    #s.Km_decay = 10.5 #mol C / m^3 scv
    s.SourceSlope = 1000
    s.SourceSlope_up = 10
    s.CriticalPressure_up = -5
    if s.noAds:
        s.CSSmax = 0.
        s.alpha = 0.
        s.kads = 0.
        s.kdes = 1.
    else:
        s.BufferPower = 80
        s.freundlichN = 0.4 #https://pmc.ncbi.nlm.nih.gov/articles/PMC7489101/#s14
        # mg0.6 kg–1 L0.4 * (mol/mg)**0.6 * kgsoil/molsoil => mol0.6 mol-1 L0.4
        s.freundlichK = 124.8 / (s.mg_per_molP)**(1-s.freundlichN) * s.solidMolarMass # 
        
        
        s.kads = 10**2 # cm3/mol # ATT: can create 1d-3d converging error if kads is too high
        s.kdes = 1.# -
        s.Qmmax = 0.45 * 0.079 # max ratio gOC-gmineral soil, see 10.1016/j.soilbio.2020.107912
        # [g OC / g mineral soil] * [g mineral soil/ cm3 bulk soil] *[ mol C/g C]
        CSSmax_ = s.Qmmax * s.bulkMassDensity_gpercm3*(1/s.molarMassP)
        s.CSSmax = CSSmax_ # mol C/cm3 bulk soil
        #s.CSSmax = s.Qmmax * s.bulkDensity_m3 / 1e6 # mol OC/mol soil * [mol soil/m3] * [m3/cm3] =  mol/cm3
        
    s.css1Function = 2 #5 # current adsorption function implemented.
    return s

def setBiochemParam(s):
    """ send the TraiRhizo biochemical parameters to dumux
        @param: the dumux soil object
    """

    s.setParameter( "Soil.css1Function", str(s.css1Function))

    s.setParameter("1.Component.LiquidDiffusionCoefficient", str(s.Ds)) #m^2/s
    s.setParameter("2.Component.LiquidDiffusionCoefficient", str(s.Dl)) #m^2/s


    s.setParameter( "Soil.doSimpleReaction", str(s.doSimpleReaction))
    
    s.setParameter( "Component.BufferPower", str(s.BufferPower))
    s.setParameter( "Component.FreundlichN", str(s.freundlichN))
    s.setParameter( "Component.FreundlichK", str(s.freundlichK))
    s.setParameter( "Soil.SourceSlope", str(s.SourceSlope ))
    s.setParameter( "Soil.SourceSlope_up", str(s.SourceSlope_up ))
    s.setParameter( "Soil.CriticalPressure_up", str(s.CriticalPressure_up ))
    if s.dimWorld == 3:
        # 1 == True
        # if we define a source or sink for the cell 
        # (from the results of the 1d models),
        # do not compute on top of biochemical reactions in dumux
        s.setParameter("Problem.reactionExclusive", "1")
    
    return s


def getCSS(s, CSW):
    """ @return concentration of adsobed carbon in the soil
        according to @param CSW [mol C / cm3 water], the concentration
        of small solutes in the soil water
        @param: s the dumux soil object
    """
    raise Exception # use the one in cpp to have the varying porosity
    return  bulkSoilDensity * s.freundlichK *(CSW ** s.freundlichN) #(s.kads * CSW * s.CSSmax)/(s.kads * CSW + s.kdes) 
    

def setIC3D(s, paramIdx, ICcc = None):
    return setIC(s, paramIdx, ICcc)

def setIC(s, paramIdx, ICcc = None):
    """ Defined the initial concentraition of the solutes
        [mol C / cm3 water] for disolved solutes and [mol C / cm3 scv]
        for solutes in the soil phase
        @param: s the dumux soil object
        @param: paramIdx index of the TraiRhizo parameter set
        @param: ICcc (optional) predefined initial conditions
    """
    if ICcc is None:
        paramSet = pd.read_csv('./output_random_rows.csv').loc[paramIdx]

        C_S_mg = 0.035 # mg P L–1, 10.1093/aob/mcaa120 ' P generally accumulates in the topsoil'
        #paramSet['CS_init'] /s.mg_per_molC## small C solutes in mol/cm3 water
        # [mg P L–1] * [mol/mg] * [L/cm3] = mol P cm^-3
        C_S = C_S_mg / s.mg_per_molP * 1e-3
        C_L = 0. #paramSet['CL_init'] /s.mg_per_molC## large C solutes in mol/cm3 water

        # concentraiton of adsobed C_S
        #s.CSS_init  = getCSS(s, C_S) #mol C/ cm3 scv
        
        
            
        unitConversion = 1.0e6 # mol/cm3  => mol/m3 
        addedVar = 1. * float(s.doSoluteFlow) # empirical factor
        
        s.ICcc = np.array([C_S *unitConversion*addedVar,
                            
                           C_L*unitConversion*addedVar, # s.CSS_init*unitConversion*addedVar,
                           0.,#s.CSS_init*unitConversion*addedVar,#C_L*unitConversion*addedVar,
                           ])# in mol/m3 water or mol/m3 scv
        if rank == 0:
            print('init s.ICcc', s.ICcc)
    else:
        s.ICcc = ICcc
    
    for i in range(s.numSoluteComp):
        #mol/m3 to mol/mol
        molarC = s.ICcc[i] / s.phaseDensity(isDissolved = (i < s.numDissolvedSoluteComp)) 
        s.setParameter( "Soil.IC.C"+str(i+1), str(molarC ))# send data to dumux
    return s


def setDefault(s):
    """ Defined some usefull default parameters
    """
    molarMassWat = 18. # [g/mol]
    densityWat = 1. #[g/cm3]
    # [mol/cm3] = [g/cm3] /  [g/mol] 
    molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 
    s.molarDensityWat = molarDensityWat

    s.setParameter("Problem.dobioChemicalReaction",str(s.doBioChemicalReaction))
    
    s.setParameter("Problem.verbose", "0")
    
    s.setParameter("Newton.Verbosity", "0") 
    # force solute mole fraction > 0 and water in possible pressure ranges
    s.setParameter("Newton.EnableChop", "true")
    
    # UpwindWeight = 1, better when we have high solute gradient.
    # UpwindWeight = 0.5, better when have high water flow and low solute gradient
    s.setParameter("Flux.UpwindWeight", "0.5")#very important because we get high solute gradient.
    
    
    s.EnableResidualCriterion = False
    s.setParameter("Newton.EnableResidualCriterion", 
                     str( s.EnableResidualCriterion ))
    s.EnableAbsoluteResidualCriterion = False
    s.setParameter("Newton.EnableAbsoluteResidualCriterion", 
                     str( s.EnableAbsoluteResidualCriterion ))
    s.SatisfyResidualAndShiftCriterion = False
    s.setParameter("Newton.SatisfyResidualAndShiftCriterion",
                     str( s.SatisfyResidualAndShiftCriterion) )  
    s.MaxTimeStepDivisions = 10
    s.setParameter("Newton.MaxTimeStepDivisions",
                     str( s.MaxTimeStepDivisions) )  
    s.MaxSteps = 30
    s.setParameter("Newton.MaxSteps",
                     str( s.MaxSteps) )  
    s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift))
    
    return s

def getSoilTextureAndShape():  
    """ soil shape and texture data
        to adapt according to the soil represented
    """
    min_b = np.array([-10, -6, -90.]) # np.array( [5, 5, 0.] )
    max_b =np.array( [10, 6, 0.]) #  np.array([-5, -5, -5.])
    cell_number = np.array( [1,1,90])#np.array( [1,1,90])#[10,6,75])# np.array( [3,12,40]) #np.array([3,4,4])# np.array( [1,1,1]) # 1cm3 #np.array([3,3,3])
    solidDensity = 2650 # [kg/m^3 solid] #taken from google docs TraiRhizo
    solidMolarMass = 60.08e-3 # [kg/mol] 
    # theta_r, theta_s, alpha, n, Ks
    # soilVG = [0.049, 0.352, 0.019, 4.887, 421.67]
    soilTextureAndShape = {'min_b' : min_b,'max_b' : max_b,
                           'cell_number':cell_number,
                           "solidDensity":solidDensity,
                        'solidMolarMass': solidMolarMass,
                           #'soilVG':soilVG
                          }
    
    return soilTextureAndShape

def setSoilParam(s):    
    """ save the soil parameters
        @param: the dumux soil object
    """
    
    s.molarMassNO3 = 62 # g/mol
    s.molarMassP = 30 # g/mol
    s.mg_per_molNO3 = s.molarMassNO3 * 1000.
    s.mg_per_molP = s.molarMassP * 1000.
    
    soilTexture = getSoilTextureAndShape()
    s.solidDensity = soilTexture['solidDensity'] #[kg/m^3 solid] 
    s.solidMolarMass = soilTexture['solidMolarMass']# [kg/mol] 

    loam = [0.08, 0.43, 0.04, 1.6, 50]
    clay = [0.1, 0.4, 0.01, 1.1, 10]
    VG_parameters_at_furrow = [
        [0.049, 0.353, 0.019, 4.887, 421.67],  # Layer 0
        [0.048, 0.355, 0.019, 4.887, 421.67],  # Layer 1
        [0.047, 0.386, 0.01910, 3.82, 529],    # Layer 2
        [0.001, 0.707, 0.04520, 1.41, 538],    # Layer 3
        [0.021, 0.305, 0.0159, 1.89, 127.67],  # Layer 4
        [0.028, 0.242, 0.0124, 1.999, 7.89]    # Layer 5
    ]
    #VG_parameters_at_furrow = [
    #    [0.049, 0.353, 0.019, 4.887, 421.67],  # Layer 0
    #]
    #VG_parameters_at_furrow=[loam,clay,loam,clay,loam,clay]
    soil =VG_parameters_at_furrow #  [[0.028, 0.242, 0.0124, 1.999, 7.89]]#Combine the hydraulic conductivity vectors from all soil layers to define soil type for simulation  
    #soil = [l1, l2, l3, l4]
    #layers_ID = [4, 4, 3, 3, 2, 2, 1, 1]  
    #layers_pos = [-120., -57., -57., -33., -33, -20, -20, 0] 
    layers_pos  = [-200, -75, -75, -60, -60, -45, -45, -30, -30, -15, -15, 0]
    layers_ID = [6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1]
 
    #layers_pos  = [-200, 0]
    #layers_ID = [1, 1]

    
    s.setVGParameters(soil)
    s.setLayersZ(layers_ID, layers_pos)

    nitrate_z = [0.,-100.]  # top soil layer of 30 cm
    nitrate_initial_values = np.array([5.e-3,1.e-3]) / 0.43 / 1000 /s.molarMassNO3 #  [kg/m3] -> [mol/cm3]
    s.setICZ_solute(nitrate_initial_values[::-1], nitrate_z[::-1], 2)  # step-wise function, ascending order

    P_z = [0., -100.]  # top soil layer of 30 cm
    P_initial_values = np.array([0.5e-4, 0.035e-4])  /s.molarMassP  #  [g/cm3] -> [mol/cm3], 10 times lower than the actual total P because majority is adsorbed
    s.setICZ_solute(P_initial_values[::-1], P_z[::-1], 1)  # step-wise function, ascending order

    #s.setLayersZ(s.points_at_furrow, s.indices_at_furrow)#points_and_indices_at_furrow[:,2], points_and_indices_at_furrow[:,:2])
    #s.setVGParameters(s.soil)
    s.vg_soil = [vg.Parameters(ss) for ss in soil]
    
    # [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
    s.solidMolDensity = s.solidDensity/s.solidMolarMass
    # [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
    s.bulkDensity_m3 = s.solidMolDensity*(1.- s.vg_soil[0].theta_S)
    s.bulkMassDensity_gpercm3 = s.solidDensity*(1.- s.vg_soil[0].theta_S)*1000/1e6

    s.setParameter( "Soil.MolarMass", str(s.solidMolarMass))
    s.setParameter( "Soil.solidDensity", str(s.solidDensity))
    #s.setVGParameters([s.soil])
    
    return s


def create_soil_model3D( usemoles, results_dir ,
                        p_mean_ = -100,paramIndx =0,
                     noAds = False, ICcc = None, doSoluteFlow = True):
    return create_soil_model( usemoles , results_dir,
                        p_mean_,paramIndx ,
                     noAds , ICcc , doSoluteFlow)

def create_soil_model( usemoles, results_dir ,
                        p_mean_ = -100,paramIndx =0,
                     noAds = False, ICcc = None, doSoluteFlow = True,
                       doBioChemicalReaction=False,
                     MaxRelativeShift = 1e-8):
    """
        Creates a soil domain from @param min_b to @param max_b with resolution @param cell_number
        homogeneous domain 
        initial potentials are linear and mean potential is @p_mean_
        @ param: noAds: turn off adsorption?
        @param: paramIndx index of the TraiRhizo parameter set to use
        @param: ICcc (optional) initial concentraiton values for the solute components
        @param: usemoles [bool] dumux uses moles (True) or grammes (False)
        returns soil_model (RichardsWrapper(RichardsSP())) and soil parameter (vg.Parameters)
    """
        
    s = RichardsWrapper(Richards4CSP())  # water and N solute          
    s.results_dir = results_dir   
    # low MaxRelativeShift == higher precision in dumux
    s.MaxRelativeShift = MaxRelativeShift
    s.MaxRelativeShift_1DS = 1e-12
    s.doBioChemicalReaction = doBioChemicalReaction
    soilTextureAndShape = getSoilTextureAndShape() 
    min_b = soilTextureAndShape['min_b']
    max_b = soilTextureAndShape['max_b']
    cell_number = soilTextureAndShape['cell_number']
    s.cell_size = np.prod((max_b - min_b) / cell_number) # cm3 
    s.setParameter( "Soil.Grid.Cells", s.dumux_str(cell_number))  # send data to dumux
    s.noAds = noAds # no adsorption?
    s.doSoluteFlow = doSoluteFlow
    
    s.initialize() 
    setDefault(s)
    
    setSoilParam(s)
    getBiochemParam(s,paramIndx)
    setBiochemParam(s)
    # setIC3D(s, paramIndx, ICcc)
    s.isPeriodic = True
    s.createGrid(min_b, max_b, cell_number, s.isPeriodic)  # [cm] 
    p_mean_ = np.linspace (-110, -20, cell_number[2])
    s = setupOther(s, p_mean_)
    
    #if rank == 0:
    #    s.base.printParams()
    
    # just print once as will not change during simulation
    write_file_array("cellVol", np.array(s.getCellVolumes()), directory_ =s.results_dir) # cm3 
    write_file_array("cellIds", np.array(s.cellIndices), directory_ =s.results_dir) # cm3
    return s
    

def setupOther(s, p_mean_):
    """ define remaining soil parameters """
    s.p_meanInit = p_mean_
    
    
    s.setParameter("SpatialParams.Temperature","293.15") # todo: redefine at each time step?
    
    # BC
    if s.dimWorld == 1: # 1d model
        s.setParameter("Soil.IC.P", s.dumux_str(p_mean_))
        s.setInnerBC("fluxCyl",  s.win)  # [cm/day]
        s.setOuterBC("fluxCyl", 0.)
    if s.dimWorld == 3:# 3d model
        s.setTopBC("noFlux")
        #in acc. with Jorda et al. (2022), however, they assume inflow if h>0
        # also easier for checking mass balance
        s.setBotBC("freeDrainage") 
    
        for i in range(1, s.numComp):# no flux
            s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(2))
            s.setParameter( "Soil.BC.Top.C"+str(i)+"Type", str(2))
            s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(0)) 
            s.setParameter( "Soil.BC.Top.C"+str(i)+"Value", str(0 ))       
        
    
    # IC
    # if p_mean_ is float or int, define before initializeProblem
    if s.dimWorld == 3:
        if isinstance(p_mean_,(int,float)):
            #print('set pmean float',p_mean_)
            s.setHomogeneousIC(p_mean_, equilibrium = True)  # cm pressure head
        elif isinstance(p_mean_,type(np.array([]))):
            pass
        else:
            print(type(p_mean_))
            raise Exception
            
    s.maxDt =  250./(3600.*24.) # [s]
    s.maxDt_1DS = s.maxDt/10. # [s], lower maxDt for 1D models
    s.initializeProblem(s.maxDt)
    s.eps_regularization = 1e-10
    s.setRegularisation(s.eps_regularization, s.eps_regularization) # needs to be low when using sand parameters. 
    
     
    # if p_mean_ is np.array, define after initializeProblem
    if isinstance(p_mean_,(int,float)):
        pass
    elif isinstance(p_mean_,type(np.array([]))):
        s.setInitialConditionHead(p_mean_)
    else:
        print(type(p_mean_))
        raise Exception
    
    # for boundary conditions constantFlow, constantFlowCyl, and atmospheric
    s.wilting_point = -15000
    s.setCriticalPressure(s.wilting_point)  
    s.ddt = 1.e-5  # [day] initial Dumux time step
    s.bulkMassErrorWater_rel = 0.
    s.bulkMassErrorWater_relLim = 0.
    #pressureinit = s.getSolutionHead()
    thetainit = s.base.getWaterContent()
    
    s.totC3dInit = sum(s.getTotCContent()) # mol
    # initial soil water and solute content
    cell_volumes = s.getCellVolumes()  # cm3
    s.buWSoilInit = sum(np.multiply(np.array(s.getWaterContent()), cell_volumes)) # cm3 water
    return s


    
def create_mapped_plant(initSim, soil_model, fname, path, 
                stochastic = False, 
                        doPhloemFlow = True, static_plant = False,
                        usemoles = True, limErr1d3d = 1e-11, spellData =None):
    """ loads a rmsl file, or creates a rootsystem opening an xml parameter set,  
        and maps it to the soil_model """
    from rhizo_modelsPlant import RhizoMappedSegments  # Helper class for cylindrical rhizosphere models
    soilTextureAndShape = getSoilTextureAndShape() 
    min_b = soilTextureAndShape['min_b']
    max_b = soilTextureAndShape['max_b']
    cell_number = soilTextureAndShape['cell_number']
    
    if fname.endswith(".rsml"):
        raise Exception # to implement
        plantModel = XylemFluxPython(fname)
        
        perirhizalModel = RhizoMappedSegments(  mode, soil_model,  usemoles, seedNum = seed, 
                                 limErr1d3dAbs = limErr1d3d)
    elif fname.endswith(".xml"):
        seed = 1
        weatherInit = weatherFunctions.weather(initSim,0)#, spellData)
        ms = pb.MappedPlant()
        perirhizalModel = RhizoMappedSegments(soilModel = soil_model, 
                             ms = pb.MappedRootSystem(),
                             limErr1d3dAbs = limErr1d3d, 
                             RichardsNCCylFoam = Richards4CCylFoam)

        perirhizalModel.ms.setSeed(seed)
        perirhizalModel.ms.readParameters(path + fname)
        if soil_model.isPeriodic:
            perirhizalModel.ms.setGeometry(pb.SDF_PlantBox(np.inf, np.inf, max_b[2]-min_b[2]))
        else:
            perirhizalModel.ms.setGeometry(pb.SDF_PlantBox( max_b[0]-min_b[0],  max_b[1]-min_b[1], max_b[2]-min_b[2]))
        #perirhizalModel.ms.setGeometry(pb.SDF_PlantBox( max_b[0]-min_b[0],  max_b[1]-min_b[1], max_b[2]-min_b[2]))
        perirhizalModel.ms.initialize(verbose = False)
        perirhizalModel.ms.simulate(initSim,verbose= False)
        if doPhloemFlow:
            plantModel = PhloemFluxPython(perirhizalModel.ms,psiXylInit = -659.8 - min_b[2],ciInit = weatherInit["cs"]*0.5) 
        else:
            plantModel = XylemFluxPython(perirhizalModel.ms)
            
    perirhizalModel.ms.constantLoc = True # segments remain in the voxel they first appeared in
    plantModel.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), 
                            pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), 
                            cut = False, # do not cut plant segments according to the voxel size
                            noChanges = True) # segments remain in the voxel they first appeared in
    
    #  function that return the index of a given position in the soil grid 
    picker = lambda x, y, z: soil_model.pick_([x,y, z])  
    # maps segments, maps root segements and soil grid indices to each other in both directions
    plantModel.rs.setSoilGrid(picker)
    
    plantModel.wilting_point = -15000.
    # set kr and kx for root system or plant
    wheat1997_conductivities.init_conductivities(r = plantModel, age_dependent = False)#not static_plant)
    if doPhloemFlow:   
        plantModel = phloemParam(plantModel, weatherInit)
    perirhizalModel.set_plantModel(plantModel)
    perirhizalModel.rhizoMassWError_rel = 0.
    perirhizalModel.rhizoMassWError_relLim = 0.
    perirhizalModel.new_soil_solute = np.array([0.])
    
    
    # to get time spent in each part of the code. TODO: currently, does not work well
    plantModel.time_start_global = timeit.default_timer()
    plantModel.time_rhizo_cumul = 0
    plantModel.time_3ds_cumul = 0
    plantModel.time_plant_cumulW = 0
    plantModel.time_plant_cumulS = 0
    plantModel.time_rhizo_i = 0
    plantModel.time_3ds_i = 0
    plantModel.time_plant_cumul = 0
    # cumulative transpiration
    plantModel.TranspirationCumul = 0 # real cumulative transpiration
    plantModel.TranspirationCumul_eval = 0 # cumulative transpiration during period with dinamic soil (for mass balance check)
    # cumulative flow    
    plantModel.seg_fluxesCumul = np.zeros((perirhizalModel.soilModel.numFluidComp, len(plantModel.plant.radii)))
     
    
    
    return perirhizalModel, plantModel
    





