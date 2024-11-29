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

# from rosi_richards2c import Richards2CSP  # C++ part (Dumux binding), macroscopic soil model
from rosi_richards10c_cyl import RichardsNCCylFoam # C++ part (Dumux binding)
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial model (a single cylindrical model is not allowed to run in parallel)
from rosi_richards10c import RichardsNCSPILU as RichardsNCSP  # C++ part (Dumux binding), macroscopic soil model
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
    # file containing the TraiRhizo parameter sets
    paramSet = pd.read_csv('./output_random_rows.csv').iloc[paramIdx].to_dict()    
    s.molarMassC = 12.011
    s.mg_per_molC = s.molarMassC * 1000.
    s.betaC = paramSet['beta_C'] # -
    s.betaO = paramSet['beta_O'] # -

    ## ATT! C_S_W_thres has to be > 0
    s.C_S_W_thresC =  paramSet['C_thres,C']/s.mg_per_molC# mgC/cm3 water => mol/cm3 water
    s.C_S_W_thresO =  paramSet['C_thres,O']/s.mg_per_molC# mgC/cm3 water => mol/cm3 water

    Ds = paramSet['DS_W'] #cm^2/d
    Dl = 0.003456 #cm^2/d
    s.Ds = Ds /(24*3600) /10000 # m^2/s
    s.Dl = Dl /(24*3600) /10000# m^2/s

    s.k_SC =  paramSet['k_C,S'] * s.mg_per_molC/paramSet['theta'] #[cm^3 water /mgC/d] => [cm^3 bulk soil/mol/d ]
    s.k_DC =  paramSet['k_d,C'] #1/d
    s.k_DO =  paramSet['k_d,O']#1/d

    s.K_L =  paramSet['K_L'] /s.mg_per_molC#[mol cm-3 soil]
    s.k_SO =  paramSet['k_O,S'] * s.mg_per_molC/paramSet['theta']
    #s.thetaInit = paramSet['theta']

    s.k_RC =  paramSet['k_r,C']#1/d
    s.k_RO =  paramSet['k_r,O']#1/d

    s.m_maxC =  paramSet['m_max,C']#1/d,  Maximum maintenance rate coefficient for the corresponding microbial groups
    s.m_maxO =  paramSet['m_max,O']#1/d,  Maximum maintenance rate coefficient for the corresponding microbial groups
    s.k_decay2 =  paramSet['p_L'] # - , Proportion of large polymers formed from dead microbial biomass due to maintenance

    s.micro_maxC =  paramSet['u_max,C']#1/d
    s.micro_maxO =  paramSet['u_max,O']#1/d

    s.v_maxL = paramSet['v_max,L'] #1/d
    s.k_decay = paramSet['Y'] # -, Maintenance yield
    s.k_growthC = paramSet['Y_C'] # -
    s.k_growthO = paramSet['Y_O'] #- 
    # can cause issue
    s.k_sorp = paramSet['k_sorp'] /s.mg_per_molC# mg C/cm 3 soil solution =>  mol C/cm3 soil solution
    
    s.alpha = 0.1 # -
    s.f_sorp = 0# 0.5
    s.k_phi = 0.1
    s.C_aOLim=1.e-10 * float(s.doSoluteFlow) # so that microbe community can always regrow
    s.C_aCLim=1.e-10 * float(s.doSoluteFlow) # so that microbe community can always regrow
    
    if s.noAds:
        s.CSSmax = 0.
        s.alpha = 0.
    else:
        s.Qmmax = 0.45 * 0.079 # max ratio gOC-gmineral soil, see 10.1016/j.soilbio.2020.107912
        # [g OC / g mineral soil] * [g mineral soil/ cm3 bulk soil] *[ mol C/g C]
        CSSmax_ = s.Qmmax * s.bulkMassDensity_gpercm3*(1/s.molarMassC)
        s.CSSmax = CSSmax_ # mol C/cm3 bulk soil
        #s.CSSmax = s.Qmmax * s.bulkDensity_m3 / 1e6 # mol OC/mol soil * [mol soil/m3] * [m3/cm3] =  mol/cm3

        
    s.css1Function = 9 # current adsorption function implemented.
    kads = 7.07e+02 # m3/kgC/yr, see 10.1016/j.soilbio.2020.107912, A.3
    yr_per_d = 1/365 # [yr/d]
    m3_per_cm3 = 1e-6; # m3/cm3
    cm3_per_m3 = 1e6; # cm3/m3
    
    # [kg/g] * [g/mol] = kg/mol
    mol_per_kgC = (1/1000) * s.molarMassC
    # [m3/kgC/yr] * [yr/d] * [cm3/m3] * [kgC/mol] = [cm3/mol/d]
    s.kads = kads * yr_per_d * cm3_per_m3 * mol_per_kgC
    
    kdes =  1.63e+03 # [1/yr] see 10.1016/j.soilbio.2020.107912, A.3
    s.kdes = kdes * yr_per_d
    return s

def setBiochemParam(s):
    """ send the TraiRhizo biochemical parameters to dumux
        @param: the dumux soil object
    """

    s.setParameter( "Soil.css1Function", str(s.css1Function))
    s.setParameter("Soil.C_aOLim", str(s.C_aOLim)) #[molC/cm3 scv]
    s.setParameter("Soil.C_aCLim", str(s.C_aCLim)) #[molC/cm3 scv]

    s.setParameter("1.Component.LiquidDiffusionCoefficient", str(s.Ds)) #m^2/s
    s.setParameter("2.Component.LiquidDiffusionCoefficient", str(s.Dl)) #m^2/s


    s.setParameter("Soil.betaC", str(s.betaC)) # -
    s.setParameter("Soil.betaO", str(s.betaO )) # -
    s.setParameter("Soil.C_S_W_thresC", str(s.C_S_W_thresC )) #mol/cm3
    s.setParameter("Soil.C_S_W_thresO", str(s.C_S_W_thresO )) #mol/cm3
    s.setParameter("Soil.k_decay", str(s.k_decay ))
    s.setParameter("Soil.k_decay2", str(s.k_decay2))
    s.setParameter("Soil.k_DC", str(s.k_DC )) # 1/d
    s.setParameter("Soil.k_DO", str(s.k_DO )) # 1/d
    s.setParameter("Soil.k_growthC", str(s.k_growthC ))
    s.setParameter("Soil.k_growthO", str(s.k_growthO ))
    s.setParameter("Soil.K_L", str(s.K_L))#[mol/cm3]
    s.setParameter("Soil.k_phi", str(s.k_phi ))
    s.setParameter("Soil.k_RC", str(s.k_RC))
    s.setParameter("Soil.k_RO", str(s.k_RO))

    s.setParameter("Soil.k_SC", str(s.k_SC )) #cm^3/mol/d
    s.setParameter("Soil.k_SO", str(s.k_SO )) #cm^3/mol/d

    s.setParameter("Soil.m_maxC", str(s.m_maxC  ))# 1/d
    s.setParameter("Soil.m_maxO", str(s.m_maxO  ))# 1/d
    s.setParameter("Soil.micro_maxC", str(s.micro_maxC ))# 1/d
    s.setParameter("Soil.micro_maxO", str(s.micro_maxO ))# 1/d
    s.setParameter("Soil.v_maxL", str(s.v_maxL))#[d-1]

    s.setParameter("Soil.k_sorp", str(s.k_sorp)) # mol / cm3 soil water or mol
    s.setParameter("Soil.f_sorp", str(s.f_sorp)) #[-]
    s.setParameter("Soil.CSSmax", str(s.CSSmax)) #[mol/cm3 scv zone 1] or mol
    s.setParameter("Soil.alpha", str(s.alpha)) #[1/d]
    s.setParameter("Soil.kads", str(s.kads)) #[cm3/mol/d] or [1/d]
    s.setParameter("Soil.kdes", str(s.kdes)) #[1/d]
    
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
    return  (s.kads * CSW * s.CSSmax)/(s.kads * CSW + s.kdes) 
    

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
        C_S = paramSet['CS_init'] /s.mg_per_molC## small C solutes in mol/cm3 water
        C_L = paramSet['CL_init'] /s.mg_per_molC## large C solutes in mol/cm3 water

        # concentraiton of adsobed C_S
        s.CSS_init  = getCSS(s, C_S) #mol C/ cm3 scv
        
            
        unitConversion = 1.0e6 # mol/cm3  => mol/m3 
        addedVar = 1. * float(s.doSoluteFlow) # empirical factor
        s.CSW_init = C_S * unitConversion
        s.ICcc = np.array([C_S *unitConversion*addedVar,
                           C_L*unitConversion*addedVar,
                            9.16666666666667e-07* unitConversion*addedVar,
                            8.33333333333333e-06* unitConversion*addedVar,
                            8.33333333333333e-07* unitConversion*addedVar,
                            8.33333333333333e-06* unitConversion*addedVar,
                            s.CSS_init*unitConversion*addedVar,
                           0.])# in mol/m3 water or mol/m3 scv
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
    s.MaxSteps = 18
    s.setParameter("Newton.MaxSteps",
                     str( s.MaxSteps) )  
    s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift))
    
    return s

def getSoilTextureAndShape():  
    """ soil shape and texture data
        to adapt according to the soil represented
    """
    min_b = np.array([-3./2, -12./2, -40.]) # np.array( [5, 5, 0.] )
    max_b =np.array( [3./2, 12./2, 0.]) #  np.array([-5, -5, -5.])
    cell_number = np.array( [3,12,40])# np.array( [3,12,40]) #np.array([3,4,4])# np.array( [1,1,1]) # 1cm3 #np.array([3,3,3])
    solidDensity = 2650 # [kg/m^3 solid] #taken from google docs TraiRhizo
    solidMolarMass = 60.08e-3 # [kg/mol] 
    # theta_r, theta_s, alpha, n, Ks
    soilVG = [0.049, 0.352, 0.019, 4.887, 421.67]
    soilTextureAndShape = {'min_b' : min_b,'max_b' : max_b,
                           'cell_number':cell_number,
                           "solidDensity":solidDensity,
                        'solidMolarMass': solidMolarMass,
                           'soilVG':soilVG}
    
    return soilTextureAndShape

def setSoilParam(s):    
    """ save the soil parameters
        @param: the dumux soil object
    """
    soilTexture = getSoilTextureAndShape()
    s.solidDensity = soilTexture['solidDensity'] #[kg/m^3 solid] 
    s.solidMolarMass = soilTexture['solidMolarMass']# [kg/mol] 
    s.soil =  soilTexture['soilVG'] 
    
    s.vg_soil = vg.Parameters(s.soil) 
    # [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
    s.solidMolDensity = s.solidDensity/s.solidMolarMass
    # [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
    s.bulkDensity_m3 = s.solidMolDensity*(1.- s.vg_soil.theta_S)
    s.bulkMassDensity_gpercm3 = s.solidDensity*(1.- s.vg_soil.theta_S)*1000/1e6

    s.setParameter( "Soil.MolarMass", str(s.solidMolarMass))
    s.setParameter( "Soil.solidDensity", str(s.solidDensity))
    s.setVGParameters([s.soil])
    
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
        
    s = RichardsWrapper(RichardsNCSP(), usemoles)  # water and N solute          
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
    setIC3D(s, paramIndx, ICcc)
    s.isPeriodic = True
    s.createGrid(min_b, max_b, cell_number, s.isPeriodic)  # [cm] 
    s = setupOther(s, p_mean_)
    
    if rank == 0:
        s.base.printParams()
    
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
        s.setBotBC("noFlux") 
    
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
    #thetainit = s.getWaterContent_()
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
        weatherInit = weatherFunctions.weather(initSim,0, spellData)
        ms = pb.MappedPlant()
        perirhizalModel = RhizoMappedSegments(soilModel = soil_model, 
                                 usemoles=usemoles,
                                              ms = ms,
                                 #seedNum = seed, 
                                 limErr1d3dAbs = limErr1d3d)

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
    # plantModel.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), 
    #                        pb.Vector3d(max_b[0], max_b[1], max_b[2]),
    #                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), 
    #                        cut = False, # do not cut plant segments according to the voxel size
    #                        noChanges = True) # segments remain in the voxel they first appeared in
    
    #  function that return the index of a given position in the soil grid 
    picker = lambda x, y, z: soil_model.pick_([x,y, z])  
    # maps segments, maps root segements and soil grid indices to each other in both directions
    plantModel.rs.setSoilGrid(picker)
    
    plantModel.wilting_point = -15000.
    # set kr and kx for root system or plant
    wheat1997_conductivities.init_conductivities(r = plantModel, age_dependent = not static_plant)
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
    plantModel.seg_fluxes0Cumul = np.array([])
    plantModel.seg_fluxes1Cumul = np.array([])
    plantModel.seg_fluxes2Cumul = np.array([])
    
    
    return perirhizalModel, plantModel
    





