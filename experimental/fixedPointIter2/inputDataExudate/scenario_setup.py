""" 
Functions to simplify setup of the scenarios for the INARI project
"""

import sys;
sys.path.append("../../../build-cmake/cpp/python_binding/");
import numpy as np
import pandas as pd
import timeit
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

from rosi_richards10c_cyl import RichardsNCCylFoam # C++ part (Dumux binding)
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial model (a single cylindrical model is not allowed to run in parallel)
from rosi_richards10c import RichardsNCSPILU as RichardsNCSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
from functional.phloem_flux import PhloemFluxPython  # root system Python hybrid solver


import plantbox as pb  # CPlantBox
import functional.van_genuchten as vg
from functional.xylem_flux import *
from datetime import *
import plantParameters
from helpfull import *
from PhloemPhotosynthesis import *
#import weatherFunctions
import evapotranspiration as evap



def getBiochemParam(s,soil_type):    
    """ define TraiRhizo biochemical parameters 
        @param: the dumux soil object
        @ param: index of the TraiRhizo parameter set to use
        
    """
    s.doSimpleReaction = 1
    
    # file containing the TraiRhizo parameter sets
    paramSet = pd.read_csv('./output_random_rows.csv').iloc[soil_type].to_dict()
    s.molarMassC = 12.011
    s.mg_per_molC = s.molarMassC * 1000.
    s.Ds = 1.e-10 #m^2/s
    s.Vmax = 7.32e-5 #mol C / m^3 scv / s
    s.Km = 10.5 #mol C / m^3 scv
    
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
    kgC_per_mol = (1/1000) * s.molarMassC
    # [m3/kgC/yr] * [yr/d] * [cm3/m3] * [kgC/mol] = [cm3/mol/d]
    s.kads = kads * yr_per_d * cm3_per_m3 * kgC_per_mol
    
    kdes =  1.63e+03 # [1/yr] see 10.1016/j.soilbio.2020.107912, A.3
    s.kdes = kdes * yr_per_d
    
    
    return s

def setBiochemParam(s):
    """ send the TraiRhizo biochemical parameters to dumux
        @param: the dumux soil object
    """
    s.setParameter( "Soil.doSimpleReaction", str(s.doSimpleReaction))
    s.setParameter( "Soil.css1Function", str(s.css1Function))

    #diffusion
    s.setParameter("1.Component.LiquidDiffusionCoefficient", str(s.Ds)) #m^2/s
    
    #decay
    s.setParameter("Soil.Vmax", str(s.Vmax)) #m^2/s
    s.setParameter("Soil.Km", str(s.Km)) #m^2/s

    #sorption
    s.setParameter("Soil.CSSmax", str(s.CSSmax)) #[mol/cm3 scv zone 1] or mol
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
    return  (s.kads * CSW * s.CSSmax)/(s.kads * CSW + s.kdes) #kd*CSW
    

def setIC3D(s, soil_type, ICcc = None):
    return setIC(s, soil_type, ICcc)

def setIC(s, soil_type, ICcc = None):
    """ Defined the initial concentraition of the solutes
        [mol C / cm3 water] for disolved solutes and [mol C / cm3 scv]
        for solutes in the soil phase
        @param: s the dumux soil object
        @param: soil_Type, 0 = loam, 1 = sand
        @param: ICcc (optional) predefined initial conditions
    """
    if ICcc is None:
        #paramSet = pd.read_csv('./output_random_rows.csv').loc[paramIdx]
        #C_S = paramSet['CS_init'] /s.mg_per_molC## small C solutes in mol/cm3 water
        #C_L = paramSet['CL_init'] /s.mg_per_molC## large C solutes in mol/cm3 water
        C_S = 0
        C_L = 0

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

    # low MaxRelativeShift == higher precision in dumux
    #s.setParameter("Problem.dobioChemicalReaction",str(s.doBioChemicalReaction))
    s.setParameter("Problem.dobioChemicalReaction",str(s.doBioChemicalReaction))
    s.setParameter("Problem.verbose", "0")
    s.setParameter("Newton.Verbosity", "0") 
    
    # force solute mole fraction > 0 and water in possible pressure ranges
    s.setParameter("Newton.EnableChop", "true")
    
    # UpwindWeight = 1, better when we have high solute gradient.
    # UpwindWeight = 0.5, better when have high water flow and low solute gradient
    s.setParameter("Flux.UpwindWeight", "1")#very important because we get high solute gradient.
    

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
    
def vg_SPP(i = int(1)):
    """ Van Genuchten parameter, called by maize()  """
        
    soil = {}
    # theta_r, theta_s, alpha, n, Ks
    soil[0] = [0.041, 0.494, 0.0256, 1.49, 245]
    soil[1] = [0.03, 0.414, 0.038, 2, 1864]
    return soil[i]

def getSoilTextureAndShape(soil_= "loam"):  
    """ soil shape and texture data
        to adapt according to the soil represented
    """
    min_b = np.array([-10., -22, -74.]) # np.array( [5, 5, 0.] )
    max_b =np.array( [10., 22, 0.]) #  np.array([-5, -5, -5.])
    cell_number = np.array([10,22,37]) # np.array( [3,12,40]) #np.array( [1,1,1]) # 1cm3
    area = 20 * 44  # cm2
    solidDensity = 2650 # [kg/m^3 solid] #taken from google docs TraiRhizo
    solidMolarMass = 60.08e-3 # [kg/mol] 
    # theta_r, theta_s, alpha, n, Ks
    
    if soil_ == "loam":
        i = 0
    else:
        i = 1
     
    soilVG = vg_SPP(i)
    
    Kc_value = np.array([1,1,1,1.2,1.2,1.2])
    Kc_days = np.array([1,42,63,98,154,288])
    
    Kc = np.zeros((Kc_days[-1]))
    dummy = 0
    for i in range(0,len(Kc)):
        if i+1 in Kc_days:
            Kc[i] = Kc_value[np.where(Kc_days==(i+1))[0]]
            dummy = dummy+1
        else:
            slope = (Kc_value[dummy]-Kc_value[dummy-1])/(Kc_days[dummy]-Kc_days[dummy-1])
            Kc[i] = Kc_value[dummy-1]+slope*((i+1)-Kc_days[dummy-1])
    
    
    soilTextureAndShape = {'min_b' : min_b,'max_b' : max_b,
                            'area':area,
                           'cell_number':cell_number,
                           "solidDensity":solidDensity,
                        'solidMolarMass': solidMolarMass,
                           'soilVG':soilVG,
                           'Kc':Kc}
    
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

def create_soil_model(simMax, usemoles, results_dir , soil_='loam',
                     noAds = False, ICcc = None, doSoluteFlow = True,
                     doBioChemicalReaction=True, 
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
    if soil_ == "loam":
        soil_type = 0
    else:
        soil_type = 1
        
    s = RichardsWrapper(RichardsNCSP(), usemoles)  # water and N solute          
    s.results_dir = results_dir   
    s.pindx = soil_type
    
    s.MaxRelativeShift = MaxRelativeShift # 1e-10
    s.MaxRelativeShift_1DS = MaxRelativeShift
    
    soilTextureAndShape = getSoilTextureAndShape(soil_type) 
    min_b = soilTextureAndShape['min_b']
    max_b = soilTextureAndShape['max_b']
    cell_number = soilTextureAndShape['cell_number']
    s.cell_size = np.prod((max_b - min_b) / cell_number) # cm3 
    s.setParameter( "Soil.Grid.Cells", s.dumux_str(cell_number))  # send data to dumux
    s.noAds = noAds # no adsorption?
    s.doSoluteFlow = doSoluteFlow
    s.doBioChemicalReaction = doBioChemicalReaction
    
    s.setParameter("Newton.Verbosity", "0") 
    s.initialize() 
    setDefault(s)
    
    setSoilParam(s)
    getBiochemParam(s,soil_type)
    setBiochemParam(s)
    setIC3D(s, soil_type, ICcc)
    s.isPeriodic = True
    s.createGrid(min_b, max_b, cell_number, s.isPeriodic)  # [cm] 
    s = setupOther(s, soil_, simMax)
    
    if rank == 0:
        s.base.printParams()
    
    # just print once as will not change during simulation
    write_file_array("cellVol", np.array(s.getCellVolumes()), directory_ =s.results_dir) # cm3 
    write_file_array("cellIds", np.array(s.cellIndices), directory_ =s.results_dir) # cm3
    
    return s
    

def setupOther(s, soil_type, simMax):
    """ define remaining soil parameters """ 
    
    # climate data 
    soilTextureAndShape = getSoilTextureAndShape(soil_type) 
    cell_number = soilTextureAndShape['cell_number']
    times, net_inf = evap.net_infiltration(soil_type, simMax, soilTextureAndShape['Kc'])
    
    s.setParameter("SpatialParams.Temperature","293.15") # todo: redefine at each time step?
    
    # BC
    p_mean_ = -100
    if s.dimWorld == 1: # 1d model
        s.setParameter("Soil.IC.P", s.dumux_str(p_mean_))
        s.setInnerBC("fluxCyl",  s.win)  # [cm/day]
        s.setOuterBC("fluxCyl", 0.)
    if s.dimWorld == 3:# 3d model

        if times is not None:
            s.setTopBC("atmospheric", 0.5, [times, net_inf])  # 0.5 is dummy value
        else:
            s.setTopBC("noFlux")

        s.setBotBC("freeDrainage") 
    
        for i in range(1, s.numComp):# no flux
            s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(2))
            s.setParameter( "Soil.BC.Top.C"+str(i)+"Type", str(2))
            s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(0)) 
            s.setParameter( "Soil.BC.Top.C"+str(i)+"Value", str(0 ))       
        
    
    # IC  
    s.maxDt =  250./(3600.*24.)
    s.maxDt_1DS = s.maxDt # [s], lower maxDt for 1D models
    s.initializeProblem(s.maxDt)
    
    s.eps_regularization = None # pcEps, krEps
    #s.setRegularisation(s.eps_regularization, s.eps_regularization) # needs to be l
     
    df = pd.read_csv("../inputDataExudate/data/init_pot_2019.csv")  # initial potential
    h  = np.flip(df[soil_type].loc[:].values) #cm
    h = np.repeat(h[:,np.newaxis],cell_number[0],axis=1) #x-axis
    h = np.repeat(h[:,:,np.newaxis],cell_number[1],axis=2) #y-axis
    h = h.flatten()
    #h = np.ones((20*45*75))*-100 #TODO
    s.setInitialConditionHead(h)  # cm
    
    # for boundary conditions constantFlow, constantFlowCyl, and atmospheric
    s.wilting_point = -15000
    s.setCriticalPressure(s.wilting_point)  
    s.ddt = 1.e-3  # [day] initial Dumux time step
    s.bulkMassErrorWater_rel = 0.
    s.bulkMassErrorWater_relLim = 0.    
    s.totC3dInit = sum(s.getTotCContent()) # mol    
    # initial soil water and solute content
    cell_volumes = s.getCellVolumes()  # cm3
    s.buWSoilInit = sum(np.multiply(np.array(s.getWaterContent()), cell_volumes)) # cm3 water
    return s


    
def create_mapped_rootsystem(initSim, simMax, soil_model, fname, path, soil_type,
                stochastic = False, 
                static_plant = False,
                usemoles = True, limErr1d3d = 1e-11):
    """ loads a rmsl file, or creates a rootsystem opening an xml parameter set,  
        and maps it to the soil_model """
    from rhizo_modelsPlant import RhizoMappedSegments  # Helper class for cylindrical rhizosphere models
    soilTextureAndShape = getSoilTextureAndShape(soil_type) 
    min_b = soilTextureAndShape['min_b']
    max_b = soilTextureAndShape['max_b']
    cell_number = soilTextureAndShape['cell_number']
    
    if fname.endswith(".rsml"):
        plantModel = XylemFluxPython(fname)
        
        perirhizalModel = RhizoMappedSegments(  mode, soil_model,  usemoles, seedNum = seed, 
                                 limErr1d3dAbs = limErr1d3d)
    elif fname.endswith(".xml"):
        seed = 1
        perirhizalModel = RhizoMappedSegments(soilModel = soil_model, 
                                 usemoles=usemoles,
                                 ms = pb.MappedRootSystem(),
                                 limErr1d3dAbs = limErr1d3d)

        perirhizalModel.ms.setSeed(seed)
        perirhizalModel.ms.readParameters(path + fname)
        
        perirhizalModel.ms.setGeometry(pb.SDF_PlantBox( max_b[0]-min_b[0],  max_b[1]-min_b[1], max_b[2]-min_b[2]))
        perirhizalModel.ms.initializeLB(5,4)
        perirhizalModel.ms.simulate(initSim,verbose= False) #initsim or 1?
        plantModel = XylemFluxPython(perirhizalModel.ms)
            
    #perirhizalModel.ms.constantLoc = True
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
    plantModel.transpiration = evap.get_transpiration(simMax, soilTextureAndShape['area'], soilTextureAndShape['Kc'], soil_type)
    # set kr and kx for root system or plant
    
    plantParameters.init_conductivities(r = plantModel)
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
    





