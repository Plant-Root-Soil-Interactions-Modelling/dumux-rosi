""" 
Functions to simplify setup of the scenarios for the INARI project
"""

import sys;
sys.path.append("../../../build-cmake/cpp/python_binding/");

import numpy as np
import pandas as pd
import timeit
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

from rosi_richards22c_cyl import RichardsNCCylFoam # C++ part (Dumux binding)
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial model (a single cylindrical model is not allowed to run in parallel)
from rosi_richards22c import RichardsNCSPILU as RichardsNCSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
from functional.phloem_flux import PhloemFluxPython  # root system Python hybrid solver


import plantbox as pb  # CPlantBox
import functional.van_genuchten as vg
from functional.xylem_flux import *
from datetime import *
import plantParameters
from helpfull import *
from PhloemPhotosynthesis import *
import weatherFunctions




def getBiochemParam(s,paramIdx):    
    """ define TraiRhizo biochemical parameters 
        @param: the dumux soil object
        @ param: index of the TraiRhizo parameter set to use
        
    """
    # file containing the TraiRhizo parameter sets
    paramSet = pd.read_csv('./output_random_rows.csv').iloc[paramIdx].to_dict() # select one specific parameter set from index
    s.molarMassC = 12.011
    s.molarMassN = 14
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
    
    s.f_Im = 4.1*1e-5*(24*3600)*1e-6 #[cm3/mol/d] or [1/d]
    s.f_Min = 2.3*1e-5*(24*3600)*1e-6 #[1/d]
    s.v_maxNH4 = 4.4*1e4*(24*3600) * 1e-3 #[cm3/mol/d] or [1/d]
    s.v_maxNO3 = 4.4*1e4*(24*3600) * 1e-3 #[1/d]
    s.F_PNU_NO3_max = 0.1333 * 24 * 3600 * 1e-3
    s.F_PNU_NH4_max = 0.1333 * 24 * 3600 * 1e-3
    
    if s.noAds:
        s.CSSmax = 0.
        s.alpha = 0.
        s.kadsN = 0.
        s.kdesN = 0.
    else:
        s.Qmmax = 0.45 * 0.079 # max ratio gOC-gmineral soil, see 10.1016/j.soilbio.2020.107912
        # [g OC / g mineral soil] * [g mineral soil/ cm3 bulk soil] *[ mol C/g C]
        CSSmax_ = s.Qmmax * s.bulkMassDensity_gpercm3*(1/s.molarMassC)
        s.CSSmax = CSSmax_ # mol C/cm3 bulk soil
        #s.CSSmax = s.Qmmax * s.bulkDensity_m3 / 1e6 # mol OC/mol soil * [mol soil/m3] * [m3/cm3] =  mol/cm3
        s.kadsN = 12.
        s.kdesN = 0.0128

        
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
    s.mucilCAffectsW = "true"
    return s

def setBiochemParam(s):
    """ send the TraiRhizo biochemical parameters to dumux
        @param: the dumux soil object
    """

    s.setParameter( "Soil.mucilCAffectsW", str(s.mucilCAffectsW))    
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
    s.setParameter("Soil.kadsN", str(s.kadsN)) #[cm3/mol/d] or [1/d]
    s.setParameter("Soil.kdesN", str(s.kdesN)) #[1/d]
    s.setParameter("Soil.f_Im", str(s.f_Im)) #[cm3/mol/d] or [1/d]
    s.setParameter("Soil.f_Min", str(s.f_Min)) #[1/d]
    s.setParameter("Soil.v_maxNH4", str(s.v_maxNH4))
    s.setParameter("Soil.v_maxNO3", str(s.v_maxNO3))
    s.setParameter("Soil.F_PNU_NO3_max", str(s.F_PNU_NO3_max))
    s.setParameter("Soil.F_PNU_NH4_max", str(s.F_PNU_NH4_max))
    
    
    if s.dimWorld == 1:
        s.setParameter("Soil.BC.dzScaling", "1")
        s.setParameter("Problem.verbose", "0")
        s.setParameter("Problem.reactionExclusive", "0")
    
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
        N_SL = C_S / 10. # random value
        N_LL = C_L / 10. # random value
        NH4 = 1e-4  * 1e-3# random value
        NO3 = 1e-4  * 1e-3# random value
        NH4S = 1e-4  * 1e-3# random value
        k_CNobj = 8.
        # concentraiton of adsobed C_S
        s.CSS_init  = getCSS(s, C_S) #mol C/ cm3 scv
        N_SS = s.CSS_init / 10.
            
        unitConversion = 1.0e6 # mol/cm3  => mol/m3 
        addedVar = 1. * float(s.doSoluteFlow) # empirical factor
        # s.CSW_init = C_S * unitConversion
        
        # TODO: update values
        s.ICcc = np.array([C_S *unitConversion*addedVar, # small solutes
                           0., # mucilage
                           C_L*unitConversion*addedVar, # heavy org mol
                           N_SL *unitConversion*addedVar, # small N solutes
                           N_LL *unitConversion*addedVar, # large N solutes
                           NH4 *unitConversion*addedVar, 
                           NO3 *unitConversion*addedVar, 
                            9.16666666666667e-07* unitConversion*addedVar, # active oligotrophes
                            8.33333333333333e-06* unitConversion*addedVar, # dormant oligo
                            8.33333333333333e-07* unitConversion*addedVar, # active copiotrophes
                            8.33333333333333e-06* unitConversion*addedVar, # dormant copio
                            s.CSS_init*unitConversion*addedVar,            # adsorbed C
                           0. ,      # CO2
                            9.16666666666667e-07* unitConversion*addedVar * k_CNobj, # active oligotrophes
                            8.33333333333333e-06* unitConversion*addedVar * k_CNobj, # dormant oligo
                            8.33333333333333e-07* unitConversion*addedVar * k_CNobj, # active copiotrophes
                            8.33333333333333e-06* unitConversion*addedVar * k_CNobj, # dormant copio
                            NH4S * unitConversion*addedVar, 
                            N_SS * unitConversion*addedVar, 
                            0. # N2
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

    # low MaxRelativeShift == higher precision in dumux
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
    s.MaxTimeStepDivisions = 100
    s.setParameter("Newton.MaxTimeStepDivisions",
                     str( s.MaxTimeStepDivisions) )  
    s.MaxSteps = 100
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
    cell_number =np.array([2,1,1])# np.array( [3,12,40]) #np.array( [1,1,1]) # 1cm3
    # #np.array([3,3,3])
    solidDensity = 2650 # [kg/m^3 solid] #taken from google docs TraiRhizo
    solidMolarMass = 60.08e-3 # [kg/mol] 
    # theta_r, theta_s, alpha, n, Ks
    soilVG = [0.08, 0.43, 0.04, 1.6, 50]
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


def create_soil_model3D(results_dir ,
                        p_mean_ = -100,paramIndx =0,
                     noAds = False, ICcc = None, doSoluteFlow = True):
    return create_soil_model( results_dir,
                        p_mean_,paramIndx ,
                     noAds , ICcc , doSoluteFlow)

def create_soil_model( results_dir ,
                        p_mean_ = -100,paramIndx =0,
                     noAds = False, ICcc = None, doSoluteFlow = True,
                       
                       doBioChemicalReaction=True,
                     MaxRelativeShift = 1e-8, doOld = False):
    """
        Creates a soil domain from @param min_b to @param max_b with resolution @param cell_number
        homogeneous domain 
        initial potentials are linear and mean potential is @p_mean_
        @ param: noAds: turn off adsorption?
        @param: paramIndx index of the TraiRhizo parameter set to use
        @param: ICcc (optional) initial concentraiton values for the solute components
        returns soil_model (RichardsWrapper(RichardsSP())) and soil parameter (vg.Parameters)
    """
        
    s = RichardsWrapper(RichardsNCSP())  # water and N solute          
    s.results_dir = results_dir   
    s.pindx = paramIndx
    
    s.MaxRelativeShift = MaxRelativeShift # 1e-10
    s.MaxRelativeShift_1DS = MaxRelativeShift
    
    soilTextureAndShape = getSoilTextureAndShape() 
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
    getBiochemParam(s,paramIndx)
    setBiochemParam(s)
    setIC3D(s, paramIndx, ICcc)
    s.isPeriodic =True 
    s.createGrid(min_b, max_b, cell_number, s.isPeriodic)  # [cm] 
    if doOld:
        s.isPeriodic =False
    s = setupOther(s, p_mean_)
    
    #if rank == 0:
    #    s.base.printParams()
    
    # just print once as will not change during simulation
    write_file_array("cellVol", np.array(s.getCellVolumes()), directory_ =s.results_dir) # cm3 
    write_file_array("cellIds", np.array(s.cellIndices), directory_ =s.results_dir) # cm3
    
    cellcenters = s.getCellCenters()
    cellcentersX = []
    cellcentersY = []
    cellcentersZ = []
    for sub_array in cellcenters:
        cellcentersX.append(sub_array[0])
        cellcentersY.append(sub_array[1])
        cellcentersZ.append(sub_array[2])

    write_file_array("cellcentersX", np.array(cellcentersX), directory_ =results_dir) # cm3
    write_file_array("cellcentersY", np.array(cellcentersY), directory_ =results_dir) # cm3
    write_file_array("cellcentersZ", np.array(cellcentersZ), directory_ =results_dir) # cm3
    
    return s
   
def setCylBC(cyl):
    # use this vector to know if it is set uptake/release, MM...
    cyl.typeBC = np.full(cyl.numSoluteComp, 3)
    cyl.typeBC[5] = 12
    cyl.typeBC[6] = 11
    
    for j in range( 1, cyl.numComp):
        cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Type", str(cyl.typeBC[j-1]))
        cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Type", str(3))
        cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Value", str(0)) 
        cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Value", str(0 )) 

def setupOther(s, p_mean_):
    """ define remaining soil parameters """
    s.p_meanInit = p_mean_
    
    # without counting the water element
    s.Cidx = np.array([1,2,3,8,9,10,11,12,13]) -1
    s.Nidx = np.array([4,5,6,7,14,15,16,17,18,19,20])-1
    s.PIFidx = np.array([1,2,4]) -1  # solutes where we have a proposed inner flux
    s.RIFidx = np.array([1,2,4,6,7]) -1 # solutes where we have a realised inner flux
    
    s.setParameter("SpatialParams.Temperature","293.15") # todo: redefine at each time step?
    
    # for boundary conditions constantFlow, constantFlowCyl, and atmospheric
    s.wilting_point = -15000 
    s.ddt = 1.e-3  # [day] initial Dumux time step
    
    # BC
    if s.dimWorld == 1: # 1d model
        s.setParameter("Soil.IC.P", s.dumux_str(p_mean_))
        s.setInnerBC("fluxCyl",  0.)  # [cm/day]
        s.setOuterBC("fluxCyl", 0.)
        
        setCylBC(s)
        s.maxDt =  250./(3600.*24.)
        
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
        
    
        if isinstance(p_mean_,(int,float)):
            #print('set pmean float',p_mean_)
            s.setHomogeneousIC(p_mean_, equilibrium = True)  # cm pressure head
        elif isinstance(p_mean_,type(np.array([]))):
            pass
        else:
            print(type(p_mean_))
            raise Exception
        s.maxDt =  250./(3600.*24.)
        #s.maxDt_1DS = s.maxDt # [s], lower maxDt for 1D models
        
    s.initializeProblem(s.maxDt)    
    s.setCriticalPressure(s.wilting_point)

    s.eps_regularization = None
    if s.eps_regularization is not None:
        s.setRegularisation(s.eps_regularization, s.eps_regularization) # needs to be low when using sand parameters. 
            

    if s.dimWorld == 3:# 3d model
    
        if isinstance(p_mean_,(int,float)):
            pass
        elif isinstance(p_mean_,type(np.array([]))):
            s.setInitialConditionHead(p_mean_)
        else:
            print(type(p_mean_))
            raise Exception
    
        s.bulkMassErrorWater_rel = 0.
        s.bulkMassErrorWater_relLim = 0.    
        totCN3dInit_ = s.getTotCNContent_each().sum(axis=1)
        s.totCN3dInit = np.array([sum(totCN3dInit_[s.Cidx]),sum(totCN3dInit_[s.Nidx])]) # sum(s.getTotCNContent()) # mol    
        # initial soil water and solute content
        cell_volumes = s.getCellVolumes()  # cm3
        s.buWSoilInit = sum(np.multiply(np.array(s.getWaterContent()), cell_volumes)) # cm3 water
    return s

        
    
def create_mapped_plant(initSim, soil_model, fname, path, 
                stochastic = False, 
                        doPhloemFlow = True,limErr1d3d = 1e-11, spellData =None):
    """ loads a rmsl file, or creates a rootsystem opening an xml parameter set,  
        and maps it to the soil_model """
    from rhizo_modelsPlant import RhizoMappedSegments  # Helper class for cylindrical rhizosphere models
    soilTextureAndShape = getSoilTextureAndShape() 
    min_b = soilTextureAndShape['min_b']
    max_b = soilTextureAndShape['max_b']
    cell_number = soilTextureAndShape['cell_number']
    
    if fname.endswith(".rsml"):
        plantModel = XylemFluxPython(fname)
        
        perirhizalModel = RhizoMappedSegments(  mode, soil_model,  seedNum = seed, 
                                 limErr1d3dAbs = limErr1d3d)
    elif fname.endswith(".xml"):
        seed = 1
        weatherInit = weatherFunctions.weather(initSim,0, spellData)
        perirhizalModel = RhizoMappedSegments(soilModel = soil_model, 
                                              ms = pb.MappedPlant(),
                                 limErr1d3dAbs = limErr1d3d,
                                 RichardsNCCylFoam = RichardsNCCylFoam)

        perirhizalModel.ms.setSeed(seed)
        perirhizalModel.ms.readParameters(path + fname)
        if soil_model.isPeriodic:
            perirhizalModel.ms.setGeometry(pb.SDF_PlantBox(np.inf, np.inf, max_b[2]-min_b[2]))
        else:
            perirhizalModel.ms.setGeometry(pb.SDF_PlantBox( max_b[0]-min_b[0],  max_b[1]-min_b[1], max_b[2]-min_b[2]))
        perirhizalModel.ms.initialize(verbose = False)
        perirhizalModel.ms.simulate(initSim,verbose= False)
        if doPhloemFlow:
            plantModel = PhloemFluxPython(perirhizalModel.ms,psiXylInit = -659.8 - min_b[2],ciInit = weatherInit["cs"]*0.5) 
        else:
            plantModel = XylemFluxPython(perirhizalModel.ms)
            
    perirhizalModel.ms.constantLoc = True
    #plantModel.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), 
    #                        pb.Vector3d(max_b[0], max_b[1], max_b[2]),
    #                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), # not used
    #                        cut = False, # do not cut plant segments according to the voxel size
    #                        noChanges = True) # segments remain in the voxel they first appeared in
    
    #  function that return the index of a given position in the soil grid 
    picker = lambda x, y, z: soil_model.pick_([x,y, z])  
    # maps segments, maps root segements and soil grid indices to each other in both directions
    plantModel.rs.setSoilGrid(picker)
    
    plantModel.wilting_point = -15000.
    # set kr and kx for root system or plant
    
    plantParameters.init_conductivities(r = plantModel)
    if doPhloemFlow:   
        plantModel = plantParameters.phloemParam(plantModel, weatherInit)
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
    plantModel.seg_fluxesCumul = np.array([np.array([]) for jj in range(len(soil_model.RIFidx)+1)])
    #plantModel.seg_fluxes1Cumul = np.array([])
    #plantModel.seg_fluxes2Cumul = np.array([])
    
    return perirhizalModel, plantModel
    





