""" 
Functions to simplify setup of the scenarios for the INARI project
"""

import sys;
sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../modules/");
sys.path.append("/data");
sys.path.append("../../../CPlantBox/src/python_modules");
sys.path.append("../../../CPlantBox/src/functional/");
sys.path.append("../../../CPlantBox/src/rsml/");
sys.path.append("../../../CPlantBox/src/visualisation/")
sys.path.append("../../../CPlantBox/src/structural/")
sys.path.append("../../../CPlantBox/src/external/")
sys.path.append("../../../CPlantBox/");

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
from functional.plant_conductivities import init_conductivities

from helpfull import *
from weather import *
from XylemPhloemPhotosynthesis import *




def getBiochemParam(s,paramIdx, noAds):    
    paramSet = pd.read_csv('./output_random_rows.csv').iloc[paramIdx].to_dict()
    
    s.mg_per_molC = 12000
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
    s.C_aOLim=1.e-10 # so that microbe community can always regrow
    s.C_aCLim=1.e-10 # so that microbe community can always regrow
    
    if noAds:
        s.CSSmax = 0.
        s.alpha = 0.
    else:
        s.Qmmax = 0.45 * 0.079 # max ratio gOC-gmineral soil, see 10.1016/j.soilbio.2020.107912
        # [g OC / g mineral soil] * [g mineral soil/ cm3 bulk soil] *[ mol C/g C]
        CSSmax_ = s.Qmmax * s.bulkMassDensity_gpercm3*(1/12.011)
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
    return s

def setIC3D(s, paramIdx, ICcc = None):
    return setIC(s, paramIdx, ICcc)

def setIC(s, paramIdx, ICcc = None):
    if ICcc is None:
        paramSet = pd.read_csv('./output_random_rows.csv').loc[paramIdx]
        C_S = paramSet['CS_init'] /s.mg_per_molC## in mol/cm3 water
        C_L = paramSet['CL_init'] /s.mg_per_molC## in mol/cm3 water

        s.CSS2_init,s.ratioInit  = s.getCSS2Init(C_S) #mol C/ cm3 scv
        if rank == 0:
            print('C_S,CSS2_init',C_S,s.CSS2_init,'CSSmax', s.CSSmax ,'ratio', s.ratioInit)
        
        if s.noAds:
            CSS2_init = 0.
            
        unitConversion = 1.0e6 # mol/cm3  => mol/m3 
        addedVar = 1.
        s.CSW_init = C_S * unitConversion
        s.ICcc = np.array([C_S *unitConversion*addedVar,
                           C_L*unitConversion*addedVar,
                            9.16666666666667e-07* unitConversion*addedVar,
                            8.33333333333333e-06* unitConversion*addedVar,
                            8.33333333333333e-07* unitConversion*addedVar,
                            8.33333333333333e-06* unitConversion*addedVar,
                            s.CSS2_init*unitConversion*addedVar,
                           0.])# in mol/m3 water or mol/m3 scv
        if rank == 0:
            print('init s.ICcc', s.ICcc,'s.k_sorp',s.k_sorp,'s.CSSmax',s.CSSmax)
            print('compute method' ,s.CSSmax, C_S , s.theta_init, s.cell_size , s.k_sorp , s.f_sorp)
    else:
        s.ICcc = ICcc
    
    for i in range(s.numSoluteComp):
        #mol/m3 to mol/mol
        molarC = s.ICcc[i] / s.phaseDensity(isDissolved = (i < s.numDissolvedSoluteComp)) 
        s.setParameter( "Soil.IC.C"+str(i+1), str(molarC ))
    return s


def setDefault(s):
    molarMassWat = 18. # [g/mol]
    densityWat = 1. #[g/cm3]
    # [mol/cm3] = [g/cm3] /  [g/mol] 
    molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 
    s.molarDensityWat = molarDensityWat

    s.MaxRelativeShift = 1e-8
    s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift))
    s.setParameter("Problem.verbose", "0")
    s.setParameter("Newton.EnableChop", "true")# force solute mole fraction > 0 and water in possible pressure ranges
    # s.setParameter("Problem.verbose_local_residual", "true")# for debug
    s.setParameter("Flux.UpwindWeight", "1")#very important because we get high solute gradient.
    
    return s

def getSoilTextureAndShape():    
    min_b = np.array([-3./2, -12./2, -41.]) # np.array( [5, 5, 0.] )
    max_b =np.array( [3./2, 12./2, 0.]) #  np.array([-5, -5, -5.])
    cell_number =np.array( [3,12,41]) # 1cm3 #np.array([3,3,3])
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

def setSoilParam(s,paramIdx):
    
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

def DtCSS2_(CSS1, CSW, CSS2):
    return 0.1/(24.*60.*60.) * (CSW - CSS2)

def create_soil_model3D( demoType, 
                    times = None, usemoles = True, dirResults = "",
                        p_mean_ = -100,paramIndx =0,
                     noAds = False, ICcc = None, DtCSS2 = DtCSS2_):
    return create_soil_model( demoType, 
                    times , usemoles , dirResults,
                        p_mean_,paramIndx ,
                     noAds , ICcc )

def create_soil_model( demoType, 
                    times = None, usemoles = True, dirResults = "",
                        p_mean_ = -100,paramIndx =0,
                     noAds = False, ICcc = None, DtCSS2 = DtCSS2_):
    """
        Creates a soil domain from @param min_b to @param max_b with resolution @param cell_number
        soil demoType is fixed and homogeneous 
        domain is periodic (if 2d or 3d)
        initial potentials are linear from @param p_top to @param p_bot
        
        returns soil_model (RichardsWrapper(RichardsSP())) and soil parameter (vg.Parameters)
    """
    if len(dirResults) == 0:
        dirResults = "./results/parallel"+str(max_rank)+"/"
    s = RichardsWrapper(RichardsNCSP(), usemoles)  # water and N solute          
    s.dirResults = dirResults
    
    m3_per_cm3 = 1e-6; # m3/cm3
    cm3_per_m3 = 1e6; # cm3/m3
    
    def getCSS2Init(CSW):
        return  (s.kads * CSW * s.CSSmax)/(s.kads * CSW + s.kdes), (s.kads * CSW)/(s.kads * CSW + s.kdes) 
    
    s.getCSS2Init = getCSS2Init
    
    
    #@see dumux-rosi\cpp\python_binding\solverbase.hh
    soilTextureAndShape = getSoilTextureAndShape() 
    min_b = soilTextureAndShape['min_b']
    max_b = soilTextureAndShape['max_b']
    cell_number = soilTextureAndShape['cell_number']
    s.cell_size = np.prod((max_b - min_b) / cell_number) # cm3 
    s.noAds = noAds
    s.initialize()
    setDefault(s)
    setSoilParam(s,paramIndx)
    s.theta_init =  vg.water_content(p_mean_, s.vg_soil)
    #assert s.theta_init == 0.4
    
    getBiochemParam(s,paramIndx,noAds )
    
    setBiochemParam(s)
    
    setIC3D(s, paramIndx, ICcc)
    s.isPeriodic = True
    s.createGrid(min_b, max_b, cell_number, s.isPeriodic)  # [cm] 
    cell_number_ = cell_number
    cell_number= s.dumux_str(cell_number)#.replace("[", "");cell_number=cell_number.replace("]", "");cell_number=cell_number.replace(",", "");
    s.setParameter( "Soil.Grid.Cells", cell_number)   
    
    s.setParameter("Problem.reactionExclusive", "1")

    s, s.vg_soil = setupOther(s, p_mean_)
    
    if False:
        DtCSS2 =  lambda CSS1, CSW, CSS2: s.kads * CSW * (s.CSSmax * cm3_per_m3 - CSS2) - s.kdes * CSS2
    else:
        # [1/d] * [d/s] * ([mol/cm3 soil scv zone 1] * [cm3/m3] * [mol/m3]/([mol/m3]+ [mol/cm3] * [cm3/m3]) - [mol/m3])
        # mol/m3/s
        def DtCSS2( CSS1, CSW, CSS2):
            return  s.alpha /(24.*60.*60.) * (s.CSSmax * cm3_per_m3 * (CSW/(CSW+s.k_sorp * cm3_per_m3))- CSS2)
        
    s.DtCSS2 = DtCSS2
    
    if rank == 0:
        s.base.printParams()
    return s, s.vg_soil
    
def setShape1D(s,r_in, r_out,length,nCells = 10, doLogarithmic=True):
    
    logbase = 0.5
    s.r_in = r_in
    s.r_out = r_out
    if doLogarithmic:
        s.points = np.logspace(np.log(r_in) / np.log(logbase), np.log(r_out) / np.log(logbase), nCells, base = logbase)
        
    else:
        s.points = [r_in + (r_out - r_in)/(nCells-1) * i for i in range(nCells) ]
        
    s.createGrid1d(s.points)
    nCells -= 1
    s.setParameter( "Soil.Grid.Cells", str(nCells))
    s.setParameter( "Problem.segLength", str(length))
    return s

def setBC1D(s):    
    s.setInnerBC("fluxCyl",  s.win)  # [cm/day] #Y 0 pressure?
    s.setOuterBC("fluxCyl",0.)
    
    s.setParameter( "Soil.BC.Bot.C1Type", str(3))
    s.setParameter( "Soil.BC.Top.C1Type", str(3))
    s.setParameter( "Soil.BC.Bot.C1Value", str(s.exuds_in)) 
    s.setParameter( "Soil.BC.Top.C1Value", str(0.)) 


    s.setParameter( "Soil.BC.Bot.C2Type", str(3))
    s.setParameter( "Soil.BC.Top.C2Type", str(3))
    s.setParameter( "Soil.BC.Bot.C2Value", str(s.exudl_in)) 
    s.setParameter( "Soil.BC.Top.C2Value", str(0.)) 

    for i in range(s.numFluidComp, s.numComp):
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(3))
        s.setParameter( "Soil.BC.Top.C"+str(i)+"Type", str(3))
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(0)) 
        s.setParameter( "Soil.BC.Top.C"+str(i)+"Value", str(0 )) 

def setupOther(s, p_mean_):
    s.p_meanInit = p_mean_
    # BC
    if s.dimWorld == 1:
        s.setParameter("Soil.IC.P", s.dumux_str(p_mean_))
        s.setInnerBC("fluxCyl",  s.win)  # [cm/day] #Y 0 pressure?
        s.setOuterBC("fluxCyl", 0.)
    if s.dimWorld == 3:
        s.setTopBC("noFlux")
        s.setBotBC("noFlux") #in acc. with Jorda et al. (2022), however, they assume inflow if h>0
        indxFluxSolute = 2
    
        for i in range(1, s.numComp):
            s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(indxFluxSolute))
            s.setParameter( "Soil.BC.Top.C"+str(i)+"Type", str(indxFluxSolute))
            s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(0)) 
            s.setParameter( "Soil.BC.Top.C"+str(i)+"Value", str(0 ))       
        
    
    # IC
    if True:
        if s.dimWorld == 3:
            if isinstance(p_mean_,(int,float)):
                #print('set pmean float',p_mean_)
                s.setHomogeneousIC(p_mean_, equilibrium = True)  # cm pressure head
            elif isinstance(p_mean_,type(np.array([]))):
                pass
            else:
                print(type(p_mean_))
                raise Exception
        


    s.initializeProblem()
    
     
    if isinstance(p_mean_,(int,float)):
        pass
    elif isinstance(p_mean_,type(np.array([]))):
        s.setInitialConditionHead(p_mean_)
    else:
        print(type(p_mean_))
        raise Exception
    
        
    s.wilting_point = -15000
    s.setCriticalPressure(s.wilting_point)  # for boundary conditions constantFlow, constantFlowCyl, and atmospheric
    s.ddt = 1.e-5  # [day] initial Dumux time step
    return s, s.vg_soil

def set_all_sd(rs, s):
    """ # sets all standard deviation to a percantage, i.e. value*s """
    for p in rs.getRootRandomParameter():
        p.a_s = p.a * s
        p.lbs = p.lb * s
        p.las = p.la * s
        p.lns = p.ln * s
        p.lmaxs = p.lmax * s
        p.rs = p.r * s
        p.thetas = p.theta * s
        p.rlts = p.rlt * s  # no used
        p.ldelays = p.ldelay * s
    seed = rs.getRootSystemParameter()  # SeedRandomParameter
    seed.firstBs = seed.firstB * s
    seed.delayBs = seed.delayB * s
    seed.maxBs = seed.maxB * s
    seed.firstSBs = seed.firstSB * s
    seed.delaySBs = seed.delaySB * s
    seed.delayRCs = seed.delayRC * s
    seed.nCs = seed.nCs * s
    seed.nzs = seed.nzs * s
    # todo seed position s


    
def create_mapped_plant( mode,initSim, soil_model, fname, path, 
                stochastic = False, mods = None, plantType = "plant",l_ks_ = "dx_2",
                        usemoles = True, limErr1d3d = 1e-11, spellData =None):
    """ loads a rmsl file, or creates a rootsystem opening an xml parameter set,  
        and maps it to the soil_model """
    
    soilTextureAndShape = getSoilTextureAndShape() 
    min_b = soilTextureAndShape['min_b']
    max_b = soilTextureAndShape['max_b']
    cell_number = soilTextureAndShape['cell_number']
    
    if fname.endswith(".rsml"):
        r = XylemFluxPython(fname)
        if plantType == "plant":
            from rhizo_modelsPlant import RhizoMappedSegments  # Helper class for cylindrical rhizosphere models
        else:
            from rhizo_modelsRS import RhizoMappedSegments  # Helper class for cylindrical rhizosphere models
        rs = RhizoMappedSegments(  mode, soil_model,  usemoles, seedNum = seed, 
                                 limErr1d3dAbs = limErr1d3d, l_ks=l_ks_)
    elif fname.endswith(".xml"):
        seed = 1
        weatherInit = weather(initSim,spellData)
        if plantType == "plant":
            from rhizo_modelsPlant import RhizoMappedSegments  # Helper class for cylindrical rhizosphere models
        else:
            from rhizo_modelsRS import RhizoMappedSegments  # Helper class for cylindrical rhizosphere models
        rs = RhizoMappedSegments(  mode, soil_model,  usemoles, seedNum = seed, 
                                 limErr1d3dAbs = limErr1d3d, l_ks=l_ks_)

        rs.setSeed(seed)
        rs.readParameters(path + fname)
        #if not stochastic:
        #    set_all_sd(rs, 0.)
            
        rs.setGeometry(pb.SDF_PlantBox( max_b[0]-min_b[0],  max_b[1]-min_b[1], max_b[2]-min_b[2]))
        rs.initialize(verbose = False)#stochastic = False)
        rs.simulate(initSim,verbose= False)
        if plantType == "plant":
            r = PhloemFluxPython(rs,psiXylInit = -659.8 - min_b[2],ciInit = weatherInit["cs"]*0.5) 
        else:
            r = XylemFluxPython(rs)
    r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), 
                            pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), 
                            cut = False, noChanges = True)
    
    picker = lambda x, y, z: soil_model.pick_([x,y, z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
    r.rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
    if plantType == "plant":    
        r = init_conductivities(r)
        r = phloemParam(r, weatherInit)
        rs.set_phloem_flux(r)
        #r.test()
        return rs, r
    else:
        r = init_conductivities_const(r)
        return rs, r
    





