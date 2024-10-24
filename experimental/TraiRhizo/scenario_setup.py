""" 
Functions to simplify setup of the scenarios for the INARI project
"""

import sys;
sys.path.append("../../build-cmake/cpp/python_binding/");
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
from rosi_richards10c import RichardsNCSPILU as RichardsNCSP 
#from rosi_richards10c import RichardsNCSPnum as RichardsNCSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
from functional.phloem_flux import PhloemFluxPython  # root system Python hybrid solver


import plantbox as pb  # CPlantBox
import van_genuchten as vg
#import evapotranspiration as evap
from functional.xylem_flux import *
from datetime import *
from functional.plant_conductivities import init_conductivities

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 164
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
oldCSV = False


def write_file_float(name, data, directory_, allranks = False):
    if (rank == 0) or allranks:
        name2 = directory_+ name+ '.txt'
        #print('write_file_float',name2, allranks)
        with open(name2, 'a') as log:
            log.write(repr( data)  +'\n')
        
def write_file_array(name, data, space =",", directory_ ="./results/", fileType = '.txt', allranks = False ):
    np.array(data).reshape(-1)
    try:
        if (rank == 0) or allranks:
            name2 = directory_+ name+ fileType
            #print('write_file_array',name2)
            with open(name2, 'a') as log:
                log.write(space.join([num for num in map(str, data)])  +'\n')
    except:
        print('write_file_array',name, data,data.shape)
        raise Exception

def vg_SPP(i = int(1)):
    """ Van Genuchten parameter, called by maize()  """
        
    soil = {}
    soil[0] = [0.045,np.nan, 0.04, 1.6, 50]# [0.041, 0.494, 0.0256, 1.49, 245]
    soil[1] = [0.03, np.nan, 0.038, 2, 1864]
    return soil[i]


def maize_SPP(soil_= "loam"):
    """ parameters for maize simulation """
    if soil_ == "loam":
        i = 0
    else:
        i = 1
    soil = vg_SPP(i)
    min_b = [-5., -5., -10.] 
    max_b = [5., 5., 0.] 
    cell_number = [4,4,5]#
    area = 20 * 44  # cm2
    
    # min_b = [-5., -5, -10.] 
    # max_b = [5., 5, 0.]
    # cell_number = [3, 3, 2]
    # area = 10*10  # cm2

    Kc_value_ = {}
    Kc_value_[0] = np.array([1,1,1,1.2,1.2,1.2])
    Kc_value_[1] = np.array([1,1,1,1.07,1.2,1.2])
    Kc_value = Kc_value_[0] #2019, 2020
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

    #plt.plot(np.linspace(0,287,288), Kc)
            
    return soil, min_b, max_b, cell_number, area, Kc

def exudation_rates():
    #cst
    kex = np.array([[1e-2],[0 ],[0]])
    #kex = np.array([[0., 2.], [0., 0.]])
    #defined per organ types
    #kex = np.array([[1e-2], [0.]])
    return kex

def init_conductivities_const(r, kr_const = 1.8e-4, kx_const = 0.1):
    """ Hydraulic conductivities  kr [1/day], kx [cm3/day] """
    r.setKr([0, kr_const, kr_const, kr_const, kr_const, kr_const])
    r.setKx([1.e3, kx_const, kx_const, kx_const, kx_const, kx_const])
    return r


def init_maize_conductivities(r, skr = 1., skx = 1.):
    """ Hydraulic conductivities for maize modified from Couvreur et al. (2012), originally from Doussan et al. (1998) """ #[age, value]

    kr00 = np.array([[0., 0.]])  # artificial shoot
    kx00 = np.array([[0., 1.e3]])  # artificial shoot

    kr0 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.000181], [8., 0.000181], [10, 0.0000648], [18, 0.0000648], [25, 0.0000173], [300, 0.0000173]])
    kr1 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.000181], [10., 0.000181], [16, 0.0000173], [300, 0.0000173]])

    kx0 = np.array([[0., 0.000864], [5., 0.00173], [12., 0.0295], [15., 0.0295], [20., 0.432], [300., 0.432]])
    kx1 = np.array([[0., 0.000864], [5., 0.0000864], [10., 0.0000864], [12., 0.0006048], [20., 0.0006048], [23., 0.00173], [300., 0.00173]])

    kr01 = np.minimum(skr * kr0[:, 1], 1.)
    kr11 = np.minimum(skr * kr1[:, 1], 1.)
    r.setKrTables([kr00[:, 1], kr01, kr11, kr11, kr01, kr01],
                  [kr00[:, 0], kr0[:, 0], kr1[:, 0], kr1[:, 0], kr0[:, 0], kr0[:, 0]])
    kx01 = np.minimum(skx * kx0[:, 1], 1.)
    kx11 = np.minimum(skx * kx1[:, 1], 1.)
    r.setKxTables([kx00[:, 1], kx01, kx11, kx11, kx01, kx01],
                  [kx00[:, 0], kx0[:, 0], kx1[:, 0], kx1[:, 0], kx0[:, 0], kx0[:, 0]])


def getBiochemParam(s,paramIdx, noAds):
    s.numFluidComp = 2
    # s.numComp = 8
    #paramSet = pd.read_csv('./TraiRhizoparam_ordered.csv').iloc[paramIdx].to_dict()
    if oldCSV:
        paramSet = pd.read_csv('./TraiRhizoparam_ordered.csv').iloc[paramIdx].to_dict()
    else:
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

    s.m_maxC =  paramSet['m_max,C']#1/d
    s.m_maxO =  paramSet['m_max,O']#1/d
    s.k_decay2 =  paramSet['p_L'] # -

    s.micro_maxC =  paramSet['u_max,C']#1/d
    s.micro_maxO =  paramSet['u_max,O']#1/d

    s.v_maxL = paramSet['v_max,L'] #1/d
    s.k_decay = paramSet['Y'] # -
    s.k_growthC = paramSet['Y_C'] # -
    s.k_growthO = paramSet['Y_O'] #- 
    # can cause issue
    s.k_sorp = paramSet['k_sorp'] /s.mg_per_molC# mg C/cm 3 soil solution =>  mol C/cm3 soil solution
    if s.css1Function == 8:
        s.k_sorp *= paramSet['theta'] * s.cell_size # mol C/cm3 soil solution => mol
    
    s.alpha = 0.1 # -
    s.f_sorp = 0# 0.5
    s.k_phi = 0.1
    s.C_aOLim=1.e-10 # so that microbe community can always regrow
    s.C_aCLim=1.e-10 # so that microbe community can always regrow
    
    if noAds:
        s.CSSmax = 0.
        s.alpha = 0.
    else:
        s.CSSmax = paramSet['CSS_max']/s.mg_per_molC # mol C/cm3 soil zone 1 to mol C/cm3 soil 
        if s.css1Function == 8: #change cssmax to content
            s.CSSmax *=  s.cell_size * s.f_sorp # mol C
        
    kads = 7.07e+02 # m3/kgC/yr, see 10.1016/j.soilbio.2020.107912, A.3
    yr_per_d = 1/365 # [yr/d]
    m3_per_cm3 = 1e-6; # m3/cm3
    cm3_per_m3 = 1e6; # cm3/m3
    
    # [kg/g] * [g/mol] = kg/mol
    mol_per_kgC = (1/1000) * s.molarMassC
    # [m3/kgC/yr] * [yr/d] * [cm3/m3] * [kgC/mol] = [cm3/mol/d]
    s.kads = kads * yr_per_d * cm3_per_m3 * mol_per_kgC
    if s.css1Function == 3: #change kads to 1/d:
        C_S = paramSet['CS_init'] /s.mg_per_molC## in mol/cm3 water
        s.kads *= C_S # [cm3/mol/d] => [1/d]
    kdes =  1.63e+03 # [1/yr] see 10.1016/j.soilbio.2020.107912, A.3
    s.kdes = kdes * yr_per_d
    return s

def setBiochemParam(s):

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
    # s.setHomogeneousIC(p_mean_, equilibrium = True) #-97.5)  # cm pressure head
    if ICcc is None:
        #paramSet = pd.read_csv('./TraiRhizoparam_ordered.csv').loc[paramIdx]
        if oldCSV:
            paramSet = pd.read_csv('./TraiRhizoparam_ordered.csv').loc[paramIdx]
        else:
            paramSet = pd.read_csv('./output_random_rows.csv').loc[paramIdx]
        C_S = paramSet['CS_init'] /s.mg_per_molC## in mol/cm3 water
        C_L = paramSet['CL_init'] /s.mg_per_molC## in mol/cm3 water

        s.CSS2_init,s.ratioInit  = s.getCSS2Init(C_S)#s.CSSmax * (C_S/(C_S+ s.k_sorp))#mol C/ cm3 scv # * (1 - s.f_sorp)
        if rank == 0:
            print('C_S,CSS2_init',C_S,s.CSS2_init,'CSSmax', s.CSSmax ,'ratio', (C_S/(C_S+ s.k_sorp)))
        
        if s.css1Function == 8: #change cssmax to content
            #mol C/ cm3 scv zone 2 = mol C * [mol C / cm3 water ] * [cm3 water / cm3 soil] * [cm3 soils]/ [mol] /[ cm3 soil ] * [cm3 soil / cm3 soil zone 1] #at init css2_zone2 == css1_zone1
            CSS2_init_zone2 = s.CSSmax * C_S * s.theta_init * s.cell_size / s.k_sorp / s.cell_size / s.f_sorp #
            #mol C/ cm3 scv zone 2 to mol C/ cm3 scv
            CSS2_init = CSS2_init_zone2 * (1 - s.f_sorp)
        if s.noAds:
            CSS2_init = 0.
        #if s.css1Function == 3: # only pde
        #    CSS2_init = C_S * (1 - s.f_sorp)
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
    else:
        s.ICcc = ICcc
    
    for i in range(s.numComp):
        #mol/m3 to mol/mol
        # print('C',i+1,s.ICcc[i], s.phaseDensity(isDissolved = (i < s.numFluidComp)),molarC )
        molarC = s.ICcc[i] / s.phaseDensity(isDissolved = (i < s.numFluidComp)) 
        # 
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
    
    s.molarMassC = 12.011
    return s

def setSoilParam(s,paramIdx):
    # paramSet = pd.read_csv('./TraiRhizoparam_ordered.csv').loc[paramIdx]
    if oldCSV:
        paramSet = pd.read_csv('./TraiRhizoparam_ordered.csv').loc[paramIdx]
    else:
        paramSet = pd.read_csv('./output_random_rows.csv').loc[paramIdx]
    s.bulkDensity =  paramSet['ro_B']*1000 #g/cm3 => kg/m3
    s.solidDensity = 2650 # [kg/m^3 solid] #taken from google docs TraiRhizo
    s.solidMolarMass = 60.08e-3 # [kg/mol] 

    # theta_r, theta_s, alpha, n, Ks
    s.soil = [0.045, np.nan, 0.04, 1.6, 50]

    s.soil[1] = 1- s.bulkDensity/s.solidDensity #== size of air volume
    s.vg_soil = vg.Parameters(s.soil) 
    # [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
    s.solidMolDensity = s.solidDensity/s.solidMolarMass
    # [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
    s.bulkDensity_m3 = s.solidMolDensity*(1.- s.vg_soil.theta_S)

    s.setParameter( "Soil.MolarMass", str(s.solidMolarMass))
    s.setParameter( "Soil.solidDensity", str(s.solidDensity))
    s.setVGParameters([s.soil])
    
    return s

def DtCSS2_(CSS1, CSW, CSS2):
    return 0.1/(24.*60.*60.) * (CSW - CSS2)

def create_soil_model3D(soil_type, year, soil_, min_b , max_b , cell_number, demoType, 
                    times = None, 
                        net_inf = None, usemoles = True, dirResults = "",
                      #"./results/parallel"+str(max_rank)+"/",
                        p_mean_ = -100,css1Function = 0,paramIndx =0,
                     noAds = False, ICcc = None, DtCSS2 = DtCSS2_):
    create_soil_model(soil_type, year, soil_, min_b , max_b , cell_number, demoType, 
                    times , net_inf , usemoles , dirResults,
                        p_mean_,css1Function,paramIndx ,
                     noAds , ICcc )
def create_soil_model(soil_type, year, soil_, min_b , max_b , cell_number, demoType, 
                    times = None, 
                        net_inf = None, usemoles = True, dirResults = "",
                      #"./results/parallel"+str(max_rank)+"/",
                        p_mean_ = -100,css1Function = 0,paramIndx =0,
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
        if rank == 0:
            print('getCSS2Init','ccsmax',s.CSSmax * cm3_per_m3, 
              'ratio',(CSW/(CSW+ s.k_sorp)),
              'css2init',s.CSSmax * (CSW/(CSW+ s.k_sorp))* cm3_per_m3)
        if s.css1Function== 3:
            return s.CSSmax * (CSW/(CSW+ s.k_sorp)), (CSW/(CSW+ s.k_sorp))
        elif s.css1Function== 9:
            return  (s.kads * CSW * s.CSSmax)/(s.kads * CSW + s.kdes), (s.kads * CSW)/(s.kads * CSW + s.kdes) 
        else:
            raise Exception
    s.getCSS2Init = getCSS2Init# lambda CSW: s.CSSmax * (CSW/(CSW+ s.k_sorp))#mol C/ cm3 scv
    #s.getCSS2Init(C_S_W)#s.CSSmax * (C_S/(C_S+ s.k_sorp))#mol C/ cm3 scv # * (1 - s.f_sorp)
    
    
    #@see dumux-rosi\cpp\python_binding\solverbase.hh
    s.cell_size = np.prod((max_b - min_b) / cell_number) # cm3 
    s.noAds = noAds
    s.initialize(doMPI_=False)
    setDefault(s)
    setSoilParam(s,paramIndx)
    s.theta_init =  vg.water_content(p_mean_, s.vg_soil)
    
    s.css1Function = css1Function
    s.setParameter( "Soil.css1Function", str(s.css1Function))
    getBiochemParam(s,paramIndx,noAds )
    
    setBiochemParam(s)
    
    setIC3D(s, paramIndx, ICcc)
    s.createGrid(min_b, max_b, cell_number, False)  # [cm] 
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
            if False :
                print('in DtCSS2', CSS1, CSW, 'CSS2',CSS2,'oldcss2',s.CSS2_init * cm3_per_m3,
                  's.CSSmax',s.CSSmax * cm3_per_m3 , 
                  'new ratio',(CSW/(CSW+s.k_sorp * cm3_per_m3)),
                  'newCSS2eq', (s.CSSmax * cm3_per_m3 * (CSW/(CSW+s.k_sorp * cm3_per_m3))),
                  'diffCSS2', (s.CSSmax * cm3_per_m3 * (CSW/(CSW+s.k_sorp * cm3_per_m3))- CSS2),
                  'ff',s.CSSmax * cm3_per_m3 * (CSW/(CSW+s.k_sorp * cm3_per_m3)) - s.CSS2_init* cm3_per_m3,
                         'old vs new',CSS2 - s.CSS2_init * cm3_per_m3,
                  CSW - s.CSW_init,
                  (CSW/(CSW+s.k_sorp * cm3_per_m3)) - s.ratioInit )
            return  s.alpha /(24.*60.*60.) * (s.CSSmax * cm3_per_m3 * (CSW/(CSW+s.k_sorp * cm3_per_m3))- CSS2)
        #DtCSS2 =  lambda CSS1, CSW, CSS2:
        
    s.DtCSS2 = DtCSS2
    #s.setComputeDtCSS2(DtCSS2)  # maps segments
    return s, s.vg_soil
    
def setShape1D(s,r_in, r_out,length,nCells = 10, doLogarithmic=True):
    
    logbase = 0.5
    #s.l = length#length in cm
    #doLogarithmic = True
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
    
    #s.setOuterBC("constantFluxCyl", 0.)  #  [cm/day]
    #s.setInnerBC("constantFluxCyl", s.win)  #  [cm/day]
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

    for i in range(s.numFluidComp + 1, s.numComp+1):
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(3))
        s.setParameter( "Soil.BC.Top.C"+str(i)+"Type", str(3))
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(0)) 
        s.setParameter( "Soil.BC.Top.C"+str(i)+"Value", str(0 )) 

def setupOther(s, p_mean_):
    s.p_meanInit = p_mean_


    #cell_number = str(cell_number) 
    # BC
    if s.dimWorld == 1:
        s.setParameter("Soil.IC.P", s.dumux_str(p_mean_))
        s.setInnerBC("fluxCyl",  s.win)  # [cm/day] #Y 0 pressure?
        s.setOuterBC("fluxCyl", 0.)
    if s.dimWorld == 3:
        s.setTopBC("noFlux")
        s.setBotBC("noFlux") #in acc. with Jorda et al. (2022), however, they assume inflow if h>0
        indxFluxSolute = 2
        # elif s.dimWorld == 1:
            # s.setInnerBC("fluxCyl", 0.)  # [cm/day] #Y 0 pressure?
            # s.setOuterBC("fluxCyl", 0.)
            # indxFluxSolute = 3
            # 
        # else:
            # raise Exception
    
        for i in range(1, s.numComp+1):
            s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(indxFluxSolute))
            s.setParameter( "Soil.BC.Top.C"+str(i)+"Type", str(indxFluxSolute))
            s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(0)) 
            s.setParameter( "Soil.BC.Top.C"+str(i)+"Value", str(0 )) 
        
        
    
    # IC
    if s.dimWorld == 3:
        if isinstance(p_mean_,(int,float)):
            #print('set pmean float',p_mean_)
            s.setHomogeneousIC(p_mean_, equilibrium = True)  # cm pressure head
        elif isinstance(p_mean_,type(np.array([]))):
            pass
        else:
            print(type(p_mean_))
            raise Exception
        
    s.initializeProblem(maxDt = 250/(3600*24))
    
        
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
    
    
    solute_conc = np.array(s.getSolution(1))
    if rank == 0:
        try:
            assert min(solute_conc) >=0
        except:
            print("soil_sol_fluxes", solute_conc)
            print("min(solute_conc)",min(solute_conc))
            raise Exception
        
    cidx = np.array(s.base.getCellIndices())
    cidx_sorted = np.sort(cidx)
    if (cidx != cidx_sorted).any():
        print('too many threads for  the number of cells: ,',cidx,cidx_sorted)
        raise Exception
    if False:
        for ncomp in range(s.numComp):
            print('after creation: s.getSolution(ncomp + 1)',
                  ncomp + 1,np.array(s.getSolution(ncomp + 1)).flatten())
        
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


def theta2H(vg,theta):#(-) to cm
    thetar =vg[0]# 0.059
    thetas = vg[1]#0.445
    alpha = vg[2]#0.00644
    n = vg[3]#1.503
    nrev = 1/(1-1/n)
    H =-(((( (thetas - thetar)/(theta - thetar))**nrev) - 1)**(1/n))/alpha
    return(H)#cm

def getCoefhours(t):
    return (np.sin(np.pi*t*2)+1)/2 #( (t%1) < 0.5)#
    
def qair2rh(qair, press,TK):
    #pressPa = press *100
    #eaPa = qair* pressPa/(0.622+0.378*qair)
    #ea = eaPa/100
    T0 = 273.16
    RH = 26.3 * press * qair  /(exp((17.67*(TK - T0))/(TK- 29.65)))
    return RH
    
    

def resistance2conductance(resistance,r, weatherX):
    resistance = resistance* (1/100) #[s/m] * [m/cm] = [s/cm]
    resistance = resistance * r.R_ph * weatherX["TairK"] / r.Patm # [s/cm] * [K] * [hPa cm3 K−1 mmol−1] * [hPa] = [s] * [cm2 mmol−1]
    resistance = resistance * (1000) * (1/10000)# [s cm2 mmol−1] * [mmol/mol] * [m2/cm2] = [s m2 mol−1]
    return 1/resistance
    
def weather(simDuration, spellData, hp:float=1):
        if simDuration == 0.:
            raise Exception
        Qnigh = 0; Qday = 960e-6 #458*2.1
        
        if  ((spellData['condition'] == "wet") or (simDuration <= spellData['spellStart']) or (simDuration > spellData['spellEnd'])):
            Tnigh = 15.8; Tday = 22
            #Tnigh = 13; Tday = 20.7
            #specificHumidity = 0.0097
            RHday = 0.6; RHnigh = 0.88
            Pair = 1010.00 #hPa
            thetaInit = 0.4#
            cs = 350e-6
        elif spellData['condition'] == "dry":
            Tnigh = 20.7; Tday = 30.27
            #Tnigh = 15.34; Tday = 23.31
            #specificHumidity = 0.0097# 0.0111
            RHday = 0.44; RHnigh = 0.78
            Pair = 1070.00 #hPa
            thetaInit = 28/100   
            #thetaInit = 10.47/100
            cs = 350e-6
        else:
            print('spellData',spellData)
            raise Exception
            
        coefhours = sinusoidal(simDuration)/2
        RH_ = RHnigh + (RHday - RHnigh) * coefhours
        TairC_ = Tnigh + (Tday - Tnigh) * coefhours
        Q_ = Qnigh + (Qday - Qnigh) * coefhours
         #co2 paartial pressure at leaf surface (mol mol-1)
        #390, 1231
        #RH = 0.5 # relative humidity
        es =  6.112 * np.exp((17.67 * TairC_)/(TairC_ + 243.5))
        ea = es*RH_#qair2ea(specificHumidity,  Pair)
        assert ea < es
        #RH = ea/es
        assert ((RH_ > 0) and(RH_ < 1))
        bl_thickness = 1/1000 #1mm * m_per_mm
        diffusivity= 2.5e-5#m2/sfor 25*C
        rbl =bl_thickness/diffusivity #s/m 13
        #cs = 350e-6
        Kcanopymean = 1e-1 # m2/s
        meanCanopyL = (2/3) * hp /2
        rcanopy = meanCanopyL/Kcanopymean
        windSpeed = 2 #m/s
        zmzh = 2 #m
        karman = 0.41 #[-]
        
        rair = 1
        if hp > 0:
            rair = np.log((zmzh - (2/3)*hp)/(0.123*hp)) * np.log((zmzh - (2/3)*hp)/(0.1*hp)) / (karman*karman*windSpeed)
            #print()
            #raise Exception
            

        # pmean = theta2H(vgSoil, thetaInit)

        weatherVar = {'TairC' : TairC_,'TairK' : TairC_ + 273.15,'Pair':Pair,"es":es,
                        'Qlight': Q_,'rbl':rbl,'rcanopy':rcanopy,'rair':rair,"ea":ea,
                        'cs':cs, 'RH':RH_, 'theta':thetaInit}
        #print("Env variables at", round(simDuration//1),"d",round((simDuration%1)*24),"hrs :\n", weatherVar)
        return weatherVar

def phloemParam(r,weatherInit ):
    r = setKrKx_phloem(r)
    r.g0 = 8e-6
    r.VcmaxrefChl1 =1.28#/2
    r.VcmaxrefChl2 = 8.33#/2
    r.a1 = 0.6/0.4#0.7/0.3#0.6/0.4 #ci/(cs - ci) for ci = 0.6*cs
    r.a3 = 1.5
    r.alpha = 0.4#0.2#/2
    r.theta = 0.6#0.9#/2
    r.k_meso = 1e-3#1e-4
    r.setKrm2([[2e-5]], False)
    r.setKrm1([[10e-2]], False)#([[2.5e-2]])
    r.setRhoSucrose([[0.51],[0.65],[0.56]], False)#0.51
    #([[14.4,9.0,0,14.4],[5.,5.],[15.]])
    rootFact = 2
    r.setRmax_st([[2.4*rootFact,1.5*rootFact,0.6*rootFact,2.4*rootFact],[2.,2.],[8.]], False)#6.0#*6 for roots, *1 for stem, *24/14*1.5 for leaves
    #r.setRmax_st([[12,9.0,6.0,12],[5.,5.],[15.]])
    r.KMrm = 0.1#VERY IMPORTANT TO KEEP IT HIGH
    #r.exud_k = np.array([2.4e-4])#*10#*(1e-1)
    #r.k_gr = 1#0
    r.sameVolume_meso_st = False
    r.sameVolume_meso_seg = True
    r.withInitVal =True
    r.initValST = 0.6#0.0
    r.initValMeso = 0.9#0.0
    r.beta_loading = 0.6
    r.Vmaxloading = 0.05 #mmol/d, needed mean loading rate:  0.3788921068507634
    r.Mloading = 0.2
    r.Gr_Y = 0.8
    r.CSTimin = 0.4
    r.surfMeso=0.0025
    r.leafGrowthZone = 2 # cm
    r.StemGrowthPerPhytomer = True # 
    r.psi_osmo_proto = -10000*1.0197 #schopfer2006
    
    r.fwr = 1e-16
    r.fw_cutoff =  0.09497583
    r.sh = 4e-4
    r.gm=0.01
    r.p_lcrit =  -15000*0.6
    
    r.limMaxErr = 1/100
    r.maxLoop = 10000
    r.minLoop=900
    
    r.C_targ = r.CSTimin
    r.C_targMesophyll = r.CSTimin
    r.k_S_ST = 5/25 *2 #daudet2002
    r.k_S_Mesophyll = 5/25*0 #daudet2002
    r.k_mucil = 1
   


    r.cs = weatherInit["cs"]

    #r.r_forPhloem(24/14*1.5, 4)
    #r.r_forPhloem(24/14, 3)
    #r.r_forPhloem(6, 2) #because roots will have high C_ST less often
    r.expression = 6
    r.update_viscosity = True
    r.solver = 1
    r.atol = 1e-12
    r.rtol = 1e-8
    #r.doNewtonRaphson = False;r.doOldEq = False
    SPAD= 41.0
    chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
    r.Chl = np.array( [chl_]) 
    return r

def setKrKx_phloem(r): #inC

    #number of vascular bundles
    VascBundle_leaf = 32
    VascBundle_stem = 52
    VascBundle_root = 1 #valid for all root type
            
    #numPerBundle
    numL = 18
    numS = 21
    numr0 = 33
    numr1 = 25
    numr2 = 25
    numr3 = 1
        
    #radius of phloem type^4 * number per bundle
    rad_s_l   = numL* (0.00025 **4)# * 2; rad_x_l_2   = (0.0005 **4) * 2   
    rad_s_s   = numS *(0.00019 **4) #* 3; rad_x_s_2   = (0.0008 **4) * 1     
    rad_s_r0  = numr0 *(0.00039 **4) #* 4    
    rad_s_r12 = numr1*(0.00035**4) #* 4; rad_x_r12_2 = (0.00087**4) * 1
    rad_s_r3  = numr3 *(0.00068**4) #* 1      

    # axial conductivity [cm^3/day] , mu is added later as it evolves with CST  
    beta = 0.9 #Thompson 2003a
    kz_l   = VascBundle_leaf * rad_s_l   * np.pi /8 * beta  
    kz_s   = VascBundle_stem * rad_s_s   * np.pi /8 * beta
    kz_r0  = VascBundle_root * rad_s_r0  * np.pi /8 * beta
    kz_r12 = VascBundle_root * rad_s_r12 * np.pi /8 * beta
    kz_r3  = VascBundle_root * rad_s_r3  * np.pi /8 * beta
    
    #print([[kz_r0,kz_r12,kz_r12,kz_r3],[kz_s,kz_s ],[kz_l]])
    #raise Exception
    #radial conductivity [1/day],
    kr_l  = 0.#3.83e-4 * hPa2cm# init: 3.83e-4 cm/d/hPa
    kr_s  = 0.#1.e-20  * hPa2cm # set to almost 0
    # kr_r0 = 1e-1
    # kr_r1 = 1e-1
    # kr_r2 = 1e-1
    # kr_r3 = 1e-1
    kr_r0 = 5e-2
    kr_r1 = 5e-2
    kr_r2 = 5e-2
    kr_r3 = 5e-2
    l_kr =  100#0.8 #cm
    
    r.setKr_st([[kr_r0,kr_r1 ,kr_r2 ,kr_r0],[kr_s,kr_s ],[kr_l]] , kr_length_= l_kr, verbose = False)
    r.setKx_st([[kz_r0,kz_r12,kz_r12,kz_r0],[kz_s,kz_s ],[kz_l]], False)
    
    a_ST = [[0.00039,0.00035,0.00035,0.00039 ],[0.00019,0.00019],[0.00025]]
    Across_s_l   = numL*VascBundle_leaf *(a_ST[2][0]**2)*np.pi# (0.00025 **2)# * 2; rad_x_l_2   = (0.0005 **4) * 2   
    Across_s_s   = numS *VascBundle_stem * (a_ST[1][0]**2)*np.pi#(0.00019 **2) #* 3; rad_x_s_2   = (0.0008 **4) * 1     
    Across_s_r0  = numr0 *VascBundle_root * (a_ST[0][0]**2)*np.pi#(0.00039 **2) #* 4    
    Across_s_r12 = numr1*VascBundle_root * (a_ST[0][1]**2)*np.pi#(0.00035**2) #* 4; rad_x_r12_2 = (0.00087**4) * 1
    Across_s_r3  =  numr3 *VascBundle_root *(a_ST[0][2]**2)*np.pi# (0.00068**2) #* 1    
    
    Perimeter_s_l   = numL*VascBundle_leaf *(a_ST[2][0])* 2 * np.pi# (0.00025 **2)# * 2; rad_x_l_2   = (0.0005 **4) * 2   
    Perimeter_s_s   = numS *VascBundle_stem * (a_ST[1][0])* 2 * np.pi#(0.00019 **2) #* 3; rad_x_s_2   = (0.0008 **4) * 1     
    Perimeter_s_r0  = numr0 *VascBundle_root * (a_ST[0][0])* 2 * np.pi#(0.00039 **2) #* 4    
    Perimeter_s_r12 = numr1*VascBundle_root * (a_ST[0][1])* 2 * np.pi#(0.00035**2) #* 4; rad_x_r12_2 = (0.00087**4) * 1
    Perimeter_s_r3  =  numr3 *VascBundle_root *(a_ST[0][2])* 2 * np.pi# (0.00068**2) #* 1  
    #print(a_ST[2][0],a_ST[1][0],a_ST[0][0],a_ST[0][1],a_ST[0][2])
    #r.a_ST = a_ST #to check for water equilibrium assumption
    #tot surface/np.pi of sieve tube  (np.pi added after)
    #r.a_ST_eqs = [[rad_s_r0,rad_s_r12,rad_s_r12,rad_s_r0],[rad_s_s,rad_s_s],[rad_s_l]]
    r.setAcross_st([[Across_s_r0,Across_s_r12,Across_s_r12,Across_s_r0],[Across_s_s,Across_s_s],[Across_s_l]], False)
    return r
    
def create_mapped_plant( nc, logbase, mode,initSim,
                min_b , max_b , cell_number, soil_model, fname, path, 
                stochastic = False, mods = None, plantType = "plant",l_ks_ = "dx_2",
                recreateComsol_ = False, usemoles = True, limErr1d3d = 1e-11, spellData =None):
    """ loads a rmsl file, or creates a rootsystem opening an xml parameter set,  
        and maps it to the soil_model """
    #global picker  # make sure it is not garbage collected away...
    
    if fname.endswith(".rsml"):
        r = XylemFluxPython(fname)
        if plantType == "plant":
            from rhizo_modelsPlant import RhizoMappedSegments  # Helper class for cylindrical rhizosphere models
        else:
            from rhizo_modelsRS import RhizoMappedSegments  # Helper class for cylindrical rhizosphere models
        rs = RhizoMappedSegments( nc, logbase, mode, soil_model, 
                                    recreateComsol_, usemoles, seedNum = seed, 
                                 limErr1d3dAbs = limErr1d3d, l_ks=l_ks_)
    elif fname.endswith(".xml"):
        seed = 1
        weatherInit = weather(initSim,spellData)
        if plantType == "plant":
            from rhizo_modelsPlant import RhizoMappedSegments  # Helper class for cylindrical rhizosphere models
        else:
            from rhizo_modelsRS import RhizoMappedSegments  # Helper class for cylindrical rhizosphere models
        rs = RhizoMappedSegments( nc, logbase, mode, soil_model, 
                                    recreateComsol_, usemoles, seedNum = seed, 
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
    #assert rs.getSeedVal() == seed
    # print('seedval at rank',rank,  rs.getSeedVal())
    # if rank == 0:
    if plantType == "plant":    
        r = init_conductivities(r)
        r = phloemParam(r, weatherInit)
        rs.set_phloem_flux(r)
        #r.test()
        return rs, r
    else:
        r = init_conductivities_const(r)
        return rs, r
    

def div0(a, b, c):   # function to avoid division by 0    
    if isinstance(a,(list, np.ndarray)):  
        return np.divide(a, b, out=np.full(len(a), c), where=b!=0)
    else:
        return div0f(a, b, c)
    
def div0f(a, b, c):    # function to avoid division by 0     
    if b != 0:
        return a/b
    else:
        return a/c

            
def write_files(file_name, psi_x, psi_i, sink, times, trans, psi_s, vol_, surf_,  depth_,  
                    dist, con, l, conc = None, c_ = None, directory ="./results/"):#krs_,
    """  saves numpy arrays ass npy files """
    
    write_file_array('psix_' + file_name, np.array(psi_x), directory_ =directory)  # xylem pressure head per segment [cm]
    write_file_array('psiinterface_' + file_name, np.array(psi_i), directory_ =directory)  # pressure head at interface per segment [cm]
    write_file_array('sink_' + file_name, -np.array(sink), directory_ =directory)  # sink per segment [cm3/day]
    write_file_array('transpiration_' + file_name, np.vstack((times, -np.array(trans))), directory_ =directory)  # time [day], transpiration [cm3/day]
    write_file_array('soil_' + file_name, np.array(psi_s), directory_ =directory)  # soil potential per cell [cm]

    write_file_array('vol_' + file_name, np.array(vol_), directory_ =directory)  # volume per subType [cm3]
    write_file_array('surf_' + file_name, np.array(surf_), directory_ =directory)  # surface per subType [cm2]
    #np.save('krs_' + file_name, np.array(krs_))  # soil potential per cell [cm2/day]
    write_file_array('depth_' + file_name, np.array(depth_), directory_ =directory)  # root system depth [cm]

    # if conc is not None:
        # write_file_array('soilc_' + file_name, np.array(conc))  # soil potential per cell [cm]
    # if c_ is not None:
        # write_file_array('carbon_' + file_name,  np.vstack((times, -np.array(c_))))  # exudation per segment  [g]

    # np.savez('solute_dist' + file_name, np.array(dist), np.array(con),  np.array(l))  # solute concentrations in cylinders



