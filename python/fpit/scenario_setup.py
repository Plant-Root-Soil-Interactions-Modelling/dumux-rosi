""" 
Functions to simplify setup of the scenarios for the INARI project
"""

import sys;
sys.path.append("../../build-cmake/cpp/python_binding/");
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
import timeit
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

# from rosi_richards2c import Richards2CSP  # C++ part (Dumux binding), macroscopic soil model
from rosi_richards10c_cyl import Richards10CCylFoam # C++ part (Dumux binding)
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial model (a single cylindrical model is not allowed to run in parallel)
from rosi_richards10c import Richards10CSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
from functional.phloem_flux import PhloemFluxPython  # root system Python hybrid solver


import plantbox as pb  # CPlantBox
import van_genuchten as vg
import evapotranspiration as evap
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


def write_file_float(name, data, directory_, allranks = False):
    if (rank == 0) or allranks:
        name2 = directory_+ name+ '.txt'
        with open(name2, 'a') as log:
            log.write(repr( data)  +'\n')
        
def write_file_array(name, data, space =",", directory_ ="./results/", fileType = '.txt', allranks = False ):
    np.array(data).reshape(-1)
    try:
        if (rank == 0) or allranks:
            name2 = directory_+ name+ fileType
            #print('write_file_array',name)
            with open(name2, 'a') as log:
                log.write(space.join([num for num in map(str, data)])  +'\n')
    except:
        print(name, data,data.shape)
        raise Exception

def vg_SPP(i = int(1)):
    """ Van Genuchten parameter, called by maize()  """
        
    soil = {}
    soil[0] = [0.045, 0.43, 0.04, 1.6, 50]# [0.041, 0.494, 0.0256, 1.49, 245]
    soil[1] = [0.03, 0.414, 0.038, 2, 1864]
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


def create_soil_model(soil_type, year, soil_, min_b , max_b , cell_number, demoType, times = None, 
                        net_inf = None, usemoles = True, dirResults = "",#"./results/parallel"+str(max_rank)+"/",
                        p_mean_ = -100,css1Function = 3):
    """
        Creates a soil domain from @param min_b to @param max_b with resolution @param cell_number
        soil demoType is fixed and homogeneous 
        domain is periodic (if 2d or 3d)
        initial potentials are linear from @param p_top to @param p_bot
        
        returns soil_model (RichardsWrapper(RichardsSP())) and soil parameter (vg.Parameters)
    """
    if len(dirResults) == 0:
        diResults = "./results/parallel"+str(max_rank)+"/"
    do1D = False
    if do1D:
        s = RichardsNoMPIWrapper(Richards10CCylFoam(), usemoles)  # water and N solute          
    else:
        s = RichardsWrapper(Richards10CSP(), usemoles)  # water and N solute          
            
    s.soil = soil_
    s.vg_soil = vg.Parameters(soil_) 
    #@see dumux-rosi\cpp\python_binding\solverbase.hh
    s.betaC = 0.001 
    s.betaO = 0.1 
    C_S = 0. # in mol/m3 water
    
    ## ATT! C_S_W_thres has to be > 0
    s.C_S_W_thresC = 0.1/1e6 # in mol/cm3 water
    s.C_S_W_thresO = 0.05/1e6
    
    
    s.k_decay = 0.2
    s.k_decay2 = 0.6 
    s.k_DC = 1. 
    s.k_DO = 1. 
    s.k_growthC = 0.2
    s.k_growthO = 0.2 
    s.K_L = 8.3
    s.k_phi = 0.1 
    s.k_RC = 0.1 
    s.k_RO = 0.1 
    
    s.k_sorp = 0.4
    s.f_sorp = 0.5
    if demoType == "dumux_10c":
        C_S = 1e-8 # in mol/m3 water
        s.CSSmax =1/1e6 # 1e-4*10000*0.
        s.alpha =0.1# 0.
        unitConversion = 1e3
        doBio = 1.
        
        CSS2_init = 0#s.CSSmax*1e6 * (C_S/(C_S+ s.k_sorp*1e6)) * (1 - s.f_sorp)#mol C/ m3 scv
    elif demoType == "dumux_3c":
        C_S = 1e-8 # in mol/m3 water
        s.CSSmax =0. # 1e-4*10000*0.
        s.alpha =0.1# 0.
        unitConversion = 1e3
        doBio = 0.
        CSS2_init = s.CSSmax*1e6 * (C_S/(C_S+ s.k_sorp*1e6)) * (1 - s.f_sorp)#mol C/ m3 scv
    elif demoType == "dumux_w":
        s.CSSmax = 0.
        s.alpha = 0.
        doBio = 0.
        unitConversion = 0.
        CSS2_init = 0.
    else:
        print('demoType',demoType)
        raise Exception
    
    
    s.css1Function = css1Function
    s.C_aOLim=1.e-10
    s.C_aCLim=1.e-10
    s.setParameter("Soil.C_aOLim", str(s.C_aOLim)) #[molC/cm3 scv]
    s.setParameter("Soil.C_aCLim", str(s.C_aCLim)) #[molC/cm3 scv]
    
    s.setParameter( "Soil.css1Function", str(s.css1Function))
    s.ICcc = np.array([C_S *unitConversion,
                       C_S/2 *unitConversion,
                        (C_S/10+s.C_aOLim)* unitConversion *doBio,
                        C_S/10* unitConversion *doBio,
                        (C_S/10+s.C_aCLim)* unitConversion *doBio,
                        C_S/10* unitConversion *doBio,
                        CSS2_init, 0.])# in mol/m3 water or mol/m3 scv
        
        
    s.initialize()
    
    if do1D:
        a_in = 0.05
        a_out = 0.5664008176051574
        cell_number = 9
        lb = 0.5
        points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), cell_number + 1, base = lb)
        s.createGrid1d(points)
    else:
        s.createGrid(min_b, max_b, cell_number, False)  # [cm] 

    #cell_number = str(cell_number)
    cell_number_ = cell_number
    cell_number= s.dumux_str(cell_number)#.replace("[", "");cell_number=cell_number.replace("]", "");cell_number=cell_number.replace(",", "");
    s.setParameter( "Soil.Grid.Cells", cell_number)    
    s.setParameter("Problem.reactionExclusive", s.dumux_str(int(not do1D)))
    
    # BC
    s.setTopBC("noFlux")
    s.setBotBC("noFlux") #in acc. with Jorda et al. (2022), however, they assume inflow if h>0
    s.solidDensity = 2700 # [kg/m^3 solid]
    s.solidMolarMass = 60.08e-3 # [kg/mol] 
    # [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
    s.solidMolDensity = s.solidDensity/s.solidMolarMass
    # [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
    s.bulkDensity_m3 = s.solidMolDensity*(1.- s.vg_soil.theta_S)

    s.setParameter( "Soil.MolarMass", str(s.solidMolarMass))
    s.setParameter( "Soil.solidDensity", str(s.solidDensity))

    s.Ds = 1e-9 # m^2/s
    s.Dl = 3e-12
    s.numComp = 8
    s.numFluidComp = 2
    s.setParameter("1.Component.LiquidDiffusionCoefficient", str(s.Ds)) #m^2/s
    s.setParameter("2.Component.LiquidDiffusionCoefficient", str(s.Dl)) #m^2/s

    s.decay = 0. #1.e-5
    
    for i in range(1, s.numComp+1):
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(2))
        s.setParameter( "Soil.BC.Top.C"+str(i)+"Type", str(2))
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(0)) 
        s.setParameter( "Soil.BC.Top.C"+str(i)+"Value", str(0 )) 
    for i in range(s.numComp):
        molarC = s.ICcc[i] / s.phaseDensity(isDissolved = (i < s.numFluidComp)) #mol/m3 to mol/mol
        s.setParameter( "Soil.IC.C"+str(i+1), str(molarC ))

    s.setParameter("Soil.betaC", str(s.betaC))
    s.setParameter("Soil.betaO", str(s.betaO ))
    s.setParameter("Soil.C_S_W_thresC", str(s.C_S_W_thresC )) #mol/cm3
    s.setParameter("Soil.C_S_W_thresO", str(s.C_S_W_thresO )) #mol/cm3
    s.setParameter("Soil.k_decay", str(s.k_decay ))
    s.setParameter("Soil.k_decay2", str(s.k_decay2))
    #s.setParameter("Soil.k_decay3", str(1 ))
    s.setParameter("Soil.k_DC", str(s.k_DC )) # 1/d
    s.setParameter("Soil.k_DO", str(s.k_DO )) # 1/d
    #s.setParameter("Soil.k_growthBis", str(1 )) #bool
    s.setParameter("Soil.k_growthC", str(s.k_growthC ))
    s.setParameter("Soil.k_growthO", str(s.k_growthO ))
    s.setParameter("Soil.K_L", str(s.K_L))#[mol/cm3]
    s.setParameter("Soil.k_phi", str(s.k_phi ))
    s.setParameter("Soil.k_RC", str(s.k_RC))
    s.setParameter("Soil.k_RO", str(s.k_RO))

    s.k_SC = 1
    s.k_SO = 10
    s.setParameter("Soil.k_SC", str(s.k_SC )) #cm^3/mol/d
    #s.setParameter("Soil.k_SCBis", str(k_SC )) #cm^3/mol/d
    s.setParameter("Soil.k_SO", str(s.k_SO )) #cm^3/mol/d
    #s.setParameter("Soil.k_SOBis", str(k_SO )) #cm^3/mol/d

    s.m_maxC = 0.1 
    s.m_maxO = 0.02 
    s.setParameter("Soil.m_maxC", str(s.m_maxC  ))# 1/d
    s.setParameter("Soil.m_maxO", str(s.m_maxO  ))# 1/d
    s.micro_maxC = 2
    s.micro_maxO = 0.01
    s.setParameter("Soil.micro_maxC", str(s.micro_maxC ))# 1/d
    s.setParameter("Soil.micro_maxO", str(s.micro_maxO ))# 1/d
    s.v_maxL = 1.5
    s.setParameter("Soil.v_maxL", str(s.v_maxL))#[d-1]

    s.setParameter("Soil.k_sorp", str(s.k_sorp)) # mol / cm3
    s.setParameter("Soil.f_sorp", str(s.f_sorp)) #[-]
    s.setParameter("Soil.CSSmax", str(s.CSSmax)) #[mol/cm3 scv]
    s.setParameter("Soil.alpha", str(s.alpha)) #[1/d]


    # Paramters
    #dumux-rosi\python\modules\richards.py
    s.setVGParameters([soil_])
    #@see dumux-rosi\cpp\python_binding\solverbase.hh
    #s.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
    s.MaxRelativeShift = 1e-8
    s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift))
    s.setParameter("Problem.verbose", "0")
    
    # IC
    if do1D:
        s.setParameter("Problem.EnableGravity", "false")
    if isinstance(p_mean_,(int,float)):
        #print('set pmean float',p_mean_)
        s.setHomogeneousIC(p_mean_, equilibrium = not do1D)  # cm pressure head
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
    
    write_file_array('getWaterContent',s.getWaterContent(), directory_ =dirResults, fileType = '.csv')
    write_file_array('getSolutionHead',s.getSolutionHead(), directory_ =dirResults, fileType = '.csv')
    
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
    
def weather(simDuration, hp:float=1):
        if simDuration == 0.:
            raise Exception
        vgSoil = [0.059, 0.45, 0.00644, 1.503, 1]
        loam = [0.08, 0.43, 0.04, 1.6, 50]
        Qnigh = 0; Qday = 960e-6 #458*2.1
        condition = "wet"
        if condition == "wet":
            Tnigh = 15.8; Tday = 22
            #Tnigh = 13; Tday = 20.7
            #specificHumidity = 0.0097
            RHday = 0.6; RHnigh = 0.88
            Pair = 1010.00 #hPa
            thetaInit = 0.4#
            cs = 350e-6
        elif condition == "dry":
            Tnigh = 20.7; Tday = 30.27
            #Tnigh = 15.34; Tday = 23.31
            #specificHumidity = 0.0097# 0.0111
            RHday = 0.44; RHnigh = 0.78
            Pair = 1070.00 #hPa
            thetaInit = 28/100   
            #thetaInit = 10.47/100
            cs = 350e-6
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
            

        pmean = theta2H(vgSoil, thetaInit)

        weatherVar = {'TairC' : TairC_,'TairK' : TairC_ + 273.15,'Pair':Pair,"es":es,
                        'Qlight': Q_,'rbl':rbl,'rcanopy':rcanopy,'rair':rair,"ea":ea,
                        'cs':cs, 'RH':RH_, 'p_mean':pmean, 'vg':loam}
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
    r.setKrm2([[2e-5]])
    r.setKrm1([[10e-2]])#([[2.5e-2]])
    r.setRhoSucrose([[0.51],[0.65],[0.56]])#0.51
    #([[14.4,9.0,0,14.4],[5.,5.],[15.]])
    rootFact = 2
    r.setRmax_st([[2.4*rootFact,1.5*rootFact,0.6*rootFact,2.4*rootFact],[2.,2.],[8.]])#6.0#*6 for roots, *1 for stem, *24/14*1.5 for leaves
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
    l_kr =  0.8 #cm
    
    r.setKr_st([[kr_r0,kr_r1 ,kr_r2 ,kr_r0],[kr_s,kr_s ],[kr_l]] , kr_length_= l_kr)
    r.setKx_st([[kz_r0,kz_r12,kz_r12,kz_r0],[kz_s,kz_s ],[kz_l]])
    
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
    r.setAcross_st([[Across_s_r0,Across_s_r12,Across_s_r12,Across_s_r0],[Across_s_s,Across_s_s],[Across_s_l]])
    return r
    
def create_mapped_plant( nc, logbase, mode,initSim,
                min_b , max_b , cell_number, soil_model, fname, path, 
                stochastic = False, mods = None, plantType = "plant",l_ks_ = "dx_2",
                recreateComsol_ = False, usemoles = True, limErr1d3d = 1e-11):
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
        weatherInit = weather(initSim)
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
        rs.initialize()#stochastic = False)
        rs.simulate(initSim, False)
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
    print('seedval at rank',rank,  rs.getSeedVal())
    assert rs.getSeedVal() == seed
    # if rank == 0:
    if plantType == "plant":    
        r = init_conductivities(r)
        r = phloemParam(r, weatherInit)
        rs.set_phloem_flux(r)
        r.test()
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



