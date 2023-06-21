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

from rosi_richards2c import Richards2CSP  # C++ part (Dumux binding), macroscopic soil model
from rosi_richardsnc import RichardsNCSP  # C++ part (Dumux binding), macroscopic soil model
from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
from functional.phloem_flux import PhloemFluxPython  # root system Python hybrid solver
from rhizo_modelsPlant import *  # Helper class for cylindrical rhizosphere models

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


def vg_SPP(i = int(1)):
    """ Van Genuchten parameter, called by maize()  """
        
    soil = {}
    soil[0] = [0.041, 0.494, 0.0256, 1.49, 245]
    soil[1] = [0.03, 0.414, 0.038, 2, 1864]
    return soil[i]


def maize_SPP(soil_= "loam"):
    """ parameters for maize simulation """
    if soil_ == "loam":
        i = 0
    else:
        i = 1
    soil = vg_SPP(i)
    min_b = [-5., -5, -20.] 
    max_b = [5., 5, 0.]
    cell_number = [5, 5, 20]
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


def create_soil_model(soil_type, year, soil_, min_b , max_b , cell_number, type, times = None, 
                        net_inf = None):
    """
        Creates a soil domain from @param min_b to @param max_b with resolution @param cell_number
        soil type is fixed and homogeneous 
        domain is periodic (if 2d or 3d)
        initial potentials are linear from @param p_top to @param p_bot
        
        returns soil_model (RichardsWrapper(RichardsSP())) and soil parameter (vg.Parameters)
    """
    soil = vg.Parameters(soil_)
    vg.create_mfp_lookup(soil, -1.e5, 1000)

    if type == "dumux":
        s = RichardsWrapper(RichardsSP())  # water only
    elif type == "dumux_dirichlet_2c":
        s = RichardsWrapper(Richards2CSP())  # water and one solute
    elif type == "dumux_dirichlet_nc":
        s = RichardsWrapper(RichardsNCSP())  # water and N solute        
    else:
        raise Exception("choose type: dumux, dumux_dirichlet_2c, dumux_dirichlet_nc")

    #@see dumux-rosi\cpp\python_binding\solverbase.hh
    s.initialize(verbose = False)
    s.createGrid(min_b, max_b, cell_number, False)  # [cm] #######################################################################

    # BC
    if times is not None:
        s.setTopBC("atmospheric", 0.5, [times, net_inf])  # 0.5 is dummy value
    else:
        s.setTopBC("noFlux")
    s.setBotBC("noFlux") #in acc. with Jorda et al. (2022), however, they assume inflow if h>0

    if (type == "dumux_dirichlet_2c") or (type == "dumux_dirichlet_nc"):  # solute BC
        s.setTopBC_solute("outflow", 0.)
        s.setBotBC_solute("outflow", 0.)

    # Paramters
    #dumux-rosi\python\modules\richards.py
    s.setVGParameters([soil_])
    #@see dumux-rosi\cpp\python_binding\solverbase.hh
    s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
    if (type == "dumux_dirichlet_2c") or (type == "dumux_dirichlet_nc"):
        s.setParameter("Component.MolarMass", "1.2e-2")  # carbon 12 g/mol
        s.setParameter("Component.LiquidDiffusionCoefficient", "0")  # 5e-10 m2 s-1 # Darrah et al. 1991
        s.setParameter("Component.BufferPower", "0")  # 5 buffer power = \rho * Kd [1]
        s.decay = 0#1.e-5
        #s.setParameter("Component.Decay", "1.e-5")  # decay [d^-1] (Awad et al. 2017) 
    s.initializeProblem()
    wilting_point = -15000
    s.setCriticalPressure(wilting_point)  # for boundary conditions constantFlow, constantFlowCyl, and atmospheric
    s.ddt = 1.e-5  # [day] initial Dumux time step

    # IC
    IgotTheFile = False
    if IgotTheFile:
        df = pd.read_csv("data_magda/init_pot_"+str(year)+".csv")  # initial potential
        h  = np.flip(df[soil_type].loc[:].values) #cm
        h = np.repeat(h[:,np.newaxis],cell_number[0],axis=1) #x-axis
        h = np.repeat(h[:,:,np.newaxis],cell_number[1],axis=2) #y-axis
        h = h.flatten()
    else:
        h = np.ones((20*45*75))*-100 #TODO
    s.setInitialConditionHead(h)  # cm

    if (type == "dumux_dirichlet_2c") or (type == "dumux_dirichlet_nc"):
        c = np.zeros((cell_number[0]*cell_number[1]*cell_number[2])) #TODO
        s.setInitialCondition(c, 1)  # kg/m3

    # plt.plot(h, np.linspace(-200., 0., h.shape[0]))
    # plt.xlabel("soil matric potential [cm]")
    # plt.ylabel("depth (cm)")
    # plt.tight_layout()
    # plt.show()
    # plt.plot(c, np.linspace(-200, 0., c.shape[0]))
    # plt.xlabel("nitrate concentration [g/cm3]")
    # plt.ylabel("depth (cm)")
    # plt.tight_layout()
    # plt.show()

    return s, soil

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
        coefhours = sinusoidal(simDuration)
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
        print("Env variables at", round(simDuration//1),"d",round((simDuration%1)*24),"hrs :\n", weatherVar)
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
    r.fwr = 0#0.1
    r.sh = 4e-4
    r.gm=0.01
    
    r.limMaxErr = 1/100
    r.maxLoop = 10000
    r.minLoop=900
    
    r.C_targ = r.CSTimin
    r.C_targMesophyll = r.CSTimin
    r.k_S_ST = 5/25 #daudet2002
    r.k_S_Mesophyll = 5/25 #daudet2002
   


    r.cs = weatherInit["cs"]

    #r.r_forPhloem(24/14*1.5, 4)
    #r.r_forPhloem(24/14, 3)
    #r.r_forPhloem(6, 2) #because roots will have high C_ST less often
    r.expression = 6
    r.update_viscosity = True
    r.solver = 1
    r.atol = 1e-10
    r.rtol = 1e-6
    #r.doNewtonRaphson = False;r.doOldEq = False
    SPAD= 41.0
    chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
    r.Chl = np.array( [chl_]) 
    r.Csoil = 1e-4
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
    l_kr = 0.8 #cm
    
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
    
def create_mapped_plant(wilting_point, nc, logbase, mode,initSim,
                min_b , max_b , cell_number, soil_model, fname, path, stochastic = False, mods = None):
    """ loads a rmsl file, or creates a rootsystem opening an xml parameter set,  
        and maps it to the soil_model """
    global picker  # make sure it is not garbage collected away...

    if fname.endswith(".rsml"):
        r = XylemFluxPython(fname)
    elif fname.endswith(".xml"):
        seed = 1
        weatherInit = weather(initSim)

        rs = RhizoMappedSegments(wilting_point, nc, logbase, mode)
        rs.setSeed(seed)
        rs.readParameters(path + fname)
        #if not stochastic:
        #    set_all_sd(rs, 0.)
            
        rs.setGeometry(pb.SDF_PlantBox( max_b[0]-min_b[0],  max_b[1]-min_b[1], max_b[2]-min_b[2]))
        rs.initialize()
        rs.simulate(initSim, False)
        r = PhloemFluxPython(rs,psiXylInit = -659.8 - min_b[2],ciInit = weatherInit["cs"]*0.5) 

    r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), 
                            pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), 
                            cut = False)

    picker = lambda x, y, z: soil_model.pick([x,y, z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
    r.rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
    # comm.barrier()
    # print("survived setSoilGrid", rank)

    # if rank == 0:
    r = init_conductivities(r)
    r = phloemParam(r, weatherInit)
    rs.set_phloem_flux(r)
    r.test()
    return rs


def write_files(file_name, psi_x, psi_i, sink, times, trans, psi_s, vol_, surf_,  depth_,  dist, con, l, conc = None, c_ = None):#krs_,
    """  saves numpy arrays ass npy files """

    np.savetxt('results/psix_' + file_name, np.array(psi_x))  # xylem pressure head per segment [cm]
    np.savetxt('results/psiinterface_' + file_name, np.array(psi_i))  # pressure head at interface per segment [cm]
    np.savetxt('results/sink_' + file_name, -np.array(sink))  # sink per segment [cm3/day]
    np.savetxt('results/transpiration_' + file_name, np.vstack((times, -np.array(trans))))  # time [day], transpiration [cm3/day]
    np.savetxt('results/soil_' + file_name, np.array(psi_s))  # soil potential per cell [cm]

    np.savetxt('results/vol_' + file_name, np.array(vol_))  # volume per subType [cm3]
    np.savetxt('results/surf_' + file_name, np.array(surf_))  # surface per subType [cm2]
    #np.save('results/krs_' + file_name, np.array(krs_))  # soil potential per cell [cm2/day]
    np.savetxt('results/depth_' + file_name, np.array(depth_))  # root system depth [cm]

    if conc is not None:
        np.savetxt('results/soilc_' + file_name, np.array(conc))  # soil potential per cell [cm]
    if c_ is not None:
        np.savetxt('results/carbon_' + file_name,  np.vstack((times, -np.array(c_))))  # exudation per segment  [g]

    np.savez('results/solute_dist' + file_name, np.array(dist), np.array(con),  np.array(l))  # solute concentrations in cylinders



def simulate_const(s, r, trans, sim_time, dt):
    """ 
        classic model:
        potential at root soil interface equals mean matric potential of surrounding finite volume
    """
    wilting_point = -15000  # cm
    skip = 6  # for output and results, skip iteration
    rs_age = 0.  # day

    start_time = timeit.default_timer()
    psi_x_, psi_s_, sink_ , x_, y_, psi_s2_ = [], [], [], [], [], []  # for post processing
    sx = s.getSolutionHead()  # inital condition, solverbase.py
    ns = len(r.rs.segments)
    if rank == 0:
        map_ = r.rs.seg2cell  # because seg2cell is a map
        mapping = np.array([map_[j] for j in range(0, ns)], dtype = np.int64)  # convert to a list

    N = int(np.ceil(sim_time / dt))

    """ simulation loop """
    for i in range(0, N):

        t = i * dt  # current simulation time

        """ 1. xylem model """
        if rank == 0:  # Root part is not parallel
            rx = r.solve(rs_age, -trans * sinusoidal2(t, dt), 0., sx, cells = True, wilting_point = wilting_point)  # xylem_flux.py
            fluxes = r.soilFluxes(rs_age, rx, sx, False)  # class XylemFlux is defined in MappedOrganism.h, approx = False
        else:
            fluxes = None
        fluxes = comm.bcast(fluxes, root = 0)  # Soil part runs parallel

        """ 2. soil model """
        s.setSource(fluxes.copy())  # richards.py
        s.solve(dt)
        sx = s.getSolutionHead()  # richards.py

        """ validity check """

        """ remember results ... """
        if rank == 0 and i % skip == 0:

            sx_ = sx[:, 0]
            psi_x_.append(rx.copy())  # cm (per root node)
            psi_s_.append(np.array([sx_[ci] for ci in mapping]))  # cm (per root segment)
            sink = np.zeros(sx_.shape)
            for k, v in fluxes.items():
                sink[k] += v
            sink_.append(sink)  # cm3/day (per soil cell)
            x_.append(t)  # day
            y_.append(np.sum(sink))  # cm3/day
            psi_s2_.append(sx_)  # cm (per soil cell)

            min_sx, min_rx, max_sx, max_rx = np.min(sx), np.min(rx), np.max(sx), np.max(rx)
            n = round(float(i) / float(N) * 100.)
            print("\n[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g}, {:g}"
                    .format(min_sx, max_sx, min_rx, max_rx, np.sum(sink), -trans * sinusoidal2(t, dt)))

    if rank == 0:
        print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_


if __name__ == '__main__':

    pass
