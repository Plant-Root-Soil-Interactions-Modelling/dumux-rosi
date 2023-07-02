
import sys; 
CPBdir = "../../../../cpb3101/CPlantBox"
sys.path.append(CPBdir+"/src");
sys.path.append(CPBdir);
sys.path.append("../../..");sys.path.append(".."); 
sys.path.append(CPBdir+"/src/python_modules");
sys.path.append("../build-cmake/cpp/python_binding/") # dumux python binding
sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../modules/") # python wrappers 
from helpuqrMasterCopy1 import *
import numpy as np

condition = "dry"
def weather(simDuration, hp):
    vgSoil = [0.059, 0.45, 0.00644, 1.503, 1]
    loam = [0.08, 0.43, 0.04, 1.6, 50]
    Qnigh = 0; Qday = 960e-6 #458*2.1
    if (condition == "wet"):
        Tnigh = 15.8; Tday = 22
        #Tnigh = 13; Tday = 20.7
        #specificHumidity = 0.0097
        RHday = 0.60; RHnigh = 0.88
        Pair = 1010.00 #hPa
        thetaInit = 40/100
        cs = 350e-6
    elif condition == "dry":
        Tnigh = 20.7; Tday = 30.27
        #Tnigh = 15.34; Tday = 23.31
        #specificHumidity = 0.0097# 0.0111
        RHday = 0.7; RHnigh =0.7#0.44; RHnigh = 0.78
        Pair = 1070.00 #hPa
        thetaInit = 28/100   
        cs = 350e-6
    else:
        print("condition",condition)
        raise Exception("condition not recognised")
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
    bl_thickness = 1/1000 #m
    diffusivity= 2.5e-5#m2/sfor 25°C
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

def resistance2conductance(resistance,r,weatherX):
    resistance = resistance* (1/100) #[s/m] * [m/cm] = [s/cm]
    resistance = resistance * r.R_ph * weatherX["TairK"] / r.Patm # [s/cm] * [K] * [hPa cm3 K−1 mmol−1] * [hPa] = [s] * [cm2 mmol−1]
    resistance = resistance * (1000) * (1/10000)# [s cm2 mmol−1] * [mmol/mol] * [m2/cm2] = [s m2 mol−1]
    return 1/resistance
def initPlant(simDuration):
    weatherInit = weather(0,0)
    #simDuration = 25 # [day] init simtime
    #spellDuration = 5
    simMax = 26#simStartSim+ spellDuration
    depth = 40
    dt = 1/24 #10min
    verbose = True

    # plant system 
    pl = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
    #pl2 = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
    path = CPBdir+"/modelparameter/plant/"
    name = "Triticum_aestivum_test_2021"#"Triticum_aestivum_adapted_2021"#

    pl.readParameters(path + name + ".xml")
    #pl2.readParameters(path + name + ".xml")



    #raise Exception
    sdf = pb.SDF_PlantBox(np.Inf, np.Inf, depth )

    pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil
    #pl2.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil


    pl.initialize(verbose = True)#, stochastic = False)
    pl.simulate(simDuration, False)#, "outputpm15.txt")

    #raise Exception
    """ Coupling to soil """



    min_b = [-3./2, -12./2, -41.]#distance between wheat plants
    max_b = [3./2, 12./2, 0.]
    rez = 0.5
    cell_number = [int(6*rez), int(24*rez), int(40*rez)]#1cm3? 
    layers = depth; soilvolume = (depth / layers) * 3 * 12
    k_soil = []
    initial = weatherInit["p_mean"]#mean matric potential [cm] pressure head

    p_mean = initial
    p_bot = p_mean + depth/2
    p_top = initial - depth/2
    sx = np.linspace(p_top, p_bot, depth)
    picker = lambda x,y,z : max(int(np.floor(-z)),-1) 
    sx_static_bu = sx    
    pl.setSoilGrid(picker)  # maps segment



    """ Parameters phloem and photosynthesis """
    r = PhloemFluxPython(pl,psiXylInit = min(sx),ciInit = weatherInit["cs"]*0.5) #XylemFluxPython(pl)#
    #r2 = PhloemFluxPython(#pl2,psiXylInit = min(sx),ciInit = weatherInit["cs"]*0.5) #XylemFluxPython(pl)#

    r = setKrKx_phloem(r)

    r.oldciEq = True

    r.Rd_ref = 0 #to avoid error (C < 0 in meso, mention this in paper)
    r.g0 = 8e-3
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
    r.initValST = 0.#0.6#0.0
    r.initValMeso = 0.#0.9#0.0
    r.beta_loading = 0.6
    r.Vmaxloading = 0.05 #mmol/d, needed mean loading rate:  0.3788921068507634
    r.Mloading = 0.2
    r.Gr_Y = 0.8
    r.CSTimin = 0.4
    r.surfMeso=0.0025
    r.leafGrowthZone = 2 # cm
    r.StemGrowthPerPhytomer = True # 
    r.psi_osmo_proto = -10000*1.0197 #schopfer2006
    r.fwr = 0

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
    
    hp = max([tempnode[2] for tempnode in r.get_nodes()]) /100 

    weatherX = weather(simDuration, hp)
    r.Patm = weatherX["Pair"]
    ##resistances
    r.g_bl = resistance2conductance(weatherX["rbl"],r,weatherX) / r.a2_bl
    r.g_canopy = resistance2conductance(weatherX["rcanopy"],r,weatherX) / r.a2_canopy
    r.g_air = resistance2conductance(weatherX["rair"],r,weatherX) / r.a2_air

    r.Qlight = weatherX["Qlight"] #; TairC = weatherX["TairC"] ; text = "night"


    r = setKrKx_xylem(weatherX["TairC"], weatherX["RH"],r)
    r.es = weatherX["es"]
    return r, weatherX
simDuration = 10
r = initPlant(simDuration)[0]
weatherX = initPlant(simDuration)[1] 
testsx =np.linspace(-100, -20000, 100)
p_errors =[]
r.usePg4Fw = True


for i, p_mean in enumerate(np.array([1])):#testsx):
    r.sh = 4e-4
    r.fwr = 0#0.001
    r.gm = 0.05
    r.loop=0
    p_mean = -10000;
    depth=60
    p_bot = p_mean + depth/2
    p_top = p_mean - depth/2
    sx = np.linspace(p_top, p_bot, depth)
    pg = np.array(r.pg)

    try:
        r.solve_photosynthesis(sim_time_ = simDuration, sxx_=sx, cells_ = True,ea_ = weatherX["ea"],
                        verbose_ = False, doLog_ = False,TairC_= weatherX["TairC"] )
    except:
        p_errors.append(int(np.mean(pg[np.where(pg<0)])))
        print(i, p_mean,np.mean(pg[np.where(pg<0)]))
        r.followTrace = True
        r.loop=0
        r.solve_photosynthesis(sim_time_ = simDuration, sxx_=sx, cells_ = True,ea_ = weatherX["ea"],
                        verbose_ = False, doLog_ = False,TairC_= weatherX["TairC"] )