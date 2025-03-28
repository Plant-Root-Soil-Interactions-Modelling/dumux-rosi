

import sys; 
CPBdir = "../../../CPlantBox"
sys.path.append(CPBdir+"/src");
sys.path.append(CPBdir);
sys.path.append("../../..");sys.path.append(".."); 
sys.path.append(CPBdir+"/src/python_modules");
sys.path.append("../build-cmake/cpp/python_binding/") # dumux python binding
sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../fixedPointIter2/modules/") # python wrappers 

#import matplotlib
#matplotlib.use('AGG') #otherwise "import matplotlib.pyplot" hangs

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
from functional.phloem_flux import PhloemFluxPython  # Python hybrid solver
#from Leuning_speedup import Leuning #about 0.7 for both
#from photosynthesis_cython import Leuning
import plantbox as pb
#from plantbox import Photosynthesis as ph
#import vtk_plot as vp
import math
import os
import numpy as np
import visualisation.vtk_plot as vp
#import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import functional.van_genuchten as vg
import pandas as pd
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()


from rosi_richards10c import RichardsNCSPILU as RichardsNCSP 

def theta2H(vg,theta):#(-) to cm
    thetar =vg[0]# 0.059
    thetas = vg[1]#0.445
    alpha = vg[2]#0.00644
    n = vg[3]#1.503
    nrev = 1/(1-1/n)
    H =-(((( (thetas - thetar)/(theta - thetar))**nrev) - 1)**(1/n))/alpha
    return(H)#cm

def sinusoidal(t):
    return (np.sin(np.pi*t*2)+1)/2 

def qair2rh(qair, press,TK):
    T0 = 273.16
    RH = 26.3 * press * qair  /(exp((17.67*(TK - T0))/(TK- 29.65)))
    return RH



def div0(a, b, c):        
    return np.divide(a, b, out=np.full(len(a), c), where=b!=0)
    
def div0f(a, b, c):    
    if b != 0:
        return a/b
    else:
        return a/c
        



def setKrKx_xylem(TairC, RH,r): #inC
    #mg/cm3
    hPa2cm = 1.0197
    dEauPure = (999.83952 + TairC * (16.952577 + TairC * 
        (- 0.0079905127 + TairC * (- 0.000046241757 + TairC * 
        (0.00000010584601 + TairC * (- 0.00000000028103006)))))) /  (1 + 0.016887236 * TairC)
    siPhi = (30 - TairC) / (91 + TairC)
    siEnne=0
    mu =  pow(10, (- 0.114 + (siPhi * (1.1 + 43.1 * pow(siEnne, 1.25) )))) 
    mu = mu /(24*60*60)/100/1000; #//mPa s to hPa d, 1.11837e-10 hPa d for pure water at 293.15K
    mu = mu * hPa2cm #hPa d to cmh2o d 

    #number of vascular bundles
    VascBundle_leaf = 32
    VascBundle_stem = 52
    VascBundle_root = 1 #valid for all root type
            
    #radius of xylem type^4 * number per bundle
    rad_x_l_1   = (0.0015 **4) * 2; rad_x_l_2   = (0.0005 **4) * 2   
    rad_x_s_1   = (0.0017 **4) * 3; rad_x_s_2   = (0.0008 **4) * 1     
    rad_x_r0_1  = (0.0015 **4) * 4    
    rad_x_r12_1 = (0.00041**4) * 4; rad_x_r12_2 = (0.00087**4) * 1
    rad_x_r3_1  = (0.00068**4) * 1      

    # axial conductivity [cm^3/day]  
    betaXylX = 1#0.1      
    kz_l  = VascBundle_leaf *(rad_x_l_1 + rad_x_l_2)    *np.pi /(mu * 8)  * betaXylX
    kz_s  = VascBundle_stem *(rad_x_s_1 + rad_x_s_2)    *np.pi /(mu * 8) * betaXylX 
    kz_r0 = VascBundle_root * rad_x_r0_1                *np.pi /(mu * 8) * betaXylX  
    kz_r1 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*np.pi /(mu * 8)  * betaXylX
    kz_r2 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*np.pi /(mu * 8)  * betaXylX 
    kz_r3 = VascBundle_root * rad_x_r3_1                *np.pi /(mu * 8) * betaXylX

    #radial conductivity [1/day],0.00014 #
    betaXyl = 1#0.1#0.1
    kr_l  = 3.83e-5 * hPa2cm * betaXyl# init: 3.83e-4 cm/d/hPa
    kr_s  = 0.#1.e-20  * hPa2cm # set to almost 0
    kr_r0 =6.37e-5 * hPa2cm * betaXyl
    kr_r1 =7.9e-5  * hPa2cm * betaXyl
    kr_r2 =7.9e-5  * hPa2cm * betaXyl
    kr_r3 =6.8e-5  * hPa2cm * betaXyl
    l_kr = 0.8 #cm
    r.setKr([[kr_r0,kr_r1,kr_r2,kr_r0],[kr_s,kr_s ],[kr_l]], kr_length_=l_kr) 
    r.setKx([[kz_r0,kz_r1,kz_r2,kz_r0],[kz_s,kz_s ],[kz_l]])
    
    
    Rgaz=8.314 #J K-1 mol-1 = cm^3*MPa/K/mol
    rho_h2o = dEauPure/1000#g/cm3
    Mh2o = 18.05 #g/mol
    MPa2hPa = 10000
    hPa2cm = 1/0.9806806
    #log(-) * (cm^3*MPa/K/mol) * (K) *(g/cm3)/ (g/mol) * (hPa/MPa) * (cm/hPa) =  cm                      
    #p_a = np.log(RH) * Rgaz * rho_h2o * (TairC + 273.15)/Mh2o * MPa2hPa * hPa2cm
    #done withint solve photosynthesis
    #r.psi_air = p_a #*MPa2hPa #used only with xylem
    return r

    

def resistance2conductance(resistance,r, weatherX):
    resistance = resistance* (1/100) #[s/m] * [m/cm] = [s/cm]
    resistance = resistance * r.R_ph * weatherX["TairK"] / r.Patm # [s/cm] * [K] * [hPa cm3 K−1 mmol−1] * [hPa] = [s] * [cm2 mmol−1]
    resistance = resistance * (1000) * (1/10000)# [s cm2 mmol−1] * [mmol/mol] * [m2/cm2] = [s m2 mol−1]
    return 1/resistance

def write_file_array(name, data, directoryN):
    name2 = directoryN+ name+ '.txt'
    with open(name2, 'a') as log:
        log.write(','.join([num for num in map(str, data)])  +'\n')

def write_file_float(name, data, directoryN):
    name2 = directoryN+  name+ '.txt'
    with open(name2, 'a') as log:
        log.write(repr( data)  +'\n')

def weather(simDuration, hp, condition, simStartSim):
    Qnigh = 0; Qday = 960e-6 #458*2.1
    if ((condition == "wet") or (simDuration <= simStartSim) or (simDuration > simStartSim +7)):
        Tnigh = 15.8; Tday = 22
        RHday = 0.6; RHnigh = 0.88
        Pair = 1010.00 #hPa
        thetaInit = 0.4#
        cs = 350e-6
    elif condition == "dry":
        Tnigh = 20.7; Tday = 30.27
        RHday = 0.44; RHnigh = 0.78
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
    es =  6.112 * np.exp((17.67 * TairC_)/(TairC_ + 243.5))
    ea = es*RH_
    assert ea < es
    
    assert ((RH_ > 0) and(RH_ < 1))
    bl_thickness = 1/1000 #1mm * m_per_mm
    diffusivity= 2.5e-5#m2/sfor 25*C
    rbl =bl_thickness/diffusivity #s/m 13
    
    Kcanopymean = 1e-1 # m2/s
    meanCanopyL = (2/3) * hp /2
    rcanopy = meanCanopyL/Kcanopymean
    windSpeed = 2 #m/s
    zmzh = 2 #m
    karman = 0.41 #[-]

    rair = 1
    if hp > 0:
        rair = np.log((zmzh - (2/3)*hp)/(0.123*hp)) * np.log((zmzh - (2/3)*hp)/(0.1*hp)) / (karman*karman*windSpeed)


    weatherVar = {'TairC' : TairC_,'TairK' : TairC_ + 273.15,'Pair':Pair,"es":es,
                    'Qlight': Q_,'rbl':rbl,'rcanopy':rcanopy,'rair':rair,"ea":ea,
                    'cs':cs, 'RH':RH_}
    print("Env variables at", round(simDuration//1),"d",round((simDuration%1)*24),"hrs :\n", weatherVar)
    return weatherVar

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
    rad_s_l   = numL* (0.00025 **4)
    rad_s_s   = numS *(0.00019 **4) 
    rad_s_r0  = numr0 *(0.00039 **4)    
    rad_s_r12 = numr1*(0.00035**4) 
    rad_s_r3  = numr3 *(0.00068**4)       

    # axial conductivity [cm^3/day] , mu is added later as it evolves with CST  
    beta = 0.9 #Thompson 2003a
    kz_l   = VascBundle_leaf * rad_s_l   * np.pi /8 * beta  
    kz_s   = VascBundle_stem * rad_s_s   * np.pi /8 * beta
    kz_r0  = VascBundle_root * rad_s_r0  * np.pi /8 * beta
    kz_r12 = VascBundle_root * rad_s_r12 * np.pi /8 * beta
    kz_r3  = VascBundle_root * rad_s_r3  * np.pi /8 * beta
    
    #radial conductivity [1/day],
    kr_l  = 0.
    kr_s  = 0.
    kr_r0 = 5e-2* 100
    kr_r1 = 5e-2* 100
    kr_r2 = 5e-2* 100
    kr_r3 = 5e-2 * 100
    l_kr =  100 #cm
       
    
    r.k_mucil = 1
    
    r.setKr_mucilRootTip([r.k_mucil, r.k_mucil, r.k_mucil, r.k_mucil] )
    r.setKr_stRootTip([[kr_r0,kr_r1 ,kr_r2 ,kr_r0],[kr_s,kr_s ],[kr_l]] )
    r.setKx_st([[kz_r0,kz_r12,kz_r12,kz_r0],[kz_s,kz_s ],[kz_l]], False)
    
    a_ST = [[0.00039,0.00035,0.00035,0.00039 ],[0.00019,0.00019],[0.00025]]
    Across_s_l   = numL*VascBundle_leaf *(a_ST[2][0]**2)*np.pi
    Across_s_s   = numS *VascBundle_stem * (a_ST[1][0]**2)*np.pi
    Across_s_r0  = numr0 *VascBundle_root * (a_ST[0][0]**2)*np.pi  
    Across_s_r12 = numr1*VascBundle_root * (a_ST[0][1]**2)*np.pi
    Across_s_r3  =  numr3 *VascBundle_root *(a_ST[0][2]**2)*np.pi
    Perimeter_s_l   = numL*VascBundle_leaf *(a_ST[2][0])* 2 * np.pi
    Perimeter_s_s   = numS *VascBundle_stem * (a_ST[1][0])* 2 * np.pi
    Perimeter_s_r0  = numr0 *VascBundle_root * (a_ST[0][0])* 2 * np.pi  
    Perimeter_s_r12 = numr1*VascBundle_root * (a_ST[0][1])* 2 * np.pi
    Perimeter_s_r3  =  numr3 *VascBundle_root *(a_ST[0][2])* 2 * np.pi
    r.setAcross_st([[Across_s_r0,Across_s_r12,Across_s_r12,Across_s_r0],[Across_s_s,Across_s_s],[Across_s_l]], False)
    return r
def setFunctionalParameters(r, weatherInit):
    r = setKrKx_phloem(r)
    r.g0 = 8e-6
    r.VcmaxrefChl1 =1.28
    r.VcmaxrefChl2 = 8.33
    r.a1 = 0.6/0.4
    r.a3 = 1.5
    r.alpha = 0.4 
    r.theta = 0.6 
    r.k_meso = 1e-3 
    
    
    r.Gr4Exud = True
    r.setKrm2([[2e-5]], False)
    r.setKrm1([[10e-2]], False) 
    r.setRhoSucrose([[0.5],[0.5],[0.5]], False)
    rootFact = 2
    stemFact = 2
    leafFact = 2
    r.setRmax_st([[2.4*rootFact,1.5*rootFact,0.6*rootFact,2.4*rootFact],
                  [2.*stemFact,2.*stemFact],
                  [8.*leafFact]], False)
    
    r.KMrm = 0.1#VERY IMPORTANT TO KEEP IT HIGH
    r.sameVolume_meso_st = False
    r.sameVolume_meso_seg = True
    r.withInitVal =True
    r.initValST = 0.4 
    r.initValMeso = 0.06 
    r.beta_loading = 0.6
    r.Vmaxloading = 0.05 #mmol/d, needed mean loading rate:  0.3788921068507634
    r.Mloading = 0.2
    r.Gr_Y = 0.8
    r.CSTimin = 0.2 
    r.CSTimin_exud = 0. 
    r.surfMeso=0.0025
    r.leafGrowthZone = 2 # cm
    r.StemGrowthPerPhytomer = True # 
                                                  
    
    #r.psi_osmo_proto = -10000*1.0197 #schopfer2006
    r.psiMin = 2500*1.0197
    r.fwr = 1e-16
    r.fw_cutoff =   0.04072 # 0.09497583
    r.sh =5e-4# 4.655e-4# 4e-4
    r.gm=0.03#0.01
    r.p_lcrit = -8375#-8847# from calibration using data of corso2020 #-15000*0.6
    
    r.limMaxErr = 1/100
    r.maxLoop = 10000
    r.minLoop=900
    
    #r.C_targ = r.CSTimin
    #r.C_targMesophyll = r.CSTimin
    #r.k_S_ST = 5/25 *2  
    #r.k_S_Mesophyll = 5/25*0   
                 
    r.C_targ = r.initValST#0.4#r.CSTimin
    r.C_targMesophyll = 0.06#r.CSTimin
    r.k_S_ST = 5/25 *100#*2 #daudet2002
    r.k_S_Mesophyll = 5/25*100 #daudet2002


    r.cs = weatherInit["cs"]

    r.expression = 6
    r.update_viscosity = True
    r.solver = 1
    r.atol = 1e-12
    r.rtol = 1e-8
    SPAD= 41.0
    chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
    r.Chl = np.array( [chl_]) 
    
    
    r.limMaxErr = 1/100
    r.maxLoop = 10000
    r.minLoop=900
    
    return r

def getVTPOutput(r, C_ST, psiXyl,Q_data, vtpindx, directoryN, simStartSim, condition):
    ana = pb.SegmentAnalyser(r.plant.mappedSegments())
    cutoff = 1e-16 #is get value too small, makes paraview crash
    mmolSuc_to_mumolC = 1e3*12
    Q_Exud_i, Q_Exud = Q_data[0]*mmolSuc_to_mumolC, Q_data[1]*mmolSuc_to_mumolC
    Q_Gr_i, Q_Grmax_i = Q_data[2]*mmolSuc_to_mumolC, Q_data[3]*mmolSuc_to_mumolC
    C_ST_p = C_ST*mmolSuc_to_mumolC
    C_ST_p[abs(C_ST_p) < cutoff] = 0
    psiXyl_p = psiXyl 
    psiXyl_p[abs(psiXyl_p) < cutoff] = 0
    #Q_Exud_i = Q_Exud_i 
    ana.addData("Q_Exud_i_nc",Q_Exud_i)
    ana.addData("Q_Exud_nc",Q_Exud)
    ana.addData("Q_Gr_i_nc",Q_Gr_i)
    ana.addData("Q_Grmax_i_nc",Q_Grmax_i)
    
    Q_Exud_i[abs(Q_Exud_i) < cutoff] = 0
    #Q_Exud = Q_Exud 
    Q_Exud[abs(Q_Exud) < cutoff] = 0
    Q_Gr_i[abs(Q_Gr_i) < cutoff] = 0
    Q_Grmax_i[abs(Q_Grmax_i) < cutoff] = 0
    subType = np.array(r.plant.subTypes)
    if False:
        ana.addData("CST", C_ST_p)

        ana.addData("psi_Xyl",psiXyl_p)
    ana.addData("Q_Exud_i",Q_Exud_i)
    ana.addData("Q_Exud",Q_Exud)
    ana.addData("Q_Gr_i",Q_Gr_i)
    ana.addData("Q_Grmax_i",Q_Grmax_i)

    ana.addData("subType",subType)
    #ana.write("results"+directoryN+"plot_"+str(int(simStartSim))+str(condition)+"at"+ str(vtpindx) +".vtp", 
    #          ["organType", "subType",
    #           "CST", "psi_Xyl"]) 

    vp.plot_plant(r.plant,p_name = [#"xylem pressure (cm)",#"sucrose concentration (mmol/cm3)", 
                                    'Q_Exud_i','Q_Exud',#'Q_Gr_i','Q_Grmax_i',
                                    "subType"],
                        vals =[ # psiXyl_p, C_ST_p,
                               Q_Exud_i,Q_Exud,#Q_Gr_i,Q_Grmax_i,
                               subType], 
                        filename = directoryN +'vtpvti/' +"plotpsi"+str(int(simStartSim))+str(condition)+"at"+ str(vtpindx), 
                  range_ = [0,5000])


    
def doSoilVTPplots(vtpindx,  plantModel, s,  
                   directoryN ='/'):
    """ vtp plots for water and C components in 3d soil and at root-soil
        interface 
    """
    soilTextureAndShape = getSoilTextureAndShape()
    min_b = soilTextureAndShape['min_b']
    max_b = soilTextureAndShape['max_b']
    cell_number = soilTextureAndShape['cell_number']
    
    swc= s.getWaterContent()
    watpotcm = theta2H(s.soil,swc)
    watpotcmbis= watpotcm + s.cellCenters[:,2]
    Skonz = s.getConcentration(1)
    Sq = Skonz * swc
    extraArray_ = [s.getWaterContent(), watpotcm, watpotcmbis ,Skonz,Sq]
    extraArrayName_ = ["theta (cm3/cm3)", "matric potential (cm)","water potential (cm)",
                       "solute concentration (mol/cm3)","solute content (mol)"]
    vp.plot_soil(s,  min_b, max_b,  cell_number, filename =  directoryN+ 'vtpvti/' +"soil"+str(vtpindx), 
                    pSoils = extraArray_, pnameSoils=extraArrayName_)
    

    

def getBiochemParam(s,paramIdx):    
    """ define TraiRhizo biochemical parameters 
        @param: the dumux soil object
        @ param: index of the TraiRhizo parameter set to use
        
    """
    # file containing the TraiRhizo parameter sets
    paramSet = pd.read_csv('./../fixedPointIter2/scripts/output_random_rows.csv').iloc[paramIdx].to_dict() # select one specific parameter set from index
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
    

def setIC3D(s, ICcc = None):
    return setIC(s, ICcc)

def setIC(s, ICcc = None):
    """ Defined the initial concentraition of the solutes
        [mol C / cm3 water] for disolved solutes and [mol C / cm3 scv]
        for solutes in the soil phase
        @param: s the dumux soil object
        @param: ICcc (optional) predefined initial conditions
    """
    if ICcc is None:
        paramSet = pd.read_csv('./../fixedPointIter2/scripts/output_random_rows.csv').loc[61]
        C_S = paramSet['CS_init'] /s.mg_per_molC## small C solutes in mol/cm3 water
        C_L = paramSet['CL_init'] /s.mg_per_molC## large C solutes in mol/cm3 water

        # concentraiton of adsobed C_S
        s.CSS_init  = getCSS(s, C_S) #mol C/ cm3 scv
        
            
        unitConversion = 1.0e6 # mol/cm3  => mol/m3 
        addedVar = 1. * float(s.doSoluteFlow) # empirical factor
        s.CSW_init = C_S * unitConversion
        s.ICcc = np.array([C_S *unitConversion*addedVar,
                           C_L*unitConversion*addedVar,
                            0.,
                            0.,
                            0.,
                            0.,
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


def pa2head(pa_, pnref = 1.e5, g = 9.81):
    """ pascal to pressure head, converts a numpy array """
    h = np.zeros(len(pa_))
    for i, p in enumerate(pa_):
        h[i] = (p - pnref) * 100. / 1000. / g
    return h


def head2pa(h_, pnref = 1.e5, g = 9.81):
    """ pressure to pascal, converts a numpy array """
    pa = np.zeros(len(h_))
    for i, h in enumerate(h_):
        pa[i] = pnref + h / 100.*1000.*g
    return pa

def getSoilTextureAndShape():  
    """ soil shape and texture data
        to adapt according to the soil represented
    """
    min_b = np.array([-3./2, -12./2, -40.]) # np.array( [5, 5, 0.] )
    max_b =np.array( [3./2, 12./2, 0.]) #  np.array([-5, -5, -5.])
    rez = 4
    cell_number = np.array((max_b - min_b)*rez ,dtype =int)#np.array( [1,1,1]) # 1cm3
    # #np.array([3,3,3])
    solidDensity = 2650 # [kg/m^3 solid] #taken from google docs TraiRhizo
    solidMolarMass = 60.08e-3 # [kg/mol] 
    # theta_r, theta_s, alpha, n, Ks
    #soilVG = [0.08, 0.43, 0.04, 1.6, 50]
    #vgSoil = [0.059, 0.45, 0.00644, 1.503, 1]
    loam = [0.08, 0.43, 0.04, 1.6, 50]
    sand = [0.045, 0.43, 0.15, 3, 1000] #
    soilTextureAndShape = {'min_b' : min_b,'max_b' : max_b,
                           'cell_number':cell_number,
                           "solidDensity":solidDensity,
                        'solidMolarMass': solidMolarMass,
                           'soilVG':sand}
    
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
                       
                       doBioChemicalReaction=True,
                     MaxRelativeShift = 1e-8, doOld = False):
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
    setIC3D(s,  ICcc)
    s.isPeriodic =True 
    s.createGrid(min_b, max_b, cell_number, s.isPeriodic)  # [cm] 
    if doOld:
        s.isPeriodic =False
    s = setupOther(s, p_mean_)
    
    #if rank == 0:
    #    s.base.printParams()
    
    # just print once as will not change during simulation
    write_file_array("cellVol", np.array(s.getCellVolumes()), directoryN =s.results_dir) # cm3 
    write_file_array("cellIds", np.array(s.cellIndices), directoryN =s.results_dir) # cm3
    
    cellcenters = s.getCellCenters()
    cellcentersX = []
    cellcentersY = []
    cellcentersZ = []
    for sub_array in cellcenters:
        cellcentersX.append(sub_array[0])
        cellcentersY.append(sub_array[1])
        cellcentersZ.append(sub_array[2])

    write_file_array("cellcentersX", np.array(cellcentersX), directoryN =results_dir) # cm3
    write_file_array("cellcentersY", np.array(cellcentersY), directoryN =results_dir) # cm3
    write_file_array("cellcentersZ", np.array(cellcentersZ), directoryN =results_dir) # cm3
    
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
    s.maxDt =  250./(3600.*24.)
    s.maxDt_1DS = s.maxDt # [s], lower maxDt for 1D models
    s.initializeProblem(s.maxDt)
    
    s.eps_regularization = 1e-10
    #s.eps_regularization = None # pcEps, krEps
    s.setRegularisation(s.eps_regularization, s.eps_regularization) # needs to be l
     
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
    s.ddt = 1.e-3  # [day] initial Dumux time step
    s.bulkMassErrorWater_rel = 0.
    s.bulkMassErrorWater_relLim = 0.    
    s.totC3dInit = sum(s.getTotCContent()) # mol    
    # initial soil water and solute content
    cell_volumes = s.getCellVolumes()  # cm3
    s.buWSoilInit = sum(np.multiply(np.array(s.getWaterContent()), cell_volumes)) # cm3 water
    return s



def compute3DSource(s, dt, sourcesInit):
    """
        compute source/sink of water and solutes in the 3d soil according to the changes in the 1d models
    """
    # when setting source: we limit according to the content at the beginning 
    # of the time step + max(net incoming flow, 0.) 
    # to get the maximum potential value of content that can be
    # in the voxel if needed

    for idComp in range(len(sourcesInit)):#cm3/day, mol/day
        SSL = sourcesInit[idComp].copy()
        if rank == 0:

            # compute maximum poentially available content for sink (cm3 water or mol)
            if idComp == 0:
                maxPotentialAvailable = (s.getWaterContent() - s.vg_soil.theta_R)  * s.CellVolumes

            for cellId, value in SSL.items(): 
                SSL[cellId] = max(SSL[cellId], -maxPotentialAvailable[cellId]/dt)    
        s.setSource(SSL,idComp)

def solve3DS(s, dt):
    """
        wrapper around the solving of the 3d soil flow
        TODO: reset maxDt after creating the solving object.
    """


    k_soil_solve = 0
    redoSolve = True
    maxRelShift = s.MaxRelativeShift

    # todo: probably to adapt
    assert dt >= 10/(24*3600)
    s.ddt = min(max(s.ddt, 1/(24*3600)), dt/10.)
    #s.maxDt = max(min(s.maxDt, dt/4.), s.ddt)
    assert dt > s.ddt
    assert s.maxDt > s.ddt
    #print('troubleshoot data',
    #      'time',s.ddt, s.maxDt, dt)
    #print('params',s.getParameter("Newton.MaxRelativeShift"),
    #      s.getParameter("Newton.MaxSteps"))
    while redoSolve:
        #s.ddt =min( 1.e-5,s.ddt)#or just reset to 1e-5?
        #s.setMaxTimeStepSize(s.maxDt) # reset each time

        try:
            decreaseMaxRelShift = False
            if rank==0:
                print("solve 3d soil", flush=True)
            s.solve(dt)  # in modules/solverbase.py

            #helpfull.run_with_timeout(60.*5,s.solve,dt) # time out error after Xmn
            if rank==0:
                print("solve 3d soil finished", flush=True)

            # if we got solute content < 0, throw exception
            solComp = [s.getSolution(ncom+1) for ncom in range(s.numSoluteComp)]
            whereError = None
            if rank == 0:
                whereError = [np.where(SC <0.) for SC in solComp]
                solComp = [min(SC) for SC in solComp]
            solComp = comm.bcast(solComp, root = 0)
            whereError = comm.bcast(whereError, root = 0)
            if min(solComp) <0.:
                print("min(solComp) <0.", rank, solComp, whereError)
                decreaseMaxRelShift = True
                raise Exception

            redoSolve = False
            ## solving succeded, reset solving parameters (in case  they were changed)
            # newton parameters are re-read at each 'solve()' calls
            s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift))# reset value
            s.setParameter("Newton.EnableResidualCriterion", "false") # sometimes helps, sometimes makes things worse
            s.setParameter("Newton.EnableAbsoluteResidualCriterion", "false")
            s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "false")

            s.setParameter("Newton.MaxSteps", "100")
            s.setParameter("Newton.MaxTimeStepDivisions", "10")
            s.createNewtonSolver() # re-create Newton solver to implement the new newton parameters

        except Exception as err:
            s.setParameter("Newton.EnableResidualCriterion", "false") # sometimes helps, sometimes makes things worse
            s.setParameter("Newton.EnableAbsoluteResidualCriterion", "false")
            s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "false")

            print(rank, f"Unexpected {err=}, {type(err)=}", 'k_soil_solve',k_soil_solve)

            if k_soil_solve > 6: # to avoid going to infinit loop
                raise Exception
            if k_soil_solve == 0: 
                # for the 1st fail, simply increase number of steps allowed
                s.setParameter("Newton.MaxSteps", "200")
                s.setParameter("Newton.MaxTimeStepDivisions", "100")
            elif k_soil_solve == 1: 
                # 2nd fail: try making the computation more precise
                print(rank,
                      'soil.solve() failed. making the computation more precise')
                s.setParameter("Newton.EnableResidualCriterion", "true") # sometimes helps, sometimes makes things worse
                s.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
                s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
                s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift/10.))# reset value
            else:
                # try decreasing MaxRelShift (if got solute content < 0)
                # try increasing maxrelshift (if got a non-convergence error)
                if decreaseMaxRelShift:
                    change = 0.1
                else:
                    change = 10
                print(rank,
                      'soil.solve() failed. NewtonMaxRelativeShift from',
                      maxRelShift,'to',maxRelShift*change)
                maxRelShift *= change
                # newton parameters are re-read at each 'solve()' calls
                s.setParameter("Newton.MaxRelativeShift", str(maxRelShift))
            s.reset() # reset solution vector
            s.createNewtonSolver() # re-create Newton solver to implement the new newton parameters
            k_soil_solve += 1