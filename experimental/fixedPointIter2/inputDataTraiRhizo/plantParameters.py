


import numpy as np
import pandas as pd
import timeit
import pandas as pd
import matplotlib; matplotlib.use('agg')
import sys;
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import timeit
from helpfull import write_file_array, write_file_float, div0, div0f
from functional.xylem_flux import sinusoidal2

import functional.van_genuchten as vg

def phloemParam(r,weatherInit ):
    """ define parameters for phloem flow and photosynthesis
    """
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
    r.CSTimin_exud = 0.2 
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
    kr_r0 = 5e-2
    kr_r1 = 5e-2
    kr_r2 = 5e-2
    kr_r3 = 5e-2
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
    


def init_conductivities(r, TairC:float = 20):
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
    
    ratio_decrease = 1.5/100
    
    r.setKrTables([[[kr_r0,kr_r0,kr_r0*ratio_decrease,kr_r0*ratio_decrease],
                    [kr_r1,kr_r1,kr_r1*ratio_decrease,kr_r1*ratio_decrease],
                    [kr_r2,kr_r2,kr_r2*ratio_decrease,kr_r2*ratio_decrease],
                    [kr_r0,kr_r0,kr_r0*ratio_decrease,kr_r0*ratio_decrease]],
                    [[kr_s],[kr_s] ],[[kr_l]]],
            [[[0,0.8,1,10000],[0,0.8,1,10000],[0,0.8,1,10000],[0,0.8,1,10000]],
            [[0],[0]],[[0]]],verbose = False, ageBased = False) 
            
    r.setKx([[kz_r0,kz_r1,kz_r2,kz_r0],[kz_s,kz_s ],[kz_l]])
    
    
    Rgaz=8.314 #J K-1 mol-1 = cm^3*MPa/K/mol
    rho_h2o = dEauPure/1000#g/cm3
    Mh2o = 18.05 #g/mol
    MPa2hPa = 10000
    hPa2cm = 1/0.9806806
    return r