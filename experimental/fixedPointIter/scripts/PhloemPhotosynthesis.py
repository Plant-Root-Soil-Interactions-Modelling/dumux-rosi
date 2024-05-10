
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

class phloemDataStorage():
    def __init__(self, PerhirizalModel, phloemModel):
        self.phloemModel = phloemModel
        self.PerhirizalModel = PerhirizalModel
        self.Nt = len(PerhirizalModel.nodes)
        self.Q_ST_init = np.array([])
        self.Q_Exud    = np.zeros(self.Nt)
        self.Q_Mucil   = np.zeros(self.Nt)
        self.Q_Exudbu    = np.zeros(self.Nt)
        self.Q_Mucilbu   = np.zeros(self.Nt)
        self.Q_in  = 0
        self.Q_Exud_inflate = 0.; 
        self.Q_Mucil_inflate = 0.
        self.Q_Exud_i = None
        self.Q_Exud_i_seg = np.array([]); 
        self.Q_Mucil_i_seg = np.array([])  
        self.error_st_abs = 0;
        self.error_st_rel=0
    def setNtbu(self):
        self.Ntbu = self.Nt
    def computePhloemFlow(self,rs_age, dt,perirhizalModel):
        
        if rank == 0:
            startphloem=rs_age
            endphloem = rs_age + dt
            stepphloem = 1
            verbose_phloem = True
            filename = perirhizalModel.results_dir +"inPM"+str(rank)+".txt"

            print("startpm",rank)
            self.phloemModel.startPM(startphloem, endphloem, stepphloem, ( perirhizalModel.weatherX["TairC"]  +273.15) , verbose_phloem, filename)
            print("endtpm",rank)

            self.Nt = len(perirhizalModel.nodes)
            Nt = self.Nt
            
            if self.phloemModel.withInitVal and (len(self.Q_ST_init) ==0) :
                self.Q_ST_init = np.array(self.phloemModel.Q_init[0:Nt])/1e3
                self.Q_meso_init = np.array(self.phloemModel.Q_init[Nt:(Nt*2)])/1e3

            # the backups

            self.Q_Exudbu = self.Q_Exud            
            self.Q_Mucilbu = self.Q_Mucil

            # att: that will be the cumulative value
            #  MMOL(/cm3) => mol(/cm3)
            self.Q_ST    = np.array(self.phloemModel.Q_out[0:Nt])/1e3
            self.Q_meso  = np.array(self.phloemModel.Q_out[Nt:(Nt*2)])/1e3
            self.Q_Rm    = np.array(self.phloemModel.Q_out[(Nt*2):(Nt*3)])/1e3
            self.Q_Exud  = np.array(self.phloemModel.Q_out[(Nt*3):(Nt*4)])/1e3 
            self.Q_Gr    = np.array(self.phloemModel.Q_out[(Nt*4):(Nt*5)])/1e3
            self.Q_Rmmax       = np.array(self.phloemModel.Q_out[(Nt*5):(Nt*6)])/1e3
            self.Q_Grmax       = np.array(self.phloemModel.Q_out[(Nt*6):(Nt*7)])/1e3
            self.Q_S_meso   = np.array(self.phloemModel.Q_out[(Nt*7):(Nt*8)])/1e3
            self.Q_S_ST   = np.array(self.phloemModel.Q_out[(Nt*8):(Nt*9)])/1e3
            self.Q_Mucil  = np.array(self.phloemModel.Q_out[(Nt*9):(Nt*10)])/1e3 #mol for nodes

            self.C_ST    = np.array(self.phloemModel.C_ST)/1e3
            self.Fl      = np.array(self.phloemModel.Fl)/1e3
            self.volST   = np.array(self.phloemModel.vol_ST)
            self.volMeso   = np.array(self.phloemModel.vol_Meso)
            self.C_S_meso   = self.Q_S_meso/self.volMeso
            self.C_S_ST   = self.Q_S_ST/self.volST
            self.C_meso  = self.Q_meso/self.volMeso
            self.Q_in   += sum(np.array(self.phloemModel.AgPhl)*dt)/1e3
            # i m missing the starch
            self.Q_out   = self.Q_Rm + self.Q_Exud + self.Q_Gr + self.Q_Mucil
            self.error_st_abs   = abs(sum(self.Q_ST + self.Q_meso + self.Q_out + self.Q_S_meso + self.Q_S_ST)- self.Q_in - sum(self.Q_ST_init)  - sum(self.Q_meso_init))
            self.error_st_rel = abs(div0(self.error_st_abs,self.Q_in,1)*100)

            self.Q_Exudbu     =   np.concatenate((self.Q_Exudbu, np.full(self.Nt - self.Ntbu, 0.))) 
            self.Q_Mucilbu       =   np.concatenate((self.Q_Mucilbu, np.full(self.Nt - self.Ntbu, 0.))) 

            self.Q_Exud_i      = (self.Q_Exud    - self.Q_Exudbu)
            self.Q_Mucil_i     = (self.Q_Mucil   - self.Q_Mucilbu)


            try:
                assert  (self.error_st_rel< 1.) or abs(self.Q_in) < 1e-13
            except:    
                print('error_st_abs',rank,self.error_st_abs, self.Q_in, self.error_st_rel)
                raise Exception

            assert self.Q_Exud_i[0] == 0#no exudation in seed node 
            assert self.Q_Mucil_i[0] == 0#no exudation in seed node  
            assert np.array(self.phloemModel.Csoil_node)[0] == 0


            try:
                assert (np.array(self.phloemModel.Csoil_seg ) == np.array(self.phloemModel.Csoil_node)[1:]).all()
            except:
                print(np.array(self.phloemModel.Csoil_seg ), np.array(self.phloemModel.Csoil_node))
                print( (np.array(self.phloemModel.Csoil_seg ) == np.array(self.phloemModel.Csoil_node)[1:]),
                         (np.array(self.phloemModel.Csoil_seg ) == np.array(self.phloemModel.Csoil_node)[1:]).all())
                raise Exception

            self.Q_Exud_i_seg = np.array( self.Q_Exud_i[1:] )  #from nod to semgment
            self.Q_Mucil_i_seg = np.array(self.Q_Mucil_i[1:]) 

            airSegsId = self.PerhirizalModel.airSegs

            try:
                assert (self.Q_Exud_i_seg[airSegsId] == 0).all()
                assert (self.Q_Mucil_i_seg[airSegsId] == 0).all()
                assert (np.array(self.phloemModel.k_mucil_)[airSegsId+1] == 0).all()
                assert (np.array(self.phloemModel.Q_Exudmax)[airSegsId+1] == 0).all()
            except:
                print("Q_Exud_i_seg", self.Q_Exud_i_seg[airSegsId] )
                print("Q_Mucil_i", self.Q_Mucil_i_seg,Q_Mucil,self.Q_Mucilbu,airSegsId)
                print("Q_Mucil_i",self.Q_Mucil_i_seg[airSegsId], self.Q_Mucil[airSegsId+1], self.Q_Mucilbu[airSegsId+1])
                print("Csoil_seg", np.array(self.phloemModel.Csoil_seg)[airSegsId])
                print("k_mucil_",self.phloemModel.k_mucil_)
                print("Q_Exudmax",np.array(self.phloemModel.Q_Exudmax)[airSegsId+1])
                print("airSegsId", airSegsId, np.where(airSegsId))
                print(len(airSegsId), len(self.phloemModel.k_mucil_))
                raise Exception

            try:
                assert min(self.Q_Exud_i_seg) >= -1e-13
                self.Q_Exud_i_seg[np.where(self.Q_Exud_i_seg<0)] = 0.
            except:
                print(self.C_ST, self.phloemModel.Csoil_node, self.Q_Exud_i_seg,self.Q_Exud)
                raise Exception

            print("sum exud", sum(self.Q_Exud_i_seg), sum(self.Q_Mucil_i_seg))
        else:
            self.Q_Exud_i_seg = None
            self.Q_Mucil_i_seg = None
            self.Q_Exud = None
            self.Q_Mucil = None 
            
    def bcastData(self):
        self.Q_Exud = comm.bcast(self.Q_Exud, root = 0) 
        self.Q_Mucil = comm.bcast(self.Q_Mucil, root = 0) 
        self.Q_Exud_i_seg = comm.bcast(self.Q_Exud_i_seg, root = 0) 
        self.Q_Mucil_i_seg = comm.bcast(self.Q_Mucil_i_seg, root = 0) 
        

def phloemParam(r,weatherInit ):
    """ define parameters for phloem flow and photosynthesis
    """
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
    r.setRhoSucrose([[0.51],[0.65],[0.56]], False)
    rootFact = 2
    r.setRmax_st([[2.4*rootFact,1.5*rootFact,0.6*rootFact,2.4*rootFact],[2.,2.],[8.]], False)
    
    r.KMrm = 0.1#VERY IMPORTANT TO KEEP IT HIGH
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

    r.expression = 6
    r.update_viscosity = True
    r.solver = 1
    r.atol = 1e-12
    r.rtol = 1e-8
    SPAD= 41.0
    chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
    r.Chl = np.array( [chl_]) 
    
    
    r.minLoop = 1000
    r.maxLoop = 5000
    
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
    
    #radial conductivity [1/day],
    kr_l  = 0.
    kr_s  = 0.
    kr_r0 = 5e-2
    kr_r1 = 5e-2
    kr_r2 = 5e-2
    kr_r3 = 5e-2
    l_kr =  100 #cm
    
    r.setKr_st([[kr_r0,kr_r1 ,kr_r2 ,kr_r0],[kr_s,kr_s ],[kr_l]] , kr_length_= l_kr, verbose = False)
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
