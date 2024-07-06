
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

class phloemDataStorage():
    def __init__(self, perirhizalModel, plantModel):
        self.plantModel = plantModel
        self.perirhizalModel = perirhizalModel
        self.plantModel = plantModel
        self.Nt = len(perirhizalModel.nodes)
        self.Q_ST_init = np.array([])
        self.C_S_meso    = np.zeros(self.Nt)
        self.C_meso    = np.zeros(self.Nt)
        self.C_S_ST    = np.zeros(self.Nt)
        self.C_ST    = np.zeros(self.Nt)
        self.Q_Exud    = np.zeros(self.Nt)
        self.Q_Mucil   = np.zeros(self.Nt)
        self.Q_Gr    = np.zeros(self.Nt)
        self.Q_Rm   = np.zeros(self.Nt)
        self.Q_Exudbu    = np.zeros(self.Nt)
        self.Q_Mucilbu   = np.zeros(self.Nt)
        self.Q_Grbu    = np.zeros(self.Nt)
        self.Q_Rmbu   = np.zeros(self.Nt)
        self.Q_in  = 0
        self.Q_Exud_cumul = 0.; 
        self.Q_Mucil_cumul = 0.
        self.Q_Exud_i =  np.zeros(self.Nt)
        self.Q_Mucil_i =  np.zeros(self.Nt)
        self.Q_Gr_i =  np.zeros(self.Nt)
        self.Q_Rm_i =  np.zeros(self.Nt)
        
        self.Q_Exud_i_seg = np.array([]); 
        self.Q_Mucil_i_seg = np.array([])  
        self.error_st_abs = 0;
        self.error_st_rel=0
        #  MMOL Suc (/cm3) => mol C (/cm3)
        self.mmolSuc_to_molC = 1/1e3*12
        
    def setNtbu(self):
        self.Ntbu = self.Nt
        
    def computePhloemFlow(self,rs_age, dt):
    
        self.plantModel.Csoil_seg = self.perirhizalModel.get_inner_solutes() * 1e3 # mol/cm3 to mmol/cm3 
        
        if rank == 0:
            startphloem=rs_age
            endphloem = rs_age + dt
            stepphloem = 1
            verbose_phloem = True
            filename =  self.perirhizalModel.results_dir +"inPM"+str(rank)+".txt"

            print("startpm",rank)
            self.plantModel.startPM(startphloem, endphloem, stepphloem, (  self.perirhizalModel.weatherX["TairC"]  +273.15) , verbose_phloem, filename)
            print("endtpm",rank)

            self.Nt = len( self.perirhizalModel.nodes)
            Nt = self.Nt
            
            #  MMOL Suc (/cm3) => mol C (/cm3)
            mmolSuc_to_molC = self.mmolSuc_to_molC
            if self.plantModel.withInitVal and (len(self.Q_ST_init) ==0) :
                self.Q_ST_init = np.array(self.plantModel.Q_init[0:Nt])* mmolSuc_to_molC
                self.Q_meso_init = np.array(self.plantModel.Q_init[Nt:(Nt*2)])* mmolSuc_to_molC

            # the backups

            self.Q_Exudbu = self.Q_Exud            
            self.Q_Mucilbu = self.Q_Mucil
            self.Q_Grbu = self.Q_Gr            
            self.Q_Rmbu = self.Q_Rm

            # att: that will be the cumulative value
            self.Q_ST    = np.array(self.plantModel.Q_out[0:Nt])* mmolSuc_to_molC
            self.Q_meso  = np.array(self.plantModel.Q_out[Nt:(Nt*2)])* mmolSuc_to_molC
            self.Q_Rm    = np.array(self.plantModel.Q_out[(Nt*2):(Nt*3)])* mmolSuc_to_molC
            self.Q_Exud  = np.array(self.plantModel.Q_out[(Nt*3):(Nt*4)])* mmolSuc_to_molC
            self.Q_Gr    = np.array(self.plantModel.Q_out[(Nt*4):(Nt*5)])* mmolSuc_to_molC
            self.Q_Rmmax       = np.array(self.plantModel.Q_out[(Nt*5):(Nt*6)])* mmolSuc_to_molC
            self.Q_Grmax       = np.array(self.plantModel.Q_out[(Nt*6):(Nt*7)])* mmolSuc_to_molC
            self.Q_S_meso   = np.array(self.plantModel.Q_out[(Nt*7):(Nt*8)])* mmolSuc_to_molC
            self.Q_S_ST   = np.array(self.plantModel.Q_out[(Nt*8):(Nt*9)])* mmolSuc_to_molC
            self.Q_Mucil  = np.array(self.plantModel.Q_out[(Nt*9):(Nt*10)])* mmolSuc_to_molC #mol for nodes

            self.C_ST    = np.array(self.plantModel.C_ST)* mmolSuc_to_molC
            self.Fl      = np.array(self.plantModel.Fl)* mmolSuc_to_molC
            self.volST   = np.array(self.plantModel.vol_ST)
            self.volMeso   = np.array(self.plantModel.vol_Meso)
            self.C_S_meso   = self.Q_S_meso/self.volMeso
            self.C_S_ST   = self.Q_S_ST/self.volST
            self.C_meso  = self.Q_meso/self.volMeso
            self.Q_in   += sum(np.array(self.plantModel.AgPhl)*dt)* mmolSuc_to_molC
            
            self.Q_out   = self.Q_Rm + self.Q_Exud + self.Q_Gr + self.Q_Mucil
            Q_S_ST_init = self.Q_ST_init # starch init value = concentration in phloem
            Q_S_meso_init = self.Q_meso_init # starch init value = concentration in mesophylle
            self.error_st_abs   = abs(sum(self.Q_ST + self.Q_meso + self.Q_out + self.Q_S_meso + self.Q_S_ST)- self.Q_in - sum(self.Q_ST_init) - sum(Q_S_ST_init) - sum(self.Q_meso_init) - sum(Q_S_meso_init))
            self.error_st_rel = abs(div0(self.error_st_abs,self.Q_in,1)*100)

            self.Q_Exudbu     =   np.concatenate((self.Q_Exudbu, np.full(self.Nt - self.Ntbu, 0.))) 
            self.Q_Mucilbu       =   np.concatenate((self.Q_Mucilbu, np.full(self.Nt - self.Ntbu, 0.))) 
            self.Q_Rmbu     =   np.concatenate((self.Q_Rmbu, np.full(self.Nt - self.Ntbu, 0.))) 
            self.Q_Grbu       =   np.concatenate((self.Q_Grbu, np.full(self.Nt - self.Ntbu, 0.))) 

            self.Q_Exud_i      = (self.Q_Exud    - self.Q_Exudbu)
            self.Q_Mucil_i     = (self.Q_Mucil   - self.Q_Mucilbu)
            self.Q_Gr_i      = (self.Q_Gr    - self.Q_Grbu)
            self.Q_Rm_i      = (self.Q_Rm    - self.Q_Rmbu)

            if False:
                # get error after longer nights: we do not get loss of C once mesophyll is empty
                try:
                    assert  (self.error_st_rel< 1.) or abs(self.Q_in) < 1e-13
                except:    
                    print('error_st_abs',rank,self.error_st_abs, self.Q_in, self.error_st_rel)
                    raise Exception

            assert self.Q_Exud_i[0] == 0#no exudation in seed node 
            assert self.Q_Mucil_i[0] == 0#no exudation in seed node  
            assert np.array(self.plantModel.Csoil_node)[0] == 0


            try:
                assert (np.array(self.plantModel.Csoil_seg ) == np.array(self.plantModel.Csoil_node)[1:]).all()
            except:
                print(np.array(self.plantModel.Csoil_seg ), np.array(self.plantModel.Csoil_node))
                print( (np.array(self.plantModel.Csoil_seg ) == np.array(self.plantModel.Csoil_node)[1:]),
                         (np.array(self.plantModel.Csoil_seg ) == np.array(self.plantModel.Csoil_node)[1:]).all())
                raise Exception

            self.Q_Exud_i_seg = np.array( self.Q_Exud_i[1:] )  #from nod to semgment
            self.Q_Mucil_i_seg = np.array(self.Q_Mucil_i[1:]) 

            airSegsId = self.perirhizalModel.airSegs

            write_file_array("Q_Exud_i_real",  self.Q_Exud_i, 
                             directory_ =self.perirhizalModel.results_dir)# to see if get val < 0
            write_file_array("Q_Mucil_i_real", self.Q_Mucil_i, 
                             directory_ =self.perirhizalModel.results_dir)# to see if get val < 0
            try:
                assert (self.Q_Exud_i_seg[airSegsId] == 0).all()
                assert (self.Q_Mucil_i_seg[airSegsId] == 0).all()
                assert (np.array(self.plantModel.k_mucil_)[airSegsId+1] == 0).all()
                assert (np.array(self.plantModel.Q_Exudmax)[airSegsId+1] == 0).all()
            except:
                print("Q_Exud_i_seg", self.Q_Exud_i_seg[airSegsId] )
                print("Q_Mucil_i", self.Q_Mucil_i_seg,Q_Mucil,self.Q_Mucilbu,airSegsId)
                print("Q_Mucil_i",self.Q_Mucil_i_seg[airSegsId], self.Q_Mucil[airSegsId+1], self.Q_Mucilbu[airSegsId+1])
                print("Csoil_seg", np.array(self.plantModel.Csoil_seg)[airSegsId])
                print("k_mucil_",self.plantModel.k_mucil_)
                print("Q_Exudmax",np.array(self.plantModel.Q_Exudmax)[airSegsId+1])
                print("airSegsId", airSegsId, np.where(airSegsId))
                print(len(airSegsId), len(self.plantModel.k_mucil_))
                raise Exception

            try: # currently, only have a release of carbon
                
                if ((min(self.Q_Exud_i_seg) < 0) or (min(self.Q_Mucil_i_seg)<0)):
                    write_file_float("timeNegativeExudMucil", rs_age, 
                                     directory_ =self.perirhizalModel.results_dir)
                    
            
                #assert min(self.Q_Exud_i_seg) >= -1e-13
                self.Q_Exud_i_seg[np.where(self.Q_Exud_i_seg<0)] = 0.
                #assert min(self.Q_Mucil_i_seg) >= -1e-13
                self.Q_Mucil_i_seg[np.where(self.Q_Mucil_i_seg<0)] = 0.
            except:
                print(self.C_ST, self.plantModel.Csoil_node, self.Q_Exud_i_seg,self.Q_Exud)
                raise Exception

            print("sum exud", sum(self.Q_Exud_i_seg), sum(self.Q_Mucil_i_seg))
        else:
            self.Q_Exud_i_seg = None
            self.Q_Mucil_i_seg = None
            self.Q_Gr_i = None;
            self.Q_Rm_i = None
            self.Q_Exud = None
            self.Q_Mucil = None 
            self.Q_Gr  = None;
            self.Q_Rm  = None
            
    def bcastData(self):
        self.Q_Exud = comm.bcast(self.Q_Exud, root = 0) 
        self.Q_Mucil = comm.bcast(self.Q_Mucil, root = 0) 
        self.Q_Exud_i_seg = comm.bcast(self.Q_Exud_i_seg, root = 0) 
        self.Q_Mucil_i_seg = comm.bcast(self.Q_Mucil_i_seg, root = 0) 
        
def computePhotosynthesis(plantModel, perirhizalModel,fpit_Helper, rs_age_i_dt, soilK):

    organTypes = np.asarray(plantModel.rs.organTypes, int)
    try:                    
        plantModel.solve_photosynthesis(sim_time_ = rs_age_i_dt, 
                    sxx_=fpit_Helper.rsx_input, 
                    cells_ = False,#(i == 0),#for 1st computation, use cell data
                    ea_ = perirhizalModel.weatherX["ea"],#not used
                    es_=perirhizalModel.weatherX["es"],#not used
                    verbose_ = False, doLog_ = False,
                    TairC_= perirhizalModel.weatherX["TairC"],#not used
                                soil_k_ = soilK, # [day-1]
                    outputDir_= "./results/rhizoplantExud")
        
        if (perirhizalModel.spellData['scenario'] == 'none') or ((perirhizalModel.spellData['scenario'] != 'baseline') and (rs_age_i_dt > perirhizalModel.spellData['spellStart']) and (rs_age_i_dt <= perirhizalModel.spellData['spellEnd'])):
            seg_fluxes = np.array(plantModel.outputFlux)# [cm3/day] 
        else:
            seg_fluxes = np.full(len(np.array(plantModel.outputFlux)),0.)



        if perirhizalModel.doPhotosynthesis and (rank == 0):
            leavesSegs = np.where(organTypes ==4)
            fluxes_leaves = seg_fluxes[leavesSegs]
            if (min(plantModel.Ev) < 0) or (min(plantModel.Jw) < 0) or (min(fluxes_leaves)<-1e-15):
                print("leaf looses water", min(plantModel.Ev),min(plantModel.Jw), min(fluxes_leaves))
                print("seg_fluxes",seg_fluxes,"leavesSegs", leavesSegs)                
                raise Exception
    except:
        plantModel.minLoop = 2
        plantModel.maxLoop = 5
        plantModel.solve_photosynthesis(sim_time_ = rs_age_i_dt, 
                    sxx_=fpit_Helper.rsx_input, 
                    cells_ = False,#(i == 0),#for 1st computation, use cell data
                    ea_ = perirhizalModel.weatherX["ea"],#not used
                    es_=perirhizalModel.weatherX["es"],#not used
                    verbose_ = True, doLog_ = True,
                    TairC_= perirhizalModel.weatherX["TairC"],#not used
                                soil_k_ = soilK, # [day-1]
                    outputDir_= ".")
        raise Exception
    
    if not perirhizalModel.doMinimumPrint:
        results_dir = perirhizalModel.results_dir
        write_file_array("fpit_Ev",np.array(plantModel.Ev),directory_ =results_dir, fileType = '.csv')
        write_file_array("fpit_Jw",np.array(plantModel.Jw),directory_ =results_dir, fileType = '.csv')
        write_file_array("fpit_fw",np.array(plantModel.fw),directory_ =results_dir, fileType = '.csv')#pg
        write_file_array("fpit_pg",np.array(plantModel.pg),directory_ =results_dir, fileType = '.csv')

        write_file_array("fpit_errPhoto", np.array(plantModel.maxErr) , directory_ =results_dir, fileType = '.csv') 
        write_file_array("fpit_errPhotoAbs", np.array(plantModel.maxErrAbs) , directory_ =results_dir, fileType = '.csv') 
        write_file_array("fpit_organTypes", organTypes, directory_ =results_dir, fileType = '.csv') 
    return seg_fluxes
        



def resistance2conductance(resistance,r, weatherX):
    resistance = resistance* (1/100) #[s/m] * [m/cm] = [s/cm]
    resistance = resistance * r.R_ph * weatherX["TairK"] / r.Patm # [s/cm] * [K] * [hPa cm3 K−1 mmol−1] * [hPa] = [s] * [cm2 mmol−1]
    resistance = resistance * (1000) * (1/10000)# [s cm2 mmol−1] * [mmol/mol] * [m2/cm2] = [s m2 mol−1]
    return 1/resistance

def computeAtmosphereData(plantModel, perirhizalModel):
    """
        story for photynthesis model the data needed for computation but
        no required as input parameter of the plantModel.solve() function
    """
    # atmosphereic pressure
    plantModel.Patm = perirhizalModel.weatherX["Pair"]

    ##resistances
    plantModel.g_bl = resistance2conductance(perirhizalModel.weatherX["rbl"],
                                             plantModel, perirhizalModel.weatherX) / plantModel.a2_bl
    plantModel.g_canopy = resistance2conductance(perirhizalModel.weatherX["rcanopy"],
                                                 plantModel, perirhizalModel.weatherX) / plantModel.a2_canopy
    plantModel.g_air = resistance2conductance(perirhizalModel.weatherX["rair"],
                                              plantModel, perirhizalModel.weatherX) / plantModel.a2_air     
    plantModel.Qlight = perirhizalModel.weatherX["Qlight"]
    
def noTranspiration(perirhizalModel, rs_age_i_dt, dt):
    """
        do we have NO transpiration?
        at night when doPhotosynthesis == True
    """
    if perirhizalModel.doPhotosynthesis:
        return (perirhizalModel.weatherX["Qlight"] == 0.)
    else:
        # TODO: check. depends how transpiration is computed for root systems
        return (sinusoidal2(rs_age_i_dt, dt) == 0)
    
def computeWaterFlow( fpit_Helper, perirhizalModel, plantModel, rs_age_i_dt, dt):
    """
        
        do plant water flow
        either photosynthesis or xylem flow
        return the plant-exterior water exchange
    """
    
    s = fpit_Helper.s # soil model
    seg_fluxes = fpit_Helper.seg_fluxes
    
    
    dist_factor =  perirhizalModel.getDeltaR() # distance between center of perihirzal zone inner cells and root surface [cm]
    
    # when there is no transpiration (night time), we use the plant wat. pot.
    # at the beginning of the time step (rsx_init). Otherwise does not converge
    if (perirhizalModel.beforeAtNight and noTranspiration(perirhizalModel, rs_age_i_dt, dt) ) or (fpit_Helper.n_iter > perirhizalModel.k_iter_2initVal) :
        fpit_Helper.rsx_input = fpit_Helper.rsx_init
        
    elif perirhizalModel.rsiCompMethod <= 1:
        fpit_Helper.rsx_input = fpit_Helper.rsx_old * ( perirhizalModel.rsiCompMethod) + fpit_Helper.rsx_init * (1. - perirhizalModel.rsiCompMethod)
    else:
        if fpit_Helper.n_iter < 4:
            fpit_Helper.rsx_input = fpit_Helper.rsx_old
            
        elif  perirhizalModel.rsiCompMethod == 2: # mean of the last 2 values
            fpit_Helper.rsx_input = (fpit_Helper.rsx_olds[-1] + fpit_Helper.rsx_olds[-2])/2
            
        elif  perirhizalModel.rsiCompMethod == 3: # mean(mean(oldvals[:-1]), oldvals[-1])            
            fpit_Helper.rsx_input = (np.array(fpit_Helper.rsx_olds[:-1]).mean(0)+ fpit_Helper.rsx_olds[-1])/2
            
        elif  perirhizalModel.rsiCompMethod == 4:# mean(oldvals)            
            fpit_Helper.rsx_input = np.array(fpit_Helper.rsx_olds).mean(0)
            
        else:
            raise Exception
                           
        
    if rank == 0:
        
        soilKIn =np.divide(vg.hydraulic_conductivity(fpit_Helper.rsx_input, 
                                                     perirhizalModel.vg_soil),
                           dist_factor) # if water enters the root
        if False:
            soilKIn_old =np.divide(vg.hydraulic_conductivity(fpit_Helper.rsx_old, 
                                                         perirhizalModel.vg_soil),
                               dist_factor) # if water enters the root
                               
            soilKIn_init =np.divide(vg.hydraulic_conductivity(fpit_Helper.rsx_init, 
                                                         perirhizalModel.vg_soil),
                               dist_factor) # if water enters the root
            soilKIn = soilKIn_old * ( perirhizalModel.rsiCompMethod) + soilKIn_init * (1. - perirhizalModel.rsiCompMethod)
            
        soilKOut = s.vg_soil.Ksat  /dist_factor# if water leaves the root # [cm/d]  / [cm]  = day-1
        
        fpit_Helper.soilK =  soilKIn 
        if( len(seg_fluxes) > 0.):# and not (perirhizalModel.beforeAtNight and (perirhizalModel.weatherX["Qlight"] == 0.)):
            fpit_Helper.soilK[np.where(seg_fluxes> 0. ) ] = soilKOut[np.where(seg_fluxes > 0. ) ] # where plant releases water
            #soilK[np.where(seg_fluxes< 0. ) ] = soilKIn[np.where(seg_fluxes < 0. ) ]

        if len(perirhizalModel.airSegs) > 0:   # infinit resistance for shoot segments and roots aboveground
            fpit_Helper.soilK[perirhizalModel.airSegs] = np.Inf



    if perirhizalModel.doPhotosynthesis:
        if (rank == 0):
            seg_fluxes = computePhotosynthesis(plantModel, perirhizalModel, fpit_Helper, rs_age_i_dt, fpit_Helper.soilK)
        else:
            seg_fluxes = None
    else: # just xylem flow
        if (rank == 0):
            transpiration = plantModel.maxTranspiration *  min(rs_age_i_dt/plantModel.maxTranspirationAge,1.)  *sinusoidal2(rs_age_i_dt, dt)
            
            rx = plantModel.solve(rs_age_i_dt, 
                         [-transpiration],
                         sx= fpit_Helper.rsx_input[0], 
                         sxx= fpit_Helper.rsx_input, 
                         cells = False, 
                         wilting_point = plantModel.wilting_point,
                          soil_k = fpit_Helper.soilK)
            plantModel.psiXyl = rx
            seg_fluxes = np.array(plantModel.segFluxes(simTime = rs_age_i_dt, 
            rx = list(rx), sx= list(fpit_Helper.rsx_input), 
                                               approx=False, cells=False, #approx, cells
                                               soil_k = list(fpit_Helper.soilK)))  #    [cm3 day-1] radial volumetric flow rate

        else :
            plantModel.psiXyl = None
            seg_fluxes = None
            
    fpit_Helper.seg_fluxes = comm.bcast(seg_fluxes,root=0)  # plant-exterior water exchanges
    plantModel.psiXyl = comm.bcast(plantModel.psiXyl, root = 0) # plant water potential
