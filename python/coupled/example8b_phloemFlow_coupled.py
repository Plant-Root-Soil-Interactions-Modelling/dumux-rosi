#!/usr/bin/env python
# coding: utf-8


import os
import sys
sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/")
sys.path.append("../../build-cmake/cpp/python_binding/")


# # Coupled carbon and water flow in CPlantBox (with a dynamic soil)
# 
# ## Simulation of water and carbon movement 
# 
# 
# In the following we will show how to compute the coupled water and carbon flow in the plant. 
# 
# 
# 
# We consider a dynamic plant and a dynamic soil. 
# To compute the carbon flux, we use the software developped by Lacointe et al. (2019).
# To compute the soil water flux, we use the software developped by Koch et al. (2021).
# 
# **Reference**
# 
# A Lacointe and P. Minchin. A mechanistic model to predict distribution of carbon among multiple sinks. *Methods in molecular biology* (Clifton, N.J.) vol. 2014, 2019.
# T Koch et al. DuMux 3 – an open-source simulator for solving flow and transport problems in porous media with a focus on model coupling,. *Computers & Mathematics with Applications* 2021	

# The sucrose flow depends on several plant, soil and atmospheric variables. For clarity, the basic functions defining those variables were moved to the file "parametersSucroseFlow".



sys.path.append("../../../CPlantBox/modelparameter/functional")
import plantbox as pb
from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
import numpy as np
import visualisation.vtk_plot as vp # for quick 3d vizualisations
import matplotlib.pyplot as plt
from functional.phloem_flux import PhloemFluxPython  
from plant_photosynthesis.wheat_FcVB_Giraud2023adapted import *
from plant_hydraulics.wheat_Giraud2023adapted import *
from plant_sucrose.wheat_phloem_Giraud2023adapted import *
from climate.dummyWeather import *
import functional.van_genuchten as vg

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import numpy as np


# ## 1. Define initial conditions




#we start with a small plant to have a lower computation time
simInit = 7 # [day] init simtime
simMax = 8
dt =2./24.
dtWater = 10./60./24.
assert dt%dtWater == 0 # only do full time steps, so need to have dt/dtWater = integer
depth = 60
weatherInit = weather(simInit)
simDuration = simInit


# plant system 
pl = pb.MappedPlant(seednum = 2) #seednum: gives the option of setting a random seed to make the simulations replicable
path = "../../../CPlantBox/modelparameter/structural/plant/"
name = "Triticum_aestivum_adapted_2023"
pl.readParameters(path + name + ".xml")

sdf = pb.SDF_PlantBox(np.Inf, np.Inf, depth )
pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil


verbose = False
pl.initialize(verbose )
pl.simulate(simInit, verbose)


#for post-processing
Q_out = 0 #sucrose lost by the plant
AnSum = 0 #assimilation
EvSum = 0 #transpiration
filename = "phloemoutputs.txt" 

Q_Rmbu      = np.array([0.])
Q_Grbu      = np.array([0.])
Q_Exudbu    = np.array([0.])
Q_STbu    = np.array([0.])

Q_Rmall      = np.array([])
Q_Grall      = np.array([])
Q_Exudall    = np.array([])
lengthTotall    = np.array([])
time = np.array([])
lengthTotInit = sum(pl.segLength())
lengthTotBU = sum(pl.segLength())


# ## 2. Define dynamic soil


min_b = np.array([-3./2, -12./2, -61.])#distance between wheat plants
max_b = np.array([3./2, 12./2, 0.])
cell_number = np.array([6, 24, 61] )#soil resolution
cell_size = (max_b - min_b)/cell_number
sri_distance = np.mean(cell_size)/2. # half of the cell size, assumed mean distance between deg. of freedom of soil and root surface
s = RichardsWrapper(RichardsSP())
s.vg_soil = vg.Parameters(weatherInit['vg']) 
s.initialize()
periodic = True
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setHomogeneousIC(weatherInit["p_mean"], True)  # cm pressure head, equilibrium
s.setTopBC("noFlux") # no rain nor evaporation
s.setBotBC("fleeFlow") # free flow at the bottom
s.setVGParameters([weatherInit['vg']]) # van Genuchten parameters
s.initializeProblem()
s.setCriticalPressure(-15000)
sx = s.getSolutionHead()  # inital condition, solverbase.py
picker = lambda x, y, z: s.pick([x, y, z])    
pl.setSoilGrid(picker)  # maps segment
            
# ## 3. create object to compute carbon and water flux
# The PhloemFluxPython class containes the functionalities of PhotosynthesisPython as well as the sucrose-related functions.



#give initial guess of leaf water potential and internal CO2 partial pressure (to start computation loop)
r = PhloemFluxPython(pl,psiXylInit = min(sx),ciInit = weatherInit["cs"]*0.5)


# ## 4. set other parameters and initial variable
# We present bellow some of the main sucrose-related parameters.



r = setPhotosynthesisParameters(r,weatherInit)

r = setKrKx_phloem(r) # conductivity of the sieve tube
r.setKrm2([[2e-5]]) #effect of the sucrose content on maintenance respiration 
r.setKrm1([[10e-2]]) #effect of structural sucrose content on maintenance respiration
r.setRhoSucrose([[0.51],[0.65],[0.56]])  #sucrose density per organ type (mmol/cm3)
r.setRmax_st([[14.4,9.0,6.0,14.4],[5.,5.],[15.]]) #maximum growth rate when water and carbon limitation is activated
r.KMfu = 0.11                                     #michaelis menten coefficient for usage of sucrose
r.beta_loading = 0.6 #feedback effect of sieve tube concentraiton on loading from mesophyll
r.Vmaxloading = 0.05 #mmol/d, max loading rate from mesophyll
r.Mloading = 0.2                                      #michaelis menten coefficient for loading of sucrose
r.Gr_Y = 0.8 # efficiency of sucrose usage for growth. if <1, we have growth respiration
r.CSTimin = 0.4 #minimum sucrose concentration below which no sucrose usage occures
r.Csoil = 1e-4 #mean soil concentration in sucrose


r.update_viscosity = True #update sucrose viscosity according to concentraiton ?
r.atol = 1e-12 #max absolute error for sucrose flow solver
r.rtol = 1e-8 #max relative error for sucrose flow solver


# ## 5. launch simulation
# The first time steps tends to require longer computation time. increasing the maximum errors allowed for the sucrose computation (r.atol, r.rtol) 
# and the minium and maximum plant segment length (dxMin, dx) can help decrease the computaiton time.



while simDuration <= simMax: 
    
    Nt = len(r.plant.nodes) 
    weatherX = weather(simDuration) #update weather variables
    r.Qlight = weatherX["Qlight"] #
    r = setKrKx_xylem(weatherX["TairC"], weatherX["RH"], r) #update xylem conductivity data

    # water flow (operator spliting)
    dtWatertot = 0
    Evdt =0
    Andt = 0
    while dtWatertot < dt:
        dtWatertot += dtWater
        ''' soil conductivity '''
        try:
            soil_k = np.array([vg.hydraulic_conductivity(sx[r.rs.seg2cell[seg_id]], s.vg_soil) if r.rs.seg2cell[seg_id] >= 0. else np.Inf for seg_id in range(Nt -1)]) #         
        except:
            soil_k = np.array([vg.hydraulic_conductivity(sx[r.rs.seg2cell[seg_id]][0], s.vg_soil) if r.rs.seg2cell[seg_id] >= 0. else np.Inf for seg_id in range(Nt -1)]) #         
        soil_k[np.where(np.array(r.outputFlux)> 0. ) ] = s.vg_soil.Ksat
        soil_k /= sri_distance
        
        """ photosynthesis """
        r.solve_photosynthesis(sim_time_ = simDuration, sxx_=sx, cells_ = True,ea_ = weatherX["ea"],es_=weatherX["es"],
            verbose_ = False, doLog_ = False,TairC_= weatherX["TairC"],outputDir_= "./results",
            soil_k_ = soil_k)
        
        """ dumux """   
        fluxesSoil = r.soilFluxes(simDuration, r.psiXyl, sx, approx=False) # plant water uptake
        
        s.setSource(fluxesSoil.copy())  # richards.py # send plant water uptake to soil solver
        s.solve(dtWater)
        sx = s.getSolutionHead()  # richards.py  

        min_sx, min_rx, max_sx, max_rx = np.min(sx), np.min(r.psiXyl), np.max(sx), np.max(r.psiXyl)
        n = round((simDuration- simInit)/(simMax-simInit) * 100.)

        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm plant, {:g} cm root collar  at {:g} days"
                .format(min_sx, max_sx, min_rx, max_rx, r.psiXyl[0], s.simTime))

        Andt += np.array(r.Ag4Phloem)*dtWater # get cumulative assimilation during the time step
        Evdt += np.sum(r.Ev)*dtWater # get cumulative transpiration during the time step

    AnSum += sum(Andt)
    EvSum += Evdt
    errLeuning = sum(r.outputFlux) #should be 0 : no storage of water in the plant
    r.Ag4Phloem = Andt/dt
    
    #simulation of phloem flow
    startphloem= simDuration
    endphloem = startphloem + dt
    stepphloem = 1
    r.startPM(startphloem, endphloem, stepphloem, ( weatherX["TairC"]  +273.15) , True, filename)
        
    #get ouput of sucrose flow computation    
    Q_ST    = np.array(r.Q_out[0:Nt])          #sieve tube sucrose content
    Q_meso  = np.array(r.Q_out[Nt:(Nt*2)])     #mesophyll sucrose content
    Q_Rm    = np.array(r.Q_out[(Nt*2):(Nt*3)]) #sucrose used for maintenance respiration
    Q_Exud  = np.array(r.Q_out[(Nt*3):(Nt*4)]) #sucrose used for exudation
    Q_Gr    = np.array(r.Q_out[(Nt*4):(Nt*5)]) #sucrose used for growth and growth respiration
    
    C_ST    = np.array(r.C_ST)           #sieve tube sucrose concentraiton
    volST   = np.array(r.vol_ST)         #sieve tube volume
    volMeso   = np.array(r.vol_Meso)      #mesophyll volume     
    C_meso  = Q_meso/volMeso              #sucrose concentration in mesophyll
    Q_out   = Q_Rm + Q_Exud + Q_Gr       #total sucrose lost/used by the plant
    error   = sum(Q_ST + Q_meso + Q_out )- AnSum  #balance residual (error)
    
    lengthTot = sum(r.plant.segLength()) #total plant length 
    
    #variation of sucrose content at the last time step (mmol)
    Q_ST_i        = Q_ST      - Q_STbu #in the sieve tubes
    Q_Rm_i        = Q_Rm      - Q_Rmbu #for maintenance
    Q_Gr_i        = Q_Gr      - Q_Grbu #for growth
    Q_Exud_i      = Q_Exud    - Q_Exudbu #for exudation
    Q_out_i       = Q_Rm_i    + Q_Exud_i      + Q_Gr_i #total usage
    
    #print some outputs
    print("\n\n\n\t\tat ", int(np.floor(simDuration)),"d", int((simDuration%1)*24),"h, PAR:",  round(r.Qlight *1e6),"mumol m-2 s-1")
    print("Error in sucrose balance:\n\tabs (mmol) {:5.2e}\trel (-) {:5.2e}".format(error, div0f(error,AnSum, 1.)))
    print("Error in water balance:\n\tabs (cm3/day) {:5.2e}".format(errLeuning))
    print("water loss (cm3):\n\ttranspiration {:5.2e}".format(EvSum))
    print("assimilated sucrose (mmol)\tAn {:5.2e}".format(AnSum)) 
    print("sucrose concentration in sieve tube (mmol ml-1):\n\tmean {:.2e}\tmin  {:5.2e} at {:d} segs \tmax  {:5.2e}".format(np.mean(C_ST), min(C_ST), len(np.where(C_ST == min(C_ST) )[0]), max(C_ST)))        
    print('cumulated \tRm   {:.2e}\tGr   {:.2e}\tExud {:5.2e}'.format(sum(Q_Rm), sum(Q_Gr), sum(Q_Exud)))
    print("aggregated sink repartition at last time step (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rm_i)/sum(Q_out_i)*100, 
         sum(Q_Gr_i)/sum(Q_out_i)*100,sum(Q_Exud_i)/sum(Q_out_i)*100))
    print("total aggregated sink repartition (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rm)/sum(Q_out)*100, 
         sum(Q_Gr)/sum(Q_out)*100,sum(Q_Exud)/sum(Q_out)*100))
    print("growth rate (cm/day)\ttotal {:5.2e}\tlast time step {:5.2e}".format(lengthTot - lengthTotInit, lengthTot - lengthTotBU))      
    
    #plant growth based on Gr * Gr_Y
    r.plant.simulate(dt, verbose)
    simDuration += dt
    
    #for post processing
    Ntbu = Nt
    Nt = len(r.plant.nodes)
    lengthTotBU = lengthTot
    Q_STbu       =   np.concatenate((Q_ST, np.full(Nt - Ntbu, 0.)))
    Q_Rmbu       =   np.concatenate((Q_Rm, np.full(Nt - Ntbu, 0.)))
    Q_Grbu       =   np.concatenate((Q_Gr, np.full(Nt - Ntbu, 0.))) 
    Q_Exudbu     =   np.concatenate((Q_Exud, np.full(Nt - Ntbu, 0.))) 
    
    
    Q_Rmall    = np.append( Q_Rmall  ,sum(Q_Rm_i))
    Q_Grall    = np.append( Q_Grall  ,sum(Q_Gr_i))
    Q_Exudall  = np.append( Q_Exudall,sum(Q_Exud_i))
    lengthTotall  = np.append( lengthTotall,lengthTot)
    time       = np.append( time ,simDuration)
    


# ## 8. plot some results



fig, axs = plt.subplots(2,2)
axs[0,0].plot(time, Q_Rmall/dt)
axs[0,0].set(xlabel='day of growth', ylabel='total Rm rate (mmol/day)')
axs[1,0].plot(time, Q_Grall/dt, 'tab:red')
axs[1,0].set(xlabel='day of growth', ylabel='total Gr rate (mmol/day)')
axs[0,1].plot(time, Q_Exudall/dt , 'tab:brown')
axs[0,1].set(xlabel='day of growth', ylabel='total exudation\nrate (mmol/day)')
axs[1,1].plot(time, lengthTotall , 'tab:green')
axs[1,1].set(xlabel='day of growth', ylabel='total plant\nlength (cm)')
fig.tight_layout()
plt.show()


# ## Take away messages
# 
# * Basic idea how to use the class *PhloemFlow*
# * The plant growth follows the rate of carbon usage for growth
