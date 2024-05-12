""" 
    plant + 3d soil + perirhizal zone for water and/or solute(s)  
"""
import matplotlib; matplotlib.use('agg')
import sys;
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import timeit
import numpy as np

sys.path.append("../modules/");
sys.path.append("../../../../CPlantBox/");
sys.path.append("../../../../CPlantBox/src")

import plantbox as pb  # CPlantBox
import vtk_plot_adapted as vp
import scenario_setup
import rhizo_modelsPlant  # Helper class for cylindrical rhizosphere models
from rhizo_modelsPlant import *
import cyl3plant as cyl3
import helpfull
from helpfull import write_file_array, write_file_float, div0, div0f, suggestNumStepsChange
from helpfull import continueLoop
from weatherFunctions import *
from PhloemPhotosynthesis import *
import printData

    

def XcGrowth(initsim, simMax,paramIndx_,spellData):
    # outer time step (outside of fixed-point iteration loop)
    dt = 20/60/24
    # max number of iteration for the fixed-point iteration
    k_iter = 20
    # plant or root system
    organism = "plant"# "RS"#
    # which functional modules to implement
    doSoluteFlow = True # only water (False) or with solutes (True)
    noAds = False # stop adsorption?
    doPhloemFlow = True
    doPhotosynthesis = True # photosynthesis-Transpiration (True) or just xylem flow (False)?
    # when there is no transpiration, we use the plant wat. pot.
    # at the beginning of the time step. Otherwise does not converge
    beforeAtNight = True  
    # static or growing organism
    static_plant = False
    # print debug messages in cyl3plant faile
    mpiVerbose = False
    # print debug messages in rhizo_mappedPlant
    mpiVerboseInner = False
    # how many files are printed. use 'False' in debug mode
    # ATT: for short ismulations only
    doMinimumPrint =  True
    # use moles (mol) and not mass (g) in dumux
    usemoles = True
    
    # @see PhloemPhotosynthesis::computeWaterFlow().
    # TODO: make weather dependent?
    maxTranspiration = 6. # cm3/day, used if not doPhotosynthesis 
    
    # get initial variables and parameters for plant and soil setup
    soilTextureAndShape = scenario_setup.getSoilTextureAndShape()
    weatherInit = scenario_setup.weather(1., spellData)
       
    # directory where the results will be printed
    results_dir="./results/ongoingCleanUp/"+str(spellData['scenario'])\
    +"_"+str(int(np.prod(soilTextureAndShape['cell_number'])))\
                    +"_"+str(paramIndx_)\
                    +"_"+str(int(initsim))+"to"+str(int(simMax))\
                    +"_"+str(int(dt*24*60))+"mn_"\
                    +str(int((dt*24*60 - int(dt*24*60))*60))+"s_"\
                    +str(max_rank)+"/"
    
    # to get printing directory/simulaiton type in the slurm.out file
    if rank == 0:
        print('results_dir',results_dir, flush = True)
    
    # if results directory already exists, make it empty
    helpfull.setupDir(results_dir)
                    
    """ initialize """

    s, soil = scenario_setup.create_soil_model(usemoles = usemoles, 
                                               results_dir = results_dir, 
                                               p_mean_ = weatherInit['p_mean'], 
                                        paramIndx=paramIndx_,
                                        noAds = noAds)

    if organism == "plant":
        path = "../../../../CPlantBox/modelparameter/structural/plant/"
        xml_name = "Triticum_aestivum_test_2021.xml"  # root growth model parameter file
    elif organism == "RS":
        path = "./"
        xml_name = "data_magda/Zeamays_synMRI_modified.xml"
    elif organism == "RSstatic":
        path = "./"
        xml_name = "../../grids/RootSystem8.rsml"
    else:
        raise Exception
    


    # all thread need a plant object, but only thread 0 will make it grow
    perirhizalModel, plantModel = scenario_setup.create_mapped_plant(initsim, s, xml_name,
                                            path, plantType = organism, 
                                            usemoles = usemoles,
                                            limErr1d3d = 5e-13, spellData = spellData)  

    # store parameters
    plantModel.maxTranspiration = maxTranspiration
    perirhizalModel.doPhotosynthesis = doPhotosynthesis
    perirhizalModel.doPhloemFlow = doPhloemFlow
    perirhizalModel.doSoluteFlow = doSoluteFlow
    perirhizalModel.doMinimumPrint = doMinimumPrint
    perirhizalModel.k_iter=k_iter
    perirhizalModel.beforeAtNight = beforeAtNight
    rs_age = initsim # age before implementation of the funcional modules
    s.mpiVerbose = mpiVerbose
    s.mpiVerboseInner = mpiVerboseInner
    perirhizalModel.mpiVerbose = mpiVerbose
    perirhizalModel.mpiVerboseInner = mpiVerboseInner
    perirhizalModel.results_dir = results_dir
    
    # used to define dry spells
    perirhizalModel.spellData = spellData
    perirhizalModel.enteredSpell = (perirhizalModel.spellData['scenario'] == 'none') or (perirhizalModel.spellData['scenario'] == 'baseline')
    perirhizalModel.leftSpell = (perirhizalModel.spellData['scenario'] == 'baseline')
    

    
    """ define initial value for loop"""
    # initial time step for the fixed-point iteration loop
    perirhizalModel.dt_inner = float(dt) # [d]
    start = True # first loop
    perirhizalModel.new_soil_solute = np.array([0.])
    
    
    # inter-cell exchanges of solute and water in 3D soil
    net_sol_flux =  np.array([np.array([]),np.array([])])
    net_flux = np.array([])
    # plant-soil water exchange
    seg_Wfluxes = np.array([])

    # cumulative transpiration
    plantModel.TranspirationCumul = 0 # real cumulative transpiration
    plantModel.TranspirationCumul_eval = 0 # cumulative transpiration during period with dinamic soil (for mass balance check)
    
    # initial soil water and solute content
    cell_volumes = s.getCellVolumes()  # cm3
    s.totC3dInit = sum(s.getTotCContent()) # mol
    s.buWSoilInit = sum(np.multiply(np.array(s.getWaterContent()), cell_volumes)) # cm3 water

    # to get time spent in each part of the code. TODO: currently, does not work well
    plantModel.time_start_global = timeit.default_timer()
    plantModel.time_rhizo_cumul = 0
    plantModel.time_3ds_cumul = 0
    plantModel.time_plant_cumulW = 0
    plantModel.time_plant_cumulS = 0
    plantModel.time_rhizo_i = 0
    plantModel.time_3ds_i = 0
    plantModel.time_plant_cumul = 0
    
    
    """ phloem variable """
    if perirhizalModel.doPhloemFlow:
        phloemData = phloemDataStorage(perirhizalModel, plantModel) # to store data and run phloem simulation
    
    
    """ prints """
    printData.initialPrint(perirhizalModel)
    
        
    while rs_age < simMax:

        rs_age += dt
        
        if perirhizalModel.doPhloemFlow:
            phloemData.setNtbu()
        
        if (rank == 0) and (not static_plant) :
            
            seg2cell_old = perirhizalModel.seg2cell
            perirhizalModel.simulate(dt, verbose = False)
            helpfull.checkseg2cellMapping(seg2cell_old, perirhizalModel)
            
                
        if (rank == 0):
            # print plant shape data for post-processing
            printData.printPlantShape(perirhizalModel,plantModel)
            

        perirhizalModel.update() # update shape data in the rhizosphere model
        
        if start: # for first loop, do extra printing to have initial error
            perirhizalModel.checkMassOMoleBalance2( diff1d3dCW_abs_lim = 1e-13) # beginning: should not be any error
            printData.printTimeAndError(perirhizalModel, rs_age)
            start = False
        # print differences between 1d and 3d soil models
        # and content in 1d and 3d
        printData.printDiff1d3d(perirhizalModel, s)             

        if perirhizalModel.doPhloemFlow:
            # go from current to cumulative exudation and mucilage release
            phloemData.Q_Exud_cumul += sum(phloemData.Q_Exud_i_seg); 
            phloemData.Q_Mucil_cumul += sum(phloemData.Q_Mucil_i_seg)

        perirhizalModel.n_iter, failedLoop, keepGoing = helpfull.resetAndSaveData1(perirhizalModel)
        
        while keepGoing or failedLoop:
            
            helpfull.resetAndSaveData2(plantModel, perirhizalModel, s)
            
            net_sol_flux, net_flux, seg_Wfluxes, real_dtinner,failedLoop, n_iter_inner_max = cyl3.simulate_const(s, 
                                                    plantModel, 
                                                    sim_time= dt,dt= perirhizalModel.dt_inner, 
                                                    rs_age=rs_age, 
                                                    Q_plant=[phloemData.Q_Exud_i_seg, 
                                                             phloemData.Q_Mucil_i_seg],
                                                    perirhizalModel= perirhizalModel,plantType = organism,
                                                    outer_R_bc_sol = net_sol_flux, 
                                                    outer_R_bc_wat = net_flux,
                                                   seg_fluxes=seg_Wfluxes,
                                                     doMinimumPrint = doMinimumPrint)
            
            keepGoing = continueLoop(perirhizalModel,perirhizalModel.n_iter, 
                                     perirhizalModel.dt_inner,
                                     failedLoop,real_dtinner,name="Outer_data",
                                     isInner = False,
                             plant = plantModel)
            perirhizalModel.n_iter += 1
            
            # we leave the loop either becuse dumux threw an error or because
            # we reached dt
            try:
                assert (abs(real_dtinner - dt) < perirhizalModel.dt_inner ) or failedLoop                
            except:
                print('real_dtinner',real_dtinner ,dt, perirhizalModel.dt_inner , failedLoop)
                write_file_array("real_dtinner_error", np.array([real_dtinner ,dt,
                                                                 perirhizalModel.dt_inner , failedLoop,abs((real_dtinner - dt)/dt*100.),rs_age]), 
                                 directory_ =results_dir, fileType = '.csv') 
                raise Exception
                
                
            perirhizalModel.dt_inner = suggestNumStepsChange(dt, perirhizalModel.dt_inner, failedLoop,
                                             np.ceil(perirhizalModel.k_iter/2), results_dir)
            
            if keepGoing or failedLoop:
                
                if rank==0:
                    print("error too high, decrease dt to",perirhizalModel.dt_inner)
                # reset data to the beginning of the iteration loops
                helpfull.resetAndSaveData3(plantModel, perirhizalModel, s)
                
                
        # manually set n_iter to 0 to see if continueLoop() still yields fales.
        # if continueLoop() = True, we had a non-convergence error 
        keepGoing_ = continueLoop(perirhizalModel,0, perirhizalModel.dt_inner,
                                  failedLoop,real_dtinner,
                                  name="TestOuter_data", plant = plantModel)
        if keepGoing_:
            print("None convergence error: only left the loop because reached the max number of iteration.")
            raise Exception  
        
        printData.printTimeAndError(perirhizalModel, rs_age)
        
        helpfull.getCumulativeTranspiration(plantModel, perirhizalModel)
        
        printData.getAndPrintErrorRates(perirhizalModel, plantModel, s, phloemData)
        
        printData.printCylData(perirhizalModel)
        
                    
        plantModel.time_start_plant = timeit.default_timer()
        if ((not static_plant) or (rs_age == initsim+dt)) and doPhloemFlow:
            phloemData.computePhloemFlow(rs_age, dt, perirhizalModel)        
        plantModel.time_plant_cumulS += (timeit.default_timer() - plantModel.time_start_plant)
            
        if doPhloemFlow:
            phloemData.bcastData()


        if (rank == 0) and doPhloemFlow:
            printData.printPhloemFlowOutput(rs_age, perirhizalModel, phloemData, plantModel)
            
    """ wrap up """
    print('finished simulation :D')
    
    
    perirhizalModel.checkMassOMoleBalance2()
    
    printData.doVTPplots(perirhizalModel, s, soilTextureAndShape)
    
    
    print("fin", rank)
    # for some reason, we get an error when the program stops
    # memory leak?
    return results_dir

        

if __name__ == '__main__': #TODO. find a way to reset maxDt after creating the solving object.
    # python3 XcGrowth.py 9 10 0 customDry 9.02 0.02
    # python3 XcGrowth.py 9 10 20 lateDry
    # python3 XcGrowth.py 12 25 98 baseline
    # python3 XcGrowth.py 9 10 0 none 
    
    if rank == 0:
        print('sys.argv',sys.argv)
        
    initsim =float(sys.argv[1])
    
    simMax = initsim + 3.
    if len(sys.argv)>2:
        simMax = float(sys.argv[2])
    paramIndx_base = 0
    if len(sys.argv)>3:
        paramIndx_base = int(sys.argv[3])
        
    
    doProfile = False
    
    scenario = "none"
    if len(sys.argv)>4:
        scenario = sys.argv[4]
    
    if scenario == "none":
        spellStart = 0 #not used
        condition = "wet"
        spellDuration =  np.Inf
    elif scenario == "baseline":
        spellStart = np.Inf
        condition = "wet"
        spellDuration =   7
    elif scenario == "earlyDry":
        spellStart = 11 
        condition = "dry"
        spellDuration =   7
    elif scenario == "lateDry":
        spellStart = 18
        condition = "dry"
        spellDuration =   7
    elif scenario == "customDry":
        spellStart = float(sys.argv[5])
        condition = "dry"
        spellDuration =   float(sys.argv[6])
    elif scenario == "customWet":
        spellStart = float(sys.argv[5])
        condition = "wet"
        spellDuration =   float(sys.argv[6])
    else :
        print("scenario", scenario,"not recognised")
        raise Exception
        
    spellEnd = spellStart + spellDuration
    spellData = {'scenario':scenario,'spellStart':spellStart,'spellEnd':spellEnd, 'condition':condition}
    
    if doProfile:
        import cProfile
        import pstats, io
        pr = cProfile.Profile()
        pr.enable()
    results_dir = XcGrowth(initsim, simMax,paramIndx_base,spellData )
    if doProfile:
        pr.disable()
        filename = results_dir+'profile'+str(rank)+'.prof' 
        s = io.StringIO()
        ps = pstats.Stats(pr, stream=s).sort_stats('tottime')
        ps.dump_stats(filename)