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
sys.path.append("../inputDataTraiRhizo/");
sys.path.append("../../../../CPlantBox/");
sys.path.append("../../../../CPlantBox/src")

sys.path.append("../../../build-cmake/cpp/python_binding/");
import plantbox as pb  # CPlantBox
import vtk_plot_adapted as vp

import rhizo_modelsPlant  # Helper class for cylindrical rhizosphere models
from rhizo_modelsPlant import *
import cyl3plant as cyl3
import helpfull
from helpfull import write_file_array, write_file_float, div0, div0f, suggestNumStepsChange
from helpfull import continueLoop
import weatherFunctions 
from PhloemPhotosynthesis import *
import printData
import scenario_setup

    

def XcGrowth(initsim, simMax,paramIndx_,spellData): 
    path = "../../../../CPlantBox/modelparameter/structural/plant/"
    xml_name = "Triticum_aestivum_test_2021.xml"  # root growth model parameter fileroot growth model parameter file
    MaxRelativeShift = 1e-8 #if paramIndx_ != 44 else 1e-10
    # outer time step (outside of fixed-point iteration loop)
    dt = 20/60/24
    dt_inner_init =  dt # 1/60/24 #
    # min, max, objective number of iteration for the fixed-point iteration
    minIter = 4 # empirical minimum number of loop to reduce error
    k_iter_2initVal = 131 # max num of iteration for loops
    k_iter = 100 # max num of iteration for loops
    targetIter= 90# target n_iter for adjusting time step of inner loop
    # which functional modules to implement
    doSoluteFlow = True # only water (False) or with solutes (True)
    noAds = False # stop adsorption?
    doPhloemFlow = True
    doPhotosynthesis = True # photosynthesis-Transpiration (True) or just xylem flow (False)?
    # when there is no transpiration, we use the plant wat. pot.
    # at the beginning of the time step. Otherwise does not converge
    do1d1dFlow = False
    beforeAtNight = False #True  
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
    
    rsiCompMethod =4# 0.5 = mixed, 0 = init, 1 = end, -1: use mean of the last 2 iterations
    # 2 : last two mean method
    # 3 : mean(oldmean, lastval)
    # 4 : mean(allvals)
    # 5: mean + weighing factor to see where min(abs(, for each seg, not just the whole setinput - real). 
    # for eahc seg, not just the whole set
    # 6: keep min obtained ksoil
    # 7 : get mean/integrated krsoil from dumux. likewise, use mean/integrated psi gotten in dumux.
    # or run plant model
    # 8 call cpb in dumux ? inner fixed-point iteration?
    
    # @see PhloemPhotosynthesis::computeWaterFlow().
    # TODO: make weather dependent?
    maxTranspiration = 12#6. # cm3/day, used if not doPhotosynthesis 
    maxTranspirationAge = 25. # day, age at which trans = maxTRans
    
    # get initial variables and parameters for plant and soil setup
    soilTextureAndShape = scenario_setup.getSoilTextureAndShape()
    weatherInit = weatherFunctions.weather(1.,dt, spellData)
       
    # directory where the results will be printed
    results_dir="./results/TraiRhizo/paperSc/"+str(rsiCompMethod*10)+str(spellData['scenario'])\
    +"_"+str(int(np.prod(soilTextureAndShape['cell_number'])))\
                    +"_"+str(paramIndx_)\
                    +"_"+str(int(initsim))+"to"+str(int(simMax))\
                    +"_"+str(int(dt_inner_init*24*60))+"mn_"\
                    +str(int((dt_inner_init*24*60 - int(dt_inner_init*24*60))*60))+"s_"\
                    +str(max_rank)+"/"
    
    # to get printing directory/simulaiton type in the slurm.out file
    if rank == 0:
        print('results_dir',results_dir, flush = True)
    
    # if results directory already exists, make it empty
    helpfull.setupDir(results_dir)
                    
    """ initialize """

    s = scenario_setup.create_soil_model(usemoles = usemoles, 
                                               results_dir = results_dir, 
                                               p_mean_ = weatherInit['p_mean'], 
                                        paramIndx=paramIndx_,
                                        noAds = noAds, doSoluteFlow = doSoluteFlow,
                                        MaxRelativeShift = MaxRelativeShift)

    


    # all thread need a plant object, but only thread 0 will make it grow
    perirhizalModel, plantModel = scenario_setup.create_mapped_plant(initsim, s, xml_name,
                                            path, 
                                            doPhloemFlow = doPhloemFlow,
                                            static_plant = static_plant,
                                            usemoles = usemoles,
                                            limErr1d3d = 5e-12, spellData = spellData)  

    # store parameters
    plantModel.maxTranspiration = maxTranspiration
    plantModel.maxTranspirationAge = maxTranspirationAge
    
    perirhizalModel.do1d1dFlow = do1d1dFlow
    perirhizalModel.getSoilTextureAndShape = scenario_setup.getSoilTextureAndShape
    perirhizalModel.k_iter_2initVal = k_iter_2initVal
    perirhizalModel.rsiCompMethod = rsiCompMethod
    perirhizalModel.doPhotosynthesis = doPhotosynthesis
    perirhizalModel.doPhloemFlow = doPhloemFlow
    perirhizalModel.doSoluteFlow = doSoluteFlow
    perirhizalModel.doMinimumPrint = doMinimumPrint
    perirhizalModel.dt_inner = dt_inner_init#float(dt) # [d]
    perirhizalModel.minIter=minIter
    perirhizalModel.k_iter=k_iter
    perirhizalModel.targetIter=targetIter
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
    start = True # first loop
    
    # inter-cell exchanges of solute and water in 3D soil
    net_sol_flux =  np.array([np.array([]),np.array([])])
    net_flux = np.array([])
    # plant-soil water exchange
    seg_Wfluxes = np.array([])
    
    """ phloem variable """
    phloemData = phloemDataStorage(perirhizalModel, plantModel) # to store data and run phloem simulation
    
    
    """ prints """
    printData.initialPrint(perirhizalModel)
    
    printData.doVTPplots(0, perirhizalModel, plantModel,s, soilTextureAndShape, 
                            datas=[], datasName=[],initPrint=True)
    
    
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
            perirhizalModel.check1d3dDiff( diff1d3dCW_abs_lim = 1e-13) # beginning: should not be any error
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
            
            if perirhizalModel.doPhloemFlow:
                Q_plant_=[phloemData.Q_Exud_i_seg, 
                        phloemData.Q_Mucil_i_seg]
            else:
                Q_plant_=[[],[]]
            
            net_sol_flux, net_flux, seg_Wfluxes, real_dtinner,failedLoop, n_iter_inner_max = cyl3.simulate_const(s, 
                                                    plantModel, 
                                                    sim_time= dt,dt= perirhizalModel.dt_inner, 
                                                    rs_age=rs_age, 
                                                    Q_plant=Q_plant_,
                                                    perirhizalModel= perirhizalModel,
                                                    outer_R_bc_sol = net_sol_flux, 
                                                    outer_R_bc_wat = net_flux,
                                                   seg_fluxes=seg_Wfluxes,
                                                     doMinimumPrint = doMinimumPrint)
            
            keepGoing = continueLoop(perirhizalModel=perirhizalModel,s=s,n_iter=perirhizalModel.n_iter, 
                                     dt_inner=perirhizalModel.dt_inner,
                                     failedLoop=failedLoop,real_dtinner=real_dtinner,name="Outer_data",
                                     isInner = False,doPrint = True, 
                         fileType = '.csv' ,
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
                
                
            perirhizalModel.dt_inner = suggestNumStepsChange(dt,  failedLoop, perirhizalModel, n_iter_inner_max)
            
            if keepGoing or failedLoop:
                
                if rank==0:
                    print("error too high, decrease dt to",perirhizalModel.dt_inner)
                # reset data to the beginning of the iteration loops
                
                helpfull.resetAndSaveData3(plantModel, perirhizalModel, s)
                
        # manually set n_iter to 0 to see if continueLoop() still yields fales.
        # if continueLoop() = True, we had a non-convergence error 
        
        keepGoing_ = continueLoop(perirhizalModel=perirhizalModel,s=s,n_iter=0, 
                                     dt_inner=perirhizalModel.dt_inner,
                                     failedLoop=failedLoop,real_dtinner=real_dtinner,name="TestOuter_data",
                                     isInner = False,doPrint = True, 
                         fileType = '.csv' ,
                             plant = plantModel)
        if keepGoing_:
            print("None convergence error: only left the loop because reached the max number of iteration.")
            raise Exception  
        
        printData.printTimeAndError(perirhizalModel, rs_age)
        
        helpfull.getCumulativeTranspirationAg(plantModel, perirhizalModel, dt)
        
        printData.getAndPrintErrorRates(perirhizalModel, plantModel, s, phloemData)
        
        printData.printCylData(perirhizalModel,rs_age )
        
                    
        plantModel.time_start_plant = timeit.default_timer()
        if ((not static_plant) or (rs_age == initsim+dt)) and doPhloemFlow:
            phloemData.computePhloemFlow(rs_age, dt)        
        plantModel.time_plant_cumulS += (timeit.default_timer() - plantModel.time_start_plant)
            
        if doPhloemFlow:
            phloemData.bcastData()


        if (rank == 0):
            printData.printOutput(rs_age, perirhizalModel, phloemData, plantModel)
            
        if int(rs_age *1000)/1000-int(rs_age) == 0.5 :# midday (TODO: change it to make it work for all outer time step)
            if rank == 0:
                datas = [
                         plantModel.psiXyl, 
                         phloemData.C_ST, phloemData.C_S_ST, 
                         phloemData.C_meso, phloemData.C_S_meso, 
                         phloemData.Q_Exud_i, phloemData.Q_Mucil_i, 
                         phloemData.Q_Gr_i,phloemData.Q_Rm_i
                        ]
                datasName = [ "psiXyl",
                             "C_ST", "C_S_ST", 
                             "C_meso", "C_S_meso", 
                             "Q_Exud", "Q_Mucil",
                             "Q_Gr","Q_Rm",
                             "Q_Exud_i","Q_Mucil_i" ,
                             "Q_Gr_i","Q_Rm_i"
                            ]
            else:
                datas = []
                datasName = []

            printData.doVTPplots(int(rs_age*10), #indx/number of vtp plot
                                perirhizalModel, plantModel,s, soilTextureAndShape, 
                                datas, datasName, initPrint=False, doSolutes = perirhizalModel.doSoluteFlow)
            
    """ wrap up """
    
    
    perirhizalModel.check1d3dDiff()
    
    
    
    print("finished simulation :D", rank,"(parting is such sweet sorrow)")
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