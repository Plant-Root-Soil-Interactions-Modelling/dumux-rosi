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
import visualisation.vtk_plot as vp
import scenario_setup
import rhizo_modelsPlant  # Helper class for cylindrical rhizosphere models
from rhizo_modelsPlant import *
import cyl3plant as cyl3
import helpfull
from helpfull import write_file_array, write_file_float, div0, div0f, suggestNumStepsChange
from helpfull import continueLoop
from weather import *
from XylemFlux import *
from PhloemPhotosynthesis import *
import printData

    

def XcGrowth(initsim, mode,simMax,extraName,paramIndx_,spellData):
    # outer time step (outside of fixed-point iteration loop)
    dt = 20/60/24
    # max number of iteration for the fixed-point iteration
    k_iter = 20
    # plant or root system
    organism = "plant"# "RS"#
    # which functional modules to implement
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
    noAds = False # stop adsorption?
    # how many files are printed. use 'False' in debug mode
    # ATT: for short ismulations only
    doMinimumPrint =  True
    # use moles (mol) and not mass (g) in dumux
    usemoles = True
    
    # get initial variables and parameters for plant and soil setup
    soilTextureAndShape = scenario_setup.getSoilTextureAndShape()
    weatherInit = scenario_setup.weather(1., spellData)
       
    # directory where the results will be printed
    results_dir="./results/newMucil4p/"+extraName+str(spellData['scenario'])\
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

    s, soil = scenario_setup.create_soil_model(demoType = mode, 
                                               usemoles = usemoles, 
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
    perirhizalModel, plantModel = scenario_setup.create_mapped_plant( mode,initsim, s, xml_name,
                                            path, plantType = organism, 
                                            usemoles = usemoles,
                                            limErr1d3d = 5e-13, spellData = spellData)  

    # store parameters
    perirhizalModel.doPhotosynthesis = doPhotosynthesis
    perirhizalModel.doPhloemFlow = doPhloemFlow
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
    dt_inner = float(dt) # [d]
    start = True # first loop
    perirhizalModel.new_soil_solute = np.array([0.])
    
    
    # inter-cell exchanges of solute and water in 3D soil
    net_sol_flux =  np.array([np.array([]),np.array([])])
    net_flux = np.array([])
    # plant-soil water exchange
    seg_fluxes_ = np.array([])

    # cumulative transpiration
    plantModel.TranspirationCumul = 0 # real cumulative transpiration
    plantModel.TranspirationCumul_eval = 0 # cumulative transpiration during period with dinamic soil (for mass balance check)
    
    # initial soil water and solute content
    cell_volumes = s.getCellVolumes()  # cm3
    buTotCSoilInit = sum(s.getTotCContent()) # mol
    buWSoilInit = sum(np.multiply(np.array(s.getWaterContent()), cell_volumes)) # cm3 water

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
            perirhizalModel.checkMassOMoleBalance2(None,None, dt,
                                    seg_fluxes =None, diff1d3dCW_abs_lim = 1e-13, takeFlux = False) # very beginning: should not be any error
            printData.printTimeAndError(perirhizalModel, rs_age)
            start = False
            
        printData.printDiff1d3d(perirhizalModel, s) # print differences between 1d and 3d soil models
        
            

        if perirhizalModel.doPhloemFlow:
            # go from current to cumulative exudation and mucilage release
            phloemData.Q_Exud_inflate += sum(phloemData.Q_Exud_i_seg); 
            phloemData.Q_Mucil_inflate += sum(phloemData.Q_Mucil_i_seg)

        perirhizalModel.rhizoMassWError_abs = 1.
        perirhizalModel.rhizoMassCError_abs = 1.
        perirhizalModel.errDiffBCs = 1.
        perirhizalModel.err = 1.
        perirhizalModel.max_err = 1.
        perirhizalModel.diff1d3dCurrant_rel = 1e6
        perirhizalModel.maxdiff1d3dCurrant_rel = 1e6
        
        n_iter = 0
        failedLoop = False # had non-convergence error in dumux
        keepGoing = True # stay in the fixed point iteration
        
        while keepGoing or failedLoop:
            plantModel.TranspirationCumul_inner = 0 # reset transpiration of inner time step to 0
            s.saveManual()
            perirhizalModel.saveManual()
            perirhizalModel.leftSpellBU = perirhizalModel.leftSpell
            perirhizalModel.enteredSpellBU = perirhizalModel.enteredSpell
            net_sol_flux_, net_flux_, seg_fluxes__, real_dtinner,failedLoop, n_iter_inner_max = cyl3.simulate_const(s, 
                                                    plantModel, 
                                                    sim_time= dt,dt= dt_inner, rs_age=rs_age, 
                                                    Q_plant=[phloemData.Q_Exud_i_seg, 
                                                             phloemData.Q_Mucil_i_seg],
                                                    r= perirhizalModel,plantType = organism,
                                                    outer_R_bc_sol = net_sol_flux, 
                                                    outer_R_bc_wat = net_flux,seg_fluxes=seg_fluxes_,
                                                    outer_n_iter = n_iter,
                                                     doMinimumPrint = doMinimumPrint)
            keepGoing = continueLoop(perirhizalModel,n_iter, dt_inner,failedLoop,real_dtinner,name="Outer_data",isInner = False,
                             plant = plantModel)
            n_iter += 1
            try:
                assert (abs(real_dtinner - dt) < dt_inner ) or failedLoop                
            except:
                print('real_dtinner',real_dtinner ,dt, dt_inner , failedLoop)
                write_file_array("real_dtinner_error", np.array([real_dtinner ,dt, dt_inner , failedLoop,abs((real_dtinner - dt)/dt*100.),rs_age]), 
                                 directory_ =results_dir, fileType = '.csv') 
                raise Exception
                
                
            nOld = int(dt/dt_inner)
            if failedLoop:# if the inner computation failed, we also need to decrease the timestep
                n_iter_inner_max = perirhizalModel.k_iter
            nNew = suggestNumStepsChange(nOld, n_iter_inner_max, np.ceil(perirhizalModel.k_iter/2), results_dir)
            try:
                assert nOld == int(nOld)
                assert nNew == int(nNew)
            except:
                print('nNew iisue',nNew , int(nNew), nOld,dt,dt_inner,(nOld == int(nOld)), (nNew == int(nNew)))
                
                write_file_array("nNew_error", np.array([nNew , int(nNew), nOld,dt,dt_inner,(nOld == int(nOld)), (nNew == int(nNew)),
                                                         real_dtinner ,dt, dt_inner , failedLoop,abs((real_dtinner - dt)/dt*100.),rs_age]), 
                                 directory_ =results_dir, fileType = '.csv') 
                
            dt_inner = dt/float(nNew)
            if keepGoing or failedLoop:
                
                if mpiVerbose and rank==0:# or (max_rank == 1):
                    print(rank, "error too high, decrease N, dt from", nOld, dt/float(nOld),"to",nNew, dt_inner)
                s.resetManual()# <= that actually only works for the last solve. I need a better method to reset to the beginning of the whole stuff
                perirhizalModel.resetManual()
                perirhizalModel.leftSpell = perirhizalModel.leftSpellBU
                perirhizalModel.enteredSpell = perirhizalModel.enteredSpellBU
                
                perirhizalModel.checkMassOMoleBalance2(None,None, dt,
                                    seg_fluxes =None, diff1d3dCW_abs_lim = np.Inf, takeFlux = False)
                # to get correct error values for sumDiff1d3dCW_relOld
                
        
        keepGoing_ = continueLoop(perirhizalModel,0, dt_inner,failedLoop,real_dtinner,name="TestOuter_data", plant = plantModel)
        if keepGoing_:
            raise Exception # only left the loop because reached the max number of iteration.
        
        
        net_sol_flux = net_sol_flux_
        net_flux = net_flux_ 
        seg_fluxes_ = seg_fluxes__
                
        printData.printTimeAndError(perirhizalModel, rs_age)
        

           
        write_file_array("OuterSuccess_error", perirhizalModel.errs, directory_ =results_dir, fileType = '.csv') 
        write_file_array("OuterSuccess_data", np.array([n_iter, perirhizalModel.err, perirhizalModel.rhizoMassWError_abs,dt_inner]), directory_ =results_dir, fileType = '.csv')
        
        if organism == "plant":
            plantModel.TranspirationCumul += plantModel.TranspirationCumul_inner #sum(np.array(plantModel.Ev) * dt) #transpiration [cm3/day] * day
            if perirhizalModel.enteredSpell and (not perirhizalModel.leftSpell):
                plantModel.TranspirationCumul_eval += plantModel.TranspirationCumul_inner
            elif perirhizalModel.leftSpell:
                plantModel.TranspirationCumul_eval = 0.
                
        else:
            plantModel.TranspirationCumul += sum(plantModel.outputFlux)
        

        buTotCAfter = sum(s.getTotCContent())   #0 get stuck here

        buWAfter = sum(np.multiply(np.array(s.getWaterContent()), cell_volumes))    

        if rank == 0:
            if (mode != "dumux_w"):
                s.bulkMassErrorCumul_abs = abs((buTotCAfter - ( buTotCSoilInit + phloemData.Q_Exud_inflate + phloemData.Q_Mucil_inflate)))#so that only works if infalte
                if buTotCAfter != 0:
                    s.bulkMassErrorCumul_rel = abs(s.bulkMassErrorCumul_abs/buTotCAfter*100)
                else:
                    s.bulkMassErrorCumul_rel =np.nan
            else:
                s.bulkMassErrorCumul_abs = np.nan
                s.bulkMassErrorCumul_rel =np.nan
            # because we have a free flow BC at the bottom, that could increase the error
            # ideally, I should add the flow at the bellow BC here
            s.bulkMassErrorWaterCumul_abs = abs(buWAfter - ( buWSoilInit - plantModel.TranspirationCumul_eval))
            s.bulkMassErrorWaterCumul_rel = abs(s.bulkMassErrorWaterCumul_abs/buWAfter*100)
        else:
            s.bulkMassErrorCumul_abs = None
            s.bulkMassErrorCumul_rel = None
            s.bulkMassErrorWaterCumul_abs = None
            s.bulkMassErrorWaterCumul_rel = None
            
        if (mode != "dumux_w"):
            #absolute and relative (%) error
            write_file_array("errorsBulkSoil", np.array([s.bulkMassCErrorPlant_abs, s.bulkMassCErrorPlant_rel, #not cumulative 
                                                s.bulkMassCError1ds_abs, s.bulkMassCError1ds_rel, 
                                                s.bulkMassErrorCumul_abs,s.bulkMassErrorCumul_rel,#cumulative
                                                s.bulkMassErrorWater_abs,s.bulkMassErrorWater_rel, #not cumulative
                                                s.bulkMassErrorWaterCumul_abs,s.bulkMassErrorWaterCumul_rel]),
                             directory_ =results_dir, fileType = '.csv')#cumulative
            write_file_array("errorMassRhizo", np.array([perirhizalModel.rhizoMassCError_abs, perirhizalModel.rhizoMassCError_rel,
                                                        perirhizalModel.rhizoMassWError_abs, perirhizalModel.rhizoMassWError_rel]), directory_ =results_dir)# not cumulativecumulative (?)

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
    
    # otherwise the shoot looks weird
    periodic = False
    
    for konz in np.array([True]):#, False]):
        extraArray_ = perirhizalModel.soilModel.getWaterContent()
        if not konz:
            extraArrayName_ = "theta (cm3)"
            extraArray_ *= cell_volumes
        else:
            extraArrayName_ = "theta (cm3/cm3)"
        print("idcomp", 0, rank)
        rp = perirhizalModel.get_concentration(0, konz)
        s.results_dir = results_dir

        # adapt to have plot_plants_and soil
        vp.plot_roots_and_soil(perirhizalModel.mappedSegments(),extraArrayName_,rp, s, periodic, 
                               soilTextureAndShape['min_b'],
                               soilTextureAndShape['max_b'],
                               soilTextureAndShape['cell_number'], 
                filename="soil_rx",sol_ind =-1,extraArray = extraArray_, extraArrayName = extraArrayName_,
                interactiveImage=False)  # VTK vizualisation
        print("idcomp_done", 0, rank)
        for i in range(1, perirhizalModel.numComp):
            print("idcomp", i, rank)
            extraArray_ = perirhizalModel.soilModel.getSolution(i) * perirhizalModel.phaseDensity(i)/1e6

            if not konz:
                extraArrayName_ = "C"+str(i)+" mol"
            else:
                extraArrayName_ = "[C"+str(i)+"] (mol/cm3)"
                if i <= perirhizalModel.numDissolvedSoluteComp:
                    extraArray_ /= (perirhizalModel.soilModel.getWaterContent() *cell_volumes) #cm3 water
                else: 
                    extraArray_ /= cell_volumes #cm3


            vp.plot_roots_and_soil(perirhizalModel.mappedSegments(),
                                   extraArrayName_ ,
                                   perirhizalModel.get_concentration(i , konz), s, 
                                   periodic,
                               soilTextureAndShape['min_b'],
                               soilTextureAndShape['max_b'],
                                   soilTextureAndShape['cell_number'], 
                    filename="C"+str(i), 
                                   sol_ind =-1,
                                   extraArray = extraArray_, 
                    extraArrayName = extraArrayName_,
                interactiveImage=False)  # VTK vizualisation
            print("idcomp_done", i, rank)
    
    print("fin", rank)
    # for some reason, we get an error when the program stops
    # memory leak?
    return results_dir

        

if __name__ == '__main__':
    # python3 XcGrowth.py 9 dumux_10c 10 0 customDry noAds 9.02 0.02
    # python3 XcGrowth.py 9 dumux_10c 10 1640 lateDry
    # python3 XcGrowth.py 12 dumux_10c 25 98 baseline
    # python3 XcGrowth.py 10 dumux_10c 25 3 none
    # python3 XcGrowth.py 9 dumux_10c 9.001 5 none
    # python3 XcGrowth.py 9 dumux_10c 9.06 0 customDry xx 9.02 0.02
    # python3 XcGrowth.py 10 dumux_10c 10.06 0 customDry xx 10.02 0.02
    if rank == 0:
        print('sys.argv',sys.argv)
    initsim =float(sys.argv[1])# initsim = 9.5
    mode = sys.argv[2] #"dumux_w" "dumux_3c" "dumux_10c" 
    
    simMax = initsim + 3.
    if len(sys.argv)>3:
        simMax = float(sys.argv[3])
    paramIndx_base = 0
    if len(sys.argv)>4:
        paramIndx_base = int(sys.argv[4])
    extraName ="" # noAds?
    #if len(sys.argv)>6:
    #    extraName = sys.argv[6]
    
    doProfile = ( extraName == "cProfile")
    
    scenario = "baseline"
    if len(sys.argv)>5:
        scenario = sys.argv[5]
    
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
        spellStart = float(sys.argv[7])
        condition = "dry"
        spellDuration =   float(sys.argv[8])
    elif scenario == "customWet":
        spellStart = float(sys.argv[7])
        condition = "wet"
        spellDuration =   float(sys.argv[8])
    else :
        print("scenario", scenario,"not recognised")
        raise Exception
        
    #simInit = 10
    #simEnd = 25
    spellEnd = spellStart + spellDuration
    spellData = {'scenario':scenario,'spellStart':spellStart,'spellEnd':spellEnd, 'condition':condition}
    
    if doProfile:
        import cProfile
        import pstats, io
        pr = cProfile.Profile()
        pr.enable()
    results_dir = XcGrowth(initsim, mode,simMax,extraName,paramIndx_base,spellData )
    if doProfile:
        pr.disable()
        filename = results_dir+'profile'+str(rank)+'.prof' 
        s = io.StringIO()
        ps = pstats.Stats(pr, stream=s).sort_stats('tottime')
        ps.dump_stats(filename)