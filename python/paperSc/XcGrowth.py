""" 
    Maize using rhizosphere models  
"""
import matplotlib; matplotlib.use('agg')
import sys;
import os
absolutePath = False
if absolutePath:
    os.chdir( "/home/rbtlm640/DUMUXr/dumux-rosi/python/paperSc")
    sys.path.append( "/home/rbtlm640/DUMUXr/dumux-rosi/python/paperSc")
    sys.path.append("/home/rbtlm640/DUMUXr/dumux-rosi/python/fpit/data/");
    sys.path.append("/home/rbtlm640/DUMUXr/dumux-rosi/python/modules_fpit/");
    sys.path.append("/home/rbtlm640/DUMUXr/CPlantBox/");#absolut path work better when
    # using pycharmSoil_solute_conc
    sys.path.append("/home/rbtlm640/DUMUXr/CPlantBox/src")
else:
    sys.path.append("../modules_fpit/");
    sys.path.append("../../../CPlantBox/");
    sys.path.append("../../../CPlantBox/src")

from importlib import reload
import plantbox as pb  # CPlantBox
import visualisation.vtk_plot as vp
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import timeit
import numpy as np

import scenario_setup
scenario_setup = reload(scenario_setup)
import rhizo_modelsPlant  # Helper class for cylindrical rhizosphere models
rhizo_modelsPlant = reload(rhizo_modelsPlant)
from rhizo_modelsPlant import *
#import evapotranspiration as evap
#import cyl_exu
import cyl3plant as cyl3
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import os
from scenario_setup import write_file_array, write_file_float, div0, div0f


def getVTPOutput(r,rs,s, datas, datasName, rs_age, results_dir, min_b, max_b, cell_number, initPrint=False):
    vtpindx= int(rs_age)
    if rank == 0:
        print('doing print at', rs_age)
        ana = pb.SegmentAnalyser(r.plant.mappedSegments())
        cutoff = 1e-15 #if value too small, makes paraview crash
        datas_p = []
        for ddindx, dd in enumerate(datas):
            dd_p = np.array(dd)
            dd_p[abs(dd_p) < cutoff] = 0
            datas_p.append(dd_p)
            ana.addData(datasName[ddindx], dd_p)

        #ana.write("results"+directoryN+"plot_"+str(int(simStartSim))+str(condition)+"at"+ str(vtpindx) +".vtp", 
        #          ["organType", "subType",
        #           "CST", "psi_Xyl"]) 

        vp.plot_plant(r.plant,p_name =datasName,# ["xylem pressure (cm)","sucrose concentration (mmol/cm3)"],
                            vals =datas_p,#[  psiXyl_p, C_ST_p], 
                            filename =results_dir+"vtpvti/plantAt"+ str(vtpindx), 
                      range_ = [0,5000])
    if not initPrint:
        periodic = False
        konz =True 

        # check because will need all ranks

        sizeSoilCell = rs.soilModel.getCellVolumes() #cm3
        extraArray_ = rs.soilModel.getWaterContent()
        if not konz:
            extraArrayName_ = "theta (cm3)"
            extraArray_ *= sizeSoilCell
        else:
            extraArrayName_ = "theta (cm3/cm3)"

        rp = rs.get_concentration(0, konz)
        s.results_dir = results_dir

        # adapt to have plot_plants_and soil
        vp.plot_roots_and_soil(rs.mappedSegments(),extraArrayName_,rp, s, periodic, min_b, max_b, cell_number, 
                filename="vtpvti/soil_rx_at"+ str(vtpindx),sol_ind =-1,extraArray = extraArray_, extraArrayName = extraArrayName_,
                interactiveImage=False)  # VTK vizualisation

        for i in range(1, rs.numComp+1):

            extraArray_ = rs.soilModel.getSolution(i) * rs.phaseDensity(i)/1e6

            if not konz:
                extraArrayName_ = "C"+str(i)+" mol"
            else:
                extraArrayName_ = "[C"+str(i)+"] (mol/cm3)"
                if i <= rs.numFluidComp:
                    extraArray_ /= (rs.soilModel.getWaterContent() *sizeSoilCell) #cm3 water
                else: 
                    extraArray_ /= sizeSoilCell #cm3

            vp.plot_roots_and_soil(rs.mappedSegments(),
                                   extraArrayName_ ,
                                   rs.get_concentration(i , konz), s, 
                                   periodic, min_b, max_b, 
                                   cell_number, 
                    filename="vtpvti/C"+str(i)+'_at'+str(vtpindx), 
                                   sol_ind =-1,
                                   extraArray = extraArray_, 
                    extraArrayName = extraArrayName_,
                interactiveImage=False)  # VTK vizualisation
        

    

"""
     * \brief Suggest a new number of time steps
     *
     * The default behavior is to suggest the old time-step number
     * scaled by the ratio between the target iterations and the
     * iterations required to actually solve the last time-step.
        // be aggressive increasing the time-step number but
        // conservative when decreasing it. the rationale is
        // that we want to avoid failing in the next iteration
        // nNew > nOld ==> dtnew < dtOld
"""
def suggestNumStepsChange(nOld, numIter_, targetIter_, results_dir):# taken from dumux
    if (numIter_ > targetIter_) :
        percent = float(numIter_ - targetIter_)/float(targetIter_)
        change = 1.0 + percent
    else:
        percent = float(targetIter_ - numIter_)/float(targetIter_)
        change = 1 /(1.0 + percent/1.2)
    #if not doMinimumPrint:
    write_file_array("suggestNumStepsChange",np.array([nOld, numIter_, targetIter_, percent,change, np.ceil(nOld * change)]), directory_ =results_dir, fileType = '.csv') 
    return int(np.ceil(nOld * change))
    

def XcGrowth(initsim, mode,simMax,extraName,paramIndx_,spellData):

    #initsim =float(sys.argv[1])# initsim = 9.5
    #mode = sys.argv[2] #"dumux_w" "dumux_3c" "dumux_10c" 
    dt = 20/60/24
    #p_mean = -100
    k_iter = 6#20
    targetIter= 5
    l_ks =  "dx_2"#"root", "dx", "dx_2"
    organism = "plant"# "RS"#
    weightBefore = False
    SRIBefore = False
    beforeAtNight = True
    adaptRSI_  = False
    static_plant = False
    useOuterFluxCyl_w = False
    useOuterFluxCyl_sol = False
    css1Function_ = 9
    lightType =""#+- "nolight" # or ""
    mpiVerbose = False
    mpiVerboseInner = False
    noAds = False
    doSimple =False
    doMinimumPrint =  True
    
    doOldCell = False
    if doSimple:
        max_b =np.array( [5, 5, 0.] )# 
        min_b = np.array([-5, -5, -5.])# 
        cell_number = np.array([3,3,3])#
    elif doOldCell:
        max_b =np.array( [5, 5, 0.]) # 
        min_b = np.array([-5, -5, -10.])# 
        cell_number =np.array([5,5,20])#
    else:
        min_b = np.array([-3./2, -12./2, -40.])#41
        cell_number =np.array( [3,12,41]) # 1cm3 np.array( [3,6,40]) #
        max_b =np.array( [3./2, 12./2, 0.])
       
    #+str(int(useOuterFluxCyl_w))+str(int(useOuterFluxCyl_sol)) \
    #+lightType+l_ks+str(int(static_plant))+str(int(weightBefore))\
    #+str(int(SRIBefore))+str(int(beforeAtNight))+str(int(adaptRSI_))\
    #+organism+str(k_iter)+"k_"+str(css1Function_)
    #+str(int(mpiVerbose))+l_ks+mode
    # 1d1dFHung
    results_dir="./results/newtrans7/"+extraName+str(spellData['scenario'])\
    +"_"+str(int(np.prod(cell_number)))\
                    +"_"+str(paramIndx_)\
                    +"_"+str(int(initsim))+"to"+str(int(simMax))\
                    +"_"+str(int(dt*24*60))+"mn_"\
                    +str(int((dt*24*60 - int(dt*24*60))*60))+"s_"\
                    +str(max_rank)+"/"
    
    #comm.barrier()
    if rank == 0:
        print('results_dir','DumuxDune27',results_dir, flush = True)
    #comm.barrier()
    if rank == 0:
        for extraText in ["","cyl_val/","printData/", "vtpvti/"]:
            if not os.path.exists(results_dir+extraText):
                os.makedirs(results_dir+extraText)
            else:
                test = os.listdir(results_dir+extraText)
                for item in test:
                    try:
                        os.remove(results_dir+extraText+item)
                    except:
                        pass
    #comm.barrier()
    """
     Cylindric rhizosphere models, C exudate from finite volumes
    """

    fname = 'L_WT_phenolics'

    """scenario"""
    year = 2013
    soil_type = "loam"
    genotype = "WT"
    comp = "phenolics"
    usemoles = True
    """ parameters   """
    soil_ = scenario_setup.vg_SPP(0)
    
    #min_b = [-5., -5, -5.] 
    #max_b = [5., 5, 0.] 
    #cell_number = [5, 5, 5]
    sim_time = 1 #154   #  [day]


    """ rhizosphere model parameters """
    recreateComsol = False

    periodic = True
    nc = 10

    logbase = 0.5  # according to Mai et al. (2019)

    """ initialize """

    weatherInit = scenario_setup.weather(1.,dt, spellData)
    s, soil = scenario_setup.create_soil_model(soil_type, year, soil_,#comp, 
                min_b, max_b, cell_number, demoType = mode, times = None, net_inf = None,
                usemoles = usemoles, dirResults = results_dir, p_mean_ = weatherInit['p_mean'], 
                                         css1Function = css1Function_,
                                        paramIndx=paramIndx_,
                                        noAds = noAds)

    if organism == "plant":
        path = "../../../CPlantBox/modelparameter/structural/plant/"
        xml_name = "Triticum_aestivum_test_2021.xml"  # root growth model parameter file
    elif organism == "RS":
        path = "./"
        xml_name = "data_magda/Zeamays_synMRI_modified.xml"
    elif organism == "RSstatic":
        path = "./"
        xml_name = "../../grids/RootSystem8.rsml"
    else:
        raise Exception
    

    write_file_array("cellVol", np.array(s.getCellVolumes()), directory_ =results_dir) # cm3 
    write_file_array("cellIds", np.array(s.cellIndices), directory_ =results_dir) # cm3
    cellcenters = s.getCellCenters()
    cellcentersX = []
    cellcentersY = []
    cellcentersZ = []
    for sub_array in cellcenters:
        cellcentersX.append(sub_array[0])
        cellcentersY.append(sub_array[1])
        cellcentersZ.append(sub_array[2])

    write_file_array("cellcentersX", np.array(cellcentersX), directory_ =results_dir) # cm3
    write_file_array("cellcentersY", np.array(cellcentersY), directory_ =results_dir) # cm3
    write_file_array("cellcentersZ", np.array(cellcentersZ), directory_ =results_dir) # cm3
    

    # all thread need a plant object, but only thread 0 will make it grow
    rs, r = scenario_setup.create_mapped_plant( nc, logbase, mode,initsim,
                                            min_b, max_b, cell_number, s, xml_name,
                                            path, plantType = organism, 
                                            recreateComsol_ = recreateComsol,
                                            usemoles = usemoles,l_ks_ = l_ks,
                                            limErr1d3d = 5e-13, spellData = spellData)  # pass parameter file for dynamic growth

    rs.weightBefore = weightBefore
    rs.SRIBefore = SRIBefore
    rs.beforeAtNight = beforeAtNight
    rs_age = initsim
    rs.useOuterFluxCyl_w = useOuterFluxCyl_w
    rs.useOuterFluxCyl_sol = useOuterFluxCyl_sol
    s.mpiVerbose = mpiVerbose
    s.mpiVerboseInner = mpiVerboseInner
    rs.mpiVerbose = mpiVerbose
    rs.mpiVerboseInner = mpiVerboseInner
    rs.spellData = spellData
    rs.enteredSpell = (rs.spellData['scenario'] == 'none') or (rs.spellData['scenario'] == 'baseline')
    rs.leftSpell = (rs.spellData['scenario'] == 'baseline')
    rs.diff1d3dCurrant = np.Inf
    rs.diff1d3dCurrant_rel = np.Inf
    rs.maxdiff1d3dCurrant = np.Inf
    rs.maxdiff1d3dCurrant_rel = np.Inf


    net_sol_flux =  np.array([np.array([]),np.array([])])
    net_flux = np.array([])


    Q_ST_init = np.array([])
    Nt = len(rs.nodes)
    Q_Exud    = np.zeros(Nt)
    Q_Mucil   = np.zeros(Nt)
    Q_Exudbu    = np.zeros(Nt)
    Q_Mucilbu   = np.zeros(Nt)
    Q_Rm    = np.zeros(Nt)
    Q_Gr   = np.zeros(Nt)
    Q_Rmbu    = np.zeros(Nt)
    Q_Grbu   = np.zeros(Nt)
    Q_in  = 0
    Nt = len(rs.nodes)
    r.minLoop = 1000
    r.maxLoop = 5000
    #simMax = initsim + 3
    failedExud = False

    r.AgCumul = 0
    r.TranspirationCumul = 0
    r.TranspirationCumul_eval = 0
    cell_volumes = s.getCellVolumes()  # cm3
    buTotCSoilInit = sum(s.getTotCContent())
    s.buWSoilInit = sum(np.multiply(np.array(s.getWaterContent()), cell_volumes))

    r.time_start_global = timeit.default_timer()
    r.time_rhizo_cumul = 0
    r.time_3ds_cumul = 0
    r.time_plant_cumulW = 0
    r.time_plant_cumulS = 0

    r.time_rhizo_i = 0
    r.time_3ds_i = 0
    r.time_plant_cumul = 0
    
    
    
    Q_Exud_inflate = 0.; Q_Mucil_inflate = 0.
    rs.results_dir = results_dir
    
    if mpiVerbose and rank == 0:# or (max_rank == 1):
        print('start loop', rank)
    secId = None
    Q_Gr_i = None;Q_Rm_i = None
    Q_Exud_i = None; Q_Mucil_i = None
    Q_Exud_i_seg = np.array([]); Q_Mucil_i_seg = np.array([])
    error_st_abs = 0;error_st_rel=0
    errs = np.array(["errRxPlant", "errW1ds", "errW3ds","errC1ds", "errC3ds",
                     "max(r.SinkLim3DS)","max(r.SinkLim1DS)","max(abs(r.OutLim1DS))",
                     "max(abs(r.InOutBC_Cdiff))",
                     "max(r.maxDiff1d3dCW_abs)",
                     "errWrsi",# "maxDiff1d3dCW_absBU",
                     "bulkMassErrorWater_abs","bulkMassErrorWater_absLim",
                     "rhizoMassWError_absLim","rhizoMassWError_abs",
                     "bulkMassErrorC_abs","bulkMassCErrorPlant",
                     "rhizoMassCError_absLim","rhizoMassCError_abs",
                     "sum(abs(diffBCS1dsFluxIn))", "sum(abs(diffBCS1dsFluxOut))",
                     "sum(abs(diffouter_R_bc_wat))",
                     "sum(abs(diffBCS1dsFluxOut_sol))",
                     "sum(abs(diffBCS1dsFluxOut_mucil))","sum(abs(diffouter_R_bc_sol))",
                     "diff1d3dCurrant","diff1d3dCurrant_rel","rhizoMassWError_rel",'err'])
    write_file_array("OuterSuccess_error", errs, directory_ =results_dir, fileType = '.csv')
    write_file_array("inner_error", errs, directory_ =results_dir, fileType = '.csv')
    
    if not doMinimumPrint:
        write_file_array("N_error", errs, directory_ =results_dir, fileType = '.csv')
        write_file_array("fpit_error", errs, directory_ =results_dir, fileType = '.csv') 
    seg_fluxes_ = np.array([])
    dt_inner = float(dt)#/float(2.)
    start = True
    rs.new_soil_solute = np.array([0.])
    
    getVTPOutput(r,rs,s,[],[], rs_age*100, results_dir, min_b, max_b, cell_number,initPrint=True)#[np.array(rs.organTypes)], ["organTypes"],
    while rs_age < simMax:

        rs_age += dt
        #comm.barrier()
        if mpiVerbose and rank == 0:# or (max_rank == 1):
            print(rank,"Day", rs_age)
        #comm.barrier()
        seg2cell_old = rs.seg2cell
        Ntbu = Nt
        #comm.barrier()
        if mpiVerbose and rank == 0:# or (max_rank == 1):
            print(rank,'simulating')
        #comm.barrier()
        if (rank == 0) and (not static_plant) :

            rs.simulate(dt, verbose = False)#(rank == 0))  # because of dumux.pick(), blocks if I do not let all threads do the simulate.
            seg2cell_new = rs.seg2cell

            write_file_array('seg2cell_keys',seg2cell_new,directory_ =results_dir, 
                             fileType = '.csv')
            write_file_array('seg2cell_vals',np.array(list(seg2cell_new.values())),
                             directory_ =results_dir, fileType = '.csv')

            # make sure that, once a segment is created, it stays in the same soil voxel
            for segid in seg2cell_old.keys():
                assert seg2cell_old[segid] == seg2cell_new[segid]
            # also segs are in only one voxel
            cell2segVals = np.concatenate((list(rs.cell2seg.values()))).flatten()
            if len(cell2segVals) != len(set(cell2segVals)):#make sure all values are unique
                print(seg2cell_new)
                print(rs.cell2seg)
                print(cell2segVals)
                print(len(cell2segVals), len(set(cell2segVals)))
                raise Exception

                
        if (rank == 0):# and (not doMinimumPrint):
            write_file_array("organTypes", np.array(rs.organTypes), directory_ =results_dir)
            write_file_array("subTypes", np.array(r.rs.subTypes), directory_ =results_dir)
            write_file_array("isRootTip", r.rs.isRootTip, directory_ =results_dir)
            write_file_array("getNodeIds", r.plant.getNodeIds(-1), directory_ =results_dir)
            length_Segs = np.array(r.rs.segLength())
            write_file_array("length_Segs", length_Segs, directory_ =results_dir)
            orgs = r.plant.getOrgans(-1, False)
            parentorgid = np.array([org.getParent().getId() for org in orgs])
            parentNidx = np.array([org.getParent().getNodeId(org.parentNI) for org in orgs])
            id_orgs = np.array([org.getId() for org in orgs])
            ot_orgs = np.array([org.organType() for org in orgs])
            st_orgs = np.array([org.getParameter("subType") for org in orgs])
            lenOrg = np.array([org.getLength(False) for org in orgs])   
            write_file_array("id_orgs", id_orgs, directory_ =results_dir)
            write_file_array("ot_orgs", ot_orgs, directory_ =results_dir)
            write_file_array("st_orgs", st_orgs, directory_ =results_dir)  
            write_file_array("lenOrg", lenOrg, directory_ =results_dir) 
            write_file_array("parentNidx", parentNidx, directory_ =results_dir) 
            
            write_file_array("nodes_X",
                             np.array([tempnode[0] for tempnode in r.get_nodes()]), 
                             directory_ =results_dir)
            write_file_array("nodes_Y", np.array([tempnode[1] for tempnode in r.get_nodes()]), directory_ =results_dir)
            write_file_array("nodes_Z", np.array([tempnode[2] for tempnode in r.get_nodes()]), directory_ =results_dir)
            idPerNode = np.concatenate([
                np.full(org.getNumberOfNodes()-1,org.getId()) for org in orgs]).reshape(-1)
            globalNodeId = np.concatenate([org.getNodeIds()[1:] for org in orgs]).reshape(-1)
            write_file_array("orgidPerNode", idPerNode, directory_ =results_dir)
            write_file_array("globalNodeId", globalNodeId, directory_ =results_dir)
            volOrg = np.array([org.orgVolume(-1,False) for org in orgs]) 
            write_file_array("volOrg", volOrg, directory_ =results_dir)

        #comm.barrier()
        if mpiVerbose and (rank == 0):
            print(rank,'simulated')
        #comm.barrier()
        rs.update()
        #comm.barrier()
        if mpiVerbose and (rank == 0):
            print(rank,'updated')
        #comm.barrier()
        if start:
            nots = len(np.array(rs.organTypes))
            rs.checkMassOMoleBalance2(None,None, dt,
                                    seg_fluxes =None, diff1d3dCW_abs_lim = 1e-13, takeFlux = False) # very beginning: should not be any error
            write_file_array("time", np.array([rs_age,0.]), directory_ =results_dir)
            write_file_array("sumErrors1ds3ds", np.concatenate((rs.sumDiff1d3dCW_abs, rs.sumDiff1d3dCW_rel)), directory_ =results_dir, fileType = '.csv')
            write_file_array("maxErrors1ds3ds", np.concatenate((rs.maxDiff1d3dCW_abs, rs.maxDiff1d3dCW_rel)), directory_ =results_dir, fileType = '.csv')             
            start = False
            

            
        write_file_array("content_diff1d3dabs"+str(0), 
                         rs.allDiff1d3dCW_abs[0], directory_ =results_dir, fileType = '.csv')
        write_file_array("content_diff1d3drel"+str(0), 
                         rs.allDiff1d3dCW_rel[0], directory_ =results_dir, fileType = '.csv')
        for nc in range(rs.numComp):# normally all 0 for nc >= numFluidComp
                write_file_array("content_diff1d3dabs"+str(nc+1), 
                                 rs.allDiff1d3dCW_abs[nc+1], directory_ =results_dir, fileType = '.csv')
                write_file_array("content_diff1d3drel"+str(nc+1), 
                                 rs.allDiff1d3dCW_rel[nc+1], directory_ =results_dir, fileType = '.csv')

                write_file_array("sol_content3d"+str(nc+1), 
                                 s.getContent(nc+1, nc < s.numFluidComp), 
                                 directory_ =results_dir, fileType = '.csv')  
                write_file_array("sol_content1d"+str(nc+1), 
                                 rs.getContentCyl(idComp = nc+1, doSum = False, reOrder = True), 
                                 directory_ =results_dir, fileType = '.csv') 
            

            
        
                           

        #comm.barrier()


        Q_Exud_inflate += sum(Q_Exud_i_seg); Q_Mucil_inflate += sum(Q_Mucil_i_seg)
        if False:
            try:
                assert abs(Q_Exud_inflate- sum(Q_Exud)) < 1e-16
                assert abs(Q_Mucil_inflate -sum(Q_Mucil)) <1e-16
            except:
                print('issue Q_Exud_inflate',rank, 
                      Q_Exud_inflate, sum(Q_Exud),abs(Q_Exud_inflate- sum(Q_Exud)), 
                      Q_Mucil_inflate, sum(Q_Mucil),abs(Q_Mucil_inflate- sum(Q_Exud)) )
                raise Exception

        #comm.barrier()
        if mpiVerbose and rank==0:# or (max_rank == 1):
            print(rank, 'to cyl3.simulate_const')
        #comm.barrier()
        n_iter = 0
        rs.rhizoMassWError_abs = 1.
        rs.rhizoMassCError_abs = 1.
        rs.errDiffBCs = 1.
        rs.err = 1.
        max_err = 1.
        rs.diff1d3dCurrant_rel = 1e6
        rs.maxdiff1d3dCurrant_rel = 1e6
        def continueLoop(rs,n_iter, dt_inner: float,failedLoop: bool,
                         real_dtinner: float,name="continueLoop", isInner = False,doPrint = True, 
                         fileType = '.csv' , plant = None):
            
            plant.time_rhizo_cumul += plant.time_rhizo_i
            plant.time_3ds_cumul += plant.time_3ds_i
            plant.time_rhizo_i = 0
            plant.time_3ds_i = 0
            plant.time_plant_cumul = plant.time_plant_cumulW + plant.time_plant_cumulS 
            
            write_file_array("totalComputetime",
                             np.array([timeit.default_timer() - plant.time_start_global,
                            plant.time_plant_cumul,plant.time_rhizo_cumul ,plant.time_3ds_cumul]) , 
                             directory_ =results_dir)
            # sumDiff1d3dCW_rel = rs.sumDiff1d3dCW_rel[:(rs.numFluidComp+1)]
            # sumDiff1d3dCW_rel = np.where(np.isnan(sumDiff1d3dCW_rel),0.,sumDiff1d3dCW_rel)
            #  or (abs(rs.rhizoMassWError_abs) > 1e-13) or (abs(rs.rhizoMassCError_abs) > 1e-9) or (max(abs(rs.errDiffBCs*0)) > 1.)
            n_iter_min = 4
            cL = ((np.floor(rs.err) > max_err) or  rs.solve_gave_up or failedLoop
                    or (np.floor(rs.diff1d3dCurrant_rel*1000.)/1000.>0.001) 
                    or (np.floor(rs.maxdiff1d3dCurrant_rel*100.)/100.>0.001) 
                    or (min(rs.new_soil_solute.reshape(-1)) < 0)  
                    or ((n_iter < n_iter_min) and (isInner)))  and (n_iter < k_iter)
            #(max(abs(sumDiff1d3dCW_rel))>1)) 
            #rs.diff1d3dCurrant_rel

            #comm.barrier()
            if rank == 0:#mpiVerbose:# or (max_rank == 1):
                print('continue loop?',rank,'n_iter',n_iter,'cL',cL,'failedLoop',failedLoop, ' np.floor(rs.err)',
                np.floor(rs.err),  'solve_gave_up',rs.solve_gave_up,
                'diff1d3dCurrant_rel',rs.diff1d3dCurrant_rel, 
                'maxdiff1d3dCurrant_rel',rs.maxdiff1d3dCurrant_rel, 
                'k_iter',k_iter)
            #comm.barrier()
            cL = comm.bcast(cL,root = 0)
            failedLoop_ = np.array( comm.bcast(comm.gather(failedLoop,root = 0),root = 0))
            #comm.barrier()
            #if mpiVerbose and (max_rank > 1):
            #    print('continue loopBis?',rank,cL,failedLoop_)
            #comm.barrier()
            assert (failedLoop_ ==failedLoop_[0]).all() # all true or all false
            
            if doPrint:
                if not os.path.isfile(results_dir+name+fileType):
                    write_file_array(name, np.array(['n_iter', 'err', 
                                                     #'rhizoMassWError_abs','rhizoMassCError_abs','maxAbsErrDiffBCs',
                                                     'diff1d3dCurrant_rel','maxdiff1d3dCurrant_rel',
                                                     'solve_gave_up', 'min__soil_solute',
                                                             'dt_inner','real_dtinner','failedLoop','cL']), directory_ =results_dir, fileType = fileType)
                    if not doMinimumPrint:
                        write_file_array(name+"Bool", np.array(['n_iter',  'err', 
                                                                #'rhizoMassWError_abs','rhizoMassCError_abs','maxAbsErrDiffBCs',
                                                                'diff1d3dCurrant_rel','maxdiff1d3dCurrant_rel',
                                                     'solve_gave_up',  'min__soil_solute',
                                                                 'dt_inner','failedLoop','cL']), directory_ =results_dir, fileType = fileType)
                    
                write_file_array(name, np.array([n_iter, rs.err, 
                                                 #rs.rhizoMassWError_abs,rs.rhizoMassCError_abs,max(abs(rs.errDiffBCs)),
                                                 rs.diff1d3dCurrant_rel,rs.maxdiff1d3dCurrant_rel,
                                                 rs.solve_gave_up, min(rs.new_soil_solute.reshape(-1)),
                                                             dt_inner, real_dtinner,failedLoop,cL]), directory_ =results_dir, fileType = fileType)
                if not doMinimumPrint:
                    write_file_array(name+"2", 
                                     np.concatenate((rs.sumDiff1d3dCW_rel,rs.maxDiff1d3dCW_rel)),  
                                     directory_ =results_dir, fileType = fileType)
                    write_file_array(name+"Bool", np.array([n_iter, (np.floor(rs.err) > max_err), 
                                                            (np.floor(rs.diff1d3dCurrant_rel*10.)/10. > 0.1),
                                                            (np.floor(rs.maxdiff1d3dCurrant_rel) >1),
                                                            #(abs(rs.rhizoMassWError_abs) > 1e-13), (abs(rs.rhizoMassCError_abs) > 1e-9), 
                                                            #(max(abs(rs.errDiffBCs*0)) > 1e-5), 
                                                            rs.solve_gave_up, min(rs.new_soil_solute.reshape(-1)) < 0.,
                                                                 dt_inner,failedLoop,cL]), directory_ =results_dir, fileType = fileType)
            
            return cL
        failedLoop = False
        cL = True
        while cL or failedLoop:
            r.AgCumul_inner = 0
            r.TranspirationCumul_inner = 0 # reset transpiration of inner time step to 0
            s.base.saveManual()# <= that actually only works for the last solve. I need a better method to reset to the beginning of the whole stuff
            rs.saveManual()
            rs.leftSpellBU = rs.leftSpell
            rs.enteredSpellBU = rs.enteredSpell
            net_sol_flux_, net_flux_, seg_fluxes__, real_dtinner,failedLoop, n_iter_inner_max = cyl3.simulate_const(s, 
                                                    r, sim_time= dt,dt= dt_inner, rs_age=rs_age, 
                                                    Q_plant=[Q_Exud_i_seg, Q_Mucil_i_seg],
                                                    r= rs,plantType = organism,
                                                     adaptRSI=adaptRSI_,
                                                    outer_R_bc_sol = net_sol_flux, 
                                                    outer_R_bc_wat = net_flux,seg_fluxes=seg_fluxes_,
                                                    results_dir = results_dir,
                                                        k_iter_ = k_iter,lightType_=lightType, outer_n_iter = n_iter,
                                                                   continueLoop= continueLoop, doMinimumPrint = doMinimumPrint)
            cL = continueLoop(rs,n_iter, dt_inner,failedLoop,real_dtinner,name="Outer_data",isInner = False,
                             plant = r)
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
                n_iter_inner_max = k_iter
            nNew = suggestNumStepsChange(nOld, n_iter_inner_max, targetIter #np.ceil(k_iter/2)
                                         , results_dir)
            try:
                assert nOld == int(nOld)
                assert nNew == int(nNew)
            except:
                print('nNew iisue',nNew , int(nNew), nOld,dt,dt_inner,(nOld == int(nOld)), (nNew == int(nNew)))
                
                write_file_array("nNew_error", np.array([nNew , int(nNew), nOld,dt,dt_inner,(nOld == int(nOld)), (nNew == int(nNew)),
                                                         real_dtinner ,dt, dt_inner , failedLoop,abs((real_dtinner - dt)/dt*100.),rs_age]), 
                                 directory_ =results_dir, fileType = '.csv') 
                
            dt_inner = dt/float(nNew)
            if cL or failedLoop:
                #currentN = int(np.ceil(dt / dt_inner))
                #comm.barrier()
                
                if mpiVerbose and rank==0:# or (max_rank == 1):
                    print(rank, "error too high, decrease N, dt from", nOld, dt/float(nOld),"to",nNew, dt_inner)
                #dt_inner = dt/(float(currentN*2.))
                #dt_inner = max(1./(24.*3600.), dt_inner) # minimum: 1 second
                # print(rank, "error too high, decrease N from", dt/float(currentN),"to",min(1/(24*3600), dt_inner))
                #comm.barrier()
                r.AgCumul_inner = 0
                r.TranspirationCumul_inner = 0 # reset transpiration of inner time step to 0
                s.base.resetManual()# <= that actually only works for the last solve. I need a better method to reset to the beginning of the whole stuff
                rs.resetManual()
                rs.leftSpell = rs.leftSpellBU
                rs.enteredSpell = rs.enteredSpellBU
                
                rs.checkMassOMoleBalance2(None,None, dt,
                                    seg_fluxes =None, diff1d3dCW_abs_lim = np.Inf, takeFlux = False)
                # to get correct error values for sumDiff1d3dCW_relOld
                
        cL_ = continueLoop(rs,0, dt_inner,failedLoop,real_dtinner,name="TestOuter_data", plant = r)
        if cL_:
            raise Exception#only left the loop because reached the max number of iteration.
        
        write_file_array("OuterSuccess_error", rs.errs, directory_ =results_dir, fileType = '.csv') 
        write_file_array("OuterSuccess_data", np.array([n_iter, rs.err, rs.rhizoMassWError_abs,dt_inner]), directory_ =results_dir, fileType = '.csv')
        net_sol_flux = net_sol_flux_
        net_flux = net_flux_ 
        seg_fluxes_ = seg_fluxes__
                
                
        #comm.barrier()
        if mpiVerbose and rank==0:# or (max_rank == 1):
            print(rank, 'left cyl3.simulate_const')
        #comm.barrier()

        
        if mpiVerbose and rank==0:# or (max_rank == 1):
            print(rank,"theta")
        #comm.barrier()#
        write_file_array("pHead", np.array(s.getSolutionHead()), directory_ =results_dir, fileType = '.csv') 
        write_file_array("theta", np.array(s.getWaterContent()), directory_ =results_dir, fileType = '.csv') 
        write_file_array("watVol", np.array(s.getWaterVolumes()), directory_ =results_dir, fileType = '.csv') 
        for i in range(rs.numFluidComp):
            #comm.barrier()
            if mpiVerbose and rank==0:# or (max_rank == 1):
                print(rank,"Soil_solute_conc"+str(i+1))
            #comm.barrier()
            write_file_array("Soil_solute_conc"+str(i+1), 
                             np.array(s.getSolution(i+1)).flatten()* rs.molarDensityWat_m3/1e6, # mol C/mol w * molW/m3 W * m3 W/cm3 W
                             directory_ =results_dir, fileType = '.csv') 
        for i in range(rs.numFluidComp, rs.numComp):
            #comm.barrier()
            if mpiVerbose and rank==0:# or (max_rank == 1):
                print(rank,"Soil_solute_conc"+str(i+1))
            #comm.barrier()
            write_file_array("Soil_solute_conc"+str(i+1), 
                             np.array(s.getSolution(i+1)).flatten()* rs.bulkDensity_m3 /1e6 , 
                             directory_ =results_dir, fileType = '.csv') 
        #comm.barrier()
        if mpiVerbose and rank==0:# or (max_rank == 1):
            print(rank,"Soil_solute_conc"+str(rs.numComp+1))
        #comm.barrier()
        write_file_array("Soil_solute_conc"+str(rs.numComp+1), np.array(s.base.getCSS1_out()).flatten(),#[ mol / cm^3]
                         # [:-1]* rs.bulkDensity_m3 /1e6 ,
                         directory_ =results_dir, fileType = '.csv') 
           
        #comm.barrier() 
        if mpiVerbose and rank==0:# or (max_rank == 1):
            print(rank, 'did s.data writing')
        #comm.barrier()
        errLeuning_abs = abs(sum(r.outputFlux))
        if organism == "plant":
            r.TranspirationCumul += r.TranspirationCumul_inner #sum(np.array(r.Ev) * dt) #transpiration [cm3/day] * day
            if rs.enteredSpell and (not rs.leftSpell):
                r.TranspirationCumul_eval += r.TranspirationCumul_inner
            elif rs.leftSpell:
                r.TranspirationCumul_eval = 0.
               
            r.An =comm.bcast( r.AgCumul_inner/(dt*24*3600) , root=0)
            
        else:
            r.TranspirationCumul += sum(r.outputFlux)
        

        #comm.barrier()
        if mpiVerbose and rank==0:# or (max_rank == 1):
            print(rank, 'getTotCContent')
        buTotCAfter = sum(s.getTotCContent())   #0 get stuck here
        #comm.barrier()
        if mpiVerbose and rank==0:# or (max_rank == 1):
            print(rank, 'getWaterContent')
        #comm.barrier()
        buWAfter = sum(np.multiply(np.array(s.getWaterContent()), cell_volumes))    

        #comm.barrier()
        if mpiVerbose and rank==0:# or (max_rank == 1):
            print(rank, 'get s.errorCumul')
        #comm.barrier()
        if rank == 0:
            if (mode != "dumux_w"):
                s.bulkMassErrorCumul_abs = abs((buTotCAfter - ( buTotCSoilInit + Q_Exud_inflate + Q_Mucil_inflate)))#so that only works if infalte
                if buTotCAfter != 0:
                    s.bulkMassErrorCumul_rel = abs(s.bulkMassErrorCumul_abs/buTotCAfter*100)
                else:
                    s.bulkMassErrorCumul_rel =np.nan
            else:
                s.bulkMassErrorCumul_abs = np.nan
                s.bulkMassErrorCumul_rel =np.nan
            # because we have a free flow BC at the bottom, that could increase the error
            # ideally, I should add the flow at the bellow BC here
            s.bulkMassErrorWaterCumul_abs = abs(buWAfter - ( s.buWSoilInit - r.TranspirationCumul_eval))
            s.bulkMassErrorWaterCumul_rel = abs(s.bulkMassErrorWaterCumul_abs/buWAfter*100)
        else:
            s.bulkMassErrorCumul_abs = None
            s.bulkMassErrorCumul_rel = None
            s.bulkMassErrorWaterCumul_abs = None
            s.bulkMassErrorWaterCumul_rel = None
            
        #comm.barrier()
        if mpiVerbose and rank==0:# or (max_rank == 1):
            print(rank, 'got s.errorCumul')
        
        #comm.barrier()
        write_file_array("time", np.array([rs_age,r.Qlight]), directory_ =results_dir)
        write_file_array("sumErrors1ds3ds", np.concatenate((rs.sumDiff1d3dCW_abs, rs.sumDiff1d3dCW_rel)), directory_ =results_dir, fileType = '.csv')
        write_file_array("maxErrors1ds3ds", np.concatenate((rs.maxDiff1d3dCW_abs, rs.maxDiff1d3dCW_rel)), directory_ =results_dir, fileType = '.csv')# 
        #comm.barrier()
        if mpiVerbose and rank==0:# or (max_rank == 1):
            print(rank, 'write some otehr stuff')
        #comm.barrier()
        
        if (mode != "dumux_w"):
            #absolute and relative (%) error
            write_file_array("errorsPlant", np.array([error_st_abs,error_st_rel,#cumulative
                                                errLeuning_abs]), directory_ =results_dir, fileType = '.csv') #not cumulative
            write_file_array("errorsBulkSoil", np.array([s.bulkMassCErrorPlant_abs, s.bulkMassCErrorPlant_rel, #not cumulative 
                                                s.bulkMassCError1ds_abs, s.bulkMassCError1ds_rel, 
                                                s.bulkMassErrorCumul_abs,s.bulkMassErrorCumul_rel,#cumulative
                                                s.bulkMassErrorWater_abs,s.bulkMassErrorWater_rel, #not cumulative
                                                s.bulkMassErrorWaterCumul_abs,s.bulkMassErrorWaterCumul_rel]),
                             directory_ =results_dir, fileType = '.csv')#cumulative
            write_file_array("errorMassRhizo", np.array([rs.rhizoMassCError_abs, rs.rhizoMassCError_rel,
                                                        rs.rhizoMassWError_abs, rs.rhizoMassWError_rel]), directory_ =results_dir)# not cumulativecumulative (?)

            write_file_float("trans", r.TranspirationCumul, directory_ =results_dir)
            write_file_array("transrate",r.Jw, directory_ =results_dir, fileType = '.csv')
        #comm.barrier()
        if mpiVerbose and rank==0:# or (max_rank == 1):
            print(rank, 'finished other data writing')
        #comm.barrier()
        if False:
            try:
                assert abs(s.bulkMassCErrorPlant_abs)  < 1e-5
            except:
                print("\n\n\nissue bulk soil balance", np.array([s.bulkMassCErrorPlant_abs, s.bulkMassCErrorPlant_rel, #not cumulative 
                                                s.bulkMassErrorCumul_abs,s.bulkMassErrorCumul_rel,#cumulative
                                                s.bulkMassErrorWater_abs,s.bulkMassErrorWater_rel, #not cumulative
                                                s.bulkMassErrorWaterCumul_abs,s.bulkMassErrorWaterCumul_rel]))
                print("\n\n\n")
                raise Exception
                
        for lId, cyl in enumerate(rs.cyls):
            if not isinstance(cyl, AirSegment):
                gId = rs.eidx[lId]
                
                pHead = np.array(cyl.getSolutionHead()).flatten()
                write_file_float("Cyl_time_"+str(gId),rs_age, 
                                 directory_ =results_dir+'cyl_val/', allranks = True)
                write_file_array("Cyl_watercontent_"+str(gId),cyl.getWaterContent(), 
                                 directory_ =results_dir+'cyl_val/', allranks = True)
                write_file_array("Cyl_pressureHead_"+str(gId),pHead, 
                                 directory_ =results_dir+'cyl_val/', allranks = True)
                write_file_array("Cyl_coord_"+str(gId),
                                 cyl.getDofCoordinates().flatten(), 
                                 directory_ =results_dir+'cyl_val/', allranks = True)
                write_file_array("Cyl_cellVol_"+str(gId),
                                 cyl.getCellVolumes().flatten(), 
                                 directory_ =results_dir+'cyl_val/', allranks = True)
                for ccc in range(rs.numComp):
                    sol0 = np.array(cyl.getContent(ccc +1, ccc < 2 )).flatten()
                    write_file_array("Cyl_content"+str(ccc+1)+"_"+str(gId)+"", 
                                 sol0, 
                                 directory_ =results_dir+'cyl_val/', allranks = True)
                if max(pHead) > 0:
                    print('issue phead',gId,rank, pHead, sol0 )
                    raise Exception
        #comm.barrier()
        if mpiVerbose and rank==0:# or (max_rank == 1):
            print(rank, 'do C growth?',  (mode != "dumux_w") and (rank == 0) and ((not static_plant) or (rs_age == initsim+dt)) and (organism == "plant"))

            
        r.Csoil_seg = rs.get_inner_solutes() * 1e3 # mol/cm3 to mmol/cm3 
        
        #comm.barrier()
        if (mode != "dumux_w") and (rank == 0) and ((not static_plant) or (rs_age == initsim+dt)) and (organism == "plant"):
            assert min(r.Csoil_seg ) >= 0.
            startphloem=rs_age
            endphloem = rs_age + dt
            stepphloem = 1
            verbose_phloem = True
            filename = results_dir +"inPM"+str(rank)+".txt"
            
            print("startpm",rank)

            r.time_start_plant = timeit.default_timer()

            r.startPM(startphloem, endphloem, stepphloem, ( rs.weatherX["TairC"]  +273.15) , verbose_phloem, filename)

            r.time_plant_cumulS += (timeit.default_timer() - r.time_start_plant)


            print("endtpm",rank)
            
            Nt = len(rs.nodes)
            
            ##
            #  MMOL Suc (/cm3) => mol C (/cm3)
            mmolSuc_to_molC = 1/1e3*12
            
            if r.withInitVal and (len(Q_ST_init) ==0) :
                Q_ST_init = np.array(r.Q_init[0:Nt])* mmolSuc_to_molC
                Q_meso_init = np.array(r.Q_init[Nt:(Nt*2)])* mmolSuc_to_molC
                
            # the backups
            
            Q_Exudbu = Q_Exud            
            Q_Mucilbu = Q_Mucil
            Q_Grbu = Q_Gr
            Q_Rmbu = Q_Rm
            
            Q_ST_dot    = np.array(r.Q_out_dot[0:Nt])* mmolSuc_to_molC
            Q_meso_dot  = np.array(r.Q_out_dot[Nt:(Nt*2)])* mmolSuc_to_molC
            Q_Rm_dot    = np.array(r.Q_out_dot[(Nt*2):(Nt*3)])* mmolSuc_to_molC
            Q_Exud_dot  = np.array(r.Q_out_dot[(Nt*3):(Nt*4)])* mmolSuc_to_molC
            Q_Gr_dot    = np.array(r.Q_out_dot[(Nt*4):(Nt*5)])* mmolSuc_to_molC
            Q_Rmmax_dot       = np.array(r.Q_out_dot[(Nt*5):(Nt*6)])* mmolSuc_to_molC
            Q_Grmax_dot       = np.array(r.Q_out_dot[(Nt*6):(Nt*7)])* mmolSuc_to_molC
            Q_S_meso_dot   = np.array(r.Q_out_dot[(Nt*7):(Nt*8)])* mmolSuc_to_molC
            Q_S_ST_dot   = np.array(r.Q_out_dot[(Nt*8):(Nt*9)])* mmolSuc_to_molC
            Q_Mucil_dot  = np.array(r.Q_out_dot[(Nt*9):(Nt*10)])* mmolSuc_to_molC #mol for nodes
            ##
            # att: that will be the cumulative value
            Q_ST    = np.array(r.Q_out[0:Nt]) * mmolSuc_to_molC
            Q_meso  = np.array(r.Q_out[Nt:(Nt*2)])* mmolSuc_to_molC
            Q_Rm    = np.array(r.Q_out[(Nt*2):(Nt*3)])* mmolSuc_to_molC
            Q_Exud  = np.array(r.Q_out[(Nt*3):(Nt*4)])* mmolSuc_to_molC 
            Q_Gr    = np.array(r.Q_out[(Nt*4):(Nt*5)])* mmolSuc_to_molC
            Q_Rmmax       = np.array(r.Q_out[(Nt*5):(Nt*6)])* mmolSuc_to_molC
            Q_Grmax       = np.array(r.Q_out[(Nt*6):(Nt*7)])* mmolSuc_to_molC
            Q_S_meso   = np.array(r.Q_out[(Nt*7):(Nt*8)])* mmolSuc_to_molC
            Q_S_ST   = np.array(r.Q_out[(Nt*8):(Nt*9)])* mmolSuc_to_molC
            Q_Mucil  = np.array(r.Q_out[(Nt*9):(Nt*10)])* mmolSuc_to_molC

            C_ST    = np.array(r.C_ST)* mmolSuc_to_molC
            Fl      = np.array(r.Fl)* mmolSuc_to_molC
            volST   = np.array(r.vol_ST)
            volMeso   = np.array(r.vol_Meso)
            C_S_meso   = Q_S_meso/volMeso
            C_S_ST   = Q_S_ST/volST
            C_meso  = Q_meso/volMeso
            Q_in   += sum(np.array(r.AgPhl)*dt)* mmolSuc_to_molC
            # i m missing the starch
            Q_out   = Q_Rm + Q_Exud + Q_Gr + Q_Mucil
            error_st_abs   = abs(sum(Q_ST + Q_meso + Q_out + Q_S_meso + Q_S_ST)- Q_in - sum(Q_ST_init)*2  - sum(Q_meso_init)*2)
            error_st_rel = abs(div0(error_st_abs,Q_in,1)*100)

            Q_Exudbu     =   np.concatenate((Q_Exudbu, np.full(Nt - Ntbu, 0.))) 
            Q_Mucilbu       =   np.concatenate((Q_Mucilbu, np.full(Nt - Ntbu, 0.))) 
            Q_Grbu     =   np.concatenate((Q_Grbu, np.full(Nt - Ntbu, 0.))) 
            Q_Rmbu       =   np.concatenate((Q_Rmbu, np.full(Nt - Ntbu, 0.)))
            Q_Gr_i      = (Q_Gr    - Q_Grbu)
            Q_Rm_i      = (Q_Rm    - Q_Rmbu)
            Q_Exud_i      = (Q_Exud    - Q_Exudbu)#*inflateVal
            Q_Mucil_i     = (Q_Mucil   - Q_Mucilbu)#*inflateVal

            # can get negative respiration at night
            #try:
            #    assert Q_in > 0
            #except:
            #    print(error_st_abs, Q_in, error_st_rel, rs.weatherX, r.plant.organTypes)
            #    raise Exception


            assert Q_Exud_i[0] == 0#no exudation in seed node I guess
            assert Q_Mucil_i[0] == 0#no exudation in seed node I guess
            assert np.array(r.Csoil_node)[0] == 0


            try:
                assert (np.array(r.Csoil_seg ) == np.array(r.Csoil_node)[1:]).all()
            except:
                print(np.array(r.Csoil_seg ), np.array(r.Csoil_node))
                print( (np.array(r.Csoil_seg ) == np.array(r.Csoil_node)[1:]),
                         (np.array(r.Csoil_seg ) == np.array(r.Csoil_node)[1:]).all())
                raise Exception

            Q_Exud_i_seg = np.array( Q_Exud_i[1:] )#*100 #from nod to semgment
            Q_Mucil_i_seg = np.array(Q_Mucil_i[1:])#*0.

            
            
            airSegsId = rs.airSegs
            #np.array(list(set(np.concatenate((rs.cell2seg.get(-1),np.where(np.array(rs.organTypes) != 2)[0])) )))#aboveground
            try:
                assert (Q_Exud_i_seg[airSegsId] == 0).all()
                assert (Q_Mucil_i_seg[airSegsId] == 0).all()
                assert (np.array(r.k_mucil_)[airSegsId+1] == 0).all()
                assert (np.array(r.Q_Exudmax)[airSegsId+1] == 0).all()
            except:
                print("Q_Exud_i_seg", Q_Exud_i_seg[airSegsId] )
                print("Q_Mucil_i", Q_Mucil_i_seg,Q_Mucil,Q_Mucilbu,airSegsId)
                print("Q_Mucil_i", Q_Mucil_i_seg[airSegsId], Q_Mucil[airSegsId+1], Q_Mucilbu[airSegsId+1])
                print("Csoil_seg", np.array(r.Csoil_seg)[airSegsId])
                print("k_mucil_",np.array(r.k_mucil_)[airSegsId+1])#,np.array(r.k_mucil_).size() ,Q_Exud.size())
                print("Q_Exudmax",np.array(r.Q_Exudmax)[airSegsId+1])
                print("airSegsId", airSegsId, np.where(airSegsId))
                print(len(airSegsId), len(r.k_mucil_))
                raise Exception

            # Q_Exud_i = Q_Mucil_i
            # Q_Exud_i[np.where(Q_Exud_i > 0.)] = 1.
            #r.outputFlux = np.array(r.outputFlux)/ 10
            
            try:
                write_file_array("Q_Exud_i_real", Q_Exud_i, directory_ =results_dir)# to see if get val < 0
                write_file_array("Q_Mucil_i_real", Q_Exud_i, directory_ =results_dir)# to see if get val < 0
                if ((min(Q_Exud_i_seg) < 0) or (min(Q_Mucil_i_seg)<0)):
                    write_file_float("timeNegativeExudMucil", rs_age, directory_ =results_dir)
                    
            
                #assert (min(Q_Exud_i_seg) >= -abs(error_st_abs)) and (min(Q_Exud_i_seg) >= -r.atol)
                Q_Exud_i_seg[np.where(Q_Exud_i_seg<0)] = 0.
                
                #assert (min(Q_Mucil_i_seg) >= -abs(error_st_abs)) and (min(Q_Mucil_i_seg) >= -r.atol ) 
                Q_Mucil_i_seg[np.where(Q_Mucil_i_seg<0)] = 0.
            except:
                print('error_st_abs',error_st_abs)
                print('C_ST',repr(C_ST))
                print('Csoil_node',repr(r.Csoil_node))
                print('Q_Exud_i_seg',repr(Q_Exud_i_seg))
                print('Q_Mucil_i_seg',repr(Q_Mucil_i_seg))
                write_file_array("Csoil_node_err", r.Csoil_node, directory_ =results_dir)
                write_file_array("C_ST_err", C_ST, directory_ =results_dir)
                write_file_array("Q_Exud_i_err", Q_Exud_i_seg, directory_ =results_dir)
                write_file_array("Q_Mucil_i_err", Q_Mucil_i_seg, directory_ =results_dir)
                write_file_array("Q_Exud_tot_err", Q_Exud_inflate, directory_ =results_dir)
                write_file_array("Q_Mucil_tot_err", Q_Mucil_inflate, directory_ =results_dir)
                write_file_array("Q_Mucil_err", Q_Mucil, directory_ =results_dir)
                write_file_array("Q_Exud_err", Q_Exud, directory_ =results_dir)
                write_file_array("Q_Mucil_dot_err", Q_Mucil_dot , directory_ =results_dir)
                write_file_array("Q_Exud_dot_err", Q_Exud_dot , directory_ =results_dir)
                failedExud = True
                #raise Exception

            print("sum exud", sum(Q_Exud_i_seg), sum(Q_Mucil_i_seg))
        elif rank > 0:
            Q_Exud_i_seg = None
            Q_Mucil_i_seg = None
            Q_Exud = None
            Q_Mucil = None
            error_st_rel = 0
            error_st_abs = 0
            Q_in = None
            
        Q_in = comm.bcast(Q_in, root = 0)
        error_st_rel = comm.bcast(error_st_rel, root = 0)
        error_st_abs = comm.bcast(error_st_abs, root = 0)
        if False:
            try:
                assert  (error_st_rel< 1.) or abs(Q_in) < 1e-13
            except:    
                print('error_st_abs',rank,error_st_abs, Q_in, error_st_rel)
                raise Exception
        #print(rank, 'share Q_Exud_i_segA')
        #comm.barrier()
        if mpiVerbose and rank==0:# or (max_rank == 1):
            print(rank, 'share Q_Exud_i_seg')
        #comm.barrier()

        Q_Exud = comm.bcast(Q_Exud, root = 0) 
        Q_Mucil = comm.bcast(Q_Mucil, root = 0) 
        Q_Exud_i_seg = comm.bcast(Q_Exud_i_seg, root = 0) 
        Q_Mucil_i_seg = comm.bcast(Q_Mucil_i_seg, root = 0) 

        #comm.barrier()
        if mpiVerbose and rank==0:# or (max_rank == 1):
            print(rank, 'print data to linux')
        #comm.barrier()
        if (rank == 0) and (mode != "dumux_w")  :
            print("\n\n\n\t\tat ", int(rs_age//1),"d", int(rs_age%1*24),"h", int(rs_age%1*24%1*60),"mn",  
                  round(r.Qlight *1e6),"mumol m-2 s-1")
            print("Error in Suc_balance:\n\tabs (mmol) {:5.2e}\trel (-) {:5.2e}".format(error_st_abs, error_st_rel))
            print("Error in photos:\n\tabs (cm3/day) {:5.2e}".format(errLeuning_abs))
            print("C_ST (mol ml-1):\n\tmean {:.2e}\tmin  {:5.2e} at {:d} segs \tmax  {:5.2e}".format(np.mean(C_ST), min(C_ST),
                                                                                                      len(np.where(C_ST == min(C_ST) )[0]), max(C_ST)))        
            print("C_me (mol ml-1):\n\tmean {:.2e}\tmin  {:5.2e}\tmax  {:5.2e}".format(np.mean(C_meso), min(C_meso), max(C_meso)))        
            print('Q_X (mol Suc): \n\tST   {:.2e}\tmeso {:5.2e}\tin   {:5.2e}'.format(sum(Q_ST), sum(Q_meso), Q_in))
            print('\tRm   {:.2e}\tGr   {:.2e}\tExud {:5.2e}'.format(sum(Q_Rm), sum(Q_Gr), sum(Q_Exud)))

            if min(C_ST) < 0.0:
                print("min(C_ST) < 0.0", min(C_ST),np.mean(C_ST),max(C_ST))
                raise Exception

            tiproots, tipstems, tipleaves = r.get_organ_segments_tips()
            write_file_array("root_segments_tips",tiproots, 
                             directory_ =results_dir)
            write_file_array("rhizoSegsId", np.array(rs.rhizoSegsId), 
                             directory_ =results_dir)
            write_file_array("Q_ST", Q_ST, directory_ =results_dir)#mmol
            write_file_array("C_ST", C_ST, directory_ =results_dir)#mmol/cm3
            write_file_array("C_meso", C_meso, directory_ =results_dir)
            write_file_array("Q_meso", Q_meso, directory_ =results_dir)


            write_file_array("Q_S_ST", Q_S_ST, directory_ =results_dir)#mmol
            write_file_array("C_S_ST", C_S_ST, directory_ =results_dir)#mmol/cm3
            write_file_array("C_S_meso", C_S_meso, directory_ =results_dir)
            write_file_array("Q_S_meso", Q_S_meso, directory_ =results_dir)

            write_file_array("Q_Rm", Q_Rm, directory_ =results_dir)
            write_file_array("Q_Mucil", Q_Mucil, directory_ =results_dir)
            write_file_array("Q_Exud", Q_Exud, directory_ =results_dir)
            write_file_array("Q_Gr", Q_Gr, directory_ =results_dir)
            write_file_array("psiXyl", r.psiXyl, directory_ =results_dir)
            write_file_array("psi_guardCell", r.pg, directory_ =results_dir)
            write_file_array("gco2", r.gco2, directory_ =results_dir)
            write_file_array("Fpsi", r.Fpsi, directory_ =results_dir)
            write_file_array("fw", r.fw, directory_ =results_dir)
            write_file_array("Q_Ag_dot_inner", r.AgCumul_inner, directory_ =results_dir)
            write_file_array("Q_Ag_dot_innerbis", r.An, directory_ =results_dir)
            write_file_array("Q_Ag_dot", np.array(r.AgPhl)* mmolSuc_to_molC, directory_ =results_dir)
            write_file_float("Q_Ag", Q_in, directory_ =results_dir)
            write_file_array("C_rsi", np.array(r.Csoil_seg ), 
                             directory_ =results_dir)#mmol/cm3
            write_file_array("CSTi_delta", np.array(r.CSTi_delta ), 
                             directory_ =results_dir)
            write_file_array("CSTi_exud", np.array(r.CSTi_exud ), 
                             directory_ =results_dir)
            write_file_array("Crsi_exud", np.array(r.Crsi_exud ), 
                             directory_ =results_dir)
            write_file_array("Q_Exudmax",np.array(r.Q_Exudmax), directory_ =results_dir)
            
            
        
            write_file_array("Q_ST_dot", Q_ST_dot   , directory_ =results_dir)
            write_file_array("Q_meso_dot", Q_meso_dot  , directory_ =results_dir)
            write_file_array("Q_Rm_dot", Q_Rm_dot   , directory_ =results_dir)
            write_file_array("Q_Gr_dot", Q_Gr_dot , directory_ =results_dir)
            write_file_array("Q_Rmmax_dot", Q_Rmmax_dot , directory_ =results_dir)
            write_file_array("Q_Grmax_dot", Q_Grmax_dot  , directory_ =results_dir) 
            write_file_array("Q_S_meso_dot", Q_S_meso_dot, directory_ =results_dir)
            write_file_array("Q_S_ST_dot", Q_S_ST_dot , directory_ =results_dir)
            
            write_file_array("Q_Exud_dot", Q_Exud_dot, directory_ =results_dir)
            write_file_array("Q_Mucil_dot", Q_Mucil_dot , directory_ =results_dir)
            

            if not doMinimumPrint:
                write_file_array("TotSoilC", s.getTotCContent(), directory_ =results_dir)
            write_file_array("Q_Gr_i", Q_Gr_i, directory_ =results_dir)
            write_file_array("Q_Rm_i",Q_Rm_i, directory_ =results_dir)
            write_file_array("Q_Exud_i", Q_Exud_i, directory_ =results_dir)
            write_file_array("Q_Mucil_i",Q_Mucil_i, directory_ =results_dir)
            write_file_float("Q_Exud_tot", Q_Exud_inflate, directory_ =results_dir)
            write_file_float("Q_Mucil_tot", Q_Mucil_inflate, directory_ =results_dir)
            write_file_float("Q_Gr_tot", sum(Q_Gr), directory_ =results_dir)
            write_file_float("Q_Rm_tot", sum(Q_Rm), directory_ =results_dir)
        if rank == 0:
            datas = [np.array(r.AgPhl)* mmolSuc_to_molC,
                     C_ST, C_S_ST, C_meso, 
                     C_S_meso, Q_Gr_i,Q_Rm_i, Q_Exud_dot,Q_Mucil_dot,
                     Q_Exud,Q_Mucil,
                     r.psiXyl, Q_Exud_i, Q_Mucil_i
                     #, np.array(r.rs.isRootTip, dtype=int)
                    ]
            datasName = ["Ag",
                         "C_ST", "C_S_ST", "C_meso", 
                         "C_S_meso", "Q_Gr_i","Q_Rm_i", "Q_Exud_dot","Q_Mucil_dot",
                         "Q_Exud", "Q_Mucil",
                     "psiXyl", "Q_Exud_i","Q_Mucil_i" 
                         #,"isRootTip"
                        ]
        else:
            datas = []
            datasName = []
        if int(rs_age *1000)/1000-int(rs_age) == 0.5 :# midday
            getVTPOutput(r,rs,s,datas, datasName, rs_age*100, results_dir, min_b, max_b, cell_number)
        
        failedExud = comm.bcast(failedExud, root = 0)
        assert not failedExud


    """ output """
    if mpiVerbose and rank==0:# or (max_rank == 1):
        print('finished simulation')
        
    print("fin", rank)
    return results_dir

        

if __name__ == '__main__':
    # python3 XcGrowth.py 9 dumux_10c 10 0 customDry noAds 9.02 0.02
    # python3 XcGrowth.py 9 dumux_10c 10 17 lateDry
    # python3 XcGrowth.py 12 dumux_10c 25 98 baseline
    # python3 XcGrowth.py 10 dumux_10c 25 3 none
    # python3 XcGrowth.py 9 dumux_10c 9.001 5 none
    # python3 XcGrowth.py 9 dumux_10c 9.06 0 customDry xx 9.02 0.02
    # python3 XcGrowth.py 10 dumux_10c 10.06 0 customDry xx 10.02 0.02
    # python3 XcGrowth.py 9 dumux_10c 10 3 none
    # python3 XcGrowth.py 10 dumux_10c 14 85 none
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
