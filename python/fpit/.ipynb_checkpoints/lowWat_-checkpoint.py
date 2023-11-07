""" 
    Maize using rhizosphere models  
"""
import matplotlib; matplotlib.use('agg')
import sys;
import os
absolutePath = False
if absolutePath:
    os.chdir( "/home/rbtlm640/DUMUXr/dumux-rosi/python/fpit")
    sys.path.append( "/home/rbtlm640/DUMUXr/dumux-rosi/python/fpit")
    sys.path.append("/home/rbtlm640/DUMUXr/dumux-rosi/python/fpit/data/");
    sys.path.append("/home/rbtlm640/DUMUXr/dumux-rosi/python/modules_fpit/");
    sys.path.append("/home/rbtlm640/DUMUXr/CPlantBox/");#absolut path work better when
    # using pycharm
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

import scenario_setup as scenario
scenario = reload(scenario)
import rhizo_modelsPlant  # Helper class for cylindrical rhizosphere models
rhizo_modelsPlant = reload(rhizo_modelsPlant)
from rhizo_modelsPlant import *
import evapotranspiration as evap
#import cyl_exu
import cyl3plant as cyl3
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import os
from scenario_setup import write_file_array, write_file_float, div0, div0f


if __name__ == '__main__':

    initsim =float(sys.argv[1])# initsim = 9.5
    
    dt = 1/60/24
    p_mean = -1000
    k_iter = 20
    l_ks =  "dx"#"root", "dx", "dx_2"
    organism = "plant"# "RS"#
    weightBefore = False
    SRIBefore = False
    beforeAtNight = True
    adaptRSI_  = False
    lightType =""#+- "nolight" # or ""
    extraName = ""
    results_dir="./results/"+extraName+lightType+l_ks+str(int(weightBefore))\
                    +str(int(SRIBefore))+str(int(beforeAtNight))+str(int(adaptRSI_))\
                        +organism+str(k_iter)+"k_"+str(initsim)\
                    +"_"+str(int(dt*24*60))+"mn_"\
                    +str(int((dt*24*60 - int(dt*24*60))*60))+"s_"\
                    +str(max_rank)+"_"+str(abs(p_mean))+"/"
    print('results_dir',results_dir)
    if rank == 0:
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)
        else:
            test = os.listdir(results_dir)
            for item in test:
                try:
                    os.remove(results_dir+item)
                except:
                    pass
    comm.barrier()
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
    soil_ = scenario.vg_SPP(0)

    min_b = [-5, -5, -10.] 
    max_b = [5, 5, 0.] 
    cell_number = [5,5,20]
    #min_b = [-5., -5, -5.] 
    #max_b = [5., 5, 0.] 
    #cell_number = [5, 5, 5]
    sim_time = 1 #154   #  [day]
    # dt_inner = 360 / (24 * 3600)  # time step [day] 20


    """ rhizosphere model parameters """
    wilting_point = -15000  # cm
    recreateComsol = False

    periodic = False
    nc = 10

    logbase = 0.5  # according to Mai et al. (2019)
    mode = "dumux_w"  

    """ initialize """


    s, soil = scenario.create_soil_model(soil_type, year, soil_,#comp, 
                min_b, max_b, cell_number, demoType = mode, times = None, net_inf = None,
                usemoles = usemoles, dirResults = results_dir, p_mean_ = p_mean)

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
    


    # all thread need a plant object, but only thread 0 will make it grow
    rs, r = scenario.create_mapped_plant(wilting_point, nc, logbase, mode,initsim,
                                            min_b, max_b, cell_number, s, xml_name,
                                            path, plantType = organism, 
                                            recreateComsol_ = recreateComsol,
                                            usemoles = usemoles,l_ks_ = l_ks,
                                            limErr1d3d = 5e-8*100)  # pass parameter file for dynamic growth

    rs.weightBefore = weightBefore
    rs.SRIBefore = SRIBefore
    rs.beforeAtNight = beforeAtNight
    rs_age = initsim


    net_sol_flux =  np.array([np.array([]),np.array([])])
    net_flux = np.array([])


    Q_ST_init = np.array([])
    Nt = len(rs.nodes)
    Q_Exudbu    = np.zeros(Nt)
    Q_Mucilbu   = np.zeros(Nt)
    Q_in  = 0
    Nt = len(rs.nodes)
    r.minLoop = 1000
    r.maxLoop = 5000
    simMax = initsim + 3

    TranspirationCumul = 0
    cell_volumes = s.getCellVolumes()  # cm3
    buTotCSoilInit = sum(s.getTotCContent())
    buWSoilInit = sum(np.multiply(np.array(s.getWaterContent()), cell_volumes))

    start_time_global = timeit.default_timer()
    time_rhizo_cumul = 0
    time_3ds_cumul = 0
    r.time_plant_cumulW = 0
    r.time_plant_cumulS = 0

    r.time_rhizo_i = 0
    r.time_3ds_i = 0
    Q_Exud_inflate = 0.; Q_Mucil_inflate = 0.
    rs.results_dir = results_dir
    #testWeather = [scenario.weather(testT/100.)["Qlight"]  for testT in range(100)]
    static_plant = True
    s.buTotCAfter = 0
    s.buTotCBefore = 0
    print('start loop', rank)
    secId = None
    Q_Exud_i = None
    Q_Exud_i_seg = np.array([]); Q_Mucil_i_seg = np.array([])
    error_st_abs = 0;error_st_rel=0

    #errs =np.array(["errRxPlant", "errW1ds", "errW3ds", "max(r.SinkLim3DS)","max(r.SinkLim1DS)","max(r.maxDiff1d3dCW_abs)", 
    #            "errWrsi", "maxDiff1d3dCW_absBU", 
    #            "bulkMassErrorWater_abs","bulkMassErrorWater_absLim","rhizoMassWError_absLim","rhizoMassWError_abs",
    #            "sum(abs(diffBCS1dsFluxIn))", "sum(abs(diffBCS1dsFluxOut))","sum(abs(diffouter_R_bc_wat))",
    #            "diff1d3dCurrant","diff1d3dCurrant_rel","rhizoMassWError_rel "])
    errs = np.array(["errRxPlant", "errW1ds", "errW3ds", 
                      "max(r.SinkLim3DS)","max(r.SinkLim1DS)","max(r.maxDiff1d3dCW_abs)", 
                "errWrsi", "maxDiff1d3dCW_absBU", 
                "bulkMassErrorWater_abs","bulkMassErrorWater_absLim",
                      "rhizoMassWError_absLim","rhizoMassWError_abs",
                "sum(abs(diffBCS1dsFluxIn))", "sum(abs(diffBCS1dsFluxOut))","sum(abs(diffouter_R_bc_wat))",
                "diff1d3dCurrant","rhizoMassWError_rel",'err'])
    write_file_array("error", errs, directory_ =results_dir, fileType = '.csv')
    write_file_array("fpit_error", errs, directory_ =results_dir, fileType = '.csv') 
    seg_fluxes_ = np.array([])
    while rs_age < simMax: #for i, dt in enumerate(np.diff(times)):

        rs_age += dt
        print("Day", rs_age)
        seg2cell_old = rs.seg2cell
        Ntbu = Nt
        print(rank,'simulating')
        if (rank == 0) and (not static_plant) :

            rs.simulate(dt)#(rank == 0))  # because of dumux.pick(), blocks if I do not let all threads do the simulate.
            seg2cell_new = rs.seg2cell

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


        print(rank,'simulated')
        rs.update()
        print(rank,'updated')
        write_file_array("organTypes", np.array(rs.organTypes), directory_ =results_dir)
        
                           
        if False:   
            for lId, cyl in enumerate(rs.cyls):
                if not isinstance(cyl, AirSegment):
                    gId = rs.eidx[lId]
                    write_file_float("segIdxCyl"+str(gId),gId, directory_ =results_dir, allranks = True)
                    write_file_array("pressureHeadcyl"+str(gId),np.array(cyl.getSolutionHead()).flatten(), directory_ =results_dir, allranks = True)
                    write_file_array("coordcyl"+str(gId), cyl.getDofCoordinates().flatten(), directory_ =results_dir, allranks = True)
                    write_file_array("solute_conc_cyl"+str(gId)+"_"+str(1), np.array(cyl.getSolution(1)).flatten()* rs.molarDensityWat_m3/1e6 , 
                                     directory_ =results_dir, allranks = True) 


        comm.barrier()

        dt_inner = dt

        Q_Exud_inflate += sum(Q_Exud_i_seg); Q_Mucil_inflate += sum(Q_Mucil_i_seg)

        print(rank, 'to cyl3.simulate_const')
        net_sol_flux, net_flux, seg_fluxes_ = cyl3.simulate_const(s, 
                                                r,  dt, dt_inner, rs_age, 
                                                Q_plant=[Q_Exud_i_seg, Q_Mucil_i_seg],
                                                r= rs,plantType = organism,
                                                 adaptRSI=adaptRSI_,
                                                wilting_point = wilting_point,
                                                outer_R_bc_sol = net_sol_flux, 
                                                outer_R_bc_wat = net_flux,seg_fluxes=seg_fluxes_,
                                                results_dir = results_dir,
                                                    k_iter_ = k_iter,lightType_=lightType)
        print(rank, 'left cyl3.simulate_const')
        time_rhizo_cumul += r.time_rhizo_i
        time_3ds_cumul += r.time_3ds_i
        r.time_rhizo_i = 0
        r.time_3ds_i = 0
        time_plant_cumul = r.time_plant_cumulW + r.time_plant_cumulS 

        if True:
            write_file_array("cellVol", np.array(s.getCellVolumes()), directory_ =results_dir) # cm3 
            write_file_array("theta", np.array(s.getWaterContent()), directory_ =results_dir) 
            for i in range(rs.numFluidComp):
                write_file_array("Soil_solute_conc"+str(i+1), np.array(s.getSolution(i+1)).flatten()* rs.molarDensityWat_m3/1e6, directory_ =results_dir) 
            for i in range(rs.numFluidComp, rs.numComp):
                write_file_array("Soil_solute_conc"+str(i+1), np.array(s.getSolution(i+1)).flatten()* rs.bulkDensity_m3 /1e6 , directory_ =results_dir) 

            write_file_array("Soil_solute_conc"+str(rs.numComp+1), np.array(s.base.getCSS1_out()).flatten()[:-1]* rs.bulkDensity_m3 /1e6 , directory_ =results_dir) 
        print(rank, 'did s.data writing')
        errLeuning_abs = abs(sum(r.outputFlux))
        if organism == "plant":
            TranspirationCumul += sum(np.array(r.Ev) * dt) #transpiration [cm3/day] * day
        else:
            TranspirationCumul += sum(r.outputFlux)
        

        buTotCAfter = sum(s.getTotCContent())   
        buWAfter = sum(np.multiply(np.array(s.getWaterContent()), cell_volumes))    

        if rank == 0:
            s.bulkMassErrorCumul_abs = abs((buTotCAfter - ( buTotCSoilInit + Q_Exud_inflate + Q_Mucil_inflate)))#so that only works if infalte
            s.bulkMassErrorWaterCumul_abs = abs(buWAfter - ( buWSoilInit - TranspirationCumul))
            if buTotCAfter != 0:
                s.bulkMassErrorCumul_rel = abs(s.bulkMassErrorCumul_abs/buTotCAfter*100)
            else:
                s.bulkMassErrorCumul_rel =np.nan
            s.bulkMassErrorWaterCumul_rel = abs(s.bulkMassErrorWaterCumul_abs/buWAfter*100)

        print(rank, 'got s.errorCumul')
        write_file_array("totalComputetime",np.array([timeit.default_timer() - start_time_global,
                            time_plant_cumul,time_rhizo_cumul ,time_3ds_cumul]) , directory_ =results_dir)
        write_file_array("time", np.array([rs_age,r.Qlight]), directory_ =results_dir)
        print(rank, 'write some otehr stuff')
        if (mode == "dumux_10c"):
            write_file_array("TotSoilC", s.getTotCContent(), directory_ =results_dir)
            write_file_float("Q_Exud_i", sum(Q_Exud_i_seg), directory_ =results_dir)
            write_file_float("Q_Mucil_i", sum(Q_Mucil_i_seg), directory_ =results_dir)
            write_file_float("Q_Exud_tot", Q_Exud_inflate, directory_ =results_dir)
            write_file_float("Q_Mucil_tot", Q_Mucil_inflate, directory_ =results_dir)
        
        if rank == 0:
            #absolute and relative (%) error
            write_file_array("errorsPlant", np.array([error_st_abs,error_st_rel,#cumulative
                                                errLeuning_abs]), directory_ =results_dir) #not cumulative
            write_file_array("errorsBulkSoil", np.array([s.bulkMassErrorPlant_abs, s.bulkMassErrorPlant_rel, #not cumulative 
                                                s.bulkMassError1ds_abs, s.bulkMassError1ds_rel, 
                                                s.bulkMassErrorCumul_abs,s.bulkMassErrorCumul_rel,#cumulative
                                                s.bulkMassErrorWater_abs,s.bulkMassErrorWater_rel, #not cumulative
                                                s.bulkMassErrorWaterCumul_abs,s.bulkMassErrorWaterCumul_rel]), directory_ =results_dir)#cumulative
            write_file_array("errorMassRhizo", np.array([rs.rhizoMassCError_abs, rs.rhizoMassCError_rel,
                                                        rs.rhizoMassWError_abs, rs.rhizoMassWError_rel]), directory_ =results_dir)# not cumulative
            write_file_array("sumErrors1ds3ds", np.concatenate((rs.sumDiff1d3dCW_abs, rs.sumDiff1d3dCW_rel)), directory_ =results_dir)
            write_file_array("maxErrors1ds3ds", np.concatenate((rs.maxDiff1d3dCW_abs, rs.maxDiff1d3dCW_rel)), directory_ =results_dir)# cumulative (?)

            write_file_array("trans", r.Ev, directory_ =results_dir)
            write_file_array("transrate",r.Jw, directory_ =results_dir)
            write_file_array("transrate",r.Jw, directory_ =results_dir)
        print(rank, 'finished other data writing')
        try:
            assert abs(s.bulkMassErrorPlant_abs)  < 1e-5
        except:
            print("\n\n\nissue bulk soil balance", np.array([s.bulkMassErrorPlant_abs, s.bulkMassErrorPlant_rel, #not cumulative 
                                            s.bulkMassErrorCumul_abs,s.bulkMassErrorCumul_rel,#cumulative
                                            s.bulkMassErrorWater_abs,s.bulkMassErrorWater_rel, #not cumulative
                                            s.bulkMassErrorWaterCumul_abs,s.bulkMassErrorWaterCumul_rel]))
            print("\n\n\n")
            raise Exception

        print(rank, 'do C growth?',  (mode == "dumux_10c") and (rank == 0) and ((not static_plant) or (rs_age == initsim+dt)) and (organism == "plant"))

        if (mode == "dumux_10c") and (rank == 0) and ((not static_plant) or (rs_age == initsim+dt)) and (organism == "plant"):

            startphloem=rs_age
            endphloem = rs_age + dt
            stepphloem = 1
            verbose_phloem = True
            filename = "results/" +"inPM.txt"
            print("startpm")

            start_time_plant = timeit.default_timer()

            r.startPM(startphloem, endphloem, stepphloem, ( weatherX["TairC"]  +273.15) , verbose_phloem, filename)

            r.time_plant_cumulS += (timeit.default_timer() - start_time_plant)


            Nt = len(rs.nodes)
            if r.withInitVal and (len(Q_ST_init) ==0) :
                Q_ST_init = np.array(r.Q_init[0:Nt])/1e3
                Q_meso_init = np.array(r.Q_init[Nt:(Nt*2)])/1e3

            # att: that will be the cumulative value
            #  MMOL(/cm3) => mol(/cm3)
            inflateVal = 1#1e3
            Q_ST    = np.array(r.Q_out[0:Nt])/1e3
            Q_meso  = np.array(r.Q_out[Nt:(Nt*2)])/1e3
            Q_Rm    = np.array(r.Q_out[(Nt*2):(Nt*3)])/1e3
            Q_Exud  = np.array(r.Q_out[(Nt*3):(Nt*4)])/1e3 
            Q_Gr    = np.array(r.Q_out[(Nt*4):(Nt*5)])/1e3
            Q_Rmmax       = np.array(r.Q_out[(Nt*5):(Nt*6)])/1e3
            Q_Grmax       = np.array(r.Q_out[(Nt*6):(Nt*7)])/1e3
            Q_S_meso   = np.array(r.Q_out[(Nt*7):(Nt*8)])/1e3
            Q_S_ST   = np.array(r.Q_out[(Nt*8):(Nt*9)])/1e3
            Q_Mucil  = np.array(r.Q_out[(Nt*9):(Nt*10)])/1e3 #mol for nodes

            C_ST    = np.array(r.C_ST)/1e3
            Fl      = np.array(r.Fl)/1e3
            volST   = np.array(r.vol_ST)
            volMeso   = np.array(r.vol_Meso)
            C_S_meso   = Q_S_meso/volMeso
            C_S_ST   = Q_S_ST/volST
            C_meso  = Q_meso/volMeso
            Q_in   += sum(np.array(r.AgPhl)*dt)/1e3
            # i m missing the starch
            Q_out   = Q_Rm + Q_Exud + Q_Gr + Q_Mucil
            error_st_abs   = abs(sum(Q_ST + Q_meso + Q_out + Q_S_meso + Q_S_ST)- Q_in - sum(Q_ST_init)  - sum(Q_meso_init))
            error_st_rel = abs(div0(error_st_abs,Q_in,1)*100)

            Q_Exudbu     =   np.concatenate((Q_Exudbu, np.full(Nt - Ntbu, 0.))) 
            Q_Mucilbu       =   np.concatenate((Q_Mucilbu, np.full(Nt - Ntbu, 0.))) 

            Q_Exud_i      = (Q_Exud    - Q_Exudbu)*inflateVal
            Q_Mucil_i     = (Q_Mucil   - Q_Mucilbu)*inflateVal

            try:
                assert Q_in > 0
            except:
                print(error_st_abs, Q_in, error_st_rel, weatherX, r.plant.organTypes)
                raise Exception

            try:
                assert  (error_st_rel< 1.) or abs(Q_in) < 1e-13
            except:    
                print(error_st_abs, Q_in, error_st_rel)
                raise Exception

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

            Q_Exud_i_seg = np.array( Q_Exud_i[1:] ) #from nod to semgment
            Q_Mucil_i_seg = np.array(Q_Mucil_i[1:])

            airSegsId = np.array(list(set(np.concatenate((rs.cell2seg.get(-1),np.where(np.array(rs.organTypes) != 2)[0])) )))#aboveground
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
                print("k_mucil_",r.k_mucil_)#,np.array(r.k_mucil_).size() ,Q_Exud.size())
                print("Q_Exudmax",np.array(r.Q_Exudmax)[airSegsId+1])
                print("airSegsId", airSegsId, np.where(airSegsId))
                print(len(airSegsId), len(r.k_mucil_))
                raise Exception

            # Q_Exud_i = Q_Mucil_i
            # Q_Exud_i[np.where(Q_Exud_i > 0.)] = 1.
            #r.outputFlux = np.array(r.outputFlux)/ 10
            try:
                assert min(Q_Exud_i_seg) >= 0.
            except:
                print(C_ST, r.Csoil_node, Q_Exud_i_seg,Q_Exud)
                raise Exception

            print("sum exud", sum(Q_Exud_i_seg), sum(Q_Mucil_i_seg))
        elif rank > 0:
            Q_Exud_i_seg = None
            Q_Mucil_i_seg = None
        
        print(rank, 'share Q_Exud_i_seg')

        Q_Exud_i_seg = comm.bcast(Q_Exud_i_seg, root = 0) 
        Q_Mucil_i_seg = comm.bcast(Q_Mucil_i_seg, root = 0) 

        print(rank, 'print data to linux')
        if (rank == 0) and (mode == "dumux_10c")  :
            print("\n\n\n\t\tat ", int(np.floor(rs_age)),"d", int((rs_age%1)*24),"h",  round(r.Qlight *1e6),"mumol m-2 s-1")
            print("Error in Suc_balance:\n\tabs (mmol) {:5.2e}\trel (-) {:5.2e}".format(error_st_abs, error_st_rel))
            print("Error in photos:\n\tabs (cm3/day) {:5.2e}".format(errLeuning_abs))
            print("C_ST (mmol ml-1):\n\tmean {:.2e}\tmin  {:5.2e} at {:d} segs \tmax  {:5.2e}".format(np.mean(C_ST), min(C_ST), len(np.where(C_ST == min(C_ST) )[0]), max(C_ST)))        
            print("C_me (mmol ml-1):\n\tmean {:.2e}\tmin  {:5.2e}\tmax  {:5.2e}".format(np.mean(C_meso), min(C_meso), max(C_meso)))        
            print('Q_X (mmol Suc): \n\tST   {:.2e}\tmeso {:5.2e}\tin   {:5.2e}'.format(sum(Q_ST), sum(Q_meso), Q_in))
            print('\tRm   {:.2e}\tGr   {:.2e}\tExud {:5.2e}'.format(sum(Q_Rm), sum(Q_Gr), sum(Q_Exud)))

            if min(C_ST) < 0.0:
                print("min(C_ST) < 0.0", min(C_ST),np.mean(C_ST),max(C_ST))
                raise Exception


            write_file_array("Q_ST", Q_ST, directory_ =results_dir)#mmol
            write_file_array("C_ST", C_ST, directory_ =results_dir)#mmol/cm3
            write_file_array("C_meso", C_meso, directory_ =results_dir)
            write_file_array("Q_meso", Q_meso, directory_ =results_dir)


            write_file_array("Q_S_ST", Q_S_ST, directory_ =results_dir)#mmol
            write_file_array("C_S_ST", C_S_ST, directory_ =results_dir)#mmol/cm3
            write_file_array("C_S_meso", C_S_meso, directory_ =results_dir)
            write_file_array("Q_S_meso", Q_S_meso, directory_ =results_dir)

            write_file_array("Q_Rm", Q_Rm, directory_ =results_dir)
            write_file_array("Q_Exud", Q_Exud, directory_ =results_dir)
            write_file_array("Q_Gr", Q_Gr, directory_ =results_dir)
            write_file_array("psiXyl", r.psiXyl, directory_ =results_dir)
            write_file_array("Fpsi", r.Fpsi, directory_ =results_dir)
            write_file_array("fw", r.fw, directory_ =results_dir)
            write_file_array("gco2", r.gco2, directory_ =results_dir)
            write_file_array("Q_Ag_dot", r.AgPhl, directory_ =results_dir)
            write_file_float("Q_Ag", Q_in, directory_ =results_dir)
            write_file_array("C_rsi", np.array(r.Csoil_seg ), directory_ =results_dir)#mmol/cm3
        
        print(rank, 'print data to linux')


    """ output """
    print('finished loop, get images')
    sizeSoilCell = rs.soilModel.getCellVolumes() #cm3
    rs.checkMassOMoleBalance2( sourceWat = np.full(len(sizeSoilCell),0.), # cm3/day 
                                     sourceSol = np.full((rs.numComp, len(sizeSoilCell)),0.), # mol/day
                                     dt = 0.,        # day    
                                     seg_fluxes = 0.,# [cm3/day]
                                     doWater = True, doSolute = True, doSolid = True, diff1d3dCW_abs_lim = np.Inf)
    vp.write_soil("results/"+str(sim_time)+"_Soil", s, min_b, max_b, cell_number, ["C concentration [g/cm³]"])

    if False:

        for konz in np.array([True]):#, False]):
            extraArray_ = rs.soilModel.getWaterContent()
            if not konz:
                extraArrayName_ = "theta (cm3)"
                extraArray_ *= sizeSoilCell
            else:
                extraArrayName_ = "theta (cm3/cm3)"
            print("idcomp", 0, rank)
            rp = rs.get_concentration(0, konz)

            vp.plot_roots_and_soil(rs.mappedSegments(),extraArrayName_,rp, s, periodic, min_b, max_b, cell_number, 
                    soil_type+genotype+"_rx", sol_ind =-1,extraArray = extraArray_, extraArrayName = extraArrayName_)  # VTK vizualisation
            print("idcomp_done", 0, rank)
            for i in range(1, rs.numComp+1):
                print("idcomp", i, rank)
                extraArray_ = rs.soilModel.getSolution(i) * rs.phaseDensity(i)/1e6

                if not konz:
                    extraArrayName_ = "C"+str(i)+" mol"
                else:
                    extraArrayName_ = "[C"+str(i)+"] (mol/cm3)"
                    if i <= rs.numFluidComp:
                        extraArray_ /= (rs.soilModel.getWaterContent() *sizeSoilCell) #cm3 water
                    else: 
                        extraArray_ /= sizeSoilCell #cm3


                vp.plot_roots_and_soil(rs.mappedSegments(),extraArrayName_ ,rs.get_concentration(i , konz), s, periodic, min_b, max_b, cell_number, 
                        soil_type+genotype+"_rx", sol_ind =-1,extraArray = extraArray_, 
                        extraArrayName = extraArrayName_)  # VTK vizualisation
                print("idcomp_done", i, rank)
    # to plot: cumulative An, Exud, mucil and each element in the soil
    # Also images of the 3D soil.
    print("fin", rank)
