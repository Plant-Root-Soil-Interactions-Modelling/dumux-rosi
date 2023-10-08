""" 
    Maize using rhizosphere models  
"""
import matplotlib; matplotlib.use('agg')
import sys;
sys.path.append("/data/");
sys.path.append("../modules_mpi/");
sys.path.append("../../../CPlantBox/");
sys.path.append("../../../CPlantBox/src")
import plantbox as pb  # CPlantBox
import visualisation.vtk_plot as vp
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import timeit
import numpy as np

import scenario_setup as scenario
from rhizo_modelsPlant import *  # Helper class for cylindrical rhizosphere models
import evapotranspiration as evap
#import cyl_exu
import cyl3plant as cyl3
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import os
from scenario_setup import write_file_array, write_file_float, div0, div0f

results_dir="./results/parallel/"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    test = os.listdir(results_dir)
    for item in test:
        try:
            os.remove(results_dir+item)
        except:
            pass

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
soil_, min_b, max_b, cell_number, area, Kc = scenario.maize_SPP(soil_type)
#min_b = [-5., -5, -5.] 
#max_b = [5., 5, 0.] 
#cell_number = [5, 5, 5]
sim_time = 1 #154   #  [day]
# dt_inner = 360 / (24 * 3600)  # time step [day] 20

x_, y_, lai =[],[],[] #evap.net_infiltration(year, soil_type, genotype, sim_time, Kc)
trans_maize = []#evap.get_transpiration(year, sim_time, area, lai, Kc)


""" rhizosphere model parameters """
wilting_point = -15000  # cm
recreateComsol = False

periodic = False
if recreateComsol:
    nc = 500  # dof+1
else:
    nc = 10
    
logbase = 0.5  # according to Mai et al. (2019)
mode = "dumux_10c"  

""" initialize """

initsim = 5.5
          
s, soil = scenario.create_soil_model(soil_type, year, soil_,#comp, 
            min_b, max_b, cell_number, demoType = mode, times = x_, net_inf = y_,
            usemoles = usemoles)

          
path = "../../../CPlantBox/modelparameter/structural/plant/"
xml_name = "Triticum_aestivum_test_2021.xml"  # root growth model parameter file


# all thread need a plant object, but only thread 0 will make it grow
rs, r = scenario.create_mapped_plant(wilting_point, nc, logbase, mode,initsim,
                                        min_b, max_b, cell_number, s, xml_name,
                                        path, plantType = "plant", 
                                        recreateComsol_ = recreateComsol,
                                        usemoles = usemoles)  # pass parameter file for dynamic growth



psi_x_, psi_s_, sink_, x_, y_, psi_s2_, soil_c_, c_, mass_soil_c_, dist, conc, l, a = [], [], [], [], [], [], [], [], [], [], [], [], []

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
dt = 1/24
simMax = initsim + 1

TranspirationCumul = 0
cell_volumes = s.getCellVolumes()  # cm3
buTotCSoilInit = sum(s.getTotCContent())
buWSoilInit = sum(np.multiply(np.array(s.getWaterContent()), cell_volumes))

start_time_global = timeit.default_timer()
time_rhizo_cumul = 0
time_3ds_cumul = 0
time_plant_cumul = 0
Q_Exud_inflate = 0.; Q_Mucil_inflate = 0.
rs.results_dir = results_dir
#testWeather = [scenario.weather(testT/100.)["Qlight"]  for testT in range(100)]

print('start loop', rank)
while rs_age < simMax: #for i, dt in enumerate(np.diff(times)):


    rs_age += dt
    print("Day", rs_age)
    seg2cell_old = rs.seg2cell
    Ntbu = Nt
    
    print('simulating again' ,rank)
    if rank == 0:
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
    
    CellVolumes = s.getCellVolumes()
    if rank == 0:
        cellIds = np.fromiter(rs.cell2seg.keys(), dtype=int)
        cellIds =  np.array([x for x in cellIds if x >= 0])
        airSegs = np.array(list(set(np.concatenate((rs.cell2seg[-1],np.where(np.array(rs.organTypes) != 2)[0])) )))#aboveground
        periSegs = np.array([i for i in range(len(rs.organTypes)) if not(i in airSegs) ])
        outer_radii_all = np.array(rs.segOuterRadii(type = 0))
        outer_radii = outer_radii_all[periSegs]
        lengths = np.array(rs.segLength())[periSegs]
        radii = np.array(rs.radii)[periSegs]
        totVolR = sum(np.pi*lengths*(outer_radii**2 - radii**2))
        totVolS = sum( CellVolumes[cellIds])
        assert abs((totVolR-totVolS)/totVolS*100) < 1e-12
        print('totV rhizo and soil', totVolR,totVolS)
        for cid in cellIds:
            periSegs = np.array([i for i in rs.cell2seg[cid] if (np.array(rs.organTypes)[i] == 2) ])
            outer_radii = outer_radii_all[periSegs]
            lengths = np.array(rs.segLength())[periSegs]
            radii = np.array(rs.radii)[periSegs]
            totVolR = sum(np.pi*lengths*(outer_radii**2 - radii**2))
            totVolS =CellVolumes[cid]
            print('cid rhizo and soil', cid,totVolR,totVolS, periSegs, lengths, outer_radii, radii)
        
    comm.barrier()
    rs.update()

    comm.barrier()
    
    weatherX = scenario.weather(rs_age) 
    r.Qlight = weatherX["Qlight"]
    # send soil concentration to plant:
    print('get get_inner_heads')
    rsx = rs.get_inner_heads(weather=weatherX)  # matric potential at the root soil interface, i.e. inner values of the cylindric models (not extrapolation to the interface!) [cm]
    print('get get_inner_solutes')
    r.Csoil_seg = rs.get_inner_solutes() * 1e3 # mol/cm3 to mmol/cm3 
    
    start_time_plant = timeit.default_timer()
    comm.barrier()
    if rank == 0:
        assert min(r.Csoil_seg ) >= 0.
        try:
            r.solve_photosynthesis(sim_time_ = rs_age, 
                        sxx_=rsx, 
                        cells_ = False,#(i == 0),#for 1st computation, use cell data
                        ea_ = weatherX["ea"],#not used
                        es_=weatherX["es"],#not used
                        verbose_ = False, doLog_ = False,
                        TairC_= weatherX["TairC"],#not used
                        outputDir_= "./results/rhizoplantExud")
        except:
            r.solve_photosynthesis(sim_time_ = rs_age, 
                        sxx_=rsx, 
                        cells_ = False,#(i == 0),#for 1st computation, use cell data
                        ea_ = weatherX["ea"],#not used
                        es_=weatherX["es"],#not used
                        verbose_ = True, doLog_ = True,
                        TairC_= weatherX["TairC"],#not used
                        outputDir_= "./results/rhizoplantExud")
            
        errLeuning_abs = abs(sum(r.outputFlux))
        TranspirationCumul += sum(np.array(r.Ev) * dt) #transpiration [cm3/day] * day
        startphloem=rs_age
        endphloem = rs_age + dt
        stepphloem = 1
        verbose_phloem = True
        filename = "results/" +"inPM.txt"
        print("startpm")
        #if i == 1:
        #r.doTroubleshooting = True
        
        r.startPM(startphloem, endphloem, stepphloem, ( weatherX["TairC"]  +273.15) , verbose_phloem, filename)
        
    comm.barrier()
    time_plant_cumul += (timeit.default_timer() - start_time_plant)
    
    
    if rank == 0:
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
                
        Q_Exud_i = np.array( Q_Exud_i[1:] ) #from nod to semgment
        Q_Mucil_i = np.array(Q_Mucil_i[1:])
        
        airSegsId = np.array(np.where(np.array([isinstance(cc, AirSegment) for cc in rs.cyls]))[0])
        try:
            assert (Q_Exud_i[airSegsId] == 0).all()
            assert (Q_Mucil_i[airSegsId] == 0).all()
            assert (np.array(r.k_mucil_)[airSegsId+1] == 0).all()
            assert (np.array(r.Q_Exudmax)[airSegsId+1] == 0).all()
        except:
            print("Q_Exud_i", Q_Exud_i[airSegsId] )
            print("Q_Mucil_i", Q_Mucil_i,Q_Mucil,Q_Mucilbu,airSegsId)
            print("Q_Mucil_i", Q_Mucil_i[airSegsId], Q_Mucil[airSegsId+1], Q_Mucilbu[airSegsId+1])
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
            assert min(Q_Exud_i) >= 0.
        except:
            print(C_ST, r.Csoil_node, Q_Exud_i,Q_Exud)
            raise Exception
        
        print("sum exud", sum(Q_Exud_i), sum(Q_Mucil_i))
    else:
        Q_Exud_i = None
        Q_Mucil_i = None
    print('share plant data')
    Q_Exud_i = comm.bcast(Q_Exud_i, root = 0) 
    Q_Mucil_i = comm.bcast(Q_Mucil_i, root = 0) 
    r.psiXyl = comm.bcast(r.psiXyl, root = 0) 
    r.outputFlux = comm.bcast(r.outputFlux, root = 0) 
    
    dt_inner = dt

    # rs = RhizoMappedSegments(r, wilting_point, nc, logbase, mode)
    if mode == "dumux_dirichlet":
        cyl_type = 1
    elif mode == "dumux_dirichlet_nc":
        cyl_type = 2
    elif mode == "dumux_dirichlet_10c":
        cyl_type = 3
    elif mode == "dumux_10c":
        cyl_type = 4
    else:
        raise("unknown demoType")
    
    Q_Exud_inflate += sum(Q_Exud_i); Q_Mucil_inflate += sum(Q_Mucil_i)
    
    print('to the inner loop')
    psi_x, psi_s, sink, x, y, psi_s2, vol_, surf_,  depth_,soil_c, c, c_All,c_All1, net_sol_flux, net_flux = cyl3.simulate_const(s, 
                                            r,  dt, dt_inner, rs_age, 
                                            demoType = mode, Q_plant=[Q_Exud_i, Q_Mucil_i], plantType = "RS", r= rs,
                                            wilting_point = wilting_point,
                                            outer_R_bc_sol = net_sol_flux, 
                                            outer_R_bc_wat = net_flux)
    
    
    time_rhizo_cumul += r.time_rhizo_i
    time_3ds_cumul += r.time_3ds_i
    
    write_file_array("cellVol", np.array(s.getCellVolumes()), directory_ =results_dir) # cm3 
    write_file_array("theta", np.array(s.getWaterContent()), directory_ =results_dir) 
    for i in range(rs.numFluidComp):
        write_file_array("Soil_solute_conc"+str(i+1), np.array(s.getSolution(i+1)).flatten()* rs.molarDensityWat_m3/1e6, directory_ =results_dir) 
    for i in range(rs.numFluidComp, rs.numComp):
        write_file_array("Soil_solute_conc"+str(i+1), np.array(s.getSolution(i+1)).flatten()* rs.bulkDensity_m3 /1e6 , directory_ =results_dir) 
    for i in range(rs.numComp):
        write_file_array("Soil_old_solute_conc"+str(i+1), rs.soilContent_old[i]/(rs.soilvolumes_old[i]*1e6) , directory_ =results_dir) 
    
    write_file_array("Soil_solute_conc"+str(rs.numComp+1), np.array(s.base.getCSS1_out()).flatten()[:-1]* rs.bulkDensity_m3 /1e6 , directory_ =results_dir) 
    buTotCAfter = sum(s.getTotCContent())   
    buWAfter = sum(np.multiply(np.array(s.getWaterContent()), cell_volumes))    
    
    if rank == 0:
        s.bulkMassErrorCumul_abs = abs((buTotCAfter - ( buTotCSoilInit + Q_Exud_inflate + Q_Mucil_inflate)))#so that only works if infalte
        s.bulkMassErrorWaterCumul_abs = abs(buWAfter - ( buWSoilInit - TranspirationCumul))
        s.bulkMassErrorCumul_rel = abs(s.bulkMassErrorCumul_abs/buTotCAfter*100)
        s.bulkMassErrorWaterCumul_rel = abs(s.bulkMassErrorWaterCumul_abs/buWAfter*100)
    
            
    
    # try:
        # # maybe 0.1% of error is too large? if buTotCAfter >> plant input
        # assert s.bulkMassErrorCumul_rel < 0.1
        # assert s.bulkMassErrorWaterCumul_rel < 0.1
        # # assert check can be useless when Q_exud and Q_Mucil low
    # except:
        # print("buTotCAfter ,  buTotCSoilInit", buTotCAfter ,  buTotCSoilInit)
        # print( "sum(Q_Exud) , sum(Q_Mucil)", sum(Q_Exud) , sum(Q_Mucil))
        # print(s.bulkMassErrorWaterCumul_rel, s.bulkMassErrorCumul_rel ,s.bulkMassErrorWaterCumul_abs, s.bulkMassErrorCumul_abs )
        # raise Exception
            
    write_file_array("totalComputetime",np.array([timeit.default_timer() - start_time_global,
                        time_plant_cumul,time_rhizo_cumul ,time_3ds_cumul]) , directory_ =results_dir)#cumulative
    # write_file_float("partsComputeTime",partsComputeTime)#not cumulative
    write_file_float("time", rs_age, directory_ =results_dir)
    write_file_array("TotSoilC", s.getTotCContent(), directory_ =results_dir)
    write_file_float("Q_Exud_i", sum(Q_Exud_i), directory_ =results_dir)
    write_file_float("Q_Mucil_i", sum(Q_Mucil_i), directory_ =results_dir)
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
        
        write_file_array("sumdiffSoilData_abs", rs.sumdiffSoilData_abs, directory_ =results_dir)# cumulative (?)
        write_file_array("maxdiffSoilData_abs", rs.maxdiffSoilData_abs, directory_ =results_dir)# cumulative (?)
        write_file_array("sumdiffSoilData_rel", rs.sumdiffSoilData_rel, directory_ =results_dir)# cumulative (?)
        write_file_array("maxdiffSoilData_rel", rs.maxdiffSoilData_rel, directory_ =results_dir)# cumulative (?)
        write_file_array("diffSoilData_abs", rs.diffSoilData_abs, directory_ =results_dir)# cumulative (?)
        write_file_array("diffSoilData_rel", rs.diffSoilData_rel, directory_ =results_dir)# cumulative (?)
        
        write_file_array("trans", r.Ev, directory_ =results_dir)
        write_file_array("transrate",r.Jw, directory_ =results_dir)
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
        try:
            assert abs(s.bulkMassErrorPlant_rel)  < 1e-5
        except:
            print("\n\n\nissue bulk soil balance", np.array([s.bulkMassErrorPlant_abs, s.bulkMassErrorPlant_rel, #not cumulative 
                                            s.bulkMassErrorCumul_abs,s.bulkMassErrorCumul_rel,#cumulative
                                            s.bulkMassErrorWater_abs,s.bulkMassErrorWater_rel, #not cumulative
                                            s.bulkMassErrorWaterCumul_abs,s.bulkMassErrorWaterCumul_rel]))
            print("\n\n\n")
            raise Exception
        if True :
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
            
            
        psi_x_.extend(psi_x) #[cm]
        psi_s_.extend(psi_s) #[cm]
        sink_.extend(sink)
        x = np.array(x)
        x_.extend(x)
        y_.extend(y)
        psi_s2_.extend(psi_s2) #[cm]
        soil_c_.extend(soil_c)  # [g/cm3]
        c_.extend(c)  # [g/day]
        # mass_soil_c_.extend(mass_soil_c) #[g]
    
        dist_, conc_, l_ = rs.collect_cylinder_solute_data()
        dist.append(dist_)
        conc.append(conc_)
        l.append(l_)
        

""" output """
sizeSoilCell = rs.soilModel.getCellVolumes()
rs.checkMassOMoleBalance2( sourceWat = np.full(len(sizeSoilCell),0.), # cm3/day 
                                 sourceSol = np.full((rs.numComp, len(sizeSoilCell)),0.), # mol/day
                                 dt = 0.,        # day    
                                 seg_fluxes = 0.,# [cm3/day]
                                 doWater = True, doSolute = True, doSolid = True,
                                 useSoilData = True)
vp.write_soil("results/"+str(sim_time)+"_Soil", s, min_b, max_b, cell_number, ["C concentration [g/cm³]"])



xx = np.array(s.getSolutionHead())
for konz in np.array([True]):#, False]):
    if not konz:
        extraArrayName_ = "theta (cm3)"
        extraArray_ = rs.soilWatVol_old*1e6
    else:
        extraArrayName_ = "theta (cm3/cm3)"
        extraArray_ = rs.soilTheta_old
    print("idcomp", 0, rank)
    rp = rs.get_concentration(0, konz)
    
    vp.plot_roots_and_soil(rs.mappedSegments(),extraArrayName_,rp, s, periodic, min_b, max_b, cell_number, 
            soil_type+genotype+"_rx", sol_ind =-1,extraArray = extraArray_, extraArrayName = extraArrayName_)  # VTK vizualisation, rs.soilTheta_old
    print("idcomp_done", 0, rank)
    for i in range(1, rs.numComp+1):
        print("idcomp", i, rank)
        if not konz:
            extraArrayName_ = "C"+str(i)+" mol"
            extraArray_ = rs.soilContent_old[i -1]
        else:
            extraArrayName_ = "[C"+str(i)+"] (mol/cm3)"
            extraArray_ = rs.soilContent_old[i -1]/(rs.soilvolumes_old[i-1]*1e6)
            
        vp.plot_roots_and_soil(rs.mappedSegments(),extraArrayName_ ,rs.get_concentration(i , konz), s, periodic, min_b, max_b, cell_number, 
                soil_type+genotype+"_rx", sol_ind =-1,extraArray = extraArray_, 
                extraArrayName = extraArrayName_)  # VTK vizualisation
        print("idcomp_done", i, rank)
# to plot: cumulative An, Exud, mucil and each element in the soil
# Also images of the 3D soil.
print("fin", rank)
