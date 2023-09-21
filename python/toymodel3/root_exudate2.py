""" 
    Maize using rhizosphere models  
"""
import sys;
sys.path.append("/data/");
sys.path.append("../modules/");
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
from scenario_setup import write_file_array, write_file_float

results_dir="./results/"
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

x_, y_, lai = evap.net_infiltration(year, soil_type, genotype, sim_time, Kc)
trans_maize = evap.get_transpiration(year, sim_time, area, lai, Kc)


""" rhizosphere model parameters """
wilting_point = -15000  # cm
recreateComsol = False
if recreateComsol:
    nc = 500  # dof+1
else:
    nc = 10
    
logbase = 0.5  # according to Mai et al. (2019)
mode = "dumux_10c"  

""" initialize """
start_time = timeit.default_timer()
initsim = 9.5
s, soil = scenario.create_soil_model(soil_type, year, soil_,#comp, 
            min_b, max_b, cell_number, type = mode, times = x_, net_inf = y_,
            usemoles = usemoles)
# sri_table_lookup = cyl_exu.open_sri_lookup("data_magda/" + soil_type)
water0 = s.getWaterVolume()  # total initial water volume in domain

path = "../../../CPlantBox/modelparameter/structural/plant/"
xml_name = "Triticum_aestivum_test_2021.xml"  # root growth model parameter file
rs, r = scenario.create_mapped_plant(wilting_point, nc, logbase, mode,initsim,
                                        min_b, max_b, cell_number, s, xml_name,
                                        path, plantType = "plant", 
                                        recreateComsol_ = recreateComsol,
                                        usemoles = usemoles)  # pass parameter file for dynamic growth



if rank == 0:
    x = s.getSolutionHead_()  # initial condition of soil [cm]
    #cc = s.getSolution_(1)  # solute concentration [kg/m3].
    #x = s.getSolutionHead()[:, 0]  # "[:, 0]"  == .flatten()
x = comm.bcast(x, root = 0)  # Soil part should/will run in parallel
kexu = scenario.exudation_rates()

ns = len(r.get_segments())
dcyl = 0#int(np.floor(ns / max_rank))
repartition = np.array([dcyl for i in range(max_rank)])
repartition[max_rank -1] = 0#ns - rank * dcyl
#assert np.sum(repartition) == ns

if rank == 0:
    rs.initializeRhizo(soil_, x, np.array([i for i in range(repartition[rank])]))#, cc)
    print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, 0,repartition[rank], timeit.default_timer() - start_time))
else:
    rs.initializeRhizo(soil_, x, np.array([i for i in range(repartition[rank-1],repartition[rank])]))#, cc)
    print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, repartition[rank-1],repartition[rank], timeit.default_timer() - start_time))

if len(rs.cyls) > 1:
    cyl = rs.cyls[1] # take a random cylinder

    write_file_array("pressureHead",np.array(cyl.getSolutionHead()).flatten())
    write_file_array("coord", cyl.getDofCoordinates().flatten())
    for i in range(rs.numFluidComp):
        write_file_array("solute_conc"+str(i+1), np.array(cyl.getSolution_(i+1)).flatten()* rs.molarDensityWat_m3/1e6 ) 
    for i in range(rs.numFluidComp, rs.numComp):
        write_file_array("solute_conc"+str(i+1), np.array(cyl.getSolution_(i+1)).flatten()* rs.bulkDensity_m3 /1e6 ) 

psi_x_, psi_s_, sink_, x_, y_, psi_s2_, soil_c_, c_, mass_soil_c_, dist, conc, l, a = [], [], [], [], [], [], [], [], [], [], [], [], []

times = [0., 1./24., 2./24., 3./24.]  #h to  days
rs_age = initsim


net_sol_flux =  np.array([np.array([]),np.array([])])
net_flux = np.array([])


Q_ST_init = np.array([])
Q_Exudbu    = np.array([0.])
Q_Mucilbu   = np.array([0.])
Q_in  = 0
Nt = len(rs.nodes)
r.minLoop = 1000
r.maxLoop = 5000
dt = 1/24/2
simMax = initsim + 1
while rs_age < simMax: #for i, dt in enumerate(np.diff(times)):


    rs_age += dt
    print("Day", rs_age)
    rs.simulate(dt)  # simulate for 1 day
    
    
    start_time = timeit.default_timer()
    repartitionOld = repartition
    nsOld = sum(repartitionOld)
    ns = len(r.get_segments())
    dcyl = int(np.floor(ns / max_rank))
    repartition = np.array([dcyl for i in range(max_rank)])
    repartition[max_rank -1] = ns - rank * dcyl
    toAdd = repartition - repartitionOld
    
    if toAdd[rank] > 0: # that could be put in main file (after RS growth)

        if rank == 0:
            rs.update(np.array([i for i in range(toAdd[rank])]) + nsOld)
            print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, nsOld, toAdd[rank] + nsOld, timeit.default_timer() - start_time))
        else:
            rs.update(np.array([i for i in range(toAdd[rank-1],toAdd[rank])]) + nsOld)
            print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, toAdd[rank-1] + nsOld, toAdd[rank] + nsOld, timeit.default_timer() - start_time))
            
    
    weatherX = scenario.weather(rs_age) 
    r.Qlight = weatherX["Qlight"]
    rsx = rs.get_inner_heads(weather=weatherX)  # matric potential at the root soil interface, i.e. inner values of the cylindric models (not extrapolation to the interface!) [cm]
    try:
        r.solve_photosynthesis(sim_time_ = initsim, 
                    sxx_=rsx, 
                    cells_ = False,#(i == 0),#for 1st computation, use cell data
                    ea_ = weatherX["ea"],#not used
                    es_=weatherX["es"],#not used
                    verbose_ = False, doLog_ = False,
                    TairC_= weatherX["TairC"],#not used
                    outputDir_= "./results/rhizoplantExud")
    except:
        r.solve_photosynthesis(sim_time_ = initsim, 
                    sxx_=rsx, 
                    cells_ = False,#(i == 0),#for 1st computation, use cell data
                    ea_ = weatherX["ea"],#not used
                    es_=weatherX["es"],#not used
                    verbose_ = True, doLog_ = True,
                    TairC_= weatherX["TairC"],#not used
                    outputDir_= "./results/rhizoplantExud")
        
        
    startphloem=rs_age
    endphloem = rs_age + dt
    stepphloem = 1
    verbose_phloem = True
    filename = "results/" +"inPM.txt"
    print("startpm")
    #if i == 1:
    #    r.doTroubleshooting = True
    r.startPM(startphloem, endphloem, stepphloem, ( weatherX["TairC"]  +273.15) , verbose_phloem, filename)
    Nt = len(rs.nodes)
    if r.withInitVal and (len(Q_ST_init) ==0) :
        Q_ST_init = np.array(r.Q_init[0:Nt])
        Q_meso_init = np.array(r.Q_init[Nt:(Nt*2)])
    # att: that will be the cumulative value
    Q_ST    = np.array(r.Q_out[0:Nt])
    Q_meso  = np.array(r.Q_out[Nt:(Nt*2)])
    Q_Rm    = np.array(r.Q_out[(Nt*2):(Nt*3)])
    Q_Exud  = np.array(r.Q_out[(Nt*3):(Nt*4)])
    Q_Gr    = np.array(r.Q_out[(Nt*4):(Nt*5)])
    Q_Mucil  = np.array(r.Q_out[(Nt*9):(Nt*10)])#mmol for nodes

    C_ST    = np.array(r.C_ST)
    Fl      = np.array(r.Fl)
    volST   = np.array(r.vol_ST)
    volMeso   = np.array(r.vol_Meso)
    C_meso  = Q_meso/volMeso
    Q_in   += sum(np.array(r.AgPhl)*dt)
    Q_out   = Q_Rm + Q_Exud + Q_Gr
    error   = sum(Q_ST + Q_meso + Q_out )- Q_in - sum(Q_ST_init)  - sum(Q_meso_init)
    
    Ntbu = Nt
    Q_Exudbu     =   np.concatenate((Q_Exud, np.full(Nt - Ntbu, 0.))) 
    Q_Mucilbu       =   np.concatenate((Q_Mucilbu, np.full(Nt - Ntbu, 0.))) 
    Q_Exud_i      = Q_Exud    - Q_Exudbu
    Q_Mucil_i     = Q_Mucil   - Q_Mucilbu
    assert Q_Exud_i[0] == 0#no exudation in seed node I guess
    assert Q_Mucil_i[0] == 0#no exudation in seed node I guess
    #Q_Exud_i =np.full( Q_Exud_i[1:],1. )#np.array( Q_Exud_i[1:] )*1e3 #from nod to semgment, also fom mmol to mol
    Q_Mucil_i = np.array(Q_Mucil_i[1:])*1e3 
    Q_Exud_i = Q_Mucil_i
    Q_Exud_i[np.where(Q_Exud_i > 0.)] = 1.
    #r.outputFlux = np.array(r.outputFlux)/ 10
    assert sum(Q_Exud_i) > 0.
    
    print("sum exud", sum(Q_Exud_i), sum(Q_Mucil_i))
    
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
        raise("unknown type")
    
    psi_x, psi_s, sink, x, y, psi_s2, vol_, surf_,  depth_,soil_c, c,repartition, c_All,c_All1, net_sol_flux, net_flux = cyl3.simulate_const(s, 
                                            r,  dt, dt_inner, kexu, rs_age, repartition, 
                                            type = mode, Q_plant=[Q_Exud_i, Q_Mucil_i], plantType = "RS", r= rs,
                                            wilting_point = wilting_point, trans_maize = trans_maize,
                                            outer_R_bc_sol = net_sol_flux, 
                                            outer_R_bc_wat = net_flux)
    # raise Exception
    #(fname, s, r, sri_table_lookup, 1., dt, trans_maize, comp, rs_age, min_b, max_b, type = cyl_type)

    for i in range(rs.numFluidComp):
        write_file_array("Soil_solute_conc"+str(i+1), np.array(s.getSolution_(i+1)).flatten()* rs.molarDensityWat_m3/1e6) 
    for i in range(rs.numFluidComp, rs.numComp):
        write_file_array("Soil_solute_conc"+str(i+1), np.array(s.getSolution_(i+1)).flatten()* rs.bulkDensity_m3 /1e6 ) 
    for i in range(rs.numComp):
        write_file_array("Soil_old_solute_conc"+str(i+1), rs.soilContent_old[i]/(rs.soilvolumes_old[i]*1e6) ) 
    
    write_file_float("Q_Exud_i", sum(Q_Exud_i))
    write_file_float("Q_Mucil_i", sum(Q_Mucil_i))
    write_file_array("trans", r.Ev)
    write_file_array("transrate",r.Jw)
    write_file_array("Q_ST", Q_ST)
    write_file_array("C_ST", C_ST)
    write_file_array("C_meso", C_meso)
    write_file_array("Q_meso", Q_meso)
    write_file_array("Q_Rm", Q_Rm)
    write_file_array("Q_Exud", Q_Exud)
    write_file_array("Q_Gr", Q_Gr)
    write_file_array("psiXyl", r.psiXyl)
    write_file_array("Fpsi", r.Fpsi)

    if rank == 0:  # collect results
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
        
        if i%1==0: 
            vp.write_soil("results/vtu_vtp/Soil_day"+str(i), s, min_b, max_b, cell_number, ["C concentration [g/cm³]"])
            print('vtu written')
            ana = pb.SegmentAnalyser(r.rs)
            ana.write("results/vtu_vtp/RootSys_day"+str(i)+".vtp")
            print('vtp written')
            scenario.write_files("staticRoot_cyl3", psi_x[0], psi_s[0], sink[0], x[0], y[0], psi_s2[0],  
                                vol_[0], surf_[0],  depth_,  dist_[0], conc_[0], l_[0], soil_c[0], c[0])
            
            cyl = rs.cyls[1] # take a random cylinder (not 0 to not have air
            

""" output """
if rank == 0:

    #ana = pb.SegmentAnalyser(rs.mappedSegments())
    #ana.write("results/"+str(sim_time)+"_RootSys.vtp")
    #
    vp.write_soil("results/"+str(sim_time)+"_Soil", s, min_b, max_b, cell_number, ["C concentration [g/cm³]"])
    
    #rs.plot_cylinders()
    #rs.plot_cylinders_solute()

    print ("Overall simulation wall time", timeit.default_timer() - start_time, " s")
    sizeSoilCell = rs.soilModel.getCellVolumes_()
    rs.checkMassOMoleBalance2( sourceWat = np.full(len(sizeSoilCell),0.), # cm3/day 
                                     sourceSol = np.full((rs.numComp, len(sizeSoilCell)),0.), # mol/day
                                     dt = 0.,        # day    
                                     seg_fluxes = 0.,# [cm3/day]
                                     doWater = True, doSolute = True, doSolid = True,
                                     useSoilData = True)
    if True:
        periodic = False
        # print(c_All)
        xx = np.array(s.getSolutionHead())
        # #cc = rs.getCC(len(xx))
        # #print(xx.shape,psi_x[-1].shape,  psi_s[-1].shape, c_All[-1].shape,len(r.get_segments()) , Q_Exud.shape)
        # vp.plot_roots_and_soil(rs.mappedSegments(), "pressure head (rx)", psi_x[-1], s, periodic, min_b, max_b, cell_number, 
                # soil_type+genotype+"_rx", sol_ind =0) 
        # # not sure that works
        # #vp.plot_roots_and_soil(rs.mappedSegments(), "sucrose",c_All[-1], s, periodic, min_b, max_b, cell_number, 
        # #        soil_type+genotype+"_rx", sol_ind =-1,extraArray = cc)  # VTK vizualisation
        # vp.plot_roots_and_soil(rs.mappedSegments(), "sucrose",c_All[-1], s, periodic, min_b, max_b, cell_number, 
                # soil_type+genotype+"_rx", sol_ind =1)  # VTK vizualisation
        # c_All_ = [c_All[-1] for i in range(rs.numComp)]
        # c_All_[1] = c_All1
        
        for konz in np.array([True, False]):
            if not konz:
                extraArrayName_ = "theta (cm3)"
                extraArray_ = rs.soilWatVol_old*1e6
            else:
                extraArrayName_ = "theta (cm3/cm3)"
                extraArray_ = rs.soilTheta_old
            print("idcomp", 0)
            vp.plot_roots_and_soil(rs.mappedSegments(),extraArrayName_,rs.get_concentration(0, konz), s, periodic, min_b, max_b, cell_number, 
                    soil_type+genotype+"_rx", sol_ind =-1,extraArray = extraArray_, extraArrayName = extraArrayName_)  # VTK vizualisation, rs.soilTheta_old
            for i in range(1, rs.numComp+1):
                print("idcomp", i)
                if not konz:
                    extraArrayName_ = "C"+str(i)+" mol"
                    extraArray_ = rs.soilContent_old[i -1]
                else:
                    extraArrayName_ = "[C"+str(i)+"] (mol/cm3)"
                    extraArray_ = rs.soilContent_old[i -1]/(rs.soilvolumes_old[i-1]*1e6)
                    
                vp.plot_roots_and_soil(rs.mappedSegments(),extraArrayName_ ,rs.get_concentration(i , konz), s, periodic, min_b, max_b, cell_number, 
                        soil_type+genotype+"_rx", sol_ind =-1,extraArray = extraArray_, 
                        extraArrayName = extraArrayName_)  # VTK vizualisation
        # to plot: cumulative An, Exud, mucil and each element in the soil
    # Also images of the 3D soil.
    print("fin")
