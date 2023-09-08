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
from scenario_setup import write_file_array

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
initsim = 1
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
        write_file_array("solute_conc"+str(i+1), np.array(cyl.getSolution_(i+1)).flatten()* rs.molarDensityWat ) 
    for i in range(rs.numFluidComp, rs.numComp):
        write_file_array("solute_conc"+str(i+1), np.array(cyl.getSolution_(i+1)).flatten()* rs.bulkDensity_m3 /1e6 ) 

psi_x_, psi_s_, sink_, x_, y_, psi_s2_, soil_c_, c_, mass_soil_c_, dist, conc, l, a = [], [], [], [], [], [], [], [], [], [], [], [], []

times = [0., 1., 2.]  #h to  days
rs_age = initsim


net_sol_flux = np.array([])
net_flux = np.array([])


for i, dt in enumerate(np.diff(times)):


    rs_age += dt
    print("Day", rs_age)
    rs.simulate(dt)  # simulate for 1 day
    
    
    if i ==0: #for the initial carbon flow rate
        weatherX = scenario.weather(initsim) 
        r.minLoop = 1000
        r.maxLoop = 5000
        r.Qlight = weatherX["Qlight"]
        r.solve_photosynthesis(sim_time_ = initsim, 
                    sxx_=x, 
                    cells_ = True,#for 1st computation, use cell data
                    ea_ = weatherX["ea"],#not used
                    es_=weatherX["es"],#not used
                    verbose_ = False, doLog_ = False,
                    TairC_= weatherX["TairC"],#not used
                    outputDir_= "./results/rhizoplantExud")
        
    weatherX = scenario.weather(rs_age) 
    startphloem=rs_age
    endphloem = rs_age + dt
    stepphloem = 1
    verbose_phloem = True
    filename = "results/" +"inPM_"+str(i)+".txt"
    print("startpm")
    r.startPM(startphloem, endphloem, stepphloem, ( weatherX["TairC"]  +273.15) , verbose_phloem, filename)
    Nt = len(rs.nodes)
    QExud  = np.array(r.Q_out[(Nt*3):(Nt*4)])#mol/day for nodes
    assert QExud[0] == 0#no exudation in seed node I guess
    QExud = QExud[1:] #from nod to semgment
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
    
    psi_x, psi_s, sink, x, y, psi_s2, vol_, surf_,  depth_,soil_c, c,repartition, c_All, net_sol_flux, net_flux = cyl3.simulate_const(s, 
                                            r,  dt, dt_inner, kexu, rs_age, repartition, 
                                            type = mode, Q_Exud=QExud, plantType = "RS", r= rs,
                                            wilting_point = wilting_point, trans_maize = trans_maize,
                                            net_sol_flux = net_sol_flux, net_flux = net_flux)
    # raise Exception
    #(fname, s, r, sri_table_lookup, 1., dt, trans_maize, comp, rs_age, min_b, max_b, type = cyl_type)

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
    
    if True:
        periodic = False
        # print(c_All)
        xx = np.array(s.getSolutionHead())
        cc = rs.getCC(len(xx))
        print(xx.shape,psi_x[-1].shape,  psi_s[-1].shape, c_All[-1].shape,len(r.get_segments()) , QExud.shape)
        vp.plot_roots_and_soil(rs.mappedSegments(), "pressure head (rx)", psi_x[-1], s, periodic, min_b, max_b, cell_number, 
                soil_type+genotype+"_rx", sol_ind =0) 
        # not sure that works
        vp.plot_roots_and_soil(rs.mappedSegments(), "sucrose",c_All[-1], s, periodic, min_b, max_b, cell_number, 
                soil_type+genotype+"_rx", sol_ind =-1,extraArray = cc)  # VTK vizualisation
        vp.plot_roots_and_soil(rs.mappedSegments(), "sucrose",c_All[-1], s, periodic, min_b, max_b, cell_number, 
                soil_type+genotype+"_rx", sol_ind =1)  # VTK vizualisation
        for i in range(1, rs.numComp):
            print("idcomp", i)
            vp.plot_roots_and_soil(rs.mappedSegments(), "sucrose",c_All[-1], s, periodic, min_b, max_b, cell_number, 
                    soil_type+genotype+"_rx", sol_ind =-1,extraArray = rs.getCC(len(xx), idComp=i))  # VTK vizualisation
        
    print("fin")
