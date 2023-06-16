""" 
    Maize using rhizosphere models  
"""
import sys;
sys.path.append("/data/");
sys.path.append("../modules/");
sys.path.append("../../../CPlantBox/");
sys.path.append("../../../CPlantBox/src")
#sys.path.append("../../../CPlantBox/src/functional/");
#sys.path.append("../../../CPlantBox/src/rsml/");
#sys.path.append("../../../CPlantBox/src/visualisation/")
#sys.path.append("../../../CPlantBox/src/structural/")
#sys.path.append("../../../CPlantBox/src/external/")
#sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver
import plantbox as pb  # CPlantBox
import visualisation.vtk_plot as vp
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import timeit
import numpy as np

import scenario_setup as scenario
from functional.phloem_flux import PhloemFluxPython  # root system Python hybrid solver
from rhizo_modelsPlant import *  # Helper class for cylindrical rhizosphere models
import evapotranspiration as evap
import cyl3
import cyl_exu
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

"""
 Cylindric rhizosphere models, C exudate from finite volumes
"""

"""scenario"""
year = 2013
soil_type = "loam"
genotype = "WT"

""" parameters   """
soil_, min_b, max_b, cell_number, area, Kc = scenario.maize_SPP(soil_type)
sim_time = 1 #154 #155  #  [day]
dt = 360 / (24 * 3600)  # time step [day] 20

x_, y_, lai = evap.net_infiltration(year, soil_type, genotype, sim_time, Kc)
trans_maize = evap.get_transpiration(year, sim_time, area, lai, Kc)

""" rhizosphere model parameters """
wilting_point = -15000  # cm
nc = 10  # dof+1
logbase = 0.5  # according to Mai et al. (2019)

mode = "dumux_dirichlet_nc"  # mode = "dumux_dirichlet_nc"

""" initialize """
start_time = timeit.default_timer()
initsim = 15
s, soil = scenario.create_soil_model(soil_type, year, soil_, min_b, max_b, 
                                        cell_number, type = 2, times = x_, net_inf = y_)
water0 = s.getWaterVolume()  # total initial water volume in domain
path = "../../../CPlantBox/modelparameter/structural/plant/"
xml_name = "Triticum_aestivum_test_2021.xml"  # root growth model parameter file
rs = scenario.create_mapped_plant(wilting_point, nc, logbase, mode,initsim,
                                        min_b, max_b, cell_number, s, xml_name, path)  # pass parameter file for dynamic growth

r = rs.rs #phloemflux object
kexu = scenario.exudation_rates()


if rank == 0:
    x = s.getSolutionHead_()  # initial condition of soil [cm]
    cc = s.getSolution_(1)  # solute concentration [kg/m3].
    #x = s.getSolutionHead()[:, 0]  # "[:, 0]"  == .flatten()
x = comm.bcast(x, root = 0)  # Soil part should/will run in parallel

ns = len(rs.segments)
dcyl = int(np.floor(ns / max_rank))
repartition = np.array([dcyl for i in range(max_rank)])
repartition[max_rank -1] = ns - rank * dcyl
assert np.sum(repartition) == ns
if rank == 0:
    rs.initializeRhizo(soil_, x, np.array([i for i in range(repartition[rank])]), cc)
    print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, 0,repartition[rank], timeit.default_timer() - start_time))
else:
    rs.initializeRhizo(soil_, x, np.array([i for i in range(repartition[rank-1],repartition[rank])]), cc)
    print ("Initialized rank {:g}/{:g} [{:g}-{:g}] in {:g} s".format(rank + 1, max_rank, repartition[rank-1],repartition[rank], timeit.default_timer() - start_time))


psi_x_, psi_s_, sink_, x_, y_, psi_s2_, soil_c_, c_, dist, conc, l, a = [], [], [], [], [], [], [], [], [], [], [], []




    
ddt = 1
for i in range(0, int(sim_time)):
    rs.simulate(ddt)  # simulate for 1 day
    
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
    
    weatherX = scenario.weather(initsim+i*ddt) 
    startphloem= initsim+(i)*ddt
    endphloem = startphloem + (i+1)*ddt
    stepphloem = 1
    verbose_phloem = True
    filename = "results/" +"inPM_"+str(i)+".txt"
    print("startpm")
    r.startPM(startphloem, endphloem, stepphloem, ( weatherX["TairC"]  +273.15) , verbose_phloem, filename)
    Nt = len(rs.nodes)
    QExud  = np.array(r.Q_out[(Nt*3):(Nt*4)])#mol/day for nodes
    
    print("Day", i)
    
    if rank == 0:
        r.test()  # sanity checks
    
    rs_age = i + 1
    
    if mode == "dumux_dirichlet":
        cyl_type = 1
    elif mode == "dumux_dirichlet_nc":
        cyl_type = 2
    else:
        raise("unknown type")
    #krs_,
    psi_x, psi_s, sink, x, y, psi_s2, vol_, surf_,  depth_,soil_c, c,repartition, c_All = cyl3.simulate_const(s, rs, 1., dt,  kexu, rs_age, repartition, type = cyl_type, Q_Exud=QExud)
    
    if rank == 0:  # collect results
        psi_x_.extend(psi_x)
        psi_s_.extend(psi_s)
        sink_.extend(sink)
        x = np.array(x)
        x_.extend(x)
        y_.extend(y)
        psi_s2_.extend(psi_s2)
        soil_c_.extend(soil_c)  # [kg/m3]
        c_.extend(c)  # [cm3/day]
        vp.write_soil("results/vtu_vtp/Soil_day"+str(i), s, min_b, max_b, cell_number,["C concentration [g/cm³]"])
        print('vtu written')

        ana = pb.SegmentAnalyser(r.rs)
        ana.write("results/vtu_vtp/RootSys_day"+str(i)+".vtp")
        print('vtp written')
        dist_, conc_, l_ = rs.collect_cylinder_solute_data()
        dist.append(dist_)
        conc.append(conc_)
        l.append(l_)
        

water = s.getWaterVolume()

""" output """
if rank == 0:

    #ana = pb.SegmentAnalyser(rs.mappedSegments())
    #ana.write("results/"+str(sim_time)+"_RootSys.vtp")
    #
    vp.write_soil("results/"+str(sim_time)+"_Soil", s, min_b, max_b, cell_number, ["C concentration [g/cm³]"])
    
    #rs.plot_cylinders()
    #rs.plot_cylinders_solute()

    scenario.write_files("maize_cyl3.txt", psi_x_, psi_s_, sink_, x_, y_, psi_s2_,  vol_, surf_,  depth_,  dist, conc, l, soil_c_, c_)
    print ("Overall simulation wall time", timeit.default_timer() - start_time, " s")
    
    if True:
        periodic = False
        print(c_All)
        print(psi_s[-1].shape, c_All[-1].shape)
        vp.plot_roots_and_soil(rs.mappedSegments(), "sucrose",c_All[-1], s, periodic, min_b, max_b, cell_number, soil_type+genotype+"_rx", sol_ind =-1) 
        vp.plot_roots_and_soil(rs.mappedSegments(), "sucrose",c_All[-1], s, periodic, min_b, max_b, cell_number, soil_type+genotype+"_rx", sol_ind =1)  # VTK vizualisation
        #vp.plot_roots_and_soil(rs.mappedSegments(), "pressure head", psi_s[-1], s, periodic, min_b, max_b, cell_number, soil_type+genotype+"_rsx")  # VTK vizualisation
        #vp.plot_roots_and_soil(rs.mappedSegments(), "pressure head (rx)", psi_x[-1], s, periodic, min_b, max_b, cell_number, soil_type+genotype+"_rx")  # VTK vizualisation
        #vp.plot_roots_and_soil(rs.mappedSegments(), "realized_inner_fluxes", realized_inner_fluxes, s, periodic, min_b, max_b, cell_number, name+"_flux")  # VTK vizualisation

    print("fin")
