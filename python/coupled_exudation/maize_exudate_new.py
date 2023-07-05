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
from rhizo_models import RhizoMappedSegments
import evapotranspiration as evap
import cyl_exu

"""
 Cylindric rhizosphere models, C exudate from finite volumes
"""

fname = 'L_WT_phenolics'

"""scenario"""
year = 2019
soil_type = "loam"
genotype = "WT"
comp = "phenolics"

""" parameters   """
soil_, min_b, max_b, cell_number, area, Kc = scenario.maize_SPP(soil_type)
sim_time = 20 #154   #  [day]
dt = 360 / (24 * 3600)  # time step [day] 20

x_, y_, lai = evap.net_infiltration(year, soil_type, genotype, sim_time, Kc)
trans_maize = evap.get_transpiration(year, sim_time, area, lai, Kc)

""" rhizosphere model parameters """
wilting_point = -15000  # cm
nc = 10  # dof+1
logbase = 0.5  # according to Mai et al. (2019)
mode = "dumux_dirichlet_nc"  

""" initialize """
start_time = timeit.default_timer()
s, soil = scenario.create_soil_model(soil_type, year, soil_,comp, min_b, max_b, cell_number, type = 2, times = x_, net_inf = y_)
sri_table_lookup = cyl_exu.open_sri_lookup("data_magda/" + soil_type)
water0 = s.getWaterVolume()  # total initial water volume in domain

xml_name = "data_magda/Zea_mays_5_Leitner_optimized.xml"  # root growth model parameter file
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
scenario.init_maize_conductivities(r)

#vp.write_soil("test_vtu", s, min_b, max_b, cell_number, ["C concentration [g/cm³]"])

if rank == 0:
    r.test()  # sanity checks

psi_x_, psi_s_, sink_, x_, y_, psi_s2_, soil_c_, c_, mass_soil_c_, dist, conc, l, a = [], [], [], [], [], [], [], [], [], [], [], [], []

for i in range(0, int(sim_time)):

    print("Day", i)

    r.rs.simulate(1.)  # simulate for 1 day
    rs_age = i + 1

    rs = RhizoMappedSegments(r, wilting_point, nc, logbase, mode)
    if mode == "dumux_dirichlet":
        cyl_type = 1
    elif mode == "dumux_dirichlet_nc":
        cyl_type = 2
    else:
        raise("unknown type")

    psi_x, psi_s, sink, x, y, psi_s2, vol_, surf_, krs_, depth_,soil_c, c, mass_soil_c = cyl_exu.simulate_const(fname, s, rs, sri_table_lookup, 1., dt, trans_maize, comp, rs_age, min_b, max_b, type = cyl_type)

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
        mass_soil_c_.extend(mass_soil_c) #[g]
        if i%2==0: 
            vp.write_soil("results/vtu_vtp/Soil_day"+str(i), s, min_b, max_b, cell_number, ["C concentration [g/cm³]"])
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

    #ana = pb.SegmentAnalyser(r.rs)
    #ana.write("results/"+str(sim_time)+"_RootSys.vtp")
    #vp.write_soil("results/"+str(sim_time)+"_Soil", s, min_b, max_b, cell_number, ["C concentration [g/cm³]"])
    
    #rs.plot_cylinders()
    #rs.plot_cylinders_solute()
    psi_x_ = np.array(psi_x_, dtype=object)
    psi_s_ = np.array(psi_s_, dtype=object)
    sink_ = np.array(sink_, dtype=object)
    x_ = np.array(x_, dtype=object)
    y_ = np.array(y_, dtype=object)
    psi_s2_ = np.array(psi_s2_, dtype=object)
    vol_ = np.array(vol_, dtype=object)
    surf_ = np.array(surf_, dtype=object)
    krs_ = np.array(krs_, dtype=object)
    dist = np.array(dist, dtype=object)
    conc = np.array(conc, dtype=object)
    l = np.array(l, dtype=object)
    soil_c_ = np.array(soil_c_, dtype=object)
    c_ = np.array(c_, dtype=object)
    mass_soil_c_ = np.array(mass_soil_c_, dtype=object)
    

    scenario.write_files(fname, psi_x_, psi_s_, sink_, x_, y_, psi_s2_,  vol_, surf_, krs_, depth_,  dist, conc, l, soil_c_, c_, mass_soil_c_)
    print ("Overall simulation wall time", timeit.default_timer() - start_time, " s")
    print("fin")
