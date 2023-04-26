""" 
    Create mapped maize root system  
"""
import sys;
sys.path.append("/data/");
sys.path.append("../modules/");
sys.path.append("../../../CPlantBox/");
sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../../CPlantBox/src/functional/");
sys.path.append("../../../CPlantBox/src/rsml/");
sys.path.append("../../../CPlantBox/src/visualisation/")
sys.path.append("../../../CPlantBox/src/structural/")
sys.path.append("../../../CPlantBox/src/external/")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver
import plantbox as pb  # CPlantBox
import vtk_plot as vp
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import timeit
import numpy as np

import scenario_setup as scenario
import evapotranspiration as evap

""" parameters   """
soil_, min_b, max_b, cell_number, area, Kc = scenario.maize_SPP()

sim_time = 10 #95  #  [day]
dt = 360 / (24 * 3600)  # time step [day] 20

start_date = '2021-05-10 00:00:00'  # INARI csv data
x_, y_ = evap.net_infiltration_table_beers_csvS(start_date, sim_time, evap.lai_maize2, Kc)

s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 2, times = x_, net_inf = y_, wet = False)

kexu = scenario.exudation_rates()

xml_name = "data/Zeamays_synMRI_modified.xml"  # root growth model parameter file
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  

c_ = []
times = []
tipnum  = []
for i in range(0, int(sim_time)):

    print("Day", i)

    r.rs.simulate(1.)  # simulate for 1 day
    rs_age = i + 1
    t = i * dt  # current simulation time

    seg_sol_fluxes = np.array(r.exudate_fluxes(rs_age+1, kexu))  # [g/day]
    seg_sol_fluxes = comm.bcast(seg_sol_fluxes, root = 0)
    c_.append(np.sum(seg_sol_fluxes))  # [cm3/day]
    times.append(rs_age + t)  # day
    
    ana = pb.SegmentAnalyser(r.rs)
    #ana.write("results/test_RS/RootSys_day"+str(i)+".vtp")
    
    tipnum.append(len(r.rs.getRootTips()))

np.save('results/test_RS/carbon',  np.vstack((times, np.array(c_), np.array(tipnum))))  # exudation per segment  [g]
