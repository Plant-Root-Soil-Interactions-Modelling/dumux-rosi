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
from scenario_setup import write_file_array, write_file_float, div0, div0f

results_dir="./results/TryGrowthissueSegId/"
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

periodic = False
if recreateComsol:
    nc = 500  # dof+1
else:
    nc = 10
    
logbase = 0.5  # according to Mai et al. (2019)
mode = "dumux_10c"  

""" initialize """

initsim = 11.5
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


dt = 1/24/2
simMax = 15
rs_age = initsim
while rs_age < simMax: #for i, dt in enumerate(np.diff(times)):


    rs_age += dt
    print("Day", rs_age)
    seg2cell_old = rs.seg2cell
    rs.simulate(dt)  # simulate for 1 day
    seg2cell_new = rs.seg2cell
    
    # make sure that, once a segment is created, it stays in the same soil voxel
    for segid in seg2cell_old.keys():
        try:
            assert seg2cell_old[segid] == seg2cell_new[segid]
        except:
            print('old',seg2cell_old[segid], segid, rs.organTypes[segid])
            print('new',seg2cell_new[segid], segid, rs.organTypes[segid])
            raise Exception

print("fin")
