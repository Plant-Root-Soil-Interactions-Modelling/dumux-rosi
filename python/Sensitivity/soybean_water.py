""" 
    Soybean (water only) 
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import matplotlib.pyplot as plt

import scenario_setup as scenario
from rhizo_models import RhizoMappedSegments
from xylem_flux import sinusoidal2
import cyl
import sra
import evapotranspiration as evap

""" parameters   """
min_b = [-19, -2.5, -200.]  # Domain [38 cm Reihe, 5 cm Pflanzen]
max_b = [19, 2.5, 0.]
cell_number = [1, 1, 200]  # 1 cm3

soil0 = [0.0809, 0.52, 0.0071, 1.5734, 99.49]
soil1 = [0.0874, 0.5359, 0.0087, 1.5231, 93]
soil36 = [0.0942, 0.5569, 0.0089, 1.4974, 87.79]
soil5 = [0.0539, 0.5193, 0.024, 1.4046, 208.78]
soil59 = [0.0675, 0.5109, 0.0111, 1.4756, 107.63]
table_name = "envirotype0"
soil_ = soil0

Kc_soybean = 1.15  # book "crop evapotranspiration" Allen, et al 1998

area = (38 * 5)
trans = 0.6 * area  # cm3/day (38 * 5 = 190)

sim_time = 1. * 87.5  # [day]
dt = 360 / (24 * 3600)  # time step [day]

range_ = ['1995-03-15 00:00:00', '1995-06-10 11:00:00']
x_, y_ = evap.net_infiltration_table_beers('data/95.pkl', range_, 87.5, evap.lai_soybean, Kc_soybean)
trans_soybean = evap.get_transpiration_beers('data/95.pkl', range_, area, 87.5, evap.lai_soybean, Kc_soybean)
# plt.plot(x_, y_)
# plt.show()
# dd

""" initialize """
p_top = -330
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = p_top, p_bot = (p_top + 200), type = 1, times = x_, net_inf = y_)  # , times = x_, net_inf = y_

xml_name = "Glycine_max_Moraes2020_opt2" + "_modified" + ".xml"  # root growth model parameter file
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
# scenario.init_conductivities_const_growth(r)
scenario.init_lupine_conductivities(r)
# scenario.init_dynamic_simple_growth(r, 1.e-3, 4.e-3, 5.e-2, 2.e-3)
# r.plot_conductivities(monocot = True, plot_now = True, axes_ind = [1, 4, 5], lateral_ind = [2, 3])

# rsml_name = "results/soybean.rsml"  # created by rootsystem_soybean.py
# r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, rsml_name)

trans_f1 = lambda t, dt:-trans * sinusoidal2(t, dt) - 0.01  # Guilaumes questions
trans_f2 = lambda t, dt:-trans * sinusoidal2(t, dt) * (t / sim_time)  # growing potential transpiration

sra_table_lookup = sra.open_sra_lookup("data/" + table_name)  # make sure the soil parameters correspond to the look up table

""" sanity checks """
if rank == 0:
    r.test()  # sanity checks

""" numerical solution """
water0 = s.getWaterVolume()  # total initial water volume in domain

psi_x_, psi_s_, sink_, x_, y_, psi_s2_ = sra.simulate_dynamic(s, r, sra_table_lookup, trans, sim_time, dt, trans_soybean)
# psi_x_, psi_s_, sink_, x_, y_, psi_s2_ = sra.simulate_const(s, r, sra_table_lookup, trans, sim_time, dt)

water = s.getWaterVolume()

""" output """
if rank == 0:

    scenario.write_files("soybean_sra0", psi_x_, psi_s_, sink_, x_, y_, psi_s2_)

    print("\ntotal uptake", water0 - water, "cm3")
    print("fin")
