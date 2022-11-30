""" 
    Soybean (water only) TODO
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import numpy as np

import scenario_setup as scenario
from rhizo_models import RhizoMappedSegments
from xylem_flux import sinusoidal2
import cyl
import sra

""" parameters   """
min_b = [-19, -2.5, -200.]  # Domain [38 cm Reihe, 5 cm Pflanzen]
max_b = [19, 2.5, 0.]
cell_number = [1, 1, 200]  # 1 cm3

theta_r = 0.025  # sandy loam
theta_s = 0.403
alpha = 0.0383  # (cm-1) soil
n = 1.3774
k_sat = 60.  # (cm d-1)
soil_ = [theta_r, theta_s, alpha, n, k_sat]

trans = 10 * 0.6 * (38 * 5)  # cm3/day (38 * 5 = 190)

sim_time = 0.25 * 87.5  # [day]
dt = 1800 / (24 * 3600)  # time step [day]

kr0 = 1.e-3
kx0 = 5.e-2
Nr = 10
Nx = 8
kr0_ = np.linspace(kr0, 100 * kr0, Nr)
kx0_ = np.linspace(kx0, 100 * kx0, Nx)

quick_result = np.zeros((Nr, Nx))

sra_table_lookup = sra.open_sra_lookup("../coupled/sra/table_jan_comp")  # make sure the soil parameters correspond to the look up table
trans_f2 = lambda t, dt:-trans * sinusoidal2(t, dt) * (t / sim_time)  # growing potential transpiration

for i, kr0 in enumerate(kr0_):
    for j, kx0 in enumerate(kx0_):

        print("******************************************")
        print("run", i, j)
        print("******************************************")

        kr1 = 4 * kr0
        kx1 = kx0 / 25

        """ initialize """
        s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = -330, p_bot = -180, type = 1)
        xml_name = "Glycine_max_Moraes2020_opt2" + "_modified" + ".xml"  # root growth model parameter file
        r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth

        """ vary conductivities """
        scenario.init_dynamic_simple_growth(r, kr0, kr1, kx0, kx1)

        """ numerical solution """
        water0 = s.getWaterVolume()  # total initial water volume in domain
        psi_x_, psi_s_, sink_, x_, y_, psi_s2_ = sra.simulate_dynamic(s, r, sra_table_lookup, sim_time, dt, trans_f2)
        water = s.getWaterVolume()

        """ output """
        scenario.write_files("soybean_sra_sa_{:g}{:g}".format(i, j), psi_x_, psi_s_, sink_, x_, y_, psi_s2_)

        quick_result[i, j] = water0 - water
        print("\ntotal uptake", water0 - water, "cm3")

print(quick_result[:,:])
print("fin")
