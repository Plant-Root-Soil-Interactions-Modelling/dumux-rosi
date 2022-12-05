""" 
    Maize using rhizosphere models TODO part for nitrate is unstable 
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import timeit
import numpy as np

import scenario_setup as scenario
from rhizo_models import RhizoMappedSegments
from xylem_flux import sinusoidal2
import evapotranspiration as evap
import cyl3

"""
 Cylindric rhizosphere models, nitrate from finite volumes
"""

""" parameters   """
min_b = [-37.5, -7.5, -200.]  # Domain Mais: 60 cm Reihe, 10 cm Pflanzen
max_b = [37.5, 7.5, 0.]
cell_number = [1, 1, 200]  # 1 cm3

soil0 = [0.0809, 0.52, 0.0071, 1.5734, 99.49]
soil1 = [0.0874, 0.5359, 0.0087, 1.5231, 93]
soil36 = [0.0942, 0.5569, 0.0089, 1.4974, 87.79]
soil5 = [0.0539, 0.5193, 0.024, 1.4046, 208.78]
soil59 = [0.0675, 0.5109, 0.0111, 1.4756, 107.63]
soil_ = soil0

Kc_maize = 1.2  # book "crop evapotranspiration" Allen, et al 1998

area = 75 * 15  # cm2

sim_time = 95  #  [day]
dt = 360 / (24 * 3600)  # time step [day] 20

range_ = ['1995-03-15 00:00:00', '1995-06-20 23:00:00']  # +3
x_, y_ = evap.net_infiltration_table_beers('data/95.pkl', range_, 95, evap.lai_maize, Kc_maize)
trans_maize = evap.get_transpiration_beers('data/95.pkl', range_, area, 95, evap.lai_maize, Kc_maize)

trans = 0.6 * area  # cm3/day
trans_f1 = lambda t, dt:-trans * sinusoidal2(t, dt)  # Guilaumes questions - 0.01
trans_f2 = lambda t, dt:-trans * sinusoidal2(t, dt) * (t / (.5 * 95))  # growing potential transpiration

""" rhizosphere model parameters """
wilting_point = -15000  # cm
nc = 10  # dof+1
logbase = 0.5  # according to Mai et al. (2019)

mode = "dumux_dirichlet"  # mode = "dumux_dirichlet"

""" initialize """
start_time = timeit.default_timer()
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = -330, p_bot = -130, type = 2, times = x_, net_inf = y_)  # , times = x_, net_inf = y_
water0 = s.getWaterVolume()  # total initial water volume in domain

xml_name = "Zeamays_synMRI_modified.xml"  # root growth model parameter file
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
scenario.init_maize_conductivities(r)

if rank == 0:
    r.test()  # sanity checks

psi_x_, psi_s_, sink_, x_, y_, psi_s2_, soil_c_, c_ = [], [], [], [], [], [], [], []

for i in range(0, int(sim_time)):

    print("Day", i)

    r.rs.simulate(1.)  # simulate for 1 day
    rs_age = i + 1

    rs = RhizoMappedSegments(r, wilting_point, nc, logbase, mode)
    if mode == "dumux_dirichlet":
        cyl_type = 1
    elif momde == "dumux_dirichlet_nc":
        cyl_type = 2
    else:
        raise("unknown type")

    psi_x, psi_s, sink, x, y, psi_s2, soil_c, c = cyl3.simulate_const(s, rs, 1., dt, trans_maize, rs_age, type = cyl_type)

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

water = s.getWaterVolume()

""" output """
if rank == 0:

    scenario.write_files("maize_cyl0", psi_x_, psi_s_, sink_, x_, y_, psi_s2_, soil_c_, c_)
    print ("Overall simulation wall time", timeit.default_timer() - start_time, " s")
    print("fin")
