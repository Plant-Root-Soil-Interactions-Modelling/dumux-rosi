""" 
Single root scenario: Soil depletion due to sinusoidal transpiration over 21 days

using an aggregated root system, and steady rate approach and fix point iteration (agg)
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import scenario_setup as scenario
import sra
import agg

""" parameters   """
min_b = [-6, -1.5, -150.]  # domain 12cm x 3cm x 150cm
max_b = [6, 1.5, 0.]
cell_number = [1, 1, 150]  # 1 cm3

theta_r = 0.025
theta_s = 0.403
alpha = 0.0383  # (cm-1) soil
n = 1.3774
k_sat = 60.  # (cm d-1)
soil_ = [theta_r, theta_s, alpha, n, k_sat]

trans = 0.6 * (12 * 3)  # cm3/day

sim_time = 21  #  [day]
dt = 60. / (24 * 3600)  # time step [day]

""" 
Initialize xylem model 
"""
""" initialize """
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = -330, p_bot = -180)
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, "results/wheat.rsml")  # created by rootsystem.py
r_agg = agg.create_aggregated_rs(r, 0., min_b, max_b, cell_number)
sra_table_lookup = sra.open_sra_lookup("../coupled/table_jan_comp")  # make sure the soil parameters correspond to the look up table
# sra_table_lookup = soil  # without using the lookup table.
picker = lambda x, y, z: s.pick([0., 0., z])
r_agg.rs.setSoilGrid(picker)

""" sanity checks """
nodes = r_agg.rs.nodes
segs = r_agg.rs.segments
# print(r.rs.nodes[0], r.rs.nodes[1], r.rs.nodes[2], r.rs.nodes[3], r.rs.nodes[-4], r.rs.nodes[-3], r.rs.nodes[-2], r.rs.nodes[-1])
# print(r.rs.segments[0], r.rs.segments[1], r.rs.segments[2], r.rs.segments[3], r.rs.segments[-4], r.rs.segments[-3], r.rs.segments[-2], r.rs.segments[-1])
# vp.plot_roots(pb.SegmentAnalyser(rs), "radius")
# r.test()  # sanity checks

""" numerical solution """
water0 = s.getWaterVolume()  # total initial water volume in domain

psi_x_, psi_s_, sink_, x_, y_, psi_s2_ = agg.simulate_const(s, r_agg, sra_table_lookup, trans, sim_time, dt)

scenario.write_files("rootsystem1d_agg", psi_x_, psi_s_, sink_, x_, y_, psi_s2_)

print("\ntotal uptake", water0 - s.getWaterVolume(), "cm3")
print("fin")

