""" 
Single root scenario - soil depletion due to sinusoidal transpiration over 21 days

using steady rate approach and fix point iteration (sra)
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/");  sys.path.append("../")

import numpy as np
import timeit

from xylem_flux import sinusoidal2
import scenario_setup as scenario
import sra

""" parameters   """

min_b = [-1, -1, -150.]  # domain
max_b = [1, 1, 0.]
cell_number = [1, 1, 150]

trans = 0.6 * 4  # cm3/day
wilting_point = -15000

sim_time = 21  #  [day]
dt = 60 / (24 * 3600)  # time step [day]

""" initialize """
sra_table_lookup = sra.open_sra_lookup("../coupled/sra/table_jan_comp")  # make sure soil corresponds to look up table

s, soil = scenario.create_soil_model(min_b, max_b, cell_number, p_top = -330, p_bot = -180)

r = scenario.create_mapped_singleroot(min_b, max_b, cell_number, s, ns = 100, l = 100, a = 0.05)
r.test()  # sanity checks

nodes = r.rs.nodes
segs = r.rs.segments
ns = len(segs)
mapping = np.array([r.rs.seg2cell[j] for j in range(0, ns)])
outer_r = r.rs.segOuterRadii()  # std::sqrt(targetV/(M_PI*l)+radii[i]*radii[i]) = sqrt(4/pi+0.05*0.05) = 1.1294864075
inner_r = r.rs.radii
types = r.rs.subTypes
rho_ = np.divide(outer_r, np.array(inner_r))
# print("Krs", r.get_krs(0.))

""" Numerical solution (a) """
start_time = timeit.default_timer()

psi_x_, psi_s_, sink_ , x_, y_, psi_s2_ = [], [], [], [], [], []  # for post processing

water0 = s.getWaterVolume()
sx = s.getSolutionHead()  # inital condition, solverbase.py
cell_centers = s.getCellCenters()
cell_centers_z = np.array([cell_centers[mapping[j]][2] for j in range(0, ns)])
seg_centers_z = np.array([0.5 * (nodes[s.x].z + nodes[s.y].z) for s in segs])
hsb = np.array([sx[mapping[j]][0] for j in range(0, ns)])  # soil bulk matric potential per segment
rsx = hsb.copy()  # initial values for fix point iteration

NT = int(np.ceil(sim_time / dt))  # number of iterations
skip = 1  # for output and results, skip iteration
t = 0.
rs_age = 0.
rx = r.solve(rs_age + t, -trans * sinusoidal2(t, dt), 0., rsx, False, wilting_point, soil_k = [])
rx_old = rx.copy()
kr_ = np.array([r.kr_f(rs_age + t, types[j]) for j in range(0, len(outer_r))])
inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up #  TODO move UP

for i in range(0, NT):

    t = i * dt  # current simulation time

    wall_iteration = timeit.default_timer()
    wall_fixpoint = timeit.default_timer()

    err = 1.e6  # cm
    c = 0
    while err > 1 and c < 1000:

        """ interpolation """
        wall_interpolation = timeit.default_timer()
        rx_ = rx[1:] - seg_centers_z  # from total matric potenti    al to matric potential
        hsb_ = hsb - cell_centers_z  # from total matric potential to matric potential
        rsx = sra.soil_root_interface_table(rx_ , hsb_, inner_kr_, rho_, sra_table_lookup)  # rsx = soil_root_interface(rx[1:] , hsb, inner_kr_, rho_, soil)
        rsx = rsx + seg_centers_z  # from matric potential to total matric potential
        wall_interpolation = timeit.default_timer() - wall_interpolation

        """ xylem matric potential """
        wall_xylem = timeit.default_timer()
        rx = r.solve(rs_age + t, -trans * sinusoidal2(t, dt), 0., rsx, False, wilting_point, soil_k = [])  # xylem_flux.py, cells = False
        err = np.linalg.norm(rx - rx_old)
        wall_xylem = timeit.default_timer() - wall_xylem

        rx_old = rx.copy()
        c += 1

    print(i, c, "iterations", wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem))

    wall_fixpoint = timeit.default_timer() - wall_fixpoint

    wall_soil = timeit.default_timer()
    fluxes = r.segFluxes(rs_age + t, rx, rsx, approx = False, cells = False)
    # validity check
    collar_flux = r.collar_flux(rs_age + t, rx.copy(), rsx.copy(), k_soil = [], cells = False)
    err = np.linalg.norm(np.sum(fluxes) - collar_flux)
    if err > 1.e-10:
        print("error: summed root surface fluxes and root collar flux differ" , err)
        raise
    err2 = np.linalg.norm(-trans * sinusoidal2(t, dt) - collar_flux)
    if r.last == "neumann":
        if err2 > 1.e-10:
            print("error: potential transpiration differs root collar flux in Neumann case" , err2)
            raise
    soil_fluxes = r.sumSegFluxes(fluxes)
    s.setSource(soil_fluxes.copy())  # richards.py
    s.solve(dt)
    sx = s.getSolutionHead()[:, 0]  # richards.py
    hsb = np.array([sx[mapping[j]] for j in range(0, ns)])  # soil bulk matric potential per segment
    wall_soil = timeit.default_timer() - wall_soil

    wall_iteration = timeit.default_timer() - wall_iteration

    """ remember results ... """
    if i % skip == 0:
        psi_x_.append(rx.copy())  # cm
        psi_s_.append(rsx.copy())  # cm
        sink_.append(fluxes.copy())  # cm3/day
        x_.append(t)  # day
        y_.append(np.sum(fluxes))  # cm3/day
        psi_s2_.append(np.array(sx))  # cm

print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

scenario.write_files("singleroot_sra", psi_x_, psi_s_, sink_, x_, y_, psi_s2_)

water_end = s.getWaterVolume()
print("\ntotal uptake", water0 - water_end, "cm3")
print("fin")

