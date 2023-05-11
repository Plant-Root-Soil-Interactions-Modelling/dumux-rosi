""" 
Benchmark C1.2 for a static root system in soil (1D or 3D)
with upscaled C_comp matrix with the steady rate approach and fixed-point-iteration in HESS paper notation

similar speed for 3d, great speedup for 1d

!without matrix inversion! 
"""
import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import timeit
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse

# import vtk_plot as vp
from rhizo_models import plot_transpiration
from scenario_setup import *

r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("3D")
name = "results/c12_up_sra"

# r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("3D", age_dependent = True)
# name = "results/c12b_up_sra"
# dt = 120 / (24 * 3600)  # does not work well, maybe uspacling of kr should be based on that of Kr

# r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip = set_scenario("1D")
# name = "results/c12_up_sra_1d"

max_error = 10.
max_iter = 100

ns = len(r.rs.segments)
nodes = r.get_nodes()
inner_ = r.rs.radii
outer_ = r.rs.segOuterRadii()
rho_ = np.divide(outer_, np.array(inner_))
rho_ = np.expand_dims(rho_, axis = 1)
kr_ = r.getKr(rs_age)
inner_kr_ = np.multiply(inner_, kr_)  # multiply for table look up
inner_kr_ = np.expand_dims(inner_kr_, axis = 1)

""" Doussan """
A_d, Kr, kx0 = r.doussan_system_matrix(rs_age)
Id = sparse.identity(ns).tocsc()  # identity matrix

A_n = A_d.copy()
A_n[0, 0] -= kx0
Kr_inv = sparse.linalg.inv(Kr)  # Kr is a diagonal matrix, thus Kr_inv sparse

B, soil2matrix, matrix2soil = r.get_soil_matrix()
nmax = len(matrix2soil)
Bt = B.transpose()
BBt_inv = sparse.linalg.inv(B @ Bt)  # sparse
# print(B.shape)
# print(B.tocsc()[:, 0])
# print((B @ Bt).tocsc()[0, 0])
# dd

""" for comparison """
Ainv_dirichlet = sparse.linalg.inv(A_d).todense()  # dense
Ainv_neumann = sparse.linalg.inv(A_n).todense()  # dense
C_comp_dirichlet = Kr @ (Id - Ainv_dirichlet @ Kr)  # Neumann, Hess, Eqn (24)
c_dirichlet = (Kr @ Ainv_dirichlet)[:, 0] * (-kx0)  # # Hess (25)
C_comp_neumann = Kr @ (Id - Ainv_neumann @ Kr)  # Neum
c_neumann = (Kr @ Ainv_neumann)[:, 0]  # Hess (33)
AinvKr_neumann_up = (((B @ Ainv_neumann) @ Kr) @ Bt)
Ainv_neumann_up = B @ Ainv_neumann
C_comp_neumann_up = B @ C_comp_neumann @ Bt
c_neumann_up = B @ c_neumann
AinvKr_dirichlet_up = (((B @ Ainv_dirichlet) @ Kr) @ Bt)
Ainv_dirichlet_up = B @ Ainv_dirichlet
C_comp_dirichlet_up = B @ C_comp_dirichlet @ Bt
c_dirichlet_up = B @ c_dirichlet
# print(C_comp_neumann_up.shape, type(C_comp_neumann_up))
""" --- """

A_dq = A_d @ Kr_inv
A_nq = A_n @ Kr_inv
Bd = A_d - Kr
Bn = A_n - Kr
A_dq_up = (B @ A_dq) @ Bt
A_nq_up = (B @ A_nq) @ Bt
Bd_up = (B @ Bd) @ Bt
Bn_up = (B @ Bn) @ Bt

Kr_up = (B @ Kr) @ Bt  # sparse
Kr_up_inv = sparse.linalg.inv(Kr_up)

inner_kr_up = BBt_inv.dot(B.dot(inner_kr_))
inner_kr_up = np.maximum(inner_kr_up, np.ones(inner_kr_up.shape) * 1.e-7)  ############################################ (too keep within table)
inner_kr_up = np.minimum(inner_kr_up, np.ones(inner_kr_up.shape) * 1.e-4)  ############################################ (too keep within table)

rho_up = BBt_inv.dot(B.dot(rho_))

""" Numerical solution (a) """
start_time = timeit.default_timer()
rs_age = np.max(r.get_ages())
x_, y_, w_, cpx, cps, cf = [], [], [], [], [], []
sink1d = []

N = round(sim_time / dt)
t = 0.

rx = [0]
sx = s.getSolutionHead()  # inital condition, solverbase.py
centers = s.getCellCenters()

t_pot = -trans * sinusoidal(t)

hs_ = np.zeros((nmax, 1))  # sx -> hs_ # soil cell indices to soil matrix indices
for j in soil2matrix.keys():
        hs_[soil2matrix[j]] += sx[j] + centers[j, 2]

hs = np.transpose(np.array([[sx[mapping[j]][0] for j in range(0, ns)]]))
for j in range(0, len(nodes) - 1):  # from matric to total
    hs[j, 0] += nodes[j + 1][2]

hxd = BBt_inv.dot(AinvKr_dirichlet_up.dot(hs_) + Ainv_dirichlet_up[:, 0] * kx0 * wilting_point)
q_dirichlet_up = -Kr_up.dot(hs_ - hxd)
# if np.sum(q_dirichlet_up) > t_pot:
rx3 = hxd
# else:
#     rx2 = BBt_inv.dot(AinvKr_neumann_up.dot(hs_) + Ainv_neumann_up[:, 0] * t_pot)

# code can be reduced, always neumann for t_pot = 0, i.e. L79,81, & 82

# print("hs_", hs_.shape)
# print(BBt_inv[0, 0])

hs2 = hs_.copy()
hs2[0, 0] -= kx0 * wilting_point
b = Bd_up.dot(hs2)  # hs_
# b[0, 0] -= 0.1 * kx0 * wilting_point
q_dirichlet = BBt_inv.dot(-sparse.linalg.spsolve(A_dq_up, b))  # -3.049 * np.ones(q_dirichlet_up.shape)  #

b = Bd.dot(hs)
b[0, 0] -= kx0 * wilting_point
q_dirichlet2 = -sparse.linalg.spsolve(A_dq, b)  # not upsacled

# if np.sum(q_dirichlet) > t_pot:
rx = hs_ - Kr_up_inv.dot(-q_dirichlet)  # [:, np.newaxis]
rx2 = hs - Kr_inv.dot(-q_dirichlet2[:, np.newaxis])
# else:
#     b = Bn_up.dot(hs_)
#     b[0, 0] -= t_pot
#     q_neumann = -sparse.linalg.spsolve(A_nq_up, b)
#     rx = hs_ - Kr_up_inv.dot(-q_neumann[:, np.newaxis])

print("q_dirichlet(rx)    ", np.min(q_dirichlet), np.max(q_dirichlet))  # not inverted, (not exactyl the same to rx3)
print("q_dirichlet2(rx2)  ", np.min(q_dirichlet2), np.max(q_dirichlet2))  # upscaled
print("q_dirichlet_up(rx3)", np.min(q_dirichlet_up), np.max(q_dirichlet_up))

print()
print("rx", np.min(rx), np.max(rx))
print("rx2", np.min(rx2), np.max(rx2))  # not inverted, (not exactyl the same to rx3)
print("rx3", np.min(rx3), np.max(rx3))  # upscaled

print()
print("hs_(rx)", np.min(hs_), np.max(hs_))  # not inverted, (not exactyl the same to rx3)
print("hs(rx2)", np.min(hs), np.max(hs))  # upscaled
dd

for i in range(0, N):

    t_pot = -trans * sinusoidal(t)  # potential transpiration ...
    print("t_pot", t_pot)

    hs_ = np.zeros((nmax, 1))  # sx -> hs_ # soil cell indices to soil matrix indices
    for j in soil2matrix.keys():
            hs_[soil2matrix[j]] += sx[j] + centers[j, 2]

    err = 1.e6
    c = 1
    rx_old = rx
    while err > max_error and c < max_iter:

        """ interpolation """
        wall_interpolation = timeit.default_timer()

        # print(rx.shape, type(rx), hs_.shape, type(hs_), inner_kr_.shape, type(inner_kr_), rho_.shape, type(rho_))
        for j in soil2matrix.keys():  # from total to matric
                rx[soil2matrix[j]] -= centers[j, 2]

        rsx = soil_root_interface_table(rx, hs_, inner_kr_up, rho_up, sra_table_lookup)

        for j in soil2matrix.keys():  # from matric to total
                rsx[soil2matrix[j]] += centers[j, 2]

        wall_interpolation = timeit.default_timer() - wall_interpolation

        """ xylem matric potential """
        wall_xylem = timeit.default_timer()

        b = Bd_up.dot(rsx)
        b[0, 0] -= kx0 * wilting_point
        q_dirichlet = -sparse.linalg.spsolve(A_dq_up, b)
        # q_dirichlet = (B @ Bt).dot(q_dirichlet)
        # q_dirichlet = BBt_inv.dot(q_dirichlet)
        # rx = rsx - Kr_up_inv.dot(-q_dirichlet[:, np.newaxis])
        #
        # hxd = BBt_inv.dot(AinvKr_dirichlet_up.dot(rsx) + Ainv_dirichlet_up[:, 0] * kx0 * wilting_point)
        # q_dirichlet_up = -Kr_up.dot(rsx - hxd)
        #
        # # rx = rsx - Kr_up_inv.dot(-q_dirichlet_up)
        # print("Dirichlet")
        # print("rx", np.min(rx), np.max(rx), np.min(hxd), np.max(hxd))
        # print("fluxes", np.sum(q_dirichlet), np.sum(q_dirichlet_up), np.sum(q_neumann), t_pot)
        # xx
        if np.sum(q_dirichlet) > t_pot:
            rx = rsx - Kr_up_inv.dot(-q_dirichlet[:, np.newaxis])
            print("Dirichlet rx", np.min(rx), np.max(rx), np.sum(q_dirichlet), np.sum(q_neumann), t_pot)
        else:
            b = Bn_up.dot(rsx)
            b[0, 0] -= t_pot
            q_neumann = -sparse.linalg.spsolve(A_nq_up, b)
            q_neumann = (B @ Bt).dot(q_neumann)
            rx = rsx - Kr_up_inv.dot(-q_neumann[:, np.newaxis])

            print("Neumann")
            print("rx", np.min(rx), np.max(rx))
            print("rsx", np.min(rsx), np.max(rsx))
            print("hs_", np.min(hs_), np.max(hs_))
            print("fluxes", np.sum(q_dirichlet), np.sum(q_neumann), t_pot)

            # print(rx.shape)
            # print(rsx.shape)
            # print(q_dirichlet.shape)
            # print(q_neumann.shape)

        err = np.linalg.norm(rx - rx_old)
        wall_xylem = timeit.default_timer() - wall_xylem
        rx_old = rx.copy()
        c += 1

    if np.sum(q_dirichlet) > t_pot:
        print("q_dirichlet", c, err, np.sum(q_dirichlet), t_pot)
        fluxes = {}
        for j in range(0, nmax):
            fluxes[matrix2soil[j]] = q_dirichlet[j]
    else:
        print("q_neumann", c, err, q_neumann.shape, np.sum(q_neumann), t_pot)
        fluxes = {}
        for j in range(0, nmax):
            fluxes[matrix2soil[j]] = q_neumann[j]

    water = s.getWaterVolume()
    s.setSource(fluxes.copy())  # richards.py
    s.solve(dt)
    soil_water = (s.getWaterVolume() - water) / dt

    old_sx = sx.copy()
    sx = s.getSolutionHead()  # richards.py

    if  i % skip == 0:
        x_.append(t)
        sum_flux = 0.
        for f in fluxes.values():
            sum_flux += f
        # print("sum_flux", sum_flux.shape, sum_flux)
        y_.append(soil_water)  # cm3/day
        cf.append(sum_flux)  # cm3/day
        cpx.append(rx[0])  # cm
        n = round(float(i) / float(N) * 100.)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days {:g}\n"
              .format(np.min(sx), np.max(sx), np.min(rx), np.max(rx), s.simTime, rx[0, 0]))

        # """ Additional sink plot """
        # if i % 60 == 0:  # every 6h
        #     ana = pb.SegmentAnalyser(r.rs)
        #     fluxes = r.segFluxes(rs_age + t, rx, old_sx, False, cells = True)  # cm3/day WRONG sx, should be old sx
        #     ana.addData("fluxes", fluxes)  # cut off for vizualisation
        #     flux1d = ana.distribution("fluxes", max_b[2], min_b[2], 15, False)
        #     sink1d.append(np.array(flux1d))
        #     print("\nSink integral = ", np.sum(np.array(flux1d)), "\n")

    t += dt

s.writeDumuxVTK(name)

""" Plot """
print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")
plot_transpiration(x_, y_, cf, lambda t: trans * sinusoidal(t))
np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter = ';')

# vp.plot_roots_and_soil(r.rs, "pressure head", rx, s, False, min_b, max_b, cell_number, name)  # VTK vizualisation
# sink1d = np.array(sink1d)
# np.save("sink1d", sink1d)
# print(sink1d.shape)
# print(sink1d[-1])

