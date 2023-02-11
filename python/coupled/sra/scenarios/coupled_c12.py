""" 
Jan's new scenarios with the upscaled sink (without SRA)

1d soil, 2 cm thick layers, dynamic conductivities (see root_conductivities.py)
"""
import sys; sys.path.append("../../../modules/"); sys.path.append("../../../../../CPlantBox/");  sys.path.append("../../../../../CPlantBox/src/python_modules")
sys.path.append("../../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../modules/fv/");

import timeit
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy import sparse
import scipy.sparse.linalg as LA

import van_genuchten as vg
from xylem_flux import XylemFluxPython  # Python hybrid solver
from HydraulicsDoussan import HydraulicsDoussan  # Doussan solver
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import plantbox as pb
import aggregated_rs as agg
from root_conductivities import *
import vtk_plot as vp
from rhizo_models import plot_transpiration

def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.

""" 
Benchmark M1.2 static root system in soil (with the classic sink)

also works parallel with mpiexec (slower, due to overhead?)
"""

""" Parameters """
min_b = [-4., -4., -15.]
max_b = [4., 4., 0.]
cell_number = [7, 7, 15]  #  [8, 8, 15]  # [16, 16, 30]  # [32, 32, 60]  # [8, 8, 15]
path = ""
fname = "../../../../grids/RootSystem8.rsml"

name = "classical360"
sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil = loam

initial = -659.8 + 7.5  # -659.8

trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

sim_time = 2  # [day] for task b
age_dependent = False  # conductivities
dt = 360. / (24 * 3600)  # [days] Time step must be very small
skip = 1

""" Initialize macroscopic soil model """
sp = vg.Parameters(soil)  # for debugging
s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_b, max_b, cell_number, False)  # [cm]
s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([soil])
s.setParameter("Newton.EnableChop", "True")
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Soil.SourceSlope", "1000")  # turns regularisation of the source term on, will change the shape of actual transpiration...
s.initializeProblem()
s.setCriticalPressure(wilting_point)

""" Initialize xylem model (a) or (b)"""

r = XylemFluxPython(fname)
r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)  # cutting
init_conductivities(r, age_dependent)

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])
cci = picker(r.rs.nodes[0].x, r.rs.nodes[0].y, r.rs.nodes[0].z)  # collar cell index
r.rs.setSoilGrid(picker)  # maps segment

n = len(r.rs.nodes)
ns = len(r.rs.segments)
nodes = r.get_nodes()
seg2cell = r.rs.seg2cell
mapping = np.array([seg2cell[j] for j in range(0, ns)])

""" experimental stuff """
Id = sparse.identity(ns).tocsc()  # identity matrix

kx_ = np.divide(r.getKx(sim_time), r.rs.segLength())  # / dl
Kx = sparse.diags(kx_).tocsc()
# print("Kx", Kx.shape, Kx[0, 0], Kx[1, 1])

kr_ = np.array(r.getEffKr(sim_time))  # times surface (2 a pi length)
Kr = sparse.diags(kr_).tocsc()
# print("Kr", Kr.shape, Kr[0, 0], Kr[1, 1])

IM = r.get_incidence_matrix().tocsc()
IMt = IM.transpose().tocsc()

L = IMt @ Kx @ IM  # Laplacian, Hess Eqn (10)
L = L[1:, 1:]  # == L_{N-1}
# print("L", L.shape)

print("inv start")

A_dirichlet = (L + Kr).tocsc()
Ainv_dirichlet = sparse.linalg.inv(A_dirichlet).todense()  # dense

A_neumann = A_dirichlet
A_neumann[0, 0] -= kx_[0]
Ainv_neumann = sparse.linalg.inv(A_neumann).todense()  # dense

C_comp_dirichlet = Kr @ (Id - Ainv_dirichlet @ Kr)  # Neumann, Hess, Eqn (24)
c_dirichlet = (Kr @ Ainv_dirichlet)[:, 0] * (-kx_[0])  # # Hess (25)
print("C_comp_dirichlet", type(C_comp_dirichlet), C_comp_dirichlet.shape)
print("c_dirichlet", type(c_dirichlet), c_dirichlet.shape)

C_comp_neumann = Kr @ (Id - Ainv_neumann @ Kr)  # Neumann, Hess, Eqn (32)
c_neumann = (Kr @ Ainv_neumann)[:, 0]  # Hess (33)
print("C_comp_neumann", type(C_comp_neumann), C_comp_neumann.shape)
print("c_neumann", type(c_neumann), c_neumann.shape)

print("inv stop")

""" Numerical solution (a) """
start_time = timeit.default_timer()
rs_age = np.max(r.get_ages())
x_, y_, w_, cpx, cps, cf = [], [], [], [], [], []
sink1d = []
sx = s.getSolutionHead()  # inital condition, solverbase.py

N = round(sim_time / dt)
t = 0.

rx = [0]

for i in range(0, N):

    t_pot = -trans * sinusoidal(t)  # potential transpiration ...
    print(t_pot)
    
    hs = np.transpose(np.array([[sx[mapping[j]][0] for j in range(0, ns)]]))
    for j in range(0, len(nodes) - 1):  # from matric to total
        hs[j, 0] += nodes[j + 1][2]

    print("hs", hs.shape, np.min(hs), np.max(hs), np.argmin(hs))
    print("sx", sx.shape, np.min(sx), np.max(sx), np.argmin(sx))


    # q_neumann = C_comp_neumann.dot(hs) + c_neumann * t_pot
    # q_dirichlet = -(C_comp_dirichlet.dot(hs) + c_dirichlet * -15000)
    
    print()
    hx = Ainv_neumann.dot(Kr.dot(hs)) + Ainv_neumann[:, 0] * t_pot  #   # Hess Eqn (29)
    hx0 = hx[0, 0] + t_pot / kx_[0]  # /kx*l
    print("hx0", hx0, hx[0, 0], hx[1, 0], hx[2, 0])
    hxd = Ainv_dirichlet.dot(Kr.dot(hs)) + Ainv_dirichlet[:, 0] * kx_[0] * wilting_point
    print("hxd", hxd[0, 0], hxd[1, 0], hxd[2, 0])    
    print()
    
    print()
    q_neumann = -Kr.dot(hs - hx)
    print("q_neumann", q_neumann.shape, np.min(q_neumann), np.max(q_neumann), np.argmin(q_neumann), np.sum(q_neumann))    
    q_dirichlet= -Kr.dot(hs - hxd) 
    print("q_dirichlet", q_dirichlet.shape, np.min(q_dirichlet), np.max(q_dirichlet), np.argmin(q_dirichlet), np.sum(q_dirichlet))
    # print("sum_q", np.sum(q), q[0], q[1])
    # hx = [0]
    print()
    #if hx0 < wilting_point:
    if np.sum(q_dirichlet) > t_pot:
        print("dirichlet", np.sum(q_dirichlet),t_pot)
        fluxes = r.sumSegFluxes(q_dirichlet[:, 0])
    else:
        print("neumann", np.sum(q_neumann), t_pot)
        fluxes = r.sumSegFluxes(q_neumann[:, 0])

    s.setSource(fluxes.copy())  # richards.py
    s.solve(dt)
    old_sx = sx.copy()
    sx = s.getSolutionHead()  # richards.py
    water = s.getWaterVolume()

    if  i % skip == 0:
        min_sx = np.min(sx)
        min_rx = np.min(rx)
        max_sx = np.max(sx)
        max_rx = np.max(rx)
        x_.append(t)
        sum_flux = 0.
        for f in fluxes.values():
            sum_flux += f
        # print("Summed fluxes ", sum_flux, "= collar flux", r.collar_flux(rs_age + t, rx, sx), "= prescribed", -trans * sinusoidal(t))
        y_.append(sum_flux)  # cm4/day
        w_.append(water)  # cm3
        # cf.append(float(r.collar_flux(rs_age + t, rx, sx)))  # cm3/day
        cpx.append(rx[0])  # cm
        cps.append(float(sx[cci]))  # cm
        n = round(float(i) / float(N) * 100.)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], soil [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days {:g}\n"
              .format(min_sx, max_sx, min_rx, max_rx, s.simTime, rx[0]))

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
plot_transpiration(x_, y_, y_, lambda t: trans * sinusoidal(t))
np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter = ';')

# vp.plot_roots_and_soil(r.rs, "pressure head", rx, s, False, min_b, max_b, cell_number, name)  # VTK vizualisation
# sink1d = np.array(sink1d)
# np.save("sink1d", sink1d)
# print(sink1d.shape)
# print(sink1d[-1])

