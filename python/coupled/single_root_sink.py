import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml.rsml_reader as rsml
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import functional.van_genuchten as vg

import numpy as np
from scipy import optimize
from scipy import integrate
import matplotlib.pyplot as plt
import timeit


def getInnerHead(p, rp, q_root, q_out, r_in, r_out, soil):
    """ returns the pressure head at the root surface according to Schroeder et al. """
    if rp < p:  # flux into root
        r = r_in
        rho = r_out / r_in
        mfp = vg.fast_mfp[soil](p) + (q_root * r_in - q_out * r_out) * (r ** 2 / r_in ** 2 / (2 * (1 - rho ** 2)) + rho ** 2 / (1 - rho ** 2) * (np.log(r_out / r) - 0.5)) + q_out * r_out * np.log(r / r_out)
        if mfp > 0:
            h = vg.fast_imfp[soil](mfp)
        else:
            h = -15000.
        if rp < h:  # flux into root
            return h
        else:  # flux into soil
            return rp  # no flux
    else:  # flux into soil
        return p  # don't use schröder

""" 
Mai et al (2019) scenario 1 water movement SINK only, No macroscopic flux
"""
N = 3  # number of cells in each dimension
loam = [0.08, 0.43, 0.04, 1.6, 50]
sp = vg.Parameters(loam)
vg.create_mfp_lookup(sp)
initial = -100.  # [cm] initial soil matric potential

r_root = 0.02  # [cm] root radius
kr = 10 * 2.e-13 * 1000 * 9.81  # [m / (Pa s)] -> [ 1 / s ]
kx = 5.e-17 * 1000 * 9.81  # [m^4 / (Pa s)] -> [m3 / s]
kr = kr * 24 * 3600  # [ 1 / s ] -> [1/day]
kx = kx * 1.e6 * 24 * 3600  # [ m3 / s ] -> [cm3/day]
# print(kr, kx)

NC = 10  # dof of the cylindrical problem
logbase = 1.5

q_r = 1.e-5 * 24 * 3600 * (2 * np.pi * r_root * 3)  # [cm / s] -> [cm3 / day]
sim_time = 15  # [day]
NT = 5000  # iteration

critP = -15000  # [cm]

""" Macroscopic soil model """
cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()
# s.setParameter("Problem.EnableGravity", "false")  # important in 1d axial-symmetric problem
s.setHomogeneousIC(initial)  # cm pressure head
s.setTopBC("noflux")
s.setBotBC("noflux")
s.createGrid([-1.5, -1.5, -3.], [1.5, 1.5, 0.], [N, N, N])  # [cm] 3x3x3
s.setVGParameters([loam])
s.initializeProblem()
s.ddt = 1.e-5  # [day] initial Dumux time step

""" Root model """
n, segs = [], []
for i in range(0, 4):  # nodes
    n.append(pb.Vector3d(0, 0, -float(i)))
for i in range(0, len(n) - 1):  # segments
    segs.append(pb.Vector2i(i, i + 1))
rs = pb.MappedSegments(n, segs, [r_root] * len(segs))  # a single root
rs.setRectangularGrid(pb.Vector3d(-1.5, -1.5, -3.), pb.Vector3d(1.5, 1.5, 0.), pb.Vector3d(3, 3, N), False)
r = XylemFluxPython(rs)
r.setKr([kr])
r.setKx([kx])

inner_radii = np.array(r.rs.radii)
outer_radii = r.rs.segOuterRadii()
seg_length = r.rs.segLength()
ns = len(seg_length)

""" Coupling soil and root model (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])
r.rs.setSoilGrid(picker)
cci = picker(0, 0, 0)  # collar cell index
# print("collar index", cci)
# for i in range(0, len(n) - 1):  # segments
#     print("soil index", r.rs.seg2cell[i])

""" Simulation """
sx = s.getSolutionHead()  # [cm]
rsx = np.zeros((ns,))  # xylem pressure at the root soil interface
seg_fluxes = np.zeros((ns,))
seg_outer_fluxes = np.zeros((ns,))
dt = sim_time / NT

uptake = []
min_rx = []

collar = []
p1d = []

water_cell = []
water_domain = []

cell_volumes = s.getCellVolumes()
net_flux = np.zeros(cell_volumes.shape)

for i in range(0, NT):

    sx = s.getSolutionHead()  # [cm]

    # solves root model
    rx = r.solve(0., -q_r, sx[cci], rsx, False, critP)  # [cm]

#     # fluxes per segment according to the classic sink
#     seg_fluxes = r.segFluxes(0., rx, sx, approx = False, cells = True)  # [cm3/day]

    # solutions of previous time step
    for j in range(0, ns):  # for each segment
        p = sx[r.rs.seg2cell[j]]  # soil cell matric potential [cm]
#         q_in = -seg_fluxes[j] / (2. * np.pi * inner_radii[j])  # [cm / day]
#         q_out = -seg_outer_fluxes[j] / (2. * np.pi * outer_radii[j])  # [cm / day]
        rsx[j] = p  # getInnerHead(p, rx[j], q_in, q_out, inner_radii[j], outer_radii[j], sp)  # [cm]

    # fluxes per segment according to schröder guess
    seg_fluxes = r.segFluxes(0., rx, rsx, approx = False, cells = False)  # [cm3/day]

    # water for net flux
    soil_water = np.multiply(np.array(s.getWaterContent()), cell_volumes)

    water_domain.append(np.sum(soil_water))
    water_cell.append(soil_water[cci])
    min_rx.append(np.min(np.array(rx)));
    p1d.append(np.min(np.array(rsx)));
    print("Minimal root xylem pressure", min_rx[-1])
    print("Cylindrical models at root surface", rsx, "cm")
    # print("soil water", np.sum(soil_water))

    seg_outer_fluxes = r.splitSoilFluxes(net_flux / dt)

    soil_fluxes = r.sumSegFluxes(seg_fluxes)  # [cm3/day]
    sum_flux = 0.
    for k, f in soil_fluxes.items():
        sum_flux += f
    sum_flux2 = 0.
    for f in seg_fluxes:
         sum_flux2 += f
    print("Summed fluxes {:g} == {:g},  predescribed {:g}".format(float(sum_flux), float(sum_flux2), -q_r))
    uptake.append(sum_flux)

    # run macroscopic soil model
    s.setSource(soil_fluxes.copy())  # [cm3/day], richards.py
    s.solve(dt)

    x0 = s.getSolutionHeadAt(cci)  # [cm]
    collar.append(x0)

    # calculate net fluxes
    net_flux = (np.multiply(np.array(s.getWaterContent()), cell_volumes) - soil_water)
    for k, v in soil_fluxes.items():
        net_flux[k] -= v * dt;
    # print(net_flux / dt)
    # print("sum", np.sum(net_flux))

fig, (ax1, ax2) = plt.subplots(2, 1)

ax1.set_title("Water amount")
ax1.plot(np.linspace(0, sim_time, NT), np.array(water_domain), label = "water in domain")
ax1.legend()
ax1.set_xlabel("Time (days)")
ax1.set_ylabel("(cm3)")

ax2.set_title("Pressure")
# ax2.plot(np.linspace(0, sim_time, NT), np.array(collar), label="soil voxel at root collar")
ax2.plot(np.linspace(0, sim_time, NT), np.array(min_rx), label = "minimum hx")
ax2.plot(np.linspace(0, sim_time, NT), np.array(p1d), label = "minimum hrs")
ax2.legend()
ax2.set_xlabel("Time (days)")
ax2.set_ylabel("Matric potential (cm)")
ax2.set_ylim(-15000, 0)

# ax3.set_title("Water uptake")
# ax3.plot(np.linspace(0, sim_time, NT), -np.array(uptake))
# ax3.set_xlabel("Time (days)")
# ax3.set_ylabel("Uptake (cm/day)")
#
# ax4.set_title("Water in domain")
# ax4.plot(np.linspace(0, sim_time, NT), np.array(water_domain))
# ax4.set_xlabel("Time (days)")
# ax4.set_ylabel("cm3")
plt.show()

# plt.title("Pressure")
# h = np.array(cyls[0].getSolutionHead())
# x = np.array(cyls[0].getDofCoordinates())
# plt.plot(x, h, "b*")
# plt.xlabel("x (cm)")
# plt.ylabel("Matric potential (cm)")
# plt.show()

# plt.title("Darcy velocity")
# h = np.array(cyls[0].getSolutionHead())
# x = np.array(cyls[0].getDofCoordinates())
# print(np.diff(x, axis=0))
# print(np.diff(h, axis=0))
# dhdx = np.divide(np.diff(h, axis=0), np.diff(x, axis=0))
# sp = vg.Parameters(loam)
# hc = []
# for h_ in h[1:]:
#     hc.append(vg.hydraulic_conductivity(h_, sp))  # mid heads
# velocity = np.multiply(np.array(hc), dhdx)
# plt.plot(x[0:-1], velocity / 100. / 24. / 3600., "b*")
# plt.xlabel("x (cm)")
# plt.ylabel("Velocity (m/s)")
# plt.show()

print("fin")
