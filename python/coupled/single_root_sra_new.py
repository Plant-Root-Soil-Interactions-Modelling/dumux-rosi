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

""" 
Mai et al (2019) scenario 1 water movement (Schöder approximation)

getInnerHead2, getInnerHead3 are expiremental and show unstable behaviour 

"""


def sra_nostress(r, p, q_root, q_out, r_in, r_out, soil):
    rho = r_out / r_in
    mfp = vg.fast_mfp[soil](p) + (q_root * r_in - q_out * r_out) * (r ** 2 / r_in ** 2 / (2 * (1 - rho ** 2)) + rho ** 2 / (1 - rho ** 2) * (np.log(r_out / r) - 0.5)) + q_out * r_out * np.log(r / r_out)
    return vg.fast_imfp[soil](mfp)


def sra_stress(r, p, q_out, r_in, r_out, soil):
    rho = r_out / r_in
    mfp = (vg.fast_mfp[soil](p) + q_out * r_out * np.log(1 / rho)) * ((r ** 2 / r_in ** 2 - 1 + 2 * rho ** 2 * np.log(r_in / r)) / (rho ** 2 - 1 + 2 * rho ** 2 * np.log(1 / rho))) + q_out * r_out * np.log(r / r_in)
    return vg.fast_imfp[soil](mfp)


def sra_flux(p, q_root, q_out, r_in, r_out, soil):
    r = r_in
    dx = 1.e-6  # [cm]
    h0 = sra_nostress(r, p, q_root, q_out, r_in, r_out, soil)
    h1 = sra_nostress(r + dx, p, q_root, q_out, r_in, r_out, soil)
    hc = vg.hydraulic_conductivity(h0, soil)
    f = hc * (h1 - h0) / dx
    return f


def stressed_flux(p, q_out, r_in, r_out, soil):
    r = r_in
    dx = 1.e-6  # [cm]
    h0 = -15000
    h1 = sra_stress(r + dx, p, q_out, r_in, r_out, soil)
    hc = vg.hydraulic_conductivity(h0, soil)
    f = hc * (h1 - h0) / dx
    return f


N = 5  # number of cells in z dimension
sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
sp = vg.Parameters(clay)
vg.create_mfp_lookup(sp)
initial = -100.  # [cm] initial soil matric potential

r_root = 0.02  # [cm] root radius
kr = 2.e-13 * 1000 * 9.81  # [m / (Pa s)] -> [ 1 / s ]
kx = 5.e-17 * 1000 * 9.81  # [m^4 / (Pa s)] -> [m3 / s]
kr = kr * 24 * 3600  # [ 1 / s ] -> [1/day]
kx = kx * 1.e6 * 24 * 3600  # [ m3 / s ] -> [cm3/day]

NC = 10  # dof of the cylindrical problem
logbase = 1.5

q_r = 1.e-5 * 24 * 3600 * (2 * np.pi * r_root * 3)  # [cm / s] -> [cm3 / day]
sim_time = 13  # [day]
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
s.createGrid([-1.5, -1.5, -3.], [1.5, 1.5, 0.], [3, 3, N])  # [cm] 3x3x3
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
rsx = np.ones((ns,)) * initial  # xylem pressure at the root soil interface
seg_soil_fluxes = np.zeros((ns,))
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

    # solves root model [classic]
    rx = r.solve(0., -q_r, sx[cci], sx, True, critP)  # [cm]

    # fluxes per segment according to the classic sink (including hydraulic conductivities as linmit)
    seg_root_fluxes = r.segFluxes(0., rx, sx, approx = False, cells = True)  # [cm3/day]

    # fluxes per segment according to schröders stressed sink
    for j in range(0, ns):  # for each segment
        p = sx[r.rs.seg2cell[j]]  # soil cell matric potential [cm]
        seg_soil_fluxes[j] = stressed_flux(p, 0., inner_radii[j], outer_radii[j], sp)  # [cm]
        seg_soil_fluxes[j] *= -2. * np.pi * inner_radii[j]  # * l (=1)

    print(seg_root_fluxes, p)
    print(seg_soil_fluxes, p)
    seg_fluxes = np.maximum(seg_root_fluxes, seg_soil_fluxes)

    # water for net flux
    soil_water = np.multiply(np.array(s.getWaterContent()), cell_volumes)

    water_domain.append(np.sum(soil_water))
    water_cell.append(soil_water[cci])
    min_rx.append(np.min(np.array(rx)));
    p1d.append(np.min(np.array(rsx)));

    soil_fluxes = r.sumSegFluxes(seg_fluxes)  # [cm3/day]
    sum_flux = 0.
    for k, f in soil_fluxes.items():
        sum_flux += f
    sum_flux2 = 0.
    for f in seg_fluxes:
         sum_flux2 += f
    # print("Summed fluxes {:g} == {:g},  predescribed {:g}".format(float(sum_flux), float(sum_flux2), -q_r))
    uptake.append(sum_flux)

    # run macroscopic soil model
    s.setSource(soil_fluxes.copy())  # [cm3/day], richards.py
    s.solve(dt)

    x0 = s.getSolutionHeadAt(cci)  # [cm]
    collar.append(x0)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

ax1.set_title("Water amount")
ax1.plot(np.linspace(0, sim_time, NT), np.array(water_cell), label = "water cell")
ax1.legend()
ax1.set_xlabel("Time (days)")
ax1.set_ylabel("(cm3)")

ax2.set_title("Pressure")
ax2.plot(np.linspace(0, sim_time, NT), np.array(collar), label = "soil voxel at root collar")
ax2.plot(np.linspace(0, sim_time, NT), np.array(min_rx), label = "root collar")
ax2.plot(np.linspace(0, sim_time, NT), np.array(p1d), label = "1d cylindrical model at root surface")
ax2.legend()
ax2.set_xlabel("Time (days)")
ax2.set_ylabel("Matric potential (cm)")
ax2.set_ylim(-15000, 0)

ax3.set_title("Water uptake")
ax3.plot(np.linspace(0, sim_time, NT), -np.array(uptake))
ax3.set_xlabel("Time (days)")
ax3.set_ylabel("Uptake (cm/day)")

ax4.set_title("Water in domain")
ax4.plot(np.linspace(0, sim_time, NT), np.array(water_domain))
ax4.set_xlabel("Time (days)")
ax4.set_ylabel("cm3")
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

