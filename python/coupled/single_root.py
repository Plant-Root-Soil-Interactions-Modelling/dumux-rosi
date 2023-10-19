import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml.rsml_reader as rsml
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import functional.van_genuchten as vg

from math import *
import numpy as np
import matplotlib.pyplot as plt
import timeit

""" 
Mai et al (2019) Scenario 1 water movement (Dumux solver for cylindrical model) 
"""
N = 1  # number of cells in each dimension
loam = [0.08, 0.43, 0.04, 1.6, 50]
sp = vg.Parameters(loam)
initial = -100.  # [cm] initial soil matric potential

r_root = 0.02  # [cm] root radius
kr = 2.e-13 * 1000 * 9.81  # [m / (Pa s)] -> [ 1 / s ]
kx = 5.e-17 * 1000 * 9.81  # [m^4 / (Pa s)] -> [m3 / s]
kr = kr * 24 * 3600  # [ 1 / s ] -> [1/day]
kx = kx * 1.e6 * 24 * 3600  # [ m3 / s ] -> [cm3/day]
# print(kr, kx)

NC = 10  #  NC-1 are the dof of the cylindrical problem
logbase = 1.5

q_r = 1.e-5 * 24 * 3600 * (2 * np.pi * r_root * 3)  # [cm / s] -> [cm3 / day] # transpiration
sim_time = 15  # [day]
NT = 200  # iteration

critP = -15000  # [cm]

""" Initialize macroscopic soil model """
s = RichardsWrapper(RichardsSP())
s.initialize()
s.setParameter("Problem.EnableGravity", "false")  # important in 1d axial-symmetric problem
s.setHomogeneousIC(initial)  # cm pressure head
s.setTopBC("noflux")
s.setBotBC("noflux")
s.createGrid([-1.5, -1.5, -3.], [1.5, 1.5, 0.], [3, 3, N])  # [cm] 3x3x3
s.setVGParameters([loam])
s.initializeProblem()
s.ddt = 1.e-5  # [day] initial Dumux time step

""" Initialize xylem model """
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

picker = lambda x, y, z: s.pick([x, y, z])
r.rs.setSoilGrid(picker)
cci = picker(0, 0, 0)  # collar cell index
# print("collar index", cci)
# for i in range(0, len(n) - 1):  # segments
#     print("soil index", r.rs.seg2cell[i])

r_outer = r.rs.segOuterRadii()
seg_length = r.rs.segLength()

""" Initialize local soil models (around each root segment) """
cyls = []
points = np.logspace(np.log(r_root) / np.log(logbase), np.log(r_outer[i]) / np.log(logbase), NC, base = logbase)
ns = len(seg_length)  # number of segments
for i in range(0, ns):
    cyl = RichardsWrapper(RichardsCylFoam())
    cyl.initialize()
    cyl.createGrid1d(points)
    cyl.setHomogeneousIC(initial)  # cm pressure head
    cyl.setVGParameters([loam])
    cyl.setOuterBC("fluxCyl", 0.)  # [cm/day]
    cyl.setInnerBC("fluxCyl", 0.)  # [cm/day]
    cyl.initializeProblem()
    cyl.setCriticalPressure(critP)  # cm pressure head
    cyls.append(cyl)

""" Simulation """
dt = sim_time / NT

min_rx, min_rsx, collar_sx = [], [], []  # cm
water_uptake, water_collar_cell, water_cyl, water_domain = [], [], [], []  # cm3

rsx = np.zeros((ns,))  # matric potential at the root soil interface [cm]
cell_volumes = s.getCellVolumes()
net_flux = np.zeros(cell_volumes.shape)
realized_inner_fluxes = np.zeros((len(cyls),))

for i in range(0, NT):

    """ 
    Xylem model 
    """
    csx = s.getSolutionHeadAt(cci)
    for j, cyl in enumerate(cyls):  # for each segment
        rsx[j] = cyl.getInnerHead()  # [cm]

    soil_k = vg.hydraulic_conductivity(rsx, sp) / r_root
    rx = r.solve(0., -q_r, 0., rsx, False, critP, soil_k)  # [cm]

    min_rsx.append(np.min(np.array(rsx)))
    collar_sx.append(csx)
    min_rx.append(np.min(np.array(rx)))
    print("Minimum of cylindrical model {:g} cm, minimal root xylem pressure {:g} cm".format(min_rsx[-1], min_rx[-1]))

    """
    Local soil model
    """
    proposed_inner_fluxes = r.segFluxes(0., rx, rsx, approx = False)  # [cm3/day]
    proposed_outer_fluxes = r.splitSoilFluxes(net_flux / dt)
    for j, cyl in enumerate(cyls):  # set cylindrical model fluxes
        l = seg_length[j]
        cyl.setInnerFluxCyl(proposed_inner_fluxes[j] / (2 * np.pi * r_root * l))  # [cm3/day] -> [cm /day]
        cyl.setOuterFluxCyl(proposed_outer_fluxes[j] / (2 * np.pi * r_outer[j] * l))  # [cm3/day] -> [cm /day]
        cyl.ddt = 1.e-5  # [day] initial time step
        cyl.solve(dt)
        realized_inner_fluxes[j] = -cyl.getInnerFlux() * (2 * np.pi * r_root * l) / r_root  # [cm/day] -> [cm3/day], ('/r_root' comes from cylindrical implementation)
        print("Propsed inner flux {:g} [cm3 day-1] realized {:g} [cm3 day-1]".format(proposed_inner_fluxes[j], realized_inner_fluxes[j]))

    """
    Macroscopic soil model
    """
    soil_water = np.multiply(np.array(s.getWaterContent()), cell_volumes)  # water per cell [cm3]
    soil_fluxes = r.sumSegFluxes(realized_inner_fluxes)  # [cm3/day]
    s.setSource(soil_fluxes.copy())  # [cm3/day], richards.py
    s.solve(dt)

    new_soil_water = np.multiply(np.array(s.getWaterContent()), cell_volumes)
    net_flux = new_soil_water - soil_water  # change in water per cell [cm3]
    for k, root_flux in soil_fluxes.items():
        net_flux[k] -= root_flux * dt
    print("Summed net flux {:g}, max movement {:g} cm3".format(np.sum(net_flux), np.max(net_flux)))  # summed fluxes should equal zero

    """ 
    Water (for output)
    """
    soil_water = new_soil_water
    water_domain.append(np.sum(soil_water))  # from previous time step
    sum_flux = 0.
    for k, f in soil_fluxes.items():
        sum_flux += f
    water_uptake.append(sum_flux)  # cm3

    cyl_water_content = cyls[0].getWaterContent()  # segment 0
    cyl_water = 0.
    for i, wc in enumerate(cyl_water_content):
        r1 = points[i]
        r2 = points[i + 1]
        cyl_water += np.pi * (r2 * r2 - r1 * r1) * 1 * wc
    print("Water volume cylindric", cyl_water, "soil", soil_water[cci])
    water_collar_cell.append(soil_water[cci])
    water_cyl.append(cyl_water)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1)

ax1.set_title("Water amount")
ax1.plot(np.linspace(0, sim_time, NT), np.array(water_collar_cell), label = "water in domain")
# ax1.plot(np.linspace(0, sim_time, NT), np.array(water_cyl), label="water cylindric")
ax1.legend()
ax1.set_xlabel("Time (days)")
ax1.set_ylabel("(cm3)")

ax2.set_title("Pressure")
# ax2.plot(np.linspace(0, sim_time, NT), np.array(collar_sx), label="soil at root collar")
ax2.plot(np.linspace(0, sim_time, NT), np.array(min_rx), label = "minimum hx")
ax2.plot(np.linspace(0, sim_time, NT), np.array(min_rsx), label = "minimum hrs")
ax2.legend()
ax2.set_xlabel("Time (days)")
ax2.set_ylabel("Matric potential (cm)")
# plt.ylim(-15000, 0)

ax3.set_title("Water uptake")
ax3.plot(np.linspace(0, sim_time, NT), -np.array(water_uptake))
ax3.set_xlabel("Time (days)")
ax3.set_ylabel("Uptake (cm/day)")

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
#
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
