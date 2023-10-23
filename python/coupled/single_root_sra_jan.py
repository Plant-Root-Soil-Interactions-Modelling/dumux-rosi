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

from scipy.optimize import fsolve

""" 

"""

N = 3  # number of cells in z dimension
sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
sp = vg.Parameters(clay)
initial = -100.  # [cm] initial soil matric potential

r_root = 0.02  # [cm] root radius
kr = 10 * 2.e-13 * 1000 * 9.81  # [m / (Pa s)] -> [ 1 / s ]
kx = 5.e-17 * 1000 * 9.81  # [m^4 / (Pa s)] -> [m3 / s]
kr = kr * 24 * 3600  # [ 1 / s ] -> [1/day]
kx = kx * 1.e6 * 24 * 3600  # [ m3 / s ] -> [cm3/day]

q_r = 1.e-5 * 24 * 3600 * (2 * np.pi * r_root * 3)  # [cm / s] -> [cm3 / day]
collar_p = -10000

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

k_soil = np.zeros((ns, 1))

""" initialize """
sx = s.getSolutionHead()  # [cm]
hx = r.solve(0., -q_r, 0, sx, True, critP)
hsb = np.array([sx[rs.seg2cell[j]] for j in range(0, ns)])
hrs = 0.5 * (hsb + hx[1:])

# b = 2  # r_root  # (1 - r_root) / (0.6 - r_root)  # todo (?)
rho = 0.6 / r_root
b = 2 * (rho * rho - 1) / (1 + 2 * rho * rho * (np.log(rho) - 0.5))

k_root = np.multiply(np.array([r.kr_f(0., rs.subTypes[j], rs.organTypes[j], 0) for j in range(0, ns)]), np.array(rs.radii))

for i in range(0, NT):

    sx = s.getSolutionHead()  # [cm]
    hsb = np.array([sx[rs.seg2cell[j]] for j in range(0, ns)])  # soil bulk matric potential per segment
    hx = np.array(r.solve(0., -q_r, 0, hrs, False, critP))

    hrs = 0.5 * (hsb + hx[1:])
    hrs0 = hrs.copy()  # remember

    """ slide 12 of jans presentation """

    for i in range(0, 20):

        hx = r.solve(0., -q_r, 0, hrs, False, critP)  # [cm] solve(self, sim_time:float, trans:list, sx:float, sxx, cells:bool, wilting_point:float, soil_k=[]):
        # hx = r.solve_dirichlet(0., collar_p, 0, hrs, False)

#         k_soil = np.array([vg.hydraulic_conductivity(sx[rs.seg2cell[j]][0], sp) for j in range(0, ns)]) * b  #  resulting hrs >= hbs ???
#         k_soil = np.array([vg.hydraulic_conductivity(hrs[j], sp) for j in range(0, ns)]) * b
        for j, _ in enumerate(hsb):
            d = hsb[j] - hrs[j]  # at least [1 cm]
            k_soil[j] = (vg.matric_flux_potential(hsb[j], sp) - vg.matric_flux_potential(hrs[j], sp)) / d

        for j in range(0, ns):
            hrs[j] = fsolve(lambda x:k_root[j] * (x - hx[j + 1]) + b * k_soil[j] * (x - hsb[j]), hrs[j])

        hrs = np.maximum(hrs, critP)

    print("bulk", hsb[1], "hrs0", hrs0[1], "hrs", hrs[1])
    # print(k_root[1], k_soil[1])

    # print(k_soil)
    seg_fluxes = r.segFluxes(0., hx, hrs, approx = False, cells = False)  # [cm3/day]

    # water for net flux
    soil_water = np.multiply(np.array(s.getWaterContent()), cell_volumes)

    water_domain.append(np.sum(soil_water))
    water_cell.append(soil_water[cci])
    min_rx.append(np.min(np.array(hx)))
    p1d.append(np.min(np.array(hrs)))

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

fig, (ax1, ax2) = plt.subplots(2, 1)

ax1.set_title("Water amount")
ax1.plot(np.linspace(0, sim_time, NT), np.array(water_domain), label = "water in domain")
ax1.legend()
ax1.set_xlabel("Time (days)")
ax1.set_ylabel("(cm3)")

ax2.set_title("Pressure")
# ax2.plot(np.linspace(0, sim_time, NT), np.array(collar), label="root collar soil cell")
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

print("fin")

