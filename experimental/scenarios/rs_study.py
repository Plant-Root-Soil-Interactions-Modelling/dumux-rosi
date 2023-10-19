
""" 
Jan's new scenarios 

Analyse root system: visualize & Krs; surfaces and SUF per layer
"""

import sys; sys.path.append("../../../modules/"); sys.path.append("../../../../../CPlantBox/");  sys.path.append("../../../../../CPlantBox/src/python_modules")
sys.path.append("../../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../modules/fv/");
sys.path.append("../");

import plantbox as pb  # CPlantBox
from xylem_flux import *  # root system Python hybrid solver
import aggregated_rs as agg
import vtk_plot as vp
import matplotlib.pyplot as plt
import numpy as np
from rhizo_models import *

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

"""
Parameters  
"""

""" soil """
min_b = [-7.5, -37.5 / 2, -110.]
max_b = [7.5, 37.5 / 2, 0.]
cell_number = [1, 1, 55]  # [8, 38, 55]  # 2cm3
periodic = True  # check data first
fname = "../../../../grids/RootSystem_verysimple2.rsml"

""" root system """
trans = 0.5 * 15 * 75  # * 15 * 75  # average per day [cm3 /day] (sinusoidal)
wilting_point = -15000  # [cm]
predefined_growth = False  # root growth by setting radial conductivities
rs_age = 78  # initial root system age

""" 
Initialize xylem model 
"""
rs = RhizoMappedSegments(fname, -150000, 1, 1.5, 0)
rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)
r = XylemFluxPython(rs)
# r = XylemFluxPython(fname)

# simplify root types to 0 for basal and 1 for laterals
types = r.rs.subTypes
types = (np.array(types) >= 12) * 1  # all roots type 0, only >=12 are laterals type 1
r.rs.subTypes = list(types)

""" sanity checks """
# r.test()  # sanity checks
ns = len(r.rs.segments)
print("number of segments", ns)
outer_r = r.rs.segOuterRadii()
inner_r = r.rs.radii
rho_ = np.divide(outer_r, np.array(inner_r))

""" 1. visualization """
print()
ana = pb.SegmentAnalyser(r.rs)
print("unconstrained")
print("min bounds", ana.getMinBounds())
print("max bounds", ana.getMaxBounds())
print()
# vp.plot_roots(ana, "subType")  # unconstrained vizualisation
ana.mapPeriodic(max_b[0] - min_b[0], max_b[1] - min_b[1])
print("periodic")
print("min bounds", ana.getMinBounds())
print("max bounds", ana.getMaxBounds())
print()
# vp.plot_roots(ana, "subType")  # unconstrained vizualisation

""" 2. surfaces """
fig, ax = plt.subplots(1, 1, figsize=(14, 12))
ax = [ax]
ana = pb.SegmentAnalyser(r.rs)
dz = (max_b[2] - min_b[2]) / cell_number[2]
z_ = np.linspace(-dz / 2, min_b[2] - dz / 2, cell_number[2])  # z - axis
rsd = ana.distribution("surface", 0., min_b[2], cell_number[2], False)
rsd2 = ana.distribution("surface", 0., min_b[2], cell_number[2], True)
ax[0].plot(rsd, z_, "-*", label="seg mid")
ax[0].plot(rsd2, z_, "-*", label="exact")
ax[0].set_title("Root surface per layer")
ax[0].set_xlabel("root surface (cm$^2$)")
ax[0].set_ylabel("depth (cm)")
ax[0].legend(loc="lower right")
plt.show()

fig, ax = plt.subplots(1, 1, figsize=(14, 12))
ax = [ax]
ana = pb.SegmentAnalyser(r.rs)
ana.mapPeriodic(max_b[0] - min_b[0], max_b[1] - min_b[1])
z_ = np.linspace(0, min_b[2], cell_number[2])  # z - axis
rsd = ana.distribution("surface", 0., min_b[2], cell_number[2], False)
rsd2 = ana.distribution("surface", 0., min_b[2], cell_number[2], True)
ax[0].plot(rsd, z_, "-*", label="seg mid")
ax[0].plot(rsd2, z_, "-*", label="exact")
ax[0].set_title("Root surface per layer")
ax[0].set_xlabel("root surface (cm$^2$)")
ax[0].set_ylabel("depth (cm)")
ax[0].legend(loc="lower right")
plt.show()

""" 3. Krs, SUF """


def plot_suf(r):
    fig, ax = plt.subplots(1, 1, figsize=(14, 12))
    ax = [ax]
    krs, _ = r.get_krs(0.)
    suf = r.get_suf(0.)
    ana = pb.SegmentAnalyser(r.rs)
    ana.addData("SUF", suf)
    d = np.array(ana.distribution("SUF", 0., min_b[2], cell_number[2], False))
    print()
    print("Krs", krs, "cm2/day")
    print("SUF ranges from", np.min(d[d > 0]), "to", np.max(d), "mean", np.mean(d[d > 0]), "sums up to", np.sum(d))
    ax[0].plot(d, z_, "-*", label="total")
    max_type = int(np.max(ana.data["subType"]))
    for i in range(0, max_type + 1):
        ana2 = pb.SegmentAnalyser(ana)  # copy
        ana2.filter("subType", i)
        d = ana2.distribution("SUF", 0., min_b[2], cell_number[2], False)
        ax[0].plot(d, z_, "-*", label="type {:g}".format(i))
    ax[0].set_title("SUF per layer")
    ax[0].set_xlabel("SUF (1)")
    ax[0].legend(loc="lower right")
    plt.show()


print("original conductivities")
kr_const = 0.00018100042  
kx_const = 0.173
r.setKr([kr_const])
r.setKx([kx_const])
plot_suf(r)

print("Low conductivities")
kx_const = 0.1 * 0.173
r.setKx([kx_const])
plot_suf(r)

print("High conductivities")
kx_const = 10 * 0.173
r.setKx([kx_const])
plot_suf(r)
