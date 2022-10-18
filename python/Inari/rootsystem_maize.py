
""" 
Analyse and choose (modify) maize root system parameters 
"""
import sys; sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb
import vtk_plot as vp

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

min_b = [-37.5, -7.5, -200.]  # Domain Mais: 60 cm Reihe, 10 cm Pflanzen
max_b = [37.5, 7.5, 0.]
cell_number = [75, 15, 200]  # 1 cm3
simtime = 0.1 * 95  # between 90-100 days

rs = pb.RootSystem()

# Open plant and root parameter from a file
path = "../../../CPlantBox/modelparameter/rootsystem/"
name = "Zeamays_synMRI"
# Zea_mays_1_Leitner_2010
# Zea_mays_2_Pagès_2014
# Zea_mays_3_Postma_2011,
# Zea_mays_4_Leitner_2014
# Zea_mays_5_Leitner_2014
# Zea_mays_6_Leitner_2014
# Zea_mays_Streuber_2020 (groß, viele shoot borne überirdisch)
# Zeamays_synMRI
rs.readParameters(path + name + ".xml")

srp = rs.getOrganRandomParameter(pb.OrganTypes.seed)
print(srp[0])
srp[0].delayB = 1.5
srp[0].maxB = 30

rrp = rs.getOrganRandomParameter(pb.OrganTypes.root)
for p in rrp:
    print("\nSubType", p.subType)
    print("Lmax", p.lmax)
    print("Radius", p.a)
    print("dx", p.dx, "cm")
    print("changed to 0.5 cm to be faster...")
    p.dx = 0.5  # probably enough
    p.dxMin = 0.05

rrp[2].lmax *= 2

#
# # print(rrp[0])
# rrp[1].theta = 0.8 * rrp[1].theta  # otherwise the initial peak in RLD is a bit too high
# # rrp[1].thetas = 0.1 * rrp[1].theta  # 10% std
# # print()
#
# rs.writeParameters(name + "_nobasals_modified" + ".xml")  # remember the modifications

# Initialize
print()
rs.setGeometry(pb.SDF_PlantBox(1.e6, 1.e6, np.abs(min_b[2])))
rs.initialize()

# Simulate
rs.simulate(simtime, True)

# Analyse
ana = pb.SegmentAnalyser(rs)

orders = np.array(rs.getParameter("subType"))
print("\nnumber of roots", len(rs.getRoots()))
print("types", np.sum(orders == 1), np.sum(orders == 2), np.sum(orders == 3), np.sum(orders == 4), np.sum(orders == 5))
print("number of nodes", len(ana.nodes))
print("number of segments", len(ana.segments))

print("\nunconfined")
print(ana.getMinBounds(), "-", ana.getMaxBounds())

# w = np.array(max_b) - np.array(min_b)
# ana.mapPeriodic(w[0], w[1])
# print("periodic")
# print(ana.getMinBounds(), "-", ana.getMaxBounds())

dz = 0.5
exact = True
domain_size = np.array(max_b) - np.array(min_b)
slice_volume = domain_size[0] * domain_size[1] * 1  # cm3
z_ = np.linspace(max_b[2] - dz, min_b[2] + dz, cell_number[2])
rld = np.array(ana.distribution("length", 0., min_b[2], cell_number[2], exact)) / slice_volume
rsd = np.array(ana.distribution("surface", 0., min_b[2], cell_number[2], exact)) / slice_volume

fig, ax = plt.subplots(1, 1, figsize = (10, 10))
ax = [ax]
ax[0].plot(rld, z_)
ax[0].set_xlabel("root length density [cm / cm3]")
ax[0].set_ylabel("depth [cm]")
# ax[1].plot(rsd, z_)
# ax[1].set_xlabel("root surface density [cm2 / cm3]")
# ax[1].set_ylabel("depth [cm]")
plt.tight_layout()
plt.savefig("results/" + name + "_rld" + ".png")
plt.show()

# Export final results
ana.write("results/maize.vtp")  # rs.write writes polylines, ana.write writes segments
rs.write("results/maize.rsml")

# # Plot, using vtk
# pd = vp.segs_to_polydata(ana, 1., ["creationTime", "radius", "organType"])
# plantActor, scalar_bar = vp.plot_roots(ana, "creationTime", name, False)
# iren = vp.render_window(plantActor, name, scalar_bar, pd.GetBounds())
# renWin = iren.GetRenderWindow()
# vp.write_png(renWin, name)
# print("saved", "results/" + name + "_rs" + ".png")
# iren.Start()

w = np.array(max_b) - np.array(min_b)
ana.mapPeriodic(w[0], w[1])
ana.write("results/maize_periodic.vtp")

#
# print("fin")
