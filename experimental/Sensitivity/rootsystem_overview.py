
""" 
Analyse and create the root system we want to analyse 

which root system do we pick? and whats the RLD, RSD in our simulation domain 
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


def try_rootsystem(path, name, simtime, min_b, max_b, cell_number, show):

    # Simulate
    rs = pb.RootSystem()
    rs.readParameters(path + name + ".xml")
    # rs.setGeometry(pb.SDF_PlantBox(1.e6, 1.e6, np.abs(min_b[2])))
    rs.initialize()
    rs.simulate(simtime, True)
    rs.write("results/" + name + ".vtp")
    rs.write("results/" + name + ".rsml")

    # Analyse
    width = np.array(max_b) - np.array(min_b)
    ana = pb.SegmentAnalyser(rs)
    orders = np.array(rs.getParameter("subType"))
    print("\nnumber of roots", len(rs.getRoots()))
    print("types", np.sum(orders == 1), np.sum(orders == 2), np.sum(orders == 3), np.sum(orders == 4), np.sum(orders == 5))
    print("number of nodes", len(ana.nodes))
    print("number of segments", len(ana.segments))
    print("\nunconfined")
    print(ana.getMinBounds(), "-", ana.getMaxBounds())
    ana.mapPeriodic(width[0], width[1])

    # RLD
    dz = 0.5
    exact = True
    slice_volume = width[0] * width[1] * 1  # cm3
    z_ = np.linspace(max_b[2] - dz, min_b[2] + dz, cell_number[2])
    rld = np.array(ana.distribution("length", 0., min_b[2], cell_number[2], exact)) / slice_volume
    # rsd = np.array(ana.distribution("surface", 0., min_b[2], cell_number[2], exact)) / slice_volume
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
    if show:
        plt.show()

    # save root system as png
    pd = vp.segs_to_polydata(ana, 1., ["creationTime", "radius", "subType", "organType"])
    plantActor, scalar_bar = vp.plot_roots(ana, "subType", name, False)
    iren = vp.render_window(plantActor, name, scalar_bar, pd.GetBounds())
    renWin = iren.GetRenderWindow()
    vp.write_png(renWin, name)
    print("saved", "results/" + name + "_rs" + ".png")
    if show:
        iren.Start()


path = "../../../CPlantBox/modelparameter/rootsystem/"

""" Maize parameter sets in our CPlantbox library """
# min_b = [-37.5, -7.5, -200.]  # Domain Mais: 60 cm Reihe, 10 cm Pflanzen
# max_b = [37.5, 7.5, 0.]
# cell_number = [75, 15, 200]  # 1 cm3
# simtime = 95  # between 90-100 days
# # Zea_mays_1_Leitner_2010
# # Zea_mays_2_Pagès_2014
# # Zea_mays_3_Postma_2011,
# # Zea_mays_4_Leitner_2014
# # Zea_mays_5_Leitner_2014
# # Zea_mays_6_Leitner_2014
# # Zea_mays_Streuber_2020 (groß, viele shoot borne überirdisch)
# # Zeamays_synMRI
# names = ["Zea_mays_1_Leitner_2010", "Zea_mays_2_Pagès_2014", "Zea_mays_3_Postma_2011",
#          "Zea_mays_4_Leitner_2014", "Zea_mays_5_Leitner_2014", "Zea_mays_6_Leitner_2014",
#          "Zea_mays_Streuber_2020", "Zeamays_synMRI"]
# name = "Zeamays_synMRI"

""" Soy parameter sets in our CPlantbox library """
min_b = [-19, -2.5, -200.]  # Domain [38 cm Reihe, 6 cm Pflanzen]
max_b = [19, 2.5, 0.]
cell_number = [38, 5, 200]  # 1 cm3
simtime = 87.5  # between 75-100 days
# Glycine_max_Moraes2020_opt1
# Glycine_max_Moraes2020_opt2
# Glycine_max_Moraes2020
# Glycine_max
names = ["Glycine_max_Moraes2020_opt1", "Glycine_max_Moraes2020_opt2", "Glycine_max_Moraes2020", "Glycine_max"]
name = "Glycine_max"

# try_rootsystem(path, name, simtime, min_b, max_b, cell_number, True)

for n in names:
    try_rootsystem(path, n, simtime, min_b, max_b, cell_number, False)

print("fin")

