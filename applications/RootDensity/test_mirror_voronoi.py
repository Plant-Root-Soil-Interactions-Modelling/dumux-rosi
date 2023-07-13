""" 
    tests Perirhizal.mirror()
    
    saves vtk files for paraview visualisation 
    
    (a) rootsystem, (b) nodes within a cell, (c) voronoi cells 
"""
import sys; sys.path.append("../../build-cmake/cpp/python_binding/"); sys.path.append("../../python/modules/");
sys.path.append("../../../CPlantBox"); sys.path.append("../../../CPlantBox/src")

import plantbox as pb
from functional.Perirhizal import *
import visualisation.vtk_plot as vp
import scenario_setup as scenario

import vtk
import numpy as np
import matplotlib.pyplot as plt

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


def get_outer_radii_horizons(rootsystem):

    """ setup """
    if rootsystem == "Soybean":
        soil_, table_name, min_b, max_b, _, area, Kc = scenario.soybean(0)  # 0 = envirotype
        xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file
        cell_number = [38, 2, 150]  # (2cm)^3
        simtime = 42  # 87.5
    elif rootsystem == "Maize":
        soil_, table_name, min_b, max_b, _, area, Kc = scenario.maize(0)  # 0 = envirotype
        xml_name = "data/Zeamays_synMRI_modified.xml"  # root growth model parameter file
        simtime = 56  #
        cell_number = [38, 8, 159]  # (2cm)^3
    elif rootsystem == "Spring Barley":
        soil_, table_name, min_b, max_b, _, area, Kc = scenario.springbarley(0)  # 0 = envirotype
        xml_name = "data/spring_barley_CF12.xml"  # root growth model parameter file
        simtime = 49  # 95
        # cell_number = [7, 2, 75]  # (2cm)^3
        # cell_number = [7, 2, 75]  # (2cm)^3
        cell_number = [1, 1, 1]  # (2cm)^3
    else:
        print("get_outer_radii_horizons() unknown rootsystem name", rootsystem)
        raise

    """ simulate """
    s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1)
    r_ = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
    r = r_.rs
    r.initialize(False)
    r.simulate(simtime, False)

    """ post processing ... """
    peri = PerirhizalPython(r)
    ms = peri.ms  # rename
    print(min_b, max_b)
    min_b = ms.minBound
    max_b = ms.maxBound
    print(min_b, max_b)
    width = np.array([max_b.x, max_b.y, max_b.z]) - np.array([min_b.x, min_b.y, min_b.z])

    """ root architecture outputs """
    ana = pb.SegmentAnalyser(r.mappedSegments())
    ana.write("rootsystem.vtp")
    ana.mapPeriodic(width[0], width[1])
    ana.write("rootsystem_periodic.vtp")

    """ Perirhizal analysis """
    i, j, k = 0, 0, 0  # 140  # Cell index (for vizualisation) [3, 1, 71]
    nodes_ = ms.nodes
    nodes = np.array([[n.x, n.y, n.z] for n in nodes_])  # to numpy array

    print("making periodic", width)
    print("nodes", nodes.shape)
    nodes = peri.make_periodic_(nodes, width)
    print("periodic nodes", nodes.shape)
    print()

    print("\nget_bounds")
    min_, max_ = peri.get_bounds(i, j, k)
    width_ = max_ - min_
    width_ = np.expand_dims(width_, axis = 0)
    center_ = min_ + width_ / 2.
    print("cell width", width_)
    print("cell center", center_)
    print()

    print("\nmirror")
    n, ni = peri.mirror(nodes, i, j, k)

    grid0 = peri.get_node_mesh(n)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName('mirror_nodes.vtu')  # contains the mirrored nodes
    writer.SetInputData(grid0)
    writer.Write()

    grid = peri.get_voronoi_mesh_(n, ni.shape[0])  # Voronoi
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName('mirror_test_fulldomain.vtu')
    writer.SetInputData(grid)
    writer.Write()


rootsystem = "Spring Barley"  # Maize, Soybean, Spring Barley
get_outer_radii_horizons(rootsystem)  # outer0
print("fin")
