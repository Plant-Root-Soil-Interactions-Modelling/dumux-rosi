""" figures for the icons (2D)"""
import sys; sys.path.append("../../build-cmake/cpp/python_binding/"); sys.path.append("../../python/modules/");
sys.path.append("../../../CPlantBox"); sys.path.append("../../../CPlantBox/src")

import plantbox as pb
from functional.Perirhizal import *

import numpy as np
import matplotlib.pyplot as plt
import random

from scipy.spatial import ConvexHull, Voronoi

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

""" parameters """


def generate_random_points(num_points, range_):
    points = []
    for _ in range(num_points):
        x = random.uniform(-range_, range_)
        y = random.uniform(-range_, range_)
        points.append(np.array([x, y]))
    return np.array(points)


def generate_regular_points(xnum_points, ynum_points, xrange, yrange):
    points = []
    for i in range(ynum_points):
        for j in range(xnum_points):
            x = (j + 0.5) * xrange / xnum_points - xrange / 2
            y = (i + 0.5) * yrange / ynum_points - yrange / 2
            points.append(np.array([x, y]))
    return np.array(points)


def flip_(nodes, center, axis):
    """ flips the nodes around the center according axis """
    # print(nodes.shape)
    # print(center.shape)
    # print(np.ones((nodes.shape[0], 1)).shape)
    n_ = nodes - np.ones((nodes.shape[0], 1)) @ center
    n_[:, axis] = -n_[:, axis]
    n_ = n_ + np.ones((nodes.shape[0], 1)) @ center
    return n_


def mirror(nodes):
    """ adds mirrored nodes to the 4 sides """
    width_ = np.array([[20, 20]])
    center_ = np.array([[0, 0]])
    nodes_surr = nodes
    fipped_n = [flip_(nodes, center_, i_) for i_ in range(0, 2)]  # flip them ...
    zeros = np.zeros((nodes.shape[0], 1))
    translate_ = np.ones((nodes.shape[0], 1)) @ width_
    trans0 = np.hstack((translate_[:, 0, np.newaxis], zeros))
    trans1 = np.hstack((zeros, translate_[:, 1, np.newaxis]))
    nodes_surr = np.vstack((nodes_surr, fipped_n[0] + trans0))  # add them
    nodes_surr = np.vstack((nodes_surr, fipped_n[0] - trans0))
    nodes_surr = np.vstack((nodes_surr, fipped_n[1] + trans1))
    nodes_surr = np.vstack((nodes_surr, fipped_n[1] - trans1))
    return nodes_surr


def vizualize(points):

    figure, axes = plt.subplots()
    axes.set_aspect(1)

    # Volumes
    vol = np.empty((nodes.shape[0]))
    vol[:] = np.nan
    for idx, reg_num in enumerate(vor.point_region):
        indices = vor.regions[reg_num]
        i_ = reg_num - 1
        if idx < nodes.shape[0]:
            if -1 in indices:  # some regions can be opened
                vol[idx] = np.nan
            else:
                vol[idx] = ConvexHull(vor.vertices[indices]).volume

    # Visualize regions
    lines = []
    reg_ = vor.regions
    ver_ = vor.vertices
    for reg in reg_:
        if -1 not in reg:
            for ind in range(-1, len(reg) - 1):
                lines.append([ver_[reg[ind]], ver_[reg[ind + 1]]])
    lines = np.array(lines)

    for i in range(0, lines.shape[0]):
        plt.plot([lines[i, 0, 0], lines[i, 1, 0]], [lines[i, 0, 1], lines[i, 1, 1]], "g:")

    # Visualize Perirhizal radii in 2D
    for i in range(0, points.shape[0]):
        a = np.sqrt(vol[i] / np.pi)
        # print(vol[i], a * a * np.pi)
        circ = plt.Circle((points[i, 0], points[i, 1]), a, fill = True, alpha = 0.1, linestyle = 'dotted', color = "red")
        axes.add_artist(circ)

    plt.plot(points[:, 0], points[:, 1], 'r*')
    plt.plot(ver_[:, 0], ver_[:, 1], 'b*')
    plt.xlim([-10, 10])
    plt.ylim([-10, 10])

    plt.show()


# Generate 8 random 3D points
points = generate_random_points(9, 8.)
print(points)
# [[ 4.04892037  5.56818163], [ 3.68849136  0.89770249], [ 2.40660187  1.01993938], [-1.87067736  1.70510243], [ 3.91117411 -6.93592641], [-0.9386731  -3.46114548], [ 0.79988892 -0.35288419], [ 2.62590918 -6.30607085]]

# [[-2.90345978  1.88259082]
 # [ 5.14163188 -4.20208878]
 # [-6.28865805 -6.32353276]
 # [ 6.9969697   6.5203306 ]
 # [ 6.34113581 -6.15953657]
 # [-5.17135285  7.72011243]
 # [ 4.86211936  4.5371164 ]
 # [ 6.61824728  1.09913438]
 # [-0.03276024  0.14013829]]

points = generate_regular_points(3, 3, 20, 20)
nodes = mirror(points)
vor = Voronoi(nodes)

vizualize(points)

