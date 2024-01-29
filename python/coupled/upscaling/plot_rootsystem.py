""" 
   Retrieves root system properties for maize and srpingbarley
"""
import sys; sys.path.append("../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../python/modules/");
sys.path.append("../../../../CPlantBox"); sys.path.append("../../../../CPlantBox/src")

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


def get_rootsystem(rootsystem_name, dim):
    """
        rootsystem_name     "springbarley", or "maize"   
        dim                 dimension, for correct mapping of soil centers to segments (e.g. r.get_hs(sx))
        
        @return 
        r
        mapper
        length
        root radius
        surface 
        z - coordinate
    """

    r, rho, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = scenario.set_scenario(rootsystem_name, dim, -200, "hydrus_loam", "length")

    length = r.rs.segLength()
    a = np.array(r.rs.radii)
    # z = r.rs.getSegmentZ() # segment mids
    segs = r.get_segments()
    nodes = r.get_nodes()
    z_ = np.array([nodes[s[1], 2] for s in segs ])  # segment apical node z
    surf = 2 * np.multiply(a, length) * np.pi

    return r, mapping, length, a, surf, z_


def get_outer_radius(rootsystem_name, dim, outer_method):

    r, rho, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = scenario.set_scenario(rootsystem_name, dim, -200, "hydrus_loam", outer_method)

    if outer_method == "voronoi":
        print("voronoi", outer_method)
        sys.stdout.flush()
        outer_ = PerirhizalPython(r.rs).get_outer_radii_bounded_voronoi()
        if np.sum([np.isnan(outer_)]) > 0:
            print("set_scenario(): NaNs in get_outer_radii_bounded_voronoi are replaced by mean", np.sum([np.isnan(outer_)]))
            sys.stdout.flush()
            outer_mean = np.nanmean(outer_) * np.ones(outer_.shape)
            outer_[np.isnan(outer_)] = outer_mean[np.isnan(outer_)]
    else:
        print("other:", outer_method)
        sys.stdout.flush()
        outer_ = PerirhizalPython(r.rs).get_outer_radii(outer_method)

    return outer_


if __name__ == "__main__":

    """ parameters """
    rootsystem = "springbarley"
    dim = "1D"
    r, mapper, length, a, surf, z = get_rootsystem(rootsystem, dim)

    # fig, axes = plt.subplots(3, 1, figsize = (10, 18))
    fig, axes = plt.subplots(1, 1, figsize = (9, 8))

    plt.show()

    print("fin")
