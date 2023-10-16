""" 
    sets up the parallel root system described by 
    Krs, and (per Layer) SUF, mean radii, summed length  
    
    get_aggregated_params() calculates aggregated parameters 
    create_parallel_rs() creates the parallel root system with root conductivities
"""

import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import plantbox as pb  # CPlantBox
from rhizo_models import *
from functional.PlantHydraulicParameters import PlantHydraulicParameters  # Doussan solver
from functional.PlantHydraulicModel import PlantHydraulicModel  # Doussan solver
from functional.Perirhizal import PerirhizalPython

import functional.van_genuchten as vg
import visualisation.vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt


def get_aggregated_params(r, rs_age, outer_method):
    """ returns krs, and per layer: suf_, kr_surf_, surf_, l_, a_    
        krs [cm2/day]         root system conductivity 
        suf_ [1]              standard uptake fraction
        kr_surf_ [cm2/day]    kr times surface (summed up per layer)
        surf_ [cm2]           surface (summed up per layer)
        l_ [cm]               length (summed up per layer)
        a_ [cm]               mean layer radius 
    """

    try:
        per = PerirhizalPython(r.rs.mappedSegments())  # if it is a MappedRootSystem we need the MappedSegments
    except:
        per = PerirhizalPython(r.rs)  # otherwise it is a MappedSegments instance

    krs, _ = r.get_krs(rs_age)
    print("krs", krs)

    suf = r.get_suf()  # SUF per layer
    suf_ = per.aggregate(suf)  # ana.distribution("suf", max_b[2], min_b[2], cell_number[2], False)
#     print("suf_", suf_.shape, np.sum(suf_))

    # r.test()

    kr_surf = []  # kr times surface summed up per layer
    segs = r.rs.segments
    nodeCTs = r.rs.nodeCTs
    subTypes = r.rs.subTypes
    lengths = np.array(r.rs.segLength())
    radii = np.array(r.rs.radii)
    for i, s in enumerate(segs):
        age = rs_age - nodeCTs[s.y]
        kr_surf.append(2 * radii[i] * np.pi * lengths[i] * r.params.kr_f(i, age, subTypes[i], 2))  # [cm2 / day] # kr_f(int segment_index, double age, int subType, int organType)
    kr_surf_ = per.aggregate(kr_surf)
    surf = np.multiply(2 * np.pi * radii, lengths)
    surf_ = per.aggregate(surf)
    l_ = per.aggregate(lengths)
    a_ = np.divide(surf_, 2 * np.pi * l_)
    if outer_method == "voronoi":
        outer_r = per.get_outer_radii_bounded_voronoi()
    else:
        outer_r = per.get_outer_radii(outer_method)
    outer_r_ = per.average(outer_r)
    # print("krs")
    # print(krs)
    # print("\nSUF", suf_.shape)
    # print(suf_)
    # print(list(suf_[0:100]))
    # print("\nkr*surf", kr_surf_.shape)
    # print(kr_surf_)
    # print("\nsurface", surf_.shape)
    # print(surf_)
    # print("\nlength", l_.shape)
    # print(l_)
    # print("\nradius", a_.shape)
    # print(a_)
    # print("\n\n")
    # dd

    return krs, suf_, kr_surf_, surf_, l_, a_, outer_r_


def create_parallel_rs(r, rs_age, cell_centers, min_b, max_b, cell_number, outer_method):
    """  one segment per layer connected by artificial segments"""

    r.update(rs_age)  # prepare hydraulic model

    krs, suf_, kr_surf_, surf_, l_, a_, outer_r_ = get_aggregated_params(r, rs_age, outer_method)

    n = cell_centers.shape[0]
    nodes = [pb.Vector3d(0, 0, 0), pb.Vector3d(0, 0, -0.1), pb.Vector3d(0, 0, -0.2)]  # maximal ns+1 nodes
    segs = []  # maximal ns segments
    radii = [0.1, 0.1]
    outer_r = [1. , 1.]
    segs.append(pb.Vector2i(0, 1))
    segs.append(pb.Vector2i(1, 2))

    c = 0
    for i in range(0, n):
        if kr_surf_[i] > 0:  # l_[i] > 0:
            x, y, z = cell_centers[i, 0], cell_centers[i, 1], cell_centers[i, 2]
            nodes.append(pb.Vector3d(x, y, z))  # node 2*i+3
            nodes.append(pb.Vector3d(x + l_[i], y, z))  # node 2*i+4
            segs.append(pb.Vector2i(2, 2 * c + 3))  # artificial shoot segment
            radii.append(0.)
            outer_r.append(1.)
            segs.append(pb.Vector2i(2 * c + 3, 2 * c + 4))  # normal segment
            radii.append(a_[i])
            outer_r.append(outer_r_[i])
            c += 1

    rs = pb.MappedSegments(nodes, segs, radii)
    rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), cut = False)

    params = PlantHydraulicParameters()
    r2 = PlantHydraulicModel("Doussan", rs, params)
    # r2.test()  # sanity checks

    print("cell_centers", cell_centers.shape, n)
    ll = np.abs(cell_centers[:, 2])  ###################################################
    suf_krs = suf_ * krs
    print("krs", krs)
    print("ll", ll.shape, ll)
    # print("suf_krs", suf_krs)
    # print("kr_surf_", kr_surf_)

    kr_up, kx_up = [0., 0.], [1.e1, 1.e1]
    for i in range(0, n):
        if kr_surf_[i] > 0:  # l_[i] > 0:

            kr_up.append(0.)  # artificial segment
            if surf_[i] > 0:
                kr_up.append(kr_surf_[i] / surf_[i])  # mean layer kr [1/day]
            else:  # no segments in layer
                kr_up.append(0.)  # regular segment

            if kr_surf_[i] - suf_krs[i] > 0:
                # kx_up.append(ll[i] * kx_up_[i])
                kx_up.append((ll[i] * suf_krs[i] * kr_surf_[i]) / (kr_surf_[i] - suf_krs[i]))  # artificial segment ll[i] *
                # Kxupscale=Krs*SFF*Krupscale/(Krupscale-Krs*SUF));  mit Kxupscale*(Hx-Hcollar)=Q
            else:  # no segments in layer
                print("warning", kr_surf_[i] - suf_krs[i], "<0")
                print("i", i, "kr_surf_", kr_surf_[i], "suf_krs", suf_krs[i], "l_", l_[i], "surf_", surf_[i])
                print("krs", krs, "suf", suf_[i])
                kx_up.append(1.e1)
                raise ValueError('create_parallel_rs(): could not calculate kx')
            kx_up.append(1.e2)  # regular segment

    r2.params.setKrValues(kr_up)
    r2.params.setKxValues(kx_up)
    #
    r2.update(rs_age)
    kr_ = np.array(kr_up)
    kx_ = np.array(kx_up)
    print("krs", krs)
    krs_new, _ = r2.get_krs(rs_age)
    print("new krs", krs_new)

    print("suf", np.min(suf_), np.max(suf_), np.sum(suf_))
    print("kr_up", np.min(kr_[1::2]), np.max(kr_[1::2]), np.mean(kr_[1::2]))
    print("kx_up", np.min(kx_[0::2]), np.max(kx_[0::2]), np.mean(kx_[0::2]))
    print("kx_up", kx_.shape, "kr_up", kr_.shape, "segs", len(segs))
    print("kx_up")
    # vp.plot_roots(pb.SegmentAnalyser(r2.rs), "radius")

    return r2, outer_r


if __name__ == "__main__":

    """ root system TODO set soil for aggregated params ..."""
    min_b = [-7.5, -37.5, -110.]
    max_b = [7.5, 37.5, 0.]
    cell_number = [1, 1, 55]

    fname = "../../../grids/RootSystem_verysimple2.rsml"
    params = PlantHydraulicParameters()
    r = PlantHydraulicModel("Doussan", fname, params)
    rs_age = 78  # for calculating age dependent conductivities

    types = r.rs.subTypes  # simplify root types
    types = (np.array(types) >= 12) * 1  # all roots type 0, only >=12 are laterals type 1
    r.rs.subTypes = list(types)

    kr00 = 0.  # artificial shoot
    kx00 = 1.e3  # artificial shoot
    kr_const_ = 1.73e-4  # [1/day]
    kx_const_ = 4.32e-2  # [cm3/day]
    params.setKr([[kr00, kr_const_, kr_const_, kr_const_, kr_const_, kr_const_, kr_const_]], [], 0.)
    params.setKx([[kx00, kx_const_, kx_const_, kx_const_, kx_const_, kx_const_, kx_const_]], [])

    r.update(rs_age)
    krs, suf_, kr_surf_, surf_, l_, a_, outer_r_ = get_aggregated_params(r, rs_age, "voronoi")
    z_ = np.linspace(0, -110, 55)
    plt.plot(suf_, z_)
    plt.show()

    # """ single root """
    # min_b = [-7.5, -37.5, -50.]
    # max_b = [7.5, 37.5, 0.]
    # cell_number = [1, 1, 100]
    # r = create_singleroot()
    # init_conductivities_const(r)
    #
    # krs, suf_, kr_surf_, surf_, l_, a_ = get_aggregated_params(r, 0., min_b, max_b, cell_number)
    # z_ = np.linspace(0, min_b[2], cell_number[2])
    # plt.plot(suf_, z_)
    # plt.show()
