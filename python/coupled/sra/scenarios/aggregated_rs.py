""" 
sets up the aggregated root system described by Krs, and (per Layer) SUF, mean radii, summed length  
"""

import sys; sys.path.append("../../../modules/"); sys.path.append("../../../../../CPlantBox/");  sys.path.append("../../../../../CPlantBox/src/python_modules")
sys.path.append("../");

import plantbox as pb  # CPlantBox
from xylem_flux import *  # root system Python hybrid solver
import vtk_plot as vp
import van_genuchten as vg

import numpy as np
import matplotlib.pyplot as plt


def init_conductivities(r):
    """ Hydraulic conductivities - for Jans scenarios """
    kr0 = np.array([[0., 0.000181], [8., 0.000181], [10, 0.0000648], [18, 0.0000648], [25, 0.0000173], [300, 0.0000173]])  # radial values [1/day]
    kr1 = np.array([[0., 0.000181], [10., 0.000181], [16, 0.0000173], [300, 0.0000173]])
    kr2 = np.array([[0., 0.000181], [10., 0.000181], [16, 0.0000173], [300, 0.0000173]])
    r.setKrTables([kr0[:, 1], kr1[:, 1], kr2[:, 1]], [kr0[:, 0], kr1[:, 0], kr2[:, 0]])  # values, age

    kx0 = np.array([[0., 0.0000864], [5., 0.00173], [12., 0.0295], [15., 0.0295], [20., 0.432], [300., 0.432]])  # axial values [cm3/day]
    kx1 = np.array([[0., 0.0000864], [5., 0.0000864], [10., 0.0000864], [12., 0.0006048], [20., 0.0006048], [23., 0.00173], [300., 0.00173]])
    kx2 = np.array([[0., 0.0000864], [5., 0.0000864], [10., 0.0000864], [12., 0.0006048], [20., 0.0006048], [23., 0.00173], [300., 0.00173]])
    r.setKxTables([kx0[:, 1], kx1[:, 1], kx2[:, 1]], [kx0[:, 0], kx1[:, 0], kx2[:, 0]])  # values, age


def init_conductivities_const(r):
    """ Hydraulic conductivities - for Jans scenarios, but constant """
    kr_const = 2.843148e-5 / (2 * np.pi * 0.05 * 0.5)  # in case of table look up, the values must agree, = 0.00018100042
    kr = np.array([[0., kr_const], [1e4, kr_const]])
    r.setKrTables([kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1]],
                  [kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0]])

    kx_const = 0.346 * 0.5  # [cm3/day] = 0.173
    kx = np.array([[0., kx_const], [1e4, kx_const]])
    r.setKxTables([kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1]],
                  [kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0]])  # values, age


def create_singleroot(ns = 100, l = 50 , a = 0.05):
    """ creates a single root with @param ns segments, length l, and radius a """
    radii = np.array([a] * ns)
    nodes = [pb.Vector3d(0, 0, 0)]
    segs = []
    dx = l / ns
    z_ = np.linspace(-dx, -l , ns)
    for i in range(0, ns):
        nodes.append(pb.Vector3d(0, 0, z_[i]))
        segs.append(pb.Vector2i(i, i + 1))
    rs = pb.MappedSegments(nodes, segs, radii)
    return XylemFluxPython(rs)


def get_aggregated_params(r, rs_age, min_b, max_b, cell_number):
    """ krs, and per layer: suf, mean radius, total length """
    ana = pb.SegmentAnalyser(r.rs)
    krs, _ = r.get_krs(rs_age)

    suf = r.get_suf(rs_age)  # SUF per layer
    print("suf", np.sum(suf))
    ana.addData("suf", suf)
    suf_ = ana.distribution("suf", max_b[2], min_b[2], cell_number[2], False)
    print("suf_", np.sum(suf_))

    kr_surf = []  # kr times surface summed up per layer
    segs = r.rs.segments
    nodeCTs = r.rs.nodeCTs
    subTypes = r.rs.subTypes
    lengths = r.rs.segLength()
    radii = r.rs.radii
    for i, s in enumerate(segs):
        age = rs_age - nodeCTs[s.y]
        t = subTypes[i]
        kr_surf.append(2 * radii[i] * np.pi * lengths[i] * r.kr_f(age, t))  # [cm2 / day]
    ana.addData("kr_surf", kr_surf)
    kr_surf_ = ana.distribution("kr_surf", max_b[2], min_b[2], cell_number[2], False)

    surf_ = ana.distribution("surface", max_b[2], min_b[2], cell_number[2], False)
    l_ = ana.distribution("length", max_b[2], min_b[2], cell_number[2], False)
    a_ = np.divide(surf_, l_)

    return krs, suf_, kr_surf_, surf_, l_, a_  # ALL NEEDED ????


def create_aggregated_rs(r, rs_age, min_b, max_b, cell_number):
    """  one segment per layer connected by artificial segments"""
    krs, suf_, kr_surf_, surf_, l_, a_ = get_aggregated_params(r, rs_age, min_b, max_b, cell_number)
    l = 0.5

    ns = 2 * cell_number[2]
    radii = np.array([ 0. if i % 2 == 0 else a_[int(i / 2)] for i in range(0, ns)])  # to have correct outer rhizosphere radius, artificial segments are set to 0
    nodes = [pb.Vector3d(0, 0, 0)]  # ns+1 nodes
    segs = []  # ns segments
    dx = (max_b[2] - min_b[2]) / cell_number[2]
    z_ = np.linspace(max_b[2] - dx / 2, min_b[2] + dx / 2, cell_number[2])
    # print(z_)
    for i in range(0, int(ns / 2)):
        nodes.append(pb.Vector3d(0, 0, z_[i]))  # node 2*i+1
        nodes.append(pb.Vector3d(l, 0, z_[i]))  # node 2*i+2
        segs.append(pb.Vector2i(0, 2 * i + 1))  # artificial shoot segment
        segs.append(pb.Vector2i(2 * i + 1, 2 * i + 2))  # normal segment
    rs = pb.MappedSegments(nodes, segs, radii)
    rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), cut = False)
    r2 = XylemFluxPython(rs)  # wrap the xylem    # init_conductivities_const(r)
    # r.test()  # sanity checks
    # z_ = np.linspace(0, -110, 55)
    # plt.plot(suf_, z_)
    # plt.show()

    ll = np.abs(z_)
    # print(ll)

    suf_krs = np.array(suf_) * krs
    # print(krs)
    # print(suf_krs)
    # print(kr_surf_)

    kr_up, kx_up = [], []
    for i in range(0, int(ns / 2)):
        kr_up.append(0.)  # artificial segment
        if surf_[i] > 0:
            kr_up.append(kr_surf_[i] / surf_[i])  # mean layer kr [1/day]
        else:  # no segments in layer
            kr_up.append(0.)  # regular segment
        if kr_surf_[i] - suf_krs[i] > 0:
            kx_up.append(ll[i] * (suf_krs[i] * kr_surf_[i]) / (kr_surf_[i] - suf_krs[i]))  # artificial segment
        else:  # no segments in layer
            kx_up.append(0.)
        kx_up.append(1.)  # regular segment
    r2.setKrValues(kr_up)
    r2.setKxValues(kx_up)

    kr_ = np.array(kr_up)
    kx_ = np.array(kx_up)
    # print(kr_)
    # print(kx_)
    print("krs", krs)
    print("suf", np.min(suf_), np.max(suf_), np.sum(suf_))
    print("kr_up", np.min(kr_[1::2]), np.max(kr_[1::2]), np.mean(kr_[1::2]))
    print("kx_up", np.min(kx_[0::2]), np.max(kx_[0::2]), np.mean(kx_[0::2]))

    # vp.plot_roots(pb.SegmentAnalyser(rs), "radius")
    return r2


if __name__ == "__main__":

    """ root system """
    # min_b = [-7.5, -37.5, -110.]
    # max_b = [7.5, 37.5, 0.]
    # cell_number = [1, 1, 55]
    #
    # fname = "../../../../grids/RootSystem_verysimple2.rsml"
    # r = XylemFluxPython(fname)
    # rs_age = 78  # for calculating age dependent conductivities
    #
    # types = r.rs.subTypes  # simplify root types
    # types = (np.array(types) >= 12) * 1  # all roots type 0, only >=12 are laterals type 1
    # r.rs.subTypes = list(types)
    #
    # init_conductivities(r)
    # # init_conductivities_const(r)
    # # r.test()  # sanity checks
    #
    # # krs, suf_, kr_surf_, surf_, l_, a_  = get_aggregated_params(r, rs_age, min_b, max_b, cell_number)
    # # z_ = np.linspace(0, -110, 55)
    # # plt.plot(suf_, z_)
    # # plt.show()
#     create_aggregated_rs(r, rs_age, min_b, max_b, cell_number)

    """ single root """
    min_b = [-7.5, -37.5, -50.]
    max_b = [7.5, 37.5, 0.]
    cell_number = [1, 1, 100]
    r = create_singleroot()
    init_conductivities_const(r)

    krs, suf_, kr_surf_, surf_, l_, a_ = get_aggregated_params(r, 0., min_b, max_b, cell_number)
    z_ = np.linspace(0, min_b[2], cell_number[2])
    plt.plot(suf_, z_)
    plt.show()

    r = create_aggregated_rs(r, 0., min_b, max_b, cell_number)
