"""
functions for the aggregated approach
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/");  sys.path.append("../")


# import plantbox as pb  # CPlantBox
# # from xylem_flux import *  # root system Python hybrid solver
# # from rhizo_models import *
# import vtk_plot as vp
# import van_genuchten as vg
# 
# import numpy as np
# import matplotlib.pyplot as plt


def get_aggregated_params(r, rs_age, min_b, max_b, cell_number):
    """ returns aggregated root system parameters (krs, and per layer: suf_, kr_surf_, surf_, l_, a_)
        r
        rs_age 
        min_b
        max_b
        cell_number
    
        returns:
        krs [cm2/day]         root system conductivity 
        suf_ [1]              standard uptake fraction
        kr_surf_ [cm2/day]    kr times surface (summed up per layer)
        surf_ [cm2]           surface (summed up per layer)
        l_ [cm]               length (summed up per layer)
        a_ [cm]               mean layer radius 
    """
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

    surf_ = np.array(ana.distribution("surface", max_b[2], min_b[2], cell_number[2], False))
    l_ = np.array(ana.distribution("length", max_b[2], min_b[2], cell_number[2], False))
    a_ = np.divide(surf_, 2 * np.pi * l_)
    return krs, np.array(suf_), np.array(kr_surf_), surf_, l_, a_


def create_aggregated_rs(r, rs_age, min_b, max_b, cell_number):
    """ ceates a root system with one segment per layer, 
        connected to the root collar by an artificial segments,
        with adjsuted conductivities
    
        r
        rs_age 
        min_b
        max_b
        cell_number    
    """
    krs, suf_, kr_surf_, surf_, l_, a_ = get_aggregated_params(r, rs_age, min_b, max_b, cell_number)

    print("krs")
    print(krs)
    print("\nSUF", suf_.shape)
    print(suf_)
    print(list(suf_[0:100]))
    print("\nkr*surf", kr_surf_.shape)
    print(kr_surf_)
    print("\nsurface", surf_.shape)
    print(surf_)
    print("\nlength", l_.shape)
    print(l_)
    print("\nradius", a_.shape)
    print(a_)
    print("\n\n")

    n = int(cell_number[2])
    nodes = [pb.Vector3d(0, 0, 0)]  # maximal ns+1 nodes
    segs, radii = [], []  # maximal ns segments
    dx = (max_b[2] - min_b[2]) / cell_number[2]
    z_ = np.linspace(max_b[2] - dx / 2, min_b[2] + dx / 2, n)  # cell centers
    # print(z_)

    c = 0
    for i in range(0, n):
        if l_[i] > 0:
            nodes.append(pb.Vector3d(0, 0, z_[i]))  # node 2*i+1
            nodes.append(pb.Vector3d(l_[i], 0, z_[i]))  # node 2*i+2
            segs.append(pb.Vector2i(0, 2 * c + 1))  # artificial shoot segment
            radii.append(0.)
            segs.append(pb.Vector2i(2 * c + 1, 2 * c + 2))  # normal segment
            radii.append(a_[i])
            c += 1
    # print("number of segments", len(segs))

    rs = pb.MappedSegments(nodes, segs, radii)
    rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), cut = False)
    r2 = XylemFluxPython(rs)  # wrap the xylem
    r.test()  # sanity checks
    # z_ = np.linspace(0, -150, 150)
    # plt.plot(suf_[0:100], z_[0:100])
    # plt.show()

    ll = np.abs(z_)
    # print(ll)

    suf_krs = suf_ * krs
    # print(krs)
    # print(suf_krs)
    # print(kr_surf_)

    kr_up, kx_up = [], []
    for i in range(0, n):
        if l_[i] > 0:

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
                raise ValueError('create_aggregated_rs() no segment in layer')
            kx_up.append(1.e1)  # regular segment

    r2.setKrValues(kr_up)
    r2.setKxValues(kx_up)

    kr_ = np.array(kr_up)
    kx_ = np.array(kx_up)

    # print(np.array(kr_up))
    # print(np.array(kx_up))

    print("krs", krs)
    print("suf", np.min(suf_), np.max(suf_), np.sum(suf_))
    print("kr_up", np.min(kr_[1::2]), np.max(kr_[1::2]), np.mean(kr_[1::2]))
    print("kx_up", np.min(kx_[0::2]), np.max(kx_[0::2]), np.mean(kx_[0::2]))
    print("kx_up", kx_.shape)
    print("kx_up")
    print(list(kx_[0::2]))

    print(kr_surf_[0])
    print(suf_krs[0])

    # vp.plot_roots(pb.SegmentAnalyser(rs), "radius")
    return r2




def simulate_const(s,r, sra_table_lookup, trans, sim_time, dt):
    pass









if __name__ == "__main__":

    """ root system """
    min_b = [-7.5, -37.5, -110.]
    max_b = [7.5, 37.5, 0.]
    cell_number = [1, 1, 55]

    fname = "../../../../grids/RootSystem_verysimple2.rsml"
    r = XylemFluxPython(fname)
    rs_age = 78  # for calculating age dependent conductivities

    types = r.rs.subTypes  # simplify root types
    types = (np.array(types) >= 12) * 1  # all roots type 0, only >=12 are laterals type 1
    r.rs.subTypes = list(types)

    init_conductivities(r)
    # init_conductivities_const(r)
    # r.test()  # sanity checks

    # krs, suf_, kr_surf_, surf_, l_, a_  = get_aggregated_params(r, rs_age, min_b, max_b, cell_number)
    # z_ = np.linspace(0, -110, 55)
    # plt.plot(suf_, z_)
    # plt.show()
    create_aggregated_rs(r, rs_age, min_b, max_b, cell_number)

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
    #
    # r = create_aggregated_rs(r, 0., min_b, max_b, cell_number)
