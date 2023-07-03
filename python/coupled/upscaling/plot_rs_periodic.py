""" 3d surface densities"""
import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import plantbox as pb
from functional.Perirhizal import *
import visualisation.vtk_plot as vp
import scenario_setup as scenario

import vtk
import numpy as np
import matplotlib.pyplot as plt


def field_rs(m, n, x_width, y_width, param_file, rs_age):
    """ makes an 2m+1 x 2n+1 array of root systems """
    random_seed = 1  # make them identical

    rs = pb.RootSystem()
    rs.setSeed(random_seed)
    rs.readParameters("data/" + param_name)
    params = rs.getRootRandomParameter()
    if param_name == "spring_barley_CF12.xml":  # some parameter mods
        # params[2].lmax *= 2
        params[1].theta = 1.31  # why is the tap root not always 0?
    rs.initialize()
    rs.simulate(rs_age, True)

    mid_rs_periodic = pb.SegmentAnalyser(rs)
    mid_rs_periodic.mapPeriodic(x_width, y_width)

    for j in range(0, 2 * m + 1):
        for i in range(0, 2 * n + 1):
            new_rs = pb.SegmentAnalyser(mid_rs_periodic)
            nodes = new_rs.nodes
            for n_ in nodes:
                n_.x += x_width * (i - m)
                n_.y += y_width * (j - m)
            if i == 0 and j == 0:
                all_rs = pb.SegmentAnalyser(new_rs)
            all_rs.addSegments(new_rs)
            print(len(all_rs.nodes))

    return all_rs, mid_rs_periodic


def field_bb(m, n, x_width, y_width):
    """ makes an 2m+1 x 2n+1 array of bounding boxes """
    all_ = []
    for j in range(0, 2 * m + 1):
        for i in range(0, 2 * n + 1):
            bb = pb.SDF_PlantBox(x_width, y_width, 110.)
            bbt = pb.SDF_RotateTranslate(bb, pb.Vector3d(x_width * (i - m), y_width * (j - n), 0))
            if j == m and i == n:  # mid
               mid_ = bbt
            else:
                if j == m or i == n:
                    all_.append(bbt)
    return pb.SDF_Union(all_), mid_


if __name__ == "__main__":

    plant = "springbarley"

    if plant == "maize":
        min_b, max_b, cell_number = scenario.maize_("3D")
        param_name = "Zeamays_synMRI_modified.xml"
        rs_age = 8 * 7  # 56 days
        m, n = 1, 1  # 76*16  -> (*3) -> 228 * 48
    elif plant == "soybean":
        min_b, max_b, cell_number = scenario.soybean_("3D")
        param_name = "Glycine_max_Moraes2020_opt2_modified.xml"
        rs_age = 6 * 7  # 42 days
        m, n = 2, 2
    elif plant == "springbarley":
        min_b, max_b, cell_number = scenario.springbarley_("3D")
        param_name = "spring_barley_CF12.xml"
        rs_age = 7 * 7  # 49 days
        # m, n = 8, 8  # 13*3 -> (*17) -> 221*51
        m, n = 4, 4  # 13*3 -> (*9) -> 117*27

    all_rs, mid_rs_periodic = field_rs(m, n, max_b[0] - min_b[0], max_b[1] - min_b[1], param_name, rs_age)

    vp.plot_roots(all_rs, "subType")
    # vp.plot_roots(mid_rs_periodic, "subType")

    all_rs.write("results/" + plant + "_all.vtp")
    mid_rs_periodic.write("results/" + plant + "_periodic.vtp")

    all_, mid_ = field_bb(m, n, max_b[0] - min_b[0], max_b[1] - min_b[1])
    rs = pb.RootSystem()
    rs.setGeometry(mid_)
    rs.write("results/" + plant + "_mid.py")
    rs = pb.RootSystem()
    rs.setGeometry(all_)
    rs.write("results/" + plant + "_all.py")

    print("happy paraviewing...")

#
# print("\nplant")
# print("rs_age", rs_age)
# print()
# params = r.getRootRandomParameter()
# for p in params:
#     print("SubType", p.subType)
#     print("radius", p.a)
#     print("growth rate", p.r)
#     print()
#
# print()
#
# # A_dirichlet, Kr, kx0 = r_.doussan_system_matrix(rs_age)
# # print(type(r_.rs))
# # print("segs")
# # for i in range(0, 10):
# #     print(segs[i])
# # # print("trans", trans)
# # print("nodes")
# # for i in range(0, 10):1
# #     print(nodes[i])
# # print()
# # print("kx0", kx0)
# # print()
# # print(A_dirichlet[0:8, 0:8])
# # print()
# # print(A_dirichlet[:, 0])
# # r_.test()
#
# sn = np.prod(cell_number)
# peri = PerirhizalPython(r)
#
# """ add 3d soil surface density """
# sd = peri.get_density("surface")
# sd = -np.minimum(sd, 1.1)  # limit for visualisation
# grid = vp.uniform_grid(min_b, max_b, cell_number)
# cell_sd = vtk.vtkDoubleArray()
# cell_sd.SetName("surface_density")
# cell_sd.SetNumberOfValues(sn)
# for j in range(0, sn):
#     cell_sd.SetValue(j, sd[j])
# celldata = grid.GetCellData()
# celldata.AddArray(cell_sd)
#
# # outer_radii = peri.get_outer_radii_vo git push --set-upstream origin new-features-from-sensitivityronoi()
# outer_radii = peri.get_outer_radii("surface")
#
# print("outer_radii", np.min(outer_radii), np.max(outer_radii), "median", np.median(outer_radii), "mean", np.mean(outer_radii), np.std(outer_radii))
#
# ana = pb.SegmentAnalyser(r.mappedSegments())
# # outer_radii = -np.minimum(outer_radii, 1.)  # limit for visualisation
# ana.addData("outer_r", outer_radii)
# ana.addConductivities(r_, rs_age, 0.1, 0.014)
# ana.addAge(rs_age)
# ana.addCellIds(r)
#
# # vp.plot_mesh(grid, "surface_density")
# # vp.plot_mesh_cuts(grid, "surface_density")
# print()
# print("age", rs_age, "nodes", len(nodes))
# print()
#
# # vp.plot_roots(ana, "age")
# # vp.plot_roots(ana, "kr")
# # vp.plot_roots(ana, "kx")
# vp.plot_roots(ana, "subType")
#
# ana.write("results/" + plant + ".vtp")
#
# # width = max_b - min_b
# # vp.plot_roots_and_mesh(ana, "outer_r", grid, "surface_density", True, width[0], width[1])
#
# # plt.hist(outer_radii, bins = 40, rwidth = 0.9, align = 'mid')
# # plt.show()
#
# print("fin")
