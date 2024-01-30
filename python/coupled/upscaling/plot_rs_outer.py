"""
     Outer radii visualized by 3D tubeplot
"""
import sys; sys.path.append("../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../python/modules/");
sys.path.append("../../../../CPlantBox"); sys.path.append("../../../../CPlantBox/src")

import plot_rootsystem as pr

from functional.Perirhizal import *
import visualisation.vtk_plot as vp

import matplotlib.pyplot as plt
import numpy as np


def plot_outer_vtk(method, dim, plant, soil, outer_method, plot_time = 13.):

    """ recreate root system """
    outer_radii = pr.get_outer_radius(plant, dim, outer_method)
    r, mapping, length, a, surf, z = pr.get_rootsystem(plant, dim)

    # """ load data """
    # fname = method + "_" + plant + "_" + dim + "_" + soil + "_" + outer_method
    # data = np.load("results/hsr_"+fname + ".npy")

    # """ add 3d soil surface density """
    # peri = PerirhizalPython(r)
    # sn = np.prod(cell_number)
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

    print("outer_radii min", np.min(outer_radii), "max", np.max(outer_radii),
          "median", np.median(outer_radii), "mean", np.mean(outer_radii), "std", np.std(outer_radii))

    ana = pb.SegmentAnalyser(r.rs.mappedSegments())
    outer_radii = np.minimum(outer_radii, 3.)  # limit for visualisation
    ana.addData("outer_r", outer_radii)
    # vp.plot_roots(ana, "outer_r")

    ana.addData("radius", outer_radii)
    # vp.plot_mesh(grid, "surface_density")
    # vp.plot_mesh_cuts(grid, "surface_density")
    vp.plot_roots(ana, "outer_r")
    # vp.plot_roots_and_mesh(ana, "outer_r", grid, "surface_density", True, width[0], width[1])
    # plt.hist(outer_radii, bins = 100, rwidth = 0.9, align = 'mid')
    fname = method + "_" + plant + "_" + dim + "_" + soil + "_" + outer_method
    ana.write(fname + ".vtp")


if __name__ == "__main__":

    method = "sra"
    plant = "springbarley"
    dim = "1D"
    soil = "hydrus_loam"
    outer_method = "voronoi"
    plot_outer_vtk(method, dim, plant, soil, outer_method, plot_time = 13.)  # , 4., 6.
    # plt.tight_layout()
    # plt.savefig('hsr_bin_' + plant + "_" + soil + dim + "z.png")
    # plt.show()

