"""
    prints simulation wall times  
"""
import sys; sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");

import numpy as np
# from scenario_setup import *


def print_timings(job_list):
    """ evalutate simulation wall times """
    names, times = [], []
    for job in job_list:
        method = job[0]
        plant = job[1]
        dim = job[2]
        soil = job[3]
        outer_method = job[4]
        name = "results/time_" + method + "_" + plant + "_" + dim + "_" + soil + "_" + outer_method
        names.append(name)
        try:
            time = np.load(name + ".npy")
        except:
            print("********* not found **********", name)
            time = 1.
        print(name, time)
        times.append(time)
    return names, times


def make_list():
    jobs = []

    method = ['sra', "agg", "par"]
    plant = ['springbarley', 'maize']
    dim = ["1D", "3D"]
    soil = ['hydrus_loam', 'hydrus_clay', 'hydrus_sandyloam']
    outer_radius = ['voronoi']  # 'length', 'surface', 'volume',

    print("Creating", len(method) * len(plant) * len(dim) * len(soil) * len(outer_radius), "simulations")
    print()

    for m in method:
        for p in plant:
            for d in dim:
                for s in soil:
                    for o in outer_radius:
                        jobs.append([m, p, d, s, o])
    return jobs


if __name__ == "__main__":

    # job_list = make_list()
    # names, timings = print_timings(job_list)
    # print()

    # c = 0
    # for l in range(0, 3):  # methods
    #     for i, p in enumerate(['springbarley', 'maize']):
    #         for j, d in enumerate(["1D", "3D"]):
    #             for k in range(0, 3):  # soils
    #                     base_index = 6 * i + 3 * j + k
    #                     print(names[c + base_index], ": ", 1 / (timings[c + base_index] / timings[base_index]))
    #         print()
    #     c += 12

    # for i, p in enumerate(['springbarley', 'maize']):
    #     for j, d in enumerate(["1D", "3D"]):  # dimesnions
    #         for k in range(0, 3):  # soils
    #                 base_index = 6 * i + 3 * j + k
    #                 if j == 0:
    #                     print(names[base_index], ": ", 1 / (timings[base_index] / timings[base_index + 3]))
    #     print()

    # """ 2D speed ups """
    # soils = ['hydrus_loam', 'hydrus_clay', 'hydrus_sandyloam']
    # for i in range(0, 3):
    #     n, t2d = print_timings([['sra', 'maize', '2D', soils[i], 'voronoi']])
    #     n, t3d = print_timings([['sra', 'maize', '3D', soils[i], 'voronoi']])
    #     print("speed up", t3d[0] / t2d[0])

    """ total speed ups """
    soils = ['hydrus_loam', 'hydrus_clay', 'hydrus_sandyloam']
    for i in range(0, 3):
        n, t1d = print_timings([['agg', 'maize', '1D', soils[i], 'voronoi']])
        n, t3d = print_timings([['sra', 'maize', '3D', soils[i], 'voronoi']])
        n, t2d = print_timings([['sra', 'maize', '2D', soils[i], 'voronoi']])
        print("speed up", t3d[0] / t1d[0])
        print("speed up", t3d[0] / t2d[0])
        n, t1d = print_timings([['agg', 'maize', '1D', soils[i], 'voronoi']])
        n, t3d = print_timings([['sra', 'maize', '3D', soils[i], 'voronoi']])
        print("speed up", t3d[0] / t1d[0])

    """ total speed ups """
    soils = ['hydrus_loam', 'hydrus_clay', 'hydrus_sandyloam']
    for i in range(0, 3):
        n, t1d = print_timings([['agg', 'springbarley', '1D', soils[i], 'voronoi']])
        n, t3d = print_timings([['sra', 'springbarley', '3D', soils[i], 'voronoi']])
        print("speed up", t3d[0] / t1d[0])
        n, t1d = print_timings([['agg', 'springbarley', '1D', soils[i], 'voronoi']])
        n, t3d = print_timings([['sra', 'springbarley', '3D', soils[i], 'voronoi']])
        print("speed up", t3d[0] / t1d[0])

