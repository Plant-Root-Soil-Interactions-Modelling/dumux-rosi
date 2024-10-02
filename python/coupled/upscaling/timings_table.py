"""
    prints simulation wall times  
"""
import sys; sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");

import numpy as np
# from scenario_setup import *


def get_timings(job_list):
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
        # print(name, time)
        times.append(time)
    return times[0]


def get_speedup_table(plant, dim, method, outer, dimB = "1D"):
    soils = ['hydrus_loam', 'hydrus_clay', 'hydrus_sandyloam']
    t_AAA = np.array([get_timings([[method, plant, dim, soils[i], outer]]) for i in range(0, 3)])
    dim = dimB
    t_AAB = np.array([get_timings([[method, plant, dim, soils[i], outer]]) for i in range(0, 3)])
    outer = 'length'
    t_ABB = np.array([get_timings([[method, plant, dim, soils[i], outer]]) for i in range(0, 3)])
    method = 'agg'
    t_BBB = np.array([get_timings([[method, plant, dim, soils[i], outer]]) for i in range(0, 3)])
    method = 'par'
    t_CBB = np.array([get_timings([[method, plant, dim, soils[i], outer]]) for i in range(0, 3)])
    print(plant + ": speed ups for 2 weeks of simulation time (xxB = " + dimB + ")\n")
    np.set_printoptions(precision = 0)
    print("AAA")
    print(np.divide(t_AAA, t_AAB), "AAB")
    print(np.divide(t_AAA, t_ABB), "ABB")
    print(np.divide(t_AAA, t_BBB), "BBB")
    print(np.divide(t_AAA, t_CBB), "CBB")
    row_AAA = np.array([np.divide(t_AAA, t_AAB), np.divide(t_AAA, t_ABB), np.divide(t_AAA, t_BBB), np.divide(t_AAA, t_CBB)])
    print("AAB")
    print(np.divide(t_AAB, t_ABB), "ABB")
    print(np.divide(t_AAB, t_BBB), "BBB")
    print(np.divide(t_AAB, t_CBB), "CBB")
    row_AAB = np.array([np.divide(t_AAB, t_ABB), np.divide(t_AAB, t_BBB), np.divide(t_AAB, t_CBB)])
    print("ABB")
    print(np.divide(t_ABB, t_BBB), "BBB")
    print(np.divide(t_ABB, t_CBB), "CBB")
    row_ABB = np.array([np.divide(t_ABB, t_BBB), np.divide(t_ABB, t_CBB)])
    print("BBB")
    print(np.divide(t_BBB, t_CBB), "CBB\n\n")
    row_BBB = np.array([np.divide(t_BBB, t_CBB)])
    table0 = np.zeros((5, 3 * 5))
    table0[0, 3:] = row_AAA.flat
    table0[1, 6:] = row_AAB.flat
    table0[2, 9:] = row_ABB.flat
    table0[3, 12:] = row_BBB.flat

    return table0


def table_maize(dimB):
    plant = 'maize'
    dim = '3D'
    method = 'sra'
    outer = 'voronoi'
    return get_speedup_table(plant, dim, method, outer, dimB)


def table_springbarley():
    plant = 'springbarley'
    dim = '3D'
    method = 'sra'
    outer = 'voronoi'
    return get_speedup_table(plant, dim, method, outer)


if __name__ == "__main__":
    table_maize("1D")
    table_springbarley()

