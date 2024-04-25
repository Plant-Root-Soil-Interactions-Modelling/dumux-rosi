"""
    transpiration plot (one column, number of rows as number of filenames)
    
    modify __main__ to select simualtion result    
"""
import sys; sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import numpy as np
import matplotlib.pyplot as plt

import timings_table as timings

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


def get_cumulative(method, plant, dim, soil, outer_method, label):
    """ cumulative transpiration over simulation time and half simulation time"""
    fnames = []
    for i in range(0, len(plant)):
        name = method[i] + "_" + plant[i] + "_" + dim[i] + "_" + soil[i] + "_" + outer_method[i]
        fnames.append("results/transpiration_" + method[i] + "_" + plant[i] + "_" + dim[i] + "_" + soil[i] + "_" + outer_method[i])

    n = len(fnames)
    data = [np.load(n_ + ".npy")  for n_ in fnames]
    cup_ = []
    cup2_ = []
    for i in range(0, n):
        t = data[i][0]
        y = data[i][2]
        dt = np.diff(t)
        so = np.array(y)
        cup = np.cumsum(np.multiply(so[:-1], dt))
        # print(i, "cumulative uptake " + label, cup[-1], "cm3")
        cup_.append(cup[-1])
        cup2_.append(cup[cup.shape[0] // 2])

    return np.array(cup_), np.array(cup2_)


def get_error_table(method, plant, dim, soil, outer_method, dimB = "1D"):

    cup_AAA, cupW_AAA = get_cumulative(method, plant, dim, soil, outer_method, "(AAA)")
    dim = [dimB] * 3  # xxB
    cup_AAB, cupW_AAB = get_cumulative(method, plant, dim, soil, outer_method, "(AAB)")
    outer_method = ["length"] * 3  # xBx
    cup_ABB, cupW_ABB = get_cumulative(method, plant, dim, soil, outer_method, "(ABB)")
    method = ["agg"] * 3  # Bxx
    cup_BBB, cupW_BBB = get_cumulative(method, plant, dim, soil, outer_method, "(BBB)")
    method = ["par"] * 3  # Bxx
    cup_CBB, cupW_CBB = get_cumulative(method, plant, dim, soil, outer_method, "(CBB)")

    print(plant[0] + ": percental error in cumulative uptake after 2 weeks (xxB = " + dimB + ")\n")
    np.set_printoptions(precision = 0)
    print("AAA")
    c1_aab = -100.*(np.ones(np.shape(cup_AAB)) - np.divide(cup_AAB, cup_AAA))
    print(c1_aab, "% AAB")
    c1_abb = -100.*(np.ones(np.shape(cup_ABB)) - np.divide(cup_ABB, cup_AAA))
    print(c1_abb, "% ABB")
    c1_bbb = -100.*(np.ones(np.shape(cup_BBB)) - np.divide(cup_BBB, cup_AAA))
    print(c1_bbb, "% BBB")
    c1_cbb = -100.*(np.ones(np.shape(cup_CBB)) - np.divide(cup_CBB, cup_AAA))
    print(c1_cbb, "% CBB")
    print("AAB")
    c2_abb = -100.*(np.ones(np.shape(cup_ABB)) - np.divide(cup_ABB, cup_AAB))
    print(c2_abb, "% ABB")
    c2_bbb = -100.*(np.ones(np.shape(cup_BBB)) - np.divide(cup_BBB, cup_AAB))
    print(c2_bbb, "% BBB")
    c2_cbb = -100.*(np.ones(np.shape(cup_CBB)) - np.divide(cup_CBB, cup_AAB))
    print(c2_cbb, "% CBB")
    print("ABB")
    c3_bbb = -100.*(np.ones(np.shape(cup_BBB)) - np.divide(cup_BBB, cup_ABB))
    print(c3_bbb, "% BBB")
    c3_cbb = -100.*(np.ones(np.shape(cup_CBB)) - np.divide(cup_CBB, cup_ABB))
    print(c3_cbb, "% CBB")
    print("BBB")
    c4_cbb = -100.*(np.ones(np.shape(cup_CBB)) - np.divide(cup_CBB, cup_BBB))
    print(c4_cbb, "% CBB\n\n")

    table0 = np.zeros((5, 3 * 5))
    table0[1:, 0:3] = np.array([c1_aab, c1_abb, c1_bbb, c1_cbb])
    table0[2:, 3:6] = np.array([c2_abb, c2_bbb, c2_cbb])
    table0[3:, 6:9] = np.array([c3_bbb, c3_cbb])
    table0[4:, 9:12] = np.array([c4_cbb])

    return table0


def get_error_table_abs(method, plant, dim, soil, outer_method, dimB = "1D"):

    cup_AAA, cupW_AAA = get_cumulative(method, plant, dim, soil, outer_method, "(AAA)")
    dim = [dimB] * 3  # xxB
    cup_AAB, cupW_AAB = get_cumulative(method, plant, dim, soil, outer_method, "(AAB)")
    outer_method = ["length"] * 3  # xBx
    cup_ABB, cupW_ABB = get_cumulative(method, plant, dim, soil, outer_method, "(ABB)")
    method = ["agg"] * 3  # Bxx
    cup_BBB, cupW_BBB = get_cumulative(method, plant, dim, soil, outer_method, "(BBB)")
    method = ["par"] * 3  # Bxx
    cup_CBB, cupW_CBB = get_cumulative(method, plant, dim, soil, outer_method, "(CBB)")
    print(plant[0] + ": percental error in cumulative uptake after 2 weeks (xxB = " + dimB + ")\n")
    np.set_printoptions(precision = 0)
    print("AAA")
    c1_aab = cup_AAB - cup_AAA
    print(c1_aab, "cm3 AAB")
    c1_abb = cup_ABB - cup_AAA
    print(c1_abb, "cm3 ABB")
    c1_bbb = cup_BBB - cup_AAA
    print(c1_bbb, "cm3 BBB")
    c1_cbb = cup_CBB - cup_AAA
    print(c1_cbb, "cm3 CBB")
    print("AAB")
    c2_abb = cup_ABB - cup_AAB
    print(c2_abb, "cm3 ABB")
    c2_bbb = cup_BBB - cup_AAB
    print(c2_bbb, "cm3 BBB")
    c2_cbb = cup_CBB - cup_AAB
    print(c2_cbb, "cm3 CBB")
    print("ABB")
    c3_bbb = cup_BBB - cup_ABB
    print(c3_bbb, "cm3 BBB")
    c3_cbb = cup_CBB - cup_ABB
    print(c3_cbb, "cm3 CBB")
    print("BBB")
    c4_cbb = cup_CBB - cup_BBB
    print(c4_cbb, "cm3 CBB\n\n")
    table0 = np.zeros((5, 3 * 5))
    table0[1:, 0:3] = np.array([c1_aab, c1_abb, c1_bbb, c1_cbb])
    table0[2:, 3:6] = np.array([c2_abb, c2_bbb, c2_cbb])
    table0[3:, 6:9] = np.array([c3_bbb, c3_cbb])
    table0[4:, 9:12] = np.array([c4_cbb])
    return table0


def get_error_table_maize(dimB, abs):
    plant = ["maize"] * 3
    soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    method = ["sra"] * 3  # Axx
    dim = ["3D"] * 3
    outer_method = ["voronoi"] * 3
    if abs:
        return get_error_table_abs(method, plant, dim, soil, outer_method, dimB)
    else:
        return get_error_table(method, plant, dim, soil, outer_method, dimB)


def get_error_table_springbarley(abs):
    plant = ["springbarley"] * 3
    soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    method = ["sra"] * 3  # Axx
    dim = ["3D"] * 3
    outer_method = ["voronoi"] * 3
    if abs:
        return get_error_table_abs(method, plant, dim, soil, outer_method)
    else:
        return get_error_table(method, plant, dim, soil, outer_method)


if __name__ == "__main__":

    abs = True  # absolute or relative errors

    table_error = get_error_table_maize("1D", abs)
    table_speedup = timings.table_maize("1D")
    np.savetxt("table_maize_1D.csv", table_error + table_speedup, delimiter = ";")

    table_error = get_error_table_maize("2D", abs)
    table_speedup = timings.table_maize("2D")
    np.savetxt("table_maize_2D.csv", table_error + table_speedup, delimiter = ",")

    table_error3 = get_error_table_springbarley(abs)
    table_speedup3 = timings.table_springbarley()
    np.savetxt("table_springbarley_1D.csv", table_error3 + table_speedup3, delimiter = ";")

