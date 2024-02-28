"""
    transpiration plot (one column, number of rows as number of filenames)
    
    modify __main__ to select simualtion result    
"""
import sys; sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import numpy as np
import matplotlib.pyplot as plt
from functional.xylem_flux import sinusoidal2

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


if __name__ == "__main__":

    """
    Maize
    """
    plant = ["maize"] * 3
    soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    method = ["sra"] * 3  # Axx
    dim = ["3D"] * 3
    outer_method = ["voronoi"] * 3
    cup_AAA, cupW_AAA = get_cumulative(method, plant, dim, soil, outer_method, "(AAA)")
    dim = ["1D"] * 3  # xxB
    cup_AAB, cupW_AAB = get_cumulative(method, plant, dim, soil, outer_method, "(AAB)")
    outer_method = ["length"] * 3  # xBx
    cup_ABB, cupW_ABB = get_cumulative(method, plant, dim, soil, outer_method, "(ABB)")
    method = ["agg"] * 3  # Bxx
    cup_BBB, cupW_BBB = get_cumulative(method, plant, dim, soil, outer_method, "(BBB)")
    method = ["par"] * 3  # Bxx
    cup_CBB, cupW_CBB = get_cumulative(method, plant, dim, soil, outer_method, "(CBB)")

    print("\nMaize: percental error in cumulative uptake compared to AAA after 2 weeks")
    np.set_printoptions(precision = 0)
    print("AAA")
    print(-100.*(np.ones(np.shape(cup_AAB)) - np.divide(cup_AAB, cup_AAA)), "% AAB")
    print(-100.*(np.ones(np.shape(cup_ABB)) - np.divide(cup_ABB, cup_AAA)), "% ABB")
    print(-100.*(np.ones(np.shape(cup_BBB)) - np.divide(cup_BBB, cup_AAA)), "% BBB")
    print(-100.*(np.ones(np.shape(cup_CBB)) - np.divide(cup_CBB, cup_AAA)), "% CBB")
    print()
    print("AAB")
    print(-100.*(np.ones(np.shape(cup_ABB)) - np.divide(cup_ABB, cup_AAB)), "% ABB")
    print(-100.*(np.ones(np.shape(cup_BBB)) - np.divide(cup_BBB, cup_AAB)), "% BBB")
    print(-100.*(np.ones(np.shape(cup_CBB)) - np.divide(cup_CBB, cup_AAB)), "% CBB")
    print("ABB")
    print(-100.*(np.ones(np.shape(cup_BBB)) - np.divide(cup_BBB, cup_ABB)), "% BBB")
    print(-100.*(np.ones(np.shape(cup_CBB)) - np.divide(cup_CBB, cup_ABB)), "% CBB")
    print("BBB")
    print(-100.*(np.ones(np.shape(cup_CBB)) - np.divide(cup_CBB, cup_BBB)), "% CBB")

    """
    Spring Barley
    """
    plant = ["springbarley"] * 3
    soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    method = ["sra"] * 3  # Axx
    dim = ["3D"] * 3
    outer_method = ["voronoi"] * 3
    cup_AAA, cupW_AAA = get_cumulative(method, plant, dim, soil, outer_method, "(AAA)")
    dim = ["1D"] * 3  # xxB
    cup_AAB, cupW_AAB = get_cumulative(method, plant, dim, soil, outer_method, "(AAB)")
    outer_method = ["length"] * 3  # xBx
    cup_ABB, cupW_ABB = get_cumulative(method, plant, dim, soil, outer_method, "(ABB)")
    method = ["agg"] * 3  # Bxx
    cup_BBB, cupW_BBB = get_cumulative(method, plant, dim, soil, outer_method, "(BBB)")
    method = ["par"] * 3  # Bxx
    cup_CBB, cupW_CBB = get_cumulative(method, plant, dim, soil, outer_method, "(CBB)")

    print("\n\nSpring barley: percental error in cumulative uptake compared to AAA after 2 weeks")
    np.set_printoptions(precision = 0)
    print("AAA")
    print(-100.*(np.ones(np.shape(cup_AAB)) - np.divide(cup_AAB, cup_AAA)), "% AAB")
    print(-100.*(np.ones(np.shape(cup_ABB)) - np.divide(cup_ABB, cup_AAA)), "% ABB")
    print(-100.*(np.ones(np.shape(cup_BBB)) - np.divide(cup_BBB, cup_AAA)), "% BBB")
    print(-100.*(np.ones(np.shape(cup_CBB)) - np.divide(cup_CBB, cup_AAA)), "% CBB")
    print()
    print("AAB")
    print(-100.*(np.ones(np.shape(cup_ABB)) - np.divide(cup_ABB, cup_AAB)), "% ABB")
    print(-100.*(np.ones(np.shape(cup_BBB)) - np.divide(cup_BBB, cup_AAB)), "% BBB")
    print(-100.*(np.ones(np.shape(cup_CBB)) - np.divide(cup_CBB, cup_AAB)), "% CBB")
    print("ABB")
    print(-100.*(np.ones(np.shape(cup_BBB)) - np.divide(cup_BBB, cup_ABB)), "% BBB")
    print(-100.*(np.ones(np.shape(cup_CBB)) - np.divide(cup_CBB, cup_ABB)), "% CBB")
    print("BBB")
    print(-100.*(np.ones(np.shape(cup_CBB)) - np.divide(cup_CBB, cup_BBB)), "% CBB")

