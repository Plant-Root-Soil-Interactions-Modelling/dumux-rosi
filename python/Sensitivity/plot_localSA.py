"""
plots results of local sensitivity analysis
"""
import sys; sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

import numpy as np
import matplotlib.pyplot as plt

import run_SA as sa

""" def SA """
file_name = "local_SA_test"
path = "results/"

# p = np.array([1.* 2 ** x for x in np.linspace(-2., 2., 9)])
# kr = 1.e-4
# kx = 1.e-3
# sa_lists = sa.make_lists(kr * p , kx * p , p, p, p, p, p, p, p, [2, 3, 4, 5])
# sa_len = len(sa_lists[0])  # assume they have the same size for all parameters
#
# name = ["kr", "kx", "lmax0", "lmax1", "lmax2", "theta0", "r0", "r1", "a", "src"]
# p1 = np.array([1.* 2 ** x for x in np.linspace(-1., 1., 9)])
# p2 = np.array([1.* 2 ** x for x in np.linspace(-2., 2., 9)])
# theta_ = np.linspace(-1,1, 9)
# # ranges = [kr * p , kx * p , p, p, p, p, p, p, p, [2, 3, 4, 5]]
# ranges = [p2, p2, p1, p1, p1, theta_, p1, p1, p1, [2, 3, 4, 5]]

""" font sizes """
SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title

""" make plots """
fig, ax = plt.subplots(3, 3, figsize = (16, 16))

for i in range(0, 3):
    for j in range(0, 3):
        lind = i * 3 + j
        print("\nproducing sub-plot", name[lind])
        file_start_ind = 2 + lind * sa_len
        trans = np.zeros((sa_len,))
        vol = np.zeros((sa_len,))
        krs = np.zeros((sa_len,))
        for k in range(0, sa_len):
            try:
                trans_ = np.load(path + "transpiration_" + file_name + str(file_start_ind + k) + ".npy")
                trans[k] = np.sum(trans_[1,:])
                print(trans[k])
            except:
                trans[k] = np.nan
                print("skipping file", file_name + str(file_start_ind + k))
            try:
                vol_ = np.load(path + "vol_" + file_name + str(file_start_ind + k) + ".npy")
                vol_ = vol_[:, -1]
                vol[k] = np.sum(vol_)
            except:
                vol[k] = np.nan
            try:
                krs_ = np.load(path + "krs_" + file_name + str(file_start_ind + k) + ".npy")
                krs[k] = krs_[-1]
            except:
                krs[k] = np.nan
        trans = trans / trans[sa_len // 2]  # nondimensionalize
        vol = vol / vol[sa_len // 2]  # nondimensionalize
        krs = krs / krs[sa_len // 2]  # nondimensionalize
        ax[i, j].plot(ranges[lind], trans, label = "uptake")
        ax[i, j].plot(ranges[lind], vol, '-.', label = "volume")
        ax[i, j].plot(ranges[lind], krs, ':', label = "krs")
        ax[i, j].plot([1.], [1.], 'r*')
        ax[i, j].legend()
        ax[i, j].set_title(name[lind])
        ax[i, j].set_ylim(0.5, 2)
        ax[i, j].set_yscale('log', base = 2)
        ax[i, j].set_xscale('log', base = 2)

plt.tight_layout(pad = 4.)
plt.show()
