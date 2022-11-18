"""
plots results of local sensitivity analysis
"""
import sys; sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

import numpy as np
import matplotlib.pyplot as plt

import run_SA as sa

""" def SA """
file_name = "local_SA_const"
path = "results/"
p = np.array([1.* 2 ** x for x in np.linspace(-2., 2., 9)])
kr = 1.e-4
kx = 1.e-3
sa_lists = sa.make_lists(kr * p , kx * p , p, p, p, p, p, p, p, [2, 3, 4, 5])
sa_len = len(sa_lists[0])  # assume they have the same size for all parameters
p2 = np.linspace(-2., 2., 9)

name = ["kr", "kx", "lmax0", "lmax1", "lmax2", "theta0", "r0", "r1", "a", "src"]

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
        print("producing sub-plot", name[lind])
        file_start_ind = 2 + lind * sa_len
        trans = np.zeros((sa_len,))
        vol = np.zeros((sa_len,))
        krs = np.zeros((sa_len,))
        for k in range(0, sa_len):
            trans_ = np.load(path + "transpiration_" + file_name + str(file_start_ind + k) + ".npy")
            trans[k] = np.sum(trans_[1,:])
            vol_ = np.load(path + "vol_" + file_name + str(file_start_ind + k) + ".npy")
            vol_ = vol_[:, -1]
            vol[k] = np.sum(vol_)
            krs_ = np.load(path + "krs_" + file_name + str(file_start_ind + k) + ".npy")
            krs[k] = krs_[-1]
        trans = trans / trans[sa_len // 2]  # nondimensionalize
        vol = vol / vol[sa_len // 2]  # nondimensionalize
        krs = krs / krs[sa_len // 2]  # nondimensionalize
        ax[i, j].plot(p2, trans, label = "uptake")
        ax[i, j].plot(p2, vol, label = "volume")
        ax[i, j].plot(p2, krs, label = "krs")
        ax[i, j].plot([0.], [1.], 'r*')
        ax[i, j].legend()
        ax[i, j].set_title(name[lind])

plt.show()
