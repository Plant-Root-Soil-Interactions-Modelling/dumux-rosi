"""
    plots results of local sensitivity analysis
    (see run_SA)
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import numpy as np
import matplotlib.pyplot as plt

import run_SA as sa


def start_index(ind, ranges):
    s = 0
    for i in range(0, ind):
        if len(ranges) > 1:
            s += len(ranges[i])
    return s


""" def SA """
file_name = "local_soybean_"
# file_name = "local_soybean_conductivities_"
# file_name = "local_soybean_tropisms_"
# file_name = "local_soybean_radii_"
# file_name = "local_singleroot_conductivities64_"

path = "results/"
# not_xlog = []

analysis_time = 87.5  # 87.5  # days  (day 60 is fine already)

names, ranges = sa.read_ranges(path + file_name)

print(len(names), "variables:")
print(names)
for i, r in enumerate(ranges):
    print(names[i], r)
print()

# SOYBEAN
# names[2] = "lmax, primaries"  # seminal and tap
# names[3] = "lmax, first order"
# names[4] = "lmax, second order"
# names[5] = "insertion angle, seminal roots"
# names[6] = "root radii scale"
# names[7] = "seminal root count"

alldata = np.load(path + file_name + "1" + ".npz")
trans_ = alldata["act_trans"]
times = alldata["times"]
print("Simulation time from", min(times), "to ", max(times), "days")

ind_ = np.argwhere(times > analysis_time)
if len(ind_) > 0:
    ind = ind_[0][0]
    ind += 1
    ind10 = ind // 10  # TODO check
else:
    ind = -1
    ind10 = -1

ind = -1
ind10 = -1

print("index", ind)
print("Plotting for day", times[ind])
dt_ = np.diff(times[:ind])

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
fig, ax = plt.subplots(3, 4, figsize = (16, 16))
final = []
finaln = []

ac = 0
for lind in range(0, len(names)):

    file_start_ind = 2 + start_index(lind, ranges)  # 1 is initial simulation
    sa_len = len(ranges[lind])
    print("\nproducing sub-plot", names[lind], "from", file_start_ind, "to", file_start_ind + sa_len - 1)  # to make sure it works

    if sa_len > 1:
        trans = np.zeros((sa_len,))
        vol = np.zeros((sa_len,))
        krs = np.zeros((sa_len,))
        nitrate = np.zeros((sa_len,))
        for k in range(0, sa_len):
            try:
                # trans_ = np.load(path + "transpiration_" + file_name + str(file_start_ind + k) + ".npy")
                alldata = np.load(path + file_name + str(file_start_ind + k) + ".npz")
                trans_ = alldata["act_trans"]
                trans[k] = -np.sum(np.multiply(trans_[:ind - 1], dt_))
                # print("trans[k]", trans[k])
            except:
                trans[k] = np.nan
                print("skipping file", file_name + str(file_start_ind + k))
            try:
                # alldata = np.load(path + file_name + str(file_start_ind + k) + ".npz")
                vol_ = alldata["vol"]
                # vol_ = np.load(path + "vol_" + file_name + str(file_start_ind + k) + ".npy")
                vol_ = vol_[:, ind10]
                vol[k] = np.sum(vol_)  # sum over sub-types
            except:
                vol[k] = np.nan
            try:
                # krs_ = np.load(path + "krs_" + file_name + str(file_start_ind + k) + ".npy")
                krs_ = alldata["krs"]
                krs[k] = krs_[ind10]
            except:
                krs[k] = np.nan
            try:
                n_ = np.load(path + "nitrate_" + file_name + str(file_start_ind + k) + ".npy")
                nitrate[k] = np.sum(np.multiply(n_[:ind - 1], dt_))
                # print(nitrate[k], 'g')
            except:
                nitrate[k] = np.nan

        trans = trans / trans[sa_len // 2]  # nondimensionalize
        nitrate = nitrate / nitrate[sa_len // 2]  # nondimensionalize
        vol = vol / vol[sa_len // 2]  # nondimensionalize
        krs = krs / krs[sa_len // 2]  # nondimensionalize

        ax.flat[ac].plot(np.array(ranges[lind]), trans, '*-', label = "water")  #  / ranges[lind][sa_len // 2]
        print("ranging", np.min(ranges[lind]), np.max(ranges[lind]), len(ranges[lind]))

        # ax.flat[ac].plot(ranges[lind], nitrate, '*-', label = "nitrate")
        # ax.flat[ac].plot(ranges[lind], vol, '-.', label = "volume")
        # ax.flat[ac].plot(ranges[lind], krs, ':', label = "krs")

        final.append(trans[-1] / vol[-1])
        finaln.append(nitrate[-1] / vol[-1])

        # center
        x = ranges[lind][sa_len // 2]
        ax.flat[ac].plot([x], [1.], 'r*')

        ax.flat[ac].set_title(names[lind])
        # ax.flat[ac].set_ylim(0.8, 1.2)
        # ax.flat[ac].set_yscale('log', base = 2)
        # if not lind in not_xlog:

        ax.flat[ac].set_xscale('log', base = 2)
        ac += 1

plt.tight_layout(pad = 4.)
plt.show()

print("final")
print(names)
print(final)

names = np.array(names)
final = np.array(final)
I = np.argsort(final)
print("sorted")
print(names[I])
print(final[I])
