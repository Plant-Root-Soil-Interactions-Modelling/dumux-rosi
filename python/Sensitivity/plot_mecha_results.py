import os
import shutil

import numpy as np
import matplotlib.pyplot as plt


def plot_index(ind):
    aerenchyma_percentage = [0.00, 0.07, 0.15, 0.22, 0.30]
    aerenchyma_number_of_areas = [1, 4, 8, 12]
    # 4*5
    cortex_layer = [ 3, 5, 7, 9, 11, 13, 15 ]
    cortex_diam = [ 0.0225, 0.0350, 0.0475 ]
    # 7*3 = 21
    stele_diam = [ 0.063, 0.056, 0.049 ]
    stele_celldiam = [ 0.004 ]
    # 3
    xylem_poles = [ 4 ]
    xylem_diam = [ 0.04, 0.06, 0.08 ]
    # 3
    i = ind
    ind_ap = int(ind // (3 * 3 * 21 * 4))
    ind = ind % (3 * 3 * 21 * 4)
    ind_anoa = int(ind // (3 * 3 * 21))
    ind = ind % (3 * 3 * 21)
    ind_cl = int(ind // (3 * 3 * 3))
    ind = ind % (3 * 3 * 3)
    ind_cd = int(ind // (3 * 3))
    ind = ind % (3 * 3)
    ind_sl = int(ind // 3)
    ind_xd = int(ind % 3)

    print("index ", i, ": aer p", ind_ap, "aer number", ind_anoa, "cortex layer", ind_cl, "cortex diam", ind_cd, "stele diam", ind_sl, "xylem diam", ind_xd)

    return aerenchyma_percentage[ind_ap]


def main():

    mecha_results_path = "/home/daniel/Dropbox/granar/mecha_results"

    npy_ind = []
    npy_data = []

    for file in os.listdir(mecha_results_path):
        if file.endswith(".npy"):  # Check if the file is a .npy file
            ind = file.split("shiny")[1].split(".")[0]
            file_path = os.path.join(mecha_results_path, file)
            npy_ind.append(ind)
            npy_data.append(np.load(file_path))

    print(npy_data[-1])  # kx, kr, a?
    npy_data = np.array(npy_data)
    print(npy_data.shape)

    ykx_ = npy_data[:, 0, 0]
    ykr_ = npy_data[:, 0, 1]
    ya_ = npy_data[:, 0, 2]
    kx_ = npy_data[:, 2, 0]
    kr_ = npy_data[:, 2, 1]
    a_ = npy_data[:, 2, 2]

    kx = np.array([float(x) for x in kx_])
    kr = np.array([float(x) for x in kr_])
    a = np.array([float(x) for x in a_])
    ykx = np.array([float(x) for x in ykx_])
    ykr = np.array([float(x) for x in ykr_])
    ya = np.array([float(x) for x in ya_])
    npy_ind = np.array([float(x) for x in npy_ind])
    p = []

    for i in npy_ind:
        p.append(plot_index(i))

    print("kx", np.min(kx), np.max(kx), np.mean(kx), "young", np.min(ykx), np.max(ykx), np.mean(ykx))
    print("kr", np.min(kr), np.max(kr), np.mean(kr), "young", np.min(ykr), np.max(ykr), np.mean(ykr))
    print("a", np.min(a), np.max(a), np.mean(a), "young", np.min(ya), np.max(ya), np.mean(ya))

    print("file saved")
    np.savez("mecha_results", kx = kx, kr = kr, a = a, ykx = ykx, ykr = ykr, ya = ya, aer_p = np.array(p), index = npy_ind)  # copy stable versions to data/

    plt.hist(kr)
    plt.show()


if __name__ == "__main__":
    main()
