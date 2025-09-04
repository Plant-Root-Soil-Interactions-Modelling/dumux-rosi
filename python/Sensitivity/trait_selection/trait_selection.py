import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src");
sys.path.append("../")

# import global_optimization_tools as got

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.cluster.vq import vq, whiten, kmeans
from sklearn.preprocessing import MinMaxScaler


def scatter_1D(nameX, index_y, all, nameC = None, scale = False):
    """ scatter plots one or multiple pairs of nameX vs nameY """
    if isinstance(nameX, str):
        nameX = [nameX]

    n = np.ceil(np.sqrt(len(nameX)))
    # print(n, np.sqrt(len(nameX)))
    fig, ax = plt.subplots(int(len(nameX) // n + 1 * (len(nameX) % n > 0)), int(n), figsize = (16, 16))
    if n == 1:
        ax = np.array([ax])
    ax_ = ax.flat

    slopes = []
    for i, _ in enumerate(nameX):

        x = np.array(all[:, i])
        y = np.array(all[:, index_y])
        # # print("x", len(x), type(x))
        # # print("y", len(y), type(y))
        # I = np.logical_and(~np.isnan(x), ~np.isnan(y))
        # x = x[I]
        # y = y[I]

        if scale:
            x = (np.array(x) - np.min(x)) / (np.max(x) - np.min(x))  # scale
            y = (np.array(y) - np.min(y)) / (np.max(y) - np.min(y))  # scale

        I = np.argsort(x)  # for regression

        if nameC is not None:
            ax_[i].scatter(x[I], y[I], c = nameC, alpha = 0.5)
        else:
            ax_[i].scatter(x[I], y[I])

        res = stats.linregress(x[I], y[I])
        ax_[i].plot(x[I], res.intercept + res.slope * x[I], 'k', label = 'fitted')
        # print(nameX[i], "slope", res.slope)
        slopes.append(res.slope)

        ax_[i].set_xlabel(nameX[i])
        ax_[i].set_ylabel(nameX[index_y])

    plt.tight_layout()

    return slopes


def distortion_plot(data_):
    """ kmeans distortion versus number of clusters (for kmeans) """
    plt.title("N = {:g}".format(data_.shape[0]))
    dist_ = []
    for i in range(0, 15):
        centers, dist = kmeans(data_, i + 1)
        dist_.append(dist)
    plt.plot(range(1, 16), dist_)
    plt.ylabel("distortion")
    plt.xlabel("number of clusters")
    plt.show()


file_name = "Root_Traits_Of_Interest.xlsx"
df = pd.read_excel(file_name)

#
# Figure 1
#
key_ = list(df.keys())[1:]  # [9:]  # skip the first three
print(key_)

all = []
for k in key_:
    x = df[k]
    x = (np.array(x) - np.min(x)) / (np.max(x) - np.min(x))
    all.append(x)
all = np.array(all).transpose()
all_clean = all[~np.isnan(all).any(axis = 1)]
print("data", all.shape, all_clean.shape)
distortion_plot(all_clean)

# add cluser indices as color
cluster_n = 12
centers, dist = kmeans(all_clean, cluster_n)
cluster_indices, _ = vq(all_clean, centers)

#
# Figure 2
#
key_ = list(df.keys())[1:]  # [9:]  # skip the first three
print(key_)

for i, k in enumerate(key_):

    target_names = [k]
    # key_.extend(["kr1", "kx1", "ykr1", "ykx1", "a1", "kr2", "kx2", "ykr2", "ykx2", "a2", "kr3", "kx3", "ykr3", "ykx3", "a3"])
    slopes = scatter_1D(key_, i, all_clean, nameC = cluster_indices, scale = True)  # , all, nameC = "id" # part10m

    plt.savefig('scatter_' + k + '.png')
    plt.show()

    I = np.argsort(-np.abs(slopes))
    print("\n")  # impact order
    for i in I:
        print(key_[i], slopes[i])

#
# Figure 3
#
