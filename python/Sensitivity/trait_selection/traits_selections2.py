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
from sklearn.preprocessing import StandardScaler


def scatter_1D(nameX, index_y, all, nameC = None, scale = False):
    """ scatter plots one or multiple pairs of nameX vs nameY """
    if isinstance(nameX, str):
        nameX = [nameX]
    n = np.ceil(np.sqrt(len(nameX)))
    fig, ax = plt.subplots(int(len(nameX) // n + 1 * (len(nameX) % n > 0)), int(n), figsize = (16, 16))
    if n == 1:
        ax = np.array([ax])
    ax_ = ax.flat
    slopes = []
    for i, _ in enumerate(nameX):
        x = np.array(all[:, i])
        y = np.array(all[:, index_y])
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


# Load the Excel file
data = pd.read_excel("Root_Traits_Of_Interest.xlsx")

# Remove rows with missing values
data = data.dropna()

# Separate the features and the genotype column
genotypes = data['Genotype']
features = data.drop(columns = ['Genotype'])

# Standardize the features
scaler = StandardScaler()
scaled_features = scaler.fit_transform(features)

# add cluser indices as color
cluster_n = 12
centers, dist = kmeans(scaled_features, cluster_n)
cluster_indices, _ = vq(scaled_features, centers)

key_ = list(data.keys())[1:]  # [9:]  # skip the first three

for i, k in enumerate(key_):

    target_names = [k]
    # key_.extend(["kr1", "kx1", "ykr1", "ykx1", "a1", "kr2", "kx2", "ykr2", "ykx2", "a2", "kr3", "kx3", "ykr3", "ykx3", "a3"])
    slopes = scatter_1D(range(0, 6), i, scaled_features, nameC = cluster_indices, scale = F)  # , all, nameC = "id" # part10m

    plt.savefig('scatter2_' + k + '.png')
    plt.show()

    I = np.argsort(-np.abs(slopes))
    print("\n")  # impact order
    for i in I:
        print(key_[i], slopes[i])
