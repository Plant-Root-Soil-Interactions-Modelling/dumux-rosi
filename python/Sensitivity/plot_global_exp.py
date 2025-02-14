"""
(experimental) plots data from global optimization 
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import json
import numpy as np
from scipy.cluster.vq import vq, whiten, kmeans
import matplotlib.pyplot as plt

from minisom import MiniSom

import global_optimization_tools as got

""" font sizes """
SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = BIGGER_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title


def make_som_maps(som, data):
    """ makes two maps: node2sample (dict), and sample2node """
    node2sample = {}
    sample2node = []
    for i, xx in enumerate(data):
        node = som.winner(xx)  # (j,i)
        sample2node.append(node)
        if node in node2sample:
            node2sample[node].append(i)
        else:
            node2sample[node] = [i]

    return node2sample, sample2node


def make_kmeans_maps(obs, code_block, n):
    """ makes two maps: node2sample (dict), and sample2node """
    sample2node, dist = vq(obs, code_block)
    node2sample_ = {}
    node2sample = {}
    for i, node in enumerate(sample2node):
        if node in node2sample_:
            node2sample_[node].append(i)
            node2sample[(node % n, node // n)].append(i)  # (j,i)
        else:
            node2sample_[node] = [i]
            node2sample[(node % n, node // n)] = [i]  # (j,i)

    return node2sample, sample2node, node2sample_


def distortion_plot(data):
    data_ = whiten(data)
    dist_ = []
    for i in range(0, 20):
        centers, dist = kmeans(data_, i + 1)
        dist_.append(dist)
    plt.plot(range(0, 20), dist_)
    plt.ylabel("distortion")
    plt.xlabel("number of clusters")
    plt.show()


def try_cluster(data, k):
    data_ = whiten(data)
    res = kmeans(data_, k)
    print(res)


""" """
folder_path = "results_cplantbox/"
exp_name = "soybean_all14"

""" load everything & merge npz results into input parameter json"""

all = got.load_json_files(exp_name, folder_path)
got.merge_results(folder_path, all)

target_names = ["length", "volume", "depth", "RLDz", "krs", "SUFz"]  # ["length", "surface",
data = got.fetch_features(target_names, all)

print("\ndata")
got.print_info(data, target_names)  # RAW
print(data.shape)

data = got.filter_data(data, 0, 200., 14000)  # 76 * 3  *100 * 0.6 = 13680 cm;
print("\nfiltered")
got.print_info(data, target_names)  # filtered
print(data.shape)
data2 = data.copy()
data = got.scale_data(data)

print()
print("target_names", target_names)

distortion_plot(data2)  # to determine m and n
m_neurons = 2  # 10
n_neurons = 4  # 10

# # som = MiniSom(n_neurons, m_neurons, data.shape[1], sigma = 1.5, learning_rate = .5,
# #               neighborhood_function = 'gaussian', random_seed = 0, topology = 'rectangular')
# som = MiniSom(n_neurons, m_neurons, data.shape[1], sigma = 1.5, learning_rate = .7, activation_distance = 'euclidean',
#               topology = 'hexagonal', neighborhood_function = 'gaussian', random_seed = 10)
# som.pca_weights_init(data)
# som.train(data, 1000, verbose = False)  # random training
# node2sample, sample2node = make_som_maps(som, data)

centers, dist = kmeans(data, m_neurons * n_neurons, rng = 1)
node2sample, sample2node, _ = make_kmeans_maps(data, centers, n_neurons)

ind_ = list(range(0, len(target_names)))
got.plot_targets(ind_, target_names, data2, m_neurons, n_neurons, node2sample, sample2node, "rose")  # hexagon, rose
# plot_targets(ind_, target_names, data2, m_neurons, n_neurons, node2sample, sample2node, "hexagon")  # hexagon, rose

pbounds = {
    'src_a': (3, 11),
    'src_first_a': (3, 14),
    'src_delay_a': (3, 14),
    # 'a145_a': (0.025, 0.25),
    'lmax145_a': (50, 150),
    'ln145_a': (0.5, 10.),
    'r145_a': (0.2, 7.),
    # 'tropismN145_a': (0., 7.),
    # 'tropismS145_a': (0., 1.),
    # 'hairsLength145_a': (0.1, 2.),
    # 'hairsZone145_a': (1., 10.),
    # 'hairsElongation145_a': (0.1, 2.),
    # 'a2_a': (0.01, 0.1),
    'lmax2_a': (5., 50.),
    'ln2_a': (0.5, 10.),
    'r2_a': (0.2, 7.),
    # 'tropismN2_a': (0., 7.),
    # 'tropismS2_a': (0., 1.),
    # 'hairsLength2_a': (0.1, 2.),
    # 'hairsZone2_a': (1., 10.),
    # 'hairsElongation2_a': (0.1, 2.),
    # 'a3_a': (0.01, 0.1),
    'lmax3_a': (5., 50.),
    # 'ln3_a': (0.5, 10.),
    'r3_a': (0.2, 7.),
    # 'tropismN3_a': (0., 7.),
    # 'tropismS3_a': (0., 1.),
    # 'hairsLength3_a': (0.1, 2.),
    # 'hairsZone3_a': (1., 10.),
    # 'hairsElongation3_a': (0.1, 2.),
    }
param_names = pbounds.keys()
data_params = got.fetch_features(param_names, all)

print("Parameter space", data_params.shape)
got.print_info(data_params, param_names)

node01 = (1, 0)
ind01 = np.array(node2sample[node01], dtype = np.int64)
print("\nParameter space node", node01, data_params[ind01].shape)
got.print_info(data_params[ind01], param_names)

distortion_plot(data_params[ind01])
k = 20

centers, dist = kmeans(data_params[ind01], k, rng = 1)
node2sample, sample2node, node2sample_ = make_kmeans_maps(data_params[ind01], centers, n_neurons)

for i in range(0, k):
    ind = np.array(node2sample_[i], dtype = np.int64)
    print("\nParameter space node", ind, data_params[ind].shape)
    got.print_info(data_params[ind], param_names)

#
# plot_targets(ind_, target_names, data2, m_neurons, n_neurons, node2sample, sample2node, "rose")
# plot_targets(ind_, target_names, data2, m_neurons, n_neurons, node2sample, sample2node, "hexagon")

# # DISTANCE MAP
# plt.figure(figsize = (9, 9))
# plt.title("Distance map")
# # plt.pcolor(som.distance_map().T)  # plotting the distance map as background , cmap = 'bone_r'
# plot_hexagons(som.distance_map().T, plt.gca(), cm.jet)
# # plt.colorbar()
# plt.show()
#
# # FREQUENCY MAP
# plt.figure(figsize = (9, 9))
# plt.title("Frequency map")
# frequencies = som.activation_response(data)
# # plt.pcolor(frequencies.T, cmap = 'Blues')
# plot_hexagons(frequencies.T, plt.gca(), cm.jet)
# # plt.colorbar()
# plt.show()

# TARGET PERFORMANCE
# target = 0  # "krs"
# img = np.zeros((m_neurons, n_neurons))
# for i in range(0, m_neurons):
#     for j in range(0, n_neurons):
#         node = (j, i)
#         if node in node2sample:
#             ind_ = node2sample[node]
#             # img[i, j] = len(ind_)  # FREQUENCY (by hand)
#             img[i, j] = np.sum(data[ind_, target]) / len(ind_)
#
# plt.pcolor(img)
# plt.colorbar()
# plt.show()

# max_iter = 200
# q_error = []
# t_error = []
# for i in range(max_iter):
#     rand_i = np.random.randint(len(data))
#     som.update(data[rand_i], som.winner(data[rand_i]), i, max_iter)
#     q_error.append(som.quantization_error(data))
#     t_error.append(som.topographic_error(data))
# plt.subplot(2, 1, 1)
# plt.plot(np.arange(max_iter), q_error)
# plt.ylabel('quantization error')
# plt.subplot(2, 1, 2)
# plt.plot(np.arange(max_iter), t_error)
# plt.ylabel('topographic error')
# plt.xlabel('iteration index')
# plt.tight_layout()
# plt.show()

# print("data set", len(all), "entries")
# print(all[10])

# scatter_1D("lmax1_a", "length", all)
# scatter_1D(["r1_a"], ["length"], all)
# scatter_1D(["r1_a", "ln1_a"], ["length", "volume"], all)

# scatter_1D("hairsLength_a", "krs", all,  "lmax1_a")
# scatter_1D("hairsZone_a", "krs", all,  "lmax1_a")

# """ Plots fo the Inari talk 6.2.2025 """
#
# folder_path = "results_cplantbox/"
# exp_name = "soybean_all14"
#
# all = got.load_json_files(exp_name, folder_path)
# got.merge_results(folder_path, all)
#
# target_names = ["length", "volume", "depth", "RLDz", "krs", "SUFz"]  # ["length", "surface",
# data = got.fetch_features(target_names, all)
# print("\ndata")
# got.print_info(data, target_names)  # RAW
# print(data.shape)
#
# data = got.filter_data(data, 0, 200., 14000)  # 76 * 3  *100 * 0.6 = 13680 cm;
# print("\nfiltered")
# got.print_info(data, target_names)  # filtered
# print(data.shape, "\n")
# data2 = data.copy()
# data = got.scale_data(data)
# got.print_info(data, target_names)
# print(data.shape, "\n")
#
# print("target_names", target_names)
#
# m_neurons = 4  # 10
# n_neurons = 4  # 10
#
# # som = MiniSom(n_neurons, m_neurons, data.shape[1], sigma = 1.5, learning_rate = .5,
# #               neighborhood_function = 'gaussian', random_seed = 0, topology = 'rectangular')
# som = MiniSom(n_neurons, m_neurons, data.shape[1], sigma = 1.5, learning_rate = .7, activation_distance = 'euclidean',
#               topology = 'hexagonal', neighborhood_function = 'gaussian', random_seed = 10)
#
# som.pca_weights_init(data)
# som.train(data, 1000, verbose = False)  # random training
#
# node2sample, sample2node = make_maps(som)
#
# ind_ = list(range(0, len(target_names)))
#
# plot_targets(ind_, target_names, data2, m_neurons, n_neurons, node2sample, sample2node, "rose")
# plot_targets(ind_, target_names, data2, m_neurons, n_neurons, node2sample, sample2node, "hexagon")
#
# # DISTANCE MAP
# plt.figure(figsize = (9, 9))
# plt.title("Distance map")
# plot_hexagons(som.distance_map().T, plt.gca(), cm.jet)
# plt.show()
#
# # FREQUENCY MAP
# plt.figure(figsize = (9, 9))
# plt.title("Frequency map")
# frequencies = som.activation_response(data)
# plot_hexagons(frequencies.T, plt.gca(), cm.jet)
# plt.show()

