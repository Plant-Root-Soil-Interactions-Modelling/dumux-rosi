"""
(experimental) plots data from global optimization 
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from minisom import MiniSom

import global_optimization_tools as got


def make_maps(som):
    """ makes to maps node2sample (dict), and sample2node """

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


def prepare_img(target_ind, data, m_neurons, n_neurons, node2sample, sample2node):
    """ prepares an image representing the mean target value on the nodes """
    img = np.zeros((m_neurons, n_neurons))
    for i in range(0, m_neurons):
        for j in range(0, n_neurons):
            node = (j, i)
            if node in node2sample:
                ind_ = node2sample[node]
                img[i, j] = np.mean(data[ind_, target_ind])

    return img


def plot_targets(target_indices, target_names, data, m_neurons, n_neurons, node2sample, sample2node):

    n = np.ceil(np.sqrt(len(target_indices)))
    fig, ax = plt.subplots(int(n), int(n), figsize = (16, 16))
    if n == 1:
        ax = np.array([ax])
    ax_ = ax.flat

    for i, target in enumerate(target_indices):
         img = prepare_img(target, data, m_neurons, n_neurons, node2sample, sample2node)
         ax_[i].pcolor(img)  # pcolor
         ax_[i].set_title(target_names[i])
         ax_[i].set_xticks([])
         ax_[i].set_yticks([])
         # divider = make_axes_locatable(ax_[i])
         # cax = divider.append_axes('right', size = '5%', pad = 0.05)
         # fig.colorbar(img, cax = cax, orientation = 'vertical')

    plt.show()

# """ font sizes """
# SMALL_SIZE = 12
# MEDIUM_SIZE = 16
# BIGGER_SIZE = 16
# plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
# plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
# plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
# plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
# plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
# plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
# plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title


""" """
folder_path = "results_cplantbox/"
exp_name = "soybean_all14"

""" load everything & merge npz results into input parameter json"""

all = got.load_json_files(exp_name, folder_path)
got.merge_results(folder_path, all)

target_names = ["length", "surface", "volume", "depth", "RLDmean", "RLDz", "krs", "SUFz"]
data = got.fetch_features(target_names, all)
data = got.scale_data(data)

n_neurons = 3
m_neurons = 3
som = MiniSom(n_neurons, m_neurons, data.shape[1], sigma = 1.5, learning_rate = .5,
              neighborhood_function = 'gaussian', random_seed = 0, topology = 'rectangular')

som.pca_weights_init(data)
som.train(data, 1000, verbose = False)  # random training

node2sample, sample2node = make_maps(som)

ind_ = list(range(0, len(target_names)))
plot_targets(ind_, target_names, data, m_neurons, n_neurons, node2sample, sample2node)

# # DISTANCE MAP
# plt.figure(figsize = (9, 9))
# plt.pcolor(som.distance_map().T)  # plotting the distance map as background , cmap = 'bone_r'
# plt.colorbar()
# plt.show()

# # FREQUENCY MAP
# plt.figure(figsize = (7, 7))
# frequencies = som.activation_response(data)
# plt.pcolor(frequencies.T, cmap = 'Blues')
# plt.colorbar()
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

