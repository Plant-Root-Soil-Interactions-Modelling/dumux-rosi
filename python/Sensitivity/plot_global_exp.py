"""
(experimental) plots data from global optimization 
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon, Ellipse
from matplotlib import cm, colorbar
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


def scale01(img):
    """ scales the img from 0 to 1 (e.g. for color mapping) """
    return np.divide(img - np.ones(img.shape) * np.min(img.flat), np.ones(img.shape) * (np.max(img.flat) - np.min(img.flat)))


def plot_hexagons(img, ax, f = cm.viridis):
    """ adds the img data as hexagons """
    img = scale01(img)
    xx = list(range(0, img.shape[1]))
    yy = list(range(0, img.shape[0]))
    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            wy = yy[i] * np.sqrt(3) / 2
            hex = RegularPolygon((xx[j] + 0.5 * (i % 2) , wy),
                                 numVertices = 6,
                                 radius = .95 / np.sqrt(3),
                                 facecolor = f(img[i, j]),
                                 alpha = 1.,
                                 edgecolor = 'gray')
            ax.add_patch(hex)
    ax.set_xlim((-1, img.shape[1] + 0.5))
    ax.set_ylim((-1, img.shape[0] - 0.5))


def plot_roses(target_indices, target_names, data, m_neurons, n_neurons, node2sample, sample2node):
    """ """
    n = np.ceil(np.sqrt(len(target_indices)))
    fig, ax = plt.subplots(int(n), int(n), figsize = (16, 16))
    if n == 1:
        ax = np.array([ax])
    ax_ = ax.flat

    for i, target in enumerate(target_indices):
         img = prepare_img(target, data, m_neurons, n_neurons, node2sample, sample2node)
         # ax_[i].pcolor(img, cmap = "viridis")  # pcolor
         plot_hexagons(img, ax_[i])
         ax_[i].set_title(target_names[i])
         ax_[i].set_xticks([])
         ax_[i].set_yticks([])
         # divider = make_axes_locatable(ax_[i])
         # cax = divider.append_axes('right', size = '5%', pad = 0.05)
         # fig.colorbar(img, cax = cax, orientation = 'vertical')

    plt.show()

    pass


def plot_targets(target_indices, target_names, data, m_neurons, n_neurons, node2sample, sample2node):
    """ plots multiple objectives using prepare_img """
    n = np.ceil(np.sqrt(len(target_indices)))
    fig, ax = plt.subplots(int(n), int(n), figsize = (16, 16))
    if n == 1:
        ax = np.array([ax])
    ax_ = ax.flat

    for i, target in enumerate(target_indices):
         img = prepare_img(target, data, m_neurons, n_neurons, node2sample, sample2node)
         # ax_[i].pcolor(img, cmap = "viridis")  # pcolor
         plot_hexagons(img, ax_[i])
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

target_names = ["length", "volume", "depth", "RLDz", "krs", "SUFz"]  # ["length", "surface",
data = got.fetch_features(target_names, all)
got.print_info(data, target_names)  # RAW
data = got.filter_data(data, 0, 200., 14000)  # 76 * 3  *100 * 0.6 = 13680 cm;
print("\nfiltered")
got.print_info(data, target_names)  # filtered
data = got.scale_data(data)

n_neurons = 3  # 10
m_neurons = 3  # 10
# som = MiniSom(n_neurons, m_neurons, data.shape[1], sigma = 1.5, learning_rate = .5,
#               neighborhood_function = 'gaussian', random_seed = 0, topology = 'rectangular')
som = MiniSom(n_neurons, m_neurons, data.shape[1], sigma = 1.5, learning_rate = .7, activation_distance = 'euclidean',
              topology = 'hexagonal', neighborhood_function = 'gaussian', random_seed = 10)

som.pca_weights_init(data)
som.train(data, 1000, verbose = False)  # random training

node2sample, sample2node = make_maps(som)

ind_ = list(range(0, len(target_names)))
plot_targets(ind_, target_names, data, m_neurons, n_neurons, node2sample, sample2node)

# DISTANCE MAP
plt.figure(figsize = (9, 9))
plt.title("Distance map")
# plt.pcolor(som.distance_map().T)  # plotting the distance map as background , cmap = 'bone_r'
plot_hexagons(som.distance_map().T, plt.gca(), cm.jet)
# plt.colorbar()
plt.show()

# FREQUENCY MAP
plt.figure(figsize = (9, 9))
plt.title("Frequency map")
frequencies = som.activation_response(data)
# plt.pcolor(frequencies.T, cmap = 'Blues')
plot_hexagons(frequencies.T, plt.gca(), cm.jet)
# plt.colorbar()
plt.show()

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

