"""
(experimental) plots data from global optimization 
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, RegularPolygon, Ellipse
from matplotlib import cm, colorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable

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


def make_maps(som, data):
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


def plot_rose(data, i, j, m_neurons, n_neurons, ax, f = cm.jet):
    """ """
    margin = 2. / 180. * np.pi
    theta_ = np.linspace(0, 2 * np.pi, len(data) + 1)

    for k, d in enumerate(data):
        x0 = j + 0.5 * (i % 2)
        y0 = i * np.sqrt(3) / 2
        x1 = 0.5 * d * np.cos(theta_[k] + margin) + x0
        y1 = 0.5 * d * np.sin(theta_[k] + margin) + y0
        x2 = 0.5 * d * np.cos(theta_[k + 1] - margin) + x0
        y2 = 0.5 * d * np.sin(theta_[k + 1] - margin) + y0
        tri = Polygon([[x0, y0], [x1, y1], [x2, y2], [x0, y0]],
                        facecolor = f(k / len(data)),
                        alpha = 1.,
                        edgecolor = 'gray')

        ax.add_patch(tri)
        # ax.text(x0, y0, "{:g}, {:g}".format(i, j), color = 'k', fontsize = 24)

    # ax.set_xlim((-1, n_neurons + 0.5))
    # ax.set_ylim((-1, m_neurons - 0.5))


def plot_targets(target_indices, target_names, data, m_neurons, n_neurons, node2sample, sample2node, mode = "quad"):
    """ plots multiple objectives using prepare_img """
    if mode == "quad" or mode == "hexagon":
        n = np.ceil(np.sqrt(len(target_indices)))
    else:
        n = 1
    fig, ax = plt.subplots(int(n), int(n), figsize = (16, 16))
    if n == 1:
        ax = np.array([ax])
    ax_ = ax.flat
    img_ = []
    for k, target in enumerate(target_indices):

         img = prepare_img(target, data, m_neurons, n_neurons, node2sample, sample2node)

         if mode == "quad":
             ax_[k].pcolor(img, cmap = "viridis")  # pcolor
         elif mode == "hexagon":
             plot_hexagons(img, ax_[k])
         elif mode == "rose":
             img_.append(img)
             plot_hexagons(img, ax_[0], f = lambda c: [1., 1., 1.])

         if mode == "quad" or mode == "hexagon":
             ax_[k].set_title(target_names[k])
             ax_[k].set_xticks([])
             ax_[k].set_yticks([])
         elif mode == "rose":
             ax_[0].set_xticks([])
             ax_[0].set_yticks([])

    if mode == "rose":
        img_ = np.array(img_)
        for k, target in enumerate(target_indices):
            img_[k,:,:] = img_[k,:,:] / np.max(img_[k,:,:].flat)
        for i in range(0, m_neurons):
            for j in range(0, n_neurons):
                data = img_[:, i, j]
                plot_rose(data, i, j, m_neurons, n_neurons, ax_[0])
        for k, name in enumerate(target_names):
            plt.plot(0, 0, color = cm.jet(k / len(target_names)), label = name)
        plt.legend()

    plt.show()


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

m_neurons = 4  # 10
n_neurons = 4  # 10

# som = MiniSom(n_neurons, m_neurons, data.shape[1], sigma = 1.5, learning_rate = .5,
#               neighborhood_function = 'gaussian', random_seed = 0, topology = 'rectangular')
som = MiniSom(n_neurons, m_neurons, data.shape[1], sigma = 1.5, learning_rate = .7, activation_distance = 'euclidean',
              topology = 'hexagonal', neighborhood_function = 'gaussian', random_seed = 10)
som.pca_weights_init(data)
som.train(data, 1000, verbose = False)  # random training

node2sample, sample2node = make_maps(som, data)

ind_ = list(range(0, len(target_names)))
plot_targets(ind_, target_names, data2, m_neurons, n_neurons, node2sample, sample2node, "rose")

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

