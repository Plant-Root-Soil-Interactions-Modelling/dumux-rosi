"""
(experimental) plots data from global optimization 
"""
import sys

from scipy import stats; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import json
import numpy as np
from scipy.cluster.vq import vq, whiten, kmeans
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, RegularPolygon, Ellipse
from matplotlib import cm, colorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable

from minisom import MiniSom

"""
file io                load_json_files, merge_results                
clustering             make_som_maps, make_kmeans_maps, scale_data, label_clusters
data management        fetch_features, print_info, filter_data, filter_list, pareto_data, pareto_list
plotting               plot_1D, scatter_1D, scatter_1D_cross, plot_targets, prepare_img

Naming:

the objectives of the multiple objectives optimization are also called targets 
parameters and objectives are called features

either data is represented by a list of dictionaries (all), and features are represented as list of keys (target_names)
or data is represented as numpy array (data), and features are represented as list of indices (target_ind)
"""


def load_json_files(exp_name, folder_path):
    """ opens all json files starting with exp_name in a folder and returns a list of Python dictionaries 
        (one dictionary by file) containing the input parameter sets """
    json_data_list = []
    for filename in os.listdir(folder_path):
        if filename.startswith(exp_name) and filename.endswith(".json"):
            file_path = os.path.join(folder_path, filename)
            try:
                with open(file_path, 'r', encoding = 'utf-8') as file:
                    data = json.load(file)
                    json_data_list.append(data)
            except (json.JSONDecodeError, IOError) as e:
                print(f"Error loading {filename}: {e}")

    print("global_optimization_tools.load_json_files(): opened", len(json_data_list), "files")
    return json_data_list


def merge_results(folder_path, json_data_list, i = -1):
    """ Merges the json_data (list of dict of input parameters) 
        with simulation results given within .npz with corresponding name """
    for data in json_data_list:

        filename = data["exp_name"]
        filename += ".npz"
        file_path = os.path.join(folder_path, filename)
        try:
            results = np.load(file_path)
            data["length"] = results["length"][i]
            data["surface"] = results["surface"][i]
            data["volume"] = results["volume"][i]
            data["-volume"] = -results["volume"][i]  # MINUS because we want to minimize
            data["depth"] = np.abs(results["depth"][i])  # they are stored with a (wrong) sign
            data["RLDmean"] = results["RLDmean"][i]
            data["RLDz"] = np.abs(results["RLDz"][i])  # they are stored with a (wrong) sign
            data["krs"] = results["krs"][i]
            data["SUFz"] = np.abs(results["SUFz"][i])  # they are stored with a (wrong) sign

        except (json.JSONDecodeError, IOError) as e:
            print(f"Error loading {filename}: {e}")


def make_som_maps(som, data):
    """ makes two maps from som clustering: node2sample (dict), and sample2node """
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
    """ makes two maps from kmeans clustering: node2sample (dict), and sample2node """
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


def label_clusters(all, n_neurons, m_neurons, target_names, mode):
    """ creates a som or kmeans clustering regarding the target_ind and adds "node" to each point in all"""
    som = None
    data = fetch_features(target_names, all)
    if mode == "som":
        data = scale_data(data)  # mean = 0, std = 1
        # som = MiniSom(n_neurons, m_neurons, data.shape[1], sigma = 1.5, learning_rate = .5,
        #               neighborhood_function = 'gaussian', random_seed = 0, topology = 'rectangular')
        som = MiniSom(n_neurons, m_neurons, data.shape[1], sigma = 1.5, learning_rate = .7, activation_distance = 'euclidean',
                      topology = 'hexagonal', neighborhood_function = 'gaussian', random_seed = 10)
        som.pca_weights_init(data)
        som.train(data, 1000, verbose = False)  # random training
        node2sample, sample2node = make_som_maps(som, data)
    elif mode == "kmeans":
        data = whiten(data)  # mean 0
        centers, dist = kmeans(data, m_neurons * n_neurons, rng = 1)
        node2sample, sample2node, _ = make_kmeans_maps(data, centers, n_neurons)
    else:
        raise "label_clusters(): unknown clustering method" + mode

    for k, point in enumerate(all):  # add results to all
        j, i = sample2node[k]
        # print(j, i, n_neurons, m_neurons)
        point["node"] = i * n_neurons + j

    return node2sample, sample2node, som


def get_(name, all):
    """ returns the key @param name from a @param all which is a list of dicts as np.array """
    l = np.array([a[name] for a in all])
    return l


def fetch_features(feature_names, all):
    """ returns feature_names as np.array with the shape (len(data), len(features)) """
    data = []
    for name in feature_names:
        data.append(get_(name, all))

    data = np.array(data).transpose()
    return data


def rstd_(data, feature_names):
    """ mean of the relative standard deviations, i.e. standard deviations in percent of mean value """
    rstd_ = [rstd_.append(100.* np.std(data[:, i]) / np.mean(data[:, i])) for i in range(0, data.shape[1])]
    return np.mean(rstd_)


def print_info(data, feature_names):
    """ prints features min, max, mean, and sd (parameter names labels the data columns) """
    print("name".ljust(20) + "\tmin\t\tmax\t\tmean\t\tstd\t\trstd")
    rstd_ = []
    for i, name in enumerate(feature_names):
        # print(name, "\t", np.min(data[:, i]), "\t", np.max(data[:, i]), "\t", np.mean(data[:, i]), "\t", np.std(data[:, i]))
        print("{:s}\t{:10.4f}\t{:10.4f}\t{:10.4f}\t{:10.4f}\t{:10.4f}".format(
            name.ljust(20), np.min(data[:, i]), np.max(data[:, i]), np.mean(data[:, i]),
            np.std(data[:, i]), 100.* np.std(data[:, i]) / np.mean(data[:, i])))
        rstd_.append(100.* np.std(data[:, i]) / np.mean(data[:, i]))
    print("\tmean rstd", np.mean(rstd_))


def filter_data(data, feature_ind, min_, max_):
    """ crops data (np.array) to values where feature_ind stays within (min_, max_) """
    ind_ = np.bitwise_and(data[:, feature_ind] > min_, data[:, feature_ind] < max_)
    data = data[ind_,:]
    return data


def filter_list(data, feature_str, min_, max_):
    """ crops data (list of dict) to values where feature_str stays within (min_, max_) """
    ind_ = []
    for d in data:
        if d[feature_str] > min_ and d[feature_str] < max_:
            ind_.append(d)
    return ind_


def is_pareto_efficient_(costs, return_mask = True):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :param return_mask: True to return a mask
    :return: An array of indices of pareto-efficient points.
        If return_mask is True, this will be an (n_points, ) boolean array
        Otherwise it will be a (n_efficient_points, ) integer array of indices.
    
    # from https://stackoverflow.com/questions/32791911/fast-calculation-of-pareto-front-in-python
    """
    is_efficient = np.arange(costs.shape[0])
    n_points = costs.shape[0]
    next_point_index = 0  # Next index in the is_efficient array to search for
    while next_point_index < len(costs):
        nondominated_point_mask = np.any(costs < costs[next_point_index], axis = 1)
        nondominated_point_mask[next_point_index] = True
        is_efficient = is_efficient[nondominated_point_mask]  # Remove dominated points
        costs = costs[nondominated_point_mask]
        next_point_index = np.sum(nondominated_point_mask[:next_point_index]) + 1
    if return_mask:
        is_efficient_mask = np.zeros(n_points, dtype = bool)
        is_efficient_mask[is_efficient] = True
        return is_efficient_mask
    else:
        return is_efficient


def pareto_data(data, feature_ind):
    """ crops data to its pareto set, using is_pareto_efficient """
    cost = data[:, feature_ind]
    mask = is_pareto_efficient_(-cost)
    return mask


def pareto_list(all, feature_str):
    """ crops data to its pareto set, using is_pareto_efficient """
    cost = fetch_features(feature_str, all)
    mask = is_pareto_efficient_(-cost)
    return mask


def scale_data(data):
    """ scales data to have a mean of 0 and a std of 1 """
    mean = np.mean(data, axis = 0)
    std = np.std(data, axis = 0)
    data2 = data - mean
    data2 = np.divide(data2, std)
    # print("data2 mean", np.mean(data2, axis = 0))
    # print("data2 std", np.std(data2, axis = 0))
    return data2


def plot_1D(nameX, nameY, all):
    """ plots nameX vs nameY sorting the values along the x-axis """
    x = get_(nameX, all)
    y = get_(nameY, all)
    I = np.argsort(x)
    plt.plot(x[I], y[I])
    plt.xlabel(nameX)
    plt.ylabel(nameY)
    plt.show()


def scatter_1D(nameX, nameY, all, nameC = None):
    """ scatter plots one or multiple pairs of nameX vs nameY """
    if isinstance(nameX, str):
        nameX = [nameX]
    if isinstance(nameY, str):
        nameY = [nameY]
    assert len(nameX) == len(nameY), "oh no"

    n = np.ceil(np.sqrt(len(nameX)))
    print(n, np.sqrt(len(nameX)))
    fig, ax = plt.subplots(int(len(nameX) // n) + 1, int(n), figsize = (16, 16))
    if n == 1:
        ax = np.array([ax])
    ax_ = ax.flat

    for i, _ in enumerate(nameX):
        x = get_(nameX[i], all)
        y = get_(nameY[i], all)
        if nameC is not None:
            c = get_(nameC, all)
        I = np.argsort(x)
        if nameC is not None:
            ax_[i].scatter(x[I], y[I], c = c[I], alpha = 0.5)
        else:
            ax_[i].scatter(x[I], y[I])

        res = stats.linregress(x[I], y[I])
        ax_[i].plot(x[I], res.intercept + res.slope * x[I], 'k:', label = 'fitted')

        ax_[i].set_xlabel(nameX[i])
        ax_[i].set_ylabel(nameY[i])
    plt.tight_layout()
    plt.show()


def scatter_1D_cross(names, all):
    """Generates pairs of names for 1D scatter plots and calls the scatter_1D function."""
    nameX = []
    nameY = []
    for i, nx in enumerate(names):
        for j, ny in enumerate(names):
            if i < j:
                nameX.append(nx)
                nameY.append(ny)
    scatter_1D(nameX, nameY, all)


def plot_hexagons(img, ax, f = cm.viridis):
    """ adds the img data as hexagons """
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
    """ single rose plot """
    margin = 2. / 180. * np.pi
    theta_ = np.linspace(0, 2 * np.pi, len(data) + 1)

    for k, d in enumerate(data):
        x0 = j + 0.5 * (i % 2)
        y0 = i * np.sqrt(3) / 2
        x1 = 0.5 * np.abs(d) * np.cos(theta_[k] + margin) + x0
        y1 = 0.5 * np.abs(d) * np.sin(theta_[k] + margin) + y0
        x2 = 0.5 * np.abs(d) * np.cos(theta_[k + 1] - margin) + x0
        y2 = 0.5 * np.abs(d) * np.sin(theta_[k + 1] - margin) + y0
        tri = Polygon([[x0, y0], [x1, y1], [x2, y2], [x0, y0]],
                        facecolor = f(k / len(data)),
                        alpha = 1.,
                        edgecolor = 'gray')

        ax.add_patch(tri)
        ax.text(x0, y0 + 0.5, "({:g},{:g})={:g}".format(j, i, i * n_neurons + j), size = 24)  # why j, i? (som did it)

    # ax.set_xlim((-1, n_neurons + 0.5))
    # ax.set_ylim((-1, m_neurons - 0.5))


def plot_targets(all, target_names, m_neurons, n_neurons, node2sample, sample2node, mode = "quad"):
    """ plots multiple objectives using prepare_img
    mode is "quad", "hexagon", "rose" 
    """
    data = fetch_features(target_names, all)

    if mode == "quad" or mode == "hexagon":
        n = np.ceil(np.sqrt(len(target_names)))
        fig, ax = plt.subplots(int(len(target_names) // n + 1 * (len(target_names) % n > 0)), int(n), figsize = (16, 16))
    else:
        n = 1
        fig, ax = plt.subplots(n, n, figsize = (16, 16))
    if n == 1:
        ax = np.array([ax])
    ax_ = ax.flat

    img_ = []
    for k, _ in enumerate(target_names):

         img = prepare_img(k, data, m_neurons, n_neurons, node2sample, sample2node)
         img = scale01_(img)

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
        for k, _ in enumerate(target_names):
            img_[k,:,:] = img_[k,:,:] / np.max(img_[k,:,:].flat)
        for i in range(0, m_neurons):
            for j in range(0, n_neurons):
                data = img_[:, i, j]
                plot_rose(data, i, j, m_neurons, n_neurons, ax_[0])
        for k, name in enumerate(target_names):
            plt.plot(0, 0, color = cm.jet(k / len(target_names)), label = name, linewidth = 5)
        plt.legend()

    plt.show()


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


def scale01_(img):
    """ scales the img from 0 to 1 (e.g. for color mapping) """
    img = np.abs(img)
    return np.divide(img - np.ones(img.shape) * np.min(img.flat), np.ones(img.shape) * (np.max(img.flat) - np.min(img.flat)))


if __name__ == "__main__":

    folder_path = "results_cplantbox/"
    exp_name = "soybean_test"

    # all = load_json_files(exp_name, folder_path)
    # merge_results(folder_path, all)
    #
    # # scatter_1D("lmax1_a", "length", all)
    # # scatter_1D(["r1_a"], ["length"], all)
    # # scatter_1D(["r1_a", "ln1_a"], ["length", "volume"], all)
    #
    # # scatter_1D("hairsLength_a", "krs", all,  "lmax1_a")
    # # scatter_1D("hairsZone_a", "krs", all,  "lmax1_a")
    #
    # print("data set", len(all), "entries")
    # print(all[0])
    #
    # # StandardScaler (sklearn)
    #
    # # step 1 (scale features, not targets)
    # # step 2 (pca, pca.explained_variance ratio
    #
    # feature_names = ['a1_a', 'lmax1_a', 'ln1_a', 'r1_a', 'tropismN1_a', 'hairsLength_a', 'hairsZone_a', 'hairsElongation_a']
    # data = fetch_features(feature_names, all)
    # data = scale_data(data)

