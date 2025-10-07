"""
    Macrcoscopic:

    zips the results from result_cplantbox/ folder?
    
    Daniel Leitner, 2025      
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import json
import zipfile
import numpy as np
from scipy.cluster.vq import vq, whiten, kmeans
import matplotlib.pyplot as plt
from matplotlib import cm, colorbar

from minisom import MiniSom

import global_optimization_tools as got
import run_cplantbox_exp

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


def bitwise_list_(all, ind):
    """ bitwise indexing of list all, i.e. returns all elements where ind is True """
    r = []
    for i, a in enumerate(all):
        if ind[i] == True:
            r.append(a)
    return r


def distortion_plot(all, target_names):
    """ kmeans distortion versus number of clusters (for kmeans) """
    data_ = whiten(got.fetch_features(target_names, all))  # normalize data
    plt.title("N = {:g}".format(data_.shape[0]))
    dist_ = []
    for i in range(0, 20):
        centers, dist = kmeans(data_, i + 1)
        dist_.append(dist)
    plt.plot(range(1, 21), dist_)
    plt.ylabel("distortion")
    plt.xlabel("number of clusters")
    plt.show()


def mean_of_dict(all, keys):
    """ puts the mean values of certain keys in a new dict (all is a list of dict) """
    data = got.fetch_features(keys, all)
    mean_data = np.mean(data, axis = 0)
    # print(keys)
    # print(mean_data)
    # print(mean_data.shape)
    new_dict = all[0].copy()
    for j, k in enumerate(keys):
        new_dict[k] = mean_data[j]
    return new_dict


def exemplify(all, node2sample, m_neurons, n_neurons, vis_hairs = True, N = 1):
    """ saves N png images of pareto solutions for each node"""
    for i in range(0, m_neurons):
        for j in range(0, n_neurons):
            node = (j, i)
            c = 0
            for k in range(0, N):
                ind = node2sample[node][c]
                pareto = all[ind]["pareto"]  # find next pareto solution in cluster node (j,i)
                while pareto == False and c + 1 < len(node2sample[node]):
                    c += 1
                    ind = node2sample[node][c]
                    pareto = all[ind]["pareto"]
                if pareto:
                    print(i, j, "example ind", ind, "which is choice", c)
                    run_cplantbox_exp.show_me_file(all[ind]["exp_name"], "results_cplantbox/", vis_hairs = True, mode = "png", png_name = "test({:g}, {:g})_{:g}".format(j, i, k))
                    c += 1


def dist_dict_(d1, d2, keys):
    if not d2["pareto"]:
        return 1.e6
    sum = 0.
    for k in keys:
        sum += np.square(d1[k] - d2[k])
    return np.sqrt(sum)


def exemplify_closest(all, node2sample, m_neurons, n_neurons, keys, vis_hairs = True):
    """ saves a single png image the closest pareto solution to the mean value """
    for i in range(0, m_neurons):
        for j in range(0, n_neurons):
            node = (j, i)
            ind = node2sample[node]
            data_ = got.fetch_features(keys, [all[k] for k in ind])
            data = whiten(data_)
            node_mean = np.mean(data, axis = 0)
            dists = np.sqrt(np.sum(np.square(data - node_mean), axis = 1))
            for k_, k in enumerate(ind):
                if not all[k]["pareto"]:
                    dists[k_] = 1.e6
            print("\nnode_mean", node_mean, node_mean.shape)
            print("dists", dists.shape, dists)
            min_i = np.argmin(dists)
            print(i, j, "example ind", ind[min_i], "distance", dists[min_i], "\n")
            run_cplantbox_exp.show_me_file(all[ind[min_i]]["exp_name"], "results_cplantbox/", vis_hairs = True, mode = "png", png_name = "test({:g}, {:g})_mid".format(j, i))


def write_closest(all, node2sample, m_neurons, n_neurons, keys, vis_hairs = True):
    lines = []
    for i in range(0, m_neurons):
        for j in range(0, n_neurons):
            node = (j, i)
            ind = node2sample[node]
            data_ = got.fetch_features(keys, [all[k] for k in ind])
            data = whiten(data_)
            node_mean = np.mean(data, axis = 0)
            dists = np.sqrt(np.sum(np.square(data - node_mean), axis = 1))
            for k_, k in enumerate(ind):
                if not all[k]["pareto"]:
                    dists[k_] = 1.e6
            min_i = np.argmin(dists)
            lines.append(all[ind[min_i]]["exp_name"])
    with open("my_pick.txt", 'w', encoding = 'utf-8') as file:
        for line in lines:
            file.write(line + '\n')


def exemplify_cluster_mean(all, node2sample, m_neurons, n_neurons, keys):
       for i in range(0, m_neurons):
        for j in range(0, n_neurons):
            node = (j, i)
            ind = node2sample[node]
            winner = mean_of_dict([all[i] for i in ind], keys)
            print(winner)
            print(type(winner))
            run_cplantbox_exp.show_me(winner, "results_cplantbox/", vis_hairs = True, mode = "png", png_name = "test({:g}, {:g})".format(j, i))


""" 1 load everything & merge npz results into input parameter json"""
folder_path = "results_cplantbox/"
# exp_name = "soybean_length14"
exp_name = "soybean_all14"

# all = got.load_json_files(exp_name, folder_path)  # open parameter files
# got.merge_results(folder_path, all)  # add results
# all_json = {}
# for a in all:
#     all_json[a["exp_name"]] = a
#
# # Write all JSON data to a zip file
# with zipfile.ZipFile(exp_name + ".zip", "w", zipfile.ZIP_DEFLATED) as zipf:
#     with zipf.open(exp_name + ".json", "w") as json_file:
#         json_file.write(json.dumps(all_json, indent = 4).encode("utf-8"))

