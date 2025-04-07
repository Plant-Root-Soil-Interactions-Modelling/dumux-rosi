"""
    plots results from global root structure optimization (see global_optimization_cpantbox.py)

    Clustering with SOM or Kmeans (see global_optimization_tools) 
    to experiment on data analysis, plots, etc. 
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
exp_name = "soybean_length14"
exp_name = "soybean_all14"

# all = got.load_json_files(exp_name, folder_path)  # open parameter files
# got.merge_results(folder_path, all)  # add results
#
# all_json = {}
# for a in all:
#     all_json[a["exp_name"]] = a
#
# # Write all JSON data to a zip file
# with zipfile.ZipFile(exp_name + ".zip", "w", zipfile.ZIP_DEFLATED) as zipf:
#     with zipf.open(exp_name + ".json", "w") as json_file:
#         json_file.write(json.dumps(all_json, indent = 4).encode("utf-8"))

# Read JSON data from the ZIP file
with zipfile.ZipFile(exp_name + ".zip", "r") as zipf:
    with zipf.open(exp_name + ".json", "r") as json_file:
        all = json.load(json_file)  # Deserialize JSON data
all = list(all.values())

""" 2 filter """
all = got.filter_list(all, "length", 200., 30000)  # 76 * 3  *100 * 0.6 = 13680 cm;
print("\nfiltered data")
target_names = ["surface", "volume", "depth", "RLDz", "krs", "SUFz"]
d = got.fetch_features(target_names, all)
got.print_info(d, target_names)  # filtered
print(d.shape)
# distortion_plot(all, target_names)  # to determine m and n

""" 3 cluster the targets """
cluster_targets = ["surface", "depth", "-volume", "krs", "SUFz", "RLDz"]  # nutrients, low carbon, reach water table,
# distortion_plot(all, target_names)  # to determine m and n
m_neurons = 3
n_neurons = 4
node2sample, sample2node, som = got.label_clusters(all, n_neurons, m_neurons, cluster_targets, "som")

print()
for i in range(0, m_neurons * n_neurons):
    node = (i % n_neurons, i // n_neurons)
    print("Node {:g} = ({:g},{:g})".format(i, i % n_neurons, i // n_neurons), len(node2sample[node]), "samples")

# # DISTANCE MAP
# plt.figure(figsize = (9, 9))
# plt.title("Distance map")
# got.plot_hexagons(som.distance_map().T, plt.gca(), cm.viridis)
# plt.show()
# # FREQUENCY MAP (not working?!)
# plt.figure(figsize = (9, 9))
# plt.title("Frequency map")
# frequencies = som.activation_response(got.fetch_features(target_names, all))
# T = np.zeros((m_neurons, n_neurons))
# for i in range(0, m_neurons):
#     for j in range(0, n_neurons):
#         T[i, j] = len(node2sample[(j, i)])
# print(T.shape, np.min(T), np.max(T))
# got.plot_hexagons(got.scale01_(T), plt.gca(), cm.viridis)
# plt.show()

""" 4 plot target clusters """
target_names = ["surface", "depth", "volume", "krs", "SUFz", "RLDz"]
got.plot_targets(all, target_names, m_neurons, n_neurons, node2sample, sample2node, "hexagon")  # hexagon, rose
target_names = ["surface", "depth", "volume", "krs", "SUFz", "RLDz", "area"]
got.plot_targets(all, target_names, m_neurons, n_neurons, node2sample, sample2node, "rose")  # hexagon, rose

""" 5 pareto set - adds a pareto solution key to all dicts """
target_names = ["surface", "depth", "-volume", "krs", "SUFz", "RLDz"]
ind = got.pareto_list(all, target_names)
pareto = []
for i, a in enumerate(all):
    all[i]["pareto"] = ind[i]
    if ind[i] == True:
        pareto.append(a)

print("\nNumber of pareto solutions: ", np.sum(ind), len(pareto))
pareto_nodes = []
for p in pareto:
    pareto_nodes.append(p["node"])
print("Pareto solutions are located in nodes", set(pareto_nodes))
pareto_nodes = np.array(pareto_nodes, dtype = np.int64)
for i in range(0, m_neurons * n_neurons):
    print("Node {:g} = ({:g},{:g})".format(i, i % n_neurons, i // n_neurons), "has", len(pareto_nodes[pareto_nodes == i]), "solutions")

# exemplify(all, node2sample, m_neurons, n_neurons)
# exemplify(all, node2sample, m_neurons, n_neurons, vis_hairs = True, N = 3)
# exemplify_closest(all, node2sample, m_neurons, n_neurons, cluster_targets, vis_hairs = True)
# write_closest(all, node2sample, m_neurons, n_neurons, cluster_targets, vis_hairs = True)

# """ 5 parameter space regarding nodes """
# pbounds = {
#     'src_a': (3, 11),
#     'src_delay_a': (3, 14),
#     'lmax145_a': (50, 150),
#     'ln145_a': (0.5, 10.),
#     'r145_a': (0.2, 7.),
#     'lmax2_a': (5., 50.),
#     'ln2_a': (0.5, 10.),
#     'r2_a': (0.2, 7.),
#     'lmax3_a': (5., 50.),
#     'r3_a': (0.2, 7.),
#     }
#
# param_names = pbounds.keys()
# data_params = got.fetch_features(param_names, all)
# print("\nAll parameter space", data_params.shape)
# got.print_info(data_params, param_names)
# print(data_params.shape)
#
""" 6 analyze single node """
print("\nAnalyze (1,0)")
node = (1, 0)  # interesting nodes:  (1, 0) = 1 (0,1) = 4, (0,0) = 0

for node in [(0, 1)]:  # [(1, 0), (0, 1), (0, 0)]:

    ind = np.array(node2sample[node], dtype = np.int64)

    all_node, pareto = [], []
    for i in ind:
        all_node.append(all[i])
        if all[i]["pareto"]:
            pareto.append(all[i])

    target_names = ["surface", "depth", "volume", "krs", "SUFz", "RLDz"]
    data_target = got.fetch_features(target_names, all)
    print("\nTargets node", node, data_target[ind].shape)
    got.print_info(data_target[ind], target_names)

    # write out pareto solution of node
    print("Node", node, "has", len(pareto), "Pareto optimal solutions")
    with open("my_pick.txt", 'w', encoding = 'utf-8') as file:
        for line in pareto:
            file.write(line["exp_name"] + '\n')

# node_listdata = [all[i] for i in ind]
# distortion_plot(node_listdata, param_names)
#
# k = 10
# node_data = got.fetch_features(param_names, node_listdata)
# centers, dist = kmeans(node_data, k, rng = 1)
# node2sample, sample2node, node2sample_ = got.make_kmeans_maps(node_data, centers, k)
#
# for i in range(0, k):
#     ind = np.array(node2sample_[i], dtype = np.int64)
#     node_data_k = [node_listdata[i] for i in ind]
#     print("\nParameter space node", i, len(node_data_k), ind)
#     got.print_info(got.fetch_features(param_names, node_data_k), param_names)
#
# """ fin """

