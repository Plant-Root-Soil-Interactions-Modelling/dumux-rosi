"""
    Macroscopic:

    Mean performance plot (by envirotype)
    (actual/potential transpiration) mean and sd over all files in list_filename of each envirotype 
    
    list_filename ... a list of the simulation results files to be considered
    results files are = path + name + "_" + envirotype + ".npz", for envirotypes in ["0", "1", "5", "36", "59"]    
    
    Daniel Leitner, 2025
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import numpy as np
import matplotlib.pyplot as plt

from functional.xylem_flux import sinusoidal2

import evapotranspiration as evap

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = BIGGER_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = BIGGER_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = BIGGER_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = BIGGER_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = BIGGER_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = BIGGER_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Kc_maize = 1.2
Kc_soybean = 1.15

list_filename = "data/my_pick_(1,0).txt"
path = "results (2)/"  # from results.zip

with open(list_filename, "r", encoding = "utf-8") as file:
    lines = file.readlines()
lines = [line.strip() for line in lines]

print(len(lines))
print(lines)

area = 76 * 3

pot_ = np.zeros((len(lines), len(["0", "1", "5", "36", "59"])))
act_ = np.zeros((len(lines), len(["0", "1", "5", "36", "59"])))
for i, name in enumerate(lines):
    for j, envirotype in enumerate(["0", "1", "5", "36", "59"]):
        data = np.load(path + name + "_" + envirotype + ".npz")
        times = data["times"]
        pot_trans = data["pot_trans"]
        act_trans = data["act_trans"]
        dt = np.diff(times)
        cup = -10 * np.cumsum(np.multiply(act_trans[:-1], dt)) / area
        cum_pot = -10 * np.cumsum(np.multiply(pot_trans[:-1], dt)) / area
        act_[i, j] = cup[-1]
        pot_[i, j] = cum_pot[-1]

per_ = np.divide(act_, pot_)
print("\nactual")
print(act_)
print("\npotential")
print(pot_)
print("\npercent")
print(per_)

print("\nbest performer")
opts = []
for j, envirotype in enumerate(["0", "1", "5", "36", "59"]):
    opt_ind = np.argmax(per_[:, j])
    opts.append(opt_ind)
    print(envirotype, ":", opt_ind, lines[opt_ind])

opts = np.unique(opts)  # unique

# Plot
mean_values = np.mean(per_, axis = 0)
std_values = np.std(per_, axis = 0)
extra_points = []
for opt in opts:
    extra_points.append([per_[opt, j] for j in range(0, len(["0", "1", "5", "36", "59"]))])

x = np.arange(per_.shape[1])
plt.figure(figsize = (8, 5))
plt.bar(["0", "1", "5", "36", "59"], mean_values, yerr = std_values, capsize = 5, alpha = 0.7, color = 'skyblue', edgecolor = 'black')

for i, extras in enumerate(extra_points):
    plt.scatter(["0", "1", "5", "36", "59"], extras, zorder = 10, label = f"Number {opts[i]}")

plt.xlabel('Envirotype')
plt.ylabel('Mean performance (actual/potential)')
plt.legend()
plt.grid(True, axis = 'y', linestyle = '--', alpha = 0.6)
plt.tight_layout()
plt.show()

