"""
Macroscopic:

prints actual, potential, and cumulative transpiration, and shows optimal file number indices

list_filename ... a list of the simulation results files to be considered

Daniel Leitner, 2025
"""

import evapotranspiration as evap
import figure_style
import matplotlib.pyplot as plt
import numpy as np

Kc_maize = 1.2
Kc_soybean = 1.15

list_filename = "data/my_pick_12.txt"
path = "results (2)/"  # from results.zip

with open(list_filename, "r", encoding="utf-8") as file:
    lines = file.readlines()
lines = [line.strip() for line in lines]

print(lines)
print(len(lines))

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
print("\nactual cumulative [mm]")
print(act_)
print("\npotential cumulative [mm]")
print(pot_)
print("\nperformance [1]")
print(per_)

print("\nbest performer")
for j, envirotype in enumerate(["0", "1", "5", "36", "59"]):
    print(envirotype, ":", np.argmax(per_[:, j]))

# for i in range(0, 12):
#     plt.clf()
#     plt.bar(["0", "1", "5", "36", "59"], per_[i,:])
#     plt.ylim((0.4, 1))
#     plt.savefig('plot{:g}.png'.format(i))
