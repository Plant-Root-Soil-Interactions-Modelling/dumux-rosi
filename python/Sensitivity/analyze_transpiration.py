"""
    transpiration plot (one column, number of rows as number of filenames)
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import numpy as np
import matplotlib.pyplot as plt

from functional.xylem_flux import sinusoidal2

import evapotranspiration as evap

Kc_maize = 1.2
Kc_soybean = 1.15

list_filename = "data/my_pick_12.txt"
path = "results/"

with open(list_filename, "r", encoding = "utf-8") as file:
    lines = file.readlines()
lines = [line.strip() for line in lines]

print(lines)

area = 76 * 3

path = "results/"

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

print("\nactual")
print(act_)
print("\npotential")
print(pot_)
print("\npercent")
print(np.divide(act_, pot_))
