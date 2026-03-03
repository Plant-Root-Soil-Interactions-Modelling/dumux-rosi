"""
Macroscopic:

Mean performance plot (by envirotype)
(actual/potential transpiration) mean and sd over all files in list_filename of each envirotype

list_filename ... a list of the simulation results files to be considered
results files are = path + name + "_" + envirotype + ".npz", for envirotypes in ["0", "1", "5", "36", "59"]

Daniel Leitner, 2026
"""

import figure_style
import matplotlib.pyplot as plt
import numpy as np


def analyze_transpiration_list(list_name, folder_name, bugged=False):
    """Analyze transpiration data from a list of files."""

    with open(list_filename, "r", encoding="utf-8") as file:
        lines = file.readlines()
    lines = [line.strip() for line in lines]

    area = 76 * 3

    pot_ = np.ones((len(lines), 1))
    act_ = np.zeros((len(lines), 1))
    for i, name in enumerate(lines):
        if bugged:
            name_ = name[:-4]
        else:
            name_ = name

        try:
            data = np.load(folder_name + name_ + "/" + name + ".npz")
            print("Opening file:", folder_name + name_ + "/" + name + ".npz", flush=True)
        except FileNotFoundError:
            print(f"File not found: {folder_name + name_ + '/' + name + '.npz'}")
            continue

        times = data["times"]
        pot_trans = data["pot_trans"]
        act_trans = data["act_trans"]
        dt = np.diff(times)
        cup = -10 * np.cumsum(np.multiply(act_trans[:-1], dt)) / area
        cum_pot = -10 * np.cumsum(np.multiply(pot_trans[:-1], dt)) / area
        act_[i] = cup[-1]
        pot_[i] = cum_pot[-1]

    per_ = np.divide(act_, pot_)
    return per_, act_, pot_


folder_names, per, namesx = [], [], []

namesx.append("Node 1 [et 0]]")
list_filename = "data/exp_name_list1.txt"  # list of experiment names
folder_name = "ex_name_list1_0_120/"  # node 1, envirotype 0, water table at 120 cm
bugged = True  # wrong subfolder name in previous version
per_, act_, pot_ = analyze_transpiration_list(list_filename, folder_name, bugged)
per.append(per_)

namesx.append("Node 1 free [et 0]")
list_filename = "data/exp_name_list1.txt"  # list of experiment names
folder_name = "exp_name_list_0_free/"  # node 1, envirotype 0, free drainage
bugged = True  # wrong subfolder name in previous version
per_, act_, pot_ = analyze_transpiration_list(list_filename, folder_name, bugged)
per.append(per_)

print("\nbest performer")
opt_ind = [np.argmax(per[i]) for i in range(len(per))]
print(opt_ind)
print("file  names of best performers:")
for i in opt_ind:
    print(lines[i])


# Plot
mean_values = [np.mean(per[i]) for i in range(len(per))]
std_values = [np.std(per[i]) for i in range(len(per))]
print("\nmean values:", mean_values)
print("std values:", std_values)

extra_points = []
for opt in opt_ind:
    extra_points.append([per[j][opt] for j in range(len(per))])

x = np.arange(per_.shape[1])
plt.figure(figsize=(8, 5))
plt.bar(namesx, mean_values, yerr=std_values, capsize=5, alpha=0.7, color="skyblue", edgecolor="black")

for i, extras in enumerate(extra_points):
    plt.scatter(namesx, extras, zorder=10, label=f"Number {opt_ind[i]}")

plt.ylabel("Mean performance (actual/potential)")
plt.legend()
plt.grid(True, axis="y", linestyle="--", alpha=0.6)
plt.tight_layout()
plt.show()
