"""
Macroscopic:

Mean performance plot (by envirotype)
(actual/potential transpiration) mean and sd over all files in list_filename of each envirotype


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
    filenames = []
    for i, name in enumerate(lines):
        if bugged:
            name_ = name[:-4]
        else:
            name_ = name
        file_name = folder_name + name_ + "/" + name + ".npz"
        filenames.append(file_name)
        try:
            data = np.load(file_name)
            print(f"Opening file: {file_name}", flush=True)
        except FileNotFoundError:
            print(f"File not found: {file_name}")
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
    return per_, act_, pot_, filenames


all_names, per, namesx = [], [], []

namesx.append("Node 1, ET0")
list_filename = "data/exp_name_list1.txt"  # list of experiment names
folder_name = "exp_name_list1_0_120/"  # node 1, envirotype 0, water table at 120 cm
bugged = True  # wrong subfolder name in previous version
per_, act_, pot_, filenames_ = analyze_transpiration_list(list_filename, folder_name, bugged)
per.append(per_)
all_names.append(filenames_)

namesx.append("Node 1  drainage, ET0")
list_filename = "data/exp_name_list1.txt"  # list of experiment names
folder_name = "exp_name_list1_0_free/"  # node 1, envirotype 0, free drainage
bugged = False  # wrong subfolder name in previous version
per_, act_, pot_, filenames_ = analyze_transpiration_list(list_filename, folder_name, bugged)
per.append(per_)
all_names.append(filenames_)

namesx.append("Node 2, ET0")
list_filename = "data/exp_name_list2.txt"  # list of experiment names
folder_name = "exp_name_list2_0_120/"  # node 2, envirotype 0, water table at 120 cm
bugged = True  # wrong subfolder name in previous version
per_, act_, pot_, filenames_ = analyze_transpiration_list(list_filename, folder_name, bugged)
per.append(per_)
all_names.append(filenames_)

namesx.append("Node 2, free, ET0")
list_filename = "data/exp_name_list2.txt"  # list of experiment names
folder_name = "exp_name_list2_0_free/"  # node 2, envirotype 0, free drainage
bugged = False  # wrong subfolder name in previous version
per_, act_, pot_, filenames_ = analyze_transpiration_list(list_filename, folder_name, bugged)
per.append(per_)
all_names.append(filenames_)

print("\nbest performer")
opt_ind = [np.argmax(per[i]) for i in range(len(per))]
print(opt_ind)
print("file  names of best performers:")
for ii, i in enumerate(opt_ind):
    print(all_names[ii][i])

# Plot
mean_values = [np.mean(per[i]) for i in range(len(per))]
std_values = [np.std(per[i]) for i in range(len(per))]
print("\nmean values:", mean_values)
print("std values:", std_values)


x = np.arange(per_.shape[1])
plt.figure(figsize=(8, 5))
plt.bar(namesx, mean_values, yerr=std_values, capsize=5, alpha=0.7, color="skyblue", edgecolor="black")

# extra points ...
optimals = [float(per[i][opt_ind[i]]) for i in range(len(per))]
yy = [
    [optimals[0], float(per[1][opt_ind[0]])],
    [optimals[1], float(per[0][opt_ind[1]])],
    [optimals[2], float(per[3][opt_ind[2]])],
    [optimals[3], float(per[2][opt_ind[3]])],
]
xx = [
    [namesx[0], namesx[1]],
    [namesx[1], namesx[0]],
    [namesx[2], namesx[3]],
    [namesx[3], namesx[2]],
]


for i in range(len(namesx)):
    print(xx[i], yy[i])
    plt.scatter(xx[i], yy[i], zorder=10, label=f"Number {opt_ind[i]}")

# reference points
plt.scatter(namesx, [0.7317886554421104, 0.5234493369758094, 0.7317886554421104, 0.5234493369758094], zorder=10, label="reference")

plt.ylabel("Mean performance (actual/potential)")
plt.legend()
plt.grid(True, axis="y", linestyle="--", alpha=0.6)
plt.tight_layout()
plt.show()
