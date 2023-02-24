import numpy as np
import matplotlib.pyplot as plt

name1 = "dry_3d"
name2 = "dry_3d_rhizo"  # "dry_3d_rhizo"
names = [name1, name2]  # works for 2,3,...
xrange = [-8, 2.]  # for second plot
yrange = [0., 350.]  # for first cumulative plot
title_names = ["Steady rate approach, 3d macrogrid, dry", "Cylindrical models, 3d macrogrid, dry"]

# name1 = "dry_1d"
# name2 = "dry_1d_rhizo"
# names = [name1, name2]  # works for 2,3,...
# xrange = [-8, 2.]  # for second plot
# yrange = [0., 350.]  # for first cumulative plot
# title_names = ["Steady rate approach, 1d macrogrid, dry", "Cylindrical models, 1d macrogrid, dry"]

# name1 = "wet_3d2"
# name2 = "wet_3d2"  # "wet_3d_rhizo"
# names = [name1, name2]
# xrange = [-100, 10]  # for second plot
# yrange = [0., 3500.]  # for first cumulative plot
# title_names = ["Steady rate approach, 3d macrogrid, wet", "Cylindrical models, 3d macrogrid, wet"]

# name1 = "wet_1d"
# name2 = "wet_1d_rhizo"
# names = [name1, name2]
# xrange = [-100, 10]  # for second plot
# yrange = [0., 3500.]  # for first cumulative plot
# title_names = ["Steady rate approach, 1d macrogrid, wet", "Cylindrical models, 1d macrogrid, wet"]

# name1 = "dry_1d"
# name2 = "dry_1d_rhizo"
# name3 = "dry_3d"
# name4 = "dry_3d_rhizo"
#
# title_names = ["Steady rate approach, 1d macrogrid, dry",
#                "Cylindrical models, 1d macrogrid, dry",
#                "Steady rate approach, 3d macrogrid, dry",
#                "Cylindrical models, 3d macrogrid, dry"]
#
# names = [name1, name2, name3, name4]  # works for 2,3,...
# xrange = [-8, 2.]  # for second plot
# yrange = [0., 350.]  # for first cumulative plot

n = len(names)

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


""" transpiration plot """
trans = 0.5 * 15 * 75
potential_trans = lambda t: trans * sinusoidal(t)

# load data
data = [np.loadtxt(n, delimiter=';')  for n in names]

fig, ax = plt.subplots(n, 1, figsize=(15, 10))

for i in range(0, n):
    t = data[i][0]
    y = data[i][1]
    ax[i].plot(t, potential_trans(t), 'k', label="potential transpiration")  # potential transpiration
    ax[i].plot(t, y, 'g', label=" actual transpiration")  # actual transpiration  according to soil model
    # ax[i].set_xlabel("Time [d]")
    ax[i].set_title(title_names[i])
    ax[i].set_ylabel("Transpiration [cm$^3$ d$^{-1}$]")
    ax[i].set_ylim([0., 1150.])
    ax[i].legend(loc='upper left')
    ax2 = ax[i].twinx()
    dt = np.diff(t)
    so = np.array(y)
    cup = np.cumsum(np.multiply(so[:-1], dt))
    ax2.plot(t[1:], cup, 'c--', label="cumulative transpiration")  # cumulative transpiration (neumann)
    ax2.set_ylabel("Cumulative soil uptake [cm$^3$]")
    ax2.set_ylim(yrange)
    ax2.legend(loc='lower right')
    print("cumulative uptake", cup[-1], "cm3")

plt.show()

""" 1d sink plot """
fig, axes = plt.subplots(1, n, figsize=(15, 15))

layers = 2 * 55
depth = -110
days = 7  # (measurements every 6h)

labels = ["peak", "redistribution"]
reds = [[1. / (i / 3 + 1), 0., 0.] for i in range(0, days)]
greens = [[0., 1. / (i / 3 + 1), 0.] for i in range(0, days)]

for j in range(0, n):

    data = np.load(names[j] + "_sink.npy")
    y = np.linspace(0, depth, layers + 1)
    y = 0.5 * (y[1:] + y[:-1])  # layer mids

    print(data.shape)

    for i in range(0, days):
        if i == 0:
            axes[j].plot(data[2 + i * 4,:], y, label=labels[0], color=reds[i])
            axes[j].plot(data[4 + i * 4,:], y, label=labels[1], color=greens[i])
        else:
            axes[j].plot(data[2 + i * 4,:], y, color=reds[i])
            axes[j].plot(data[i * 4,:], y, color=greens[i])

        axes[j].plot([0, 0], [-110, 0], 'k:')
        axes[j].set_title(title_names[j])
        axes[j].set_ylabel("Depth (cm)")
        axes[j].set_xlabel(r"Sink (cm$^3$/day)")
        axes[j].set_xlim(xrange)
        axes[j].legend()

plt.margins()
plt.show()

