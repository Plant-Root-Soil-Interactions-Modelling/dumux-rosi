import numpy as np
import matplotlib.pyplot as plt

# name1 = "dry_small_1d_rhizo"
# name2 = "dry_small_1d"
# name3 = "dry_small_agg"

name1 = "wet_small_1d_rhizo"
name2 = "wet_small_1d"
name3 = "wet_small_agg"

names = [name1, name2, name3]

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

data = []
for n in names:
    data.append(np.loadtxt(n, delimiter=';'))

fig, ax = plt.subplots(len(names), 1, figsize=(15, 15))

for i in range(0, len(names)):
    t = data[i][0]
    y = data[i][1]
    ax[i].plot(t, potential_trans(t), 'k', label="potential transpiration")  # potential transpiration
    ax[i].plot(t, y, 'g', label=" actual transpiration")  # actual transpiration  according to soil model
    ax[i].set_xlabel("Time [d]" + " (" + names[i] + ")")
    ax[i].set_ylabel("Transpiration $[cm^3 d^{-1}]$")
    ax2 = ax[i].twinx()
    dt = np.diff(t)
    so = np.array(y)
    ax2.plot(t[1:], np.cumsum(np.multiply(so[:-1], dt)), 'c--')  # cumulative transpiration (neumann)
    ax2.set_ylabel("Cumulative soil uptake $[cm^3]$")
    ax[i].legend(loc='upper right')

plt.show()

""" 1d sink plot """
fig, axes = plt.subplots(1, len(names), figsize=(15, 15))

layers = 55
depth = -110
days = 7  # (measurements every 6h)

for j in range(0, len(names)):

    data = np.load(names[j] + "_sink.npy")
    # print(data.shape)
    y = np.linspace(0, depth, layers + 1)
    y = 0.5 * (y[1:] + y[:-1])  # layer mids

    labels = ["peak", "redistribution"]

    reds = [[1. / (i / 3 + 1), 0., 0.] for i in range(0, days)]
    greens = [[0., 1. / (i / 3 + 1), 0.] for i in range(0, days)]
    
    axes[j].plot([0, 0], [0, -110], 'k:')

    for i in range(0, days - 1):
        if i == 0:
            axes[j].plot(data[2 + i * 4,:], y, label=labels[0], color=reds[i])
            axes[j].plot(data[4 + i * 4,:], y, label=labels[1], color=greens[i])
        else:
            axes[j].plot(data[2 + i * 4,:], y, color=reds[i])
            axes[j].plot(data[i * 4,:], y, color=greens[i])

        axes[j].set_title(names[j])
        axes[j].set_ylabel("Depth (cm)")
        axes[j].set_xlabel(r"Sink (cm$^3$/day)")
#        axes[j].set_xlim([-8, 1.5])
#        axes[j].set_ylim([-12, 0])
        axes[j].legend()

plt.show()

