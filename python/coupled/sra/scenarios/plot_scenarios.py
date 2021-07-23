import numpy as np
import matplotlib.pyplot as plt

name1 = "dry_3d"
name2 = "dry_3d_rhizo"
names = [name1, name2]

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

data1 = np.loadtxt(name1, delimiter=';')  # np.vstack((x_, -np.array(y_))),
data2 = np.loadtxt(name2, delimiter=';')
data = [data1, data2]

fig, ax = plt.subplots(2, 1, figsize=(15, 10))

for i in range(0, 2):
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
fig, axes = plt.subplots(1, 2, figsize=(15, 15))

layers = 55
depth = -110
days = 7  # (measurements every 6h)  

for j in range(0, 2):

    data = np.load(names[j] + "_sink.npy")
    y = np.linspace(0, depth, layers + 1)
    y = 0.5 * (y[1:] + y[:-1])  # layer mids
    
    print(data.shape)
    labels = ["peak", "redistribution"]
    
    reds = [[1. / (i / 3 + 1), 0., 0.] for i in range(0, days)]
    greens = [[0., 1. / (i / 3 + 1), 0.] for i in range(0, days)]
    
    for i in range(0, days):
        if i == 0:
            axes[j].plot(data[2 + i * 4,:], y, label=labels[0], color=reds[i])
            axes[j].plot(data[4 + i * 4,:], y, label=labels[1], color=greens[i])
        else: 
            axes[j].plot(data[2 + i * 4,:], y, color=reds[i])
            axes[j].plot(data[i * 4,:], y, color=greens[i])
        
        axes[j].set_title(names[j])
        axes[j].set_ylabel("Depth (cm)")
        axes[j].set_xlabel(r"Sink (cm$^3$/day)")
        axes[j].set_xlim([-8, 1.5])
#         axes[j].set_ylim([-12, 0])
        axes[j].legend()

plt.show()

