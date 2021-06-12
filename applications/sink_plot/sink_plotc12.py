import numpy as np
import matplotlib.pyplot as plt


SMALL_SIZE = 14
MEDIUM_SIZE = 14
BIGGER_SIZE = 14
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

fig, axes = plt.subplots(1, 1, figsize=(7, 7))
axes = [axes] # in case there is only 1

data = np.load("sink1d.npy")
y = np.linspace(0, -15, 16)
y = 0.5*(y[1:]+y[:-1]) # layer mids

print(data.shape)
labels = ["peak","redistribution"]

reds = [[1./(i/3+1), 0., 0.] for i in range(0,7)]
greens = [[0.,1./(i/3+1), 0.] for i in range(0,7)]

for i in range(0,7):
    if i == 0:
        axes[0].plot(data[2+i*4,:], y, label = labels[0], color = reds[i])
        axes[0].plot(data[4+i*4,:], y, label = labels[1], color = greens[i])
    else: 
        axes[0].plot(data[2+i*4,:], y, color = reds[i])
        axes[0].plot(data[i*4,:], y, color = greens[i])
    
    axes[0].set_ylabel("Depth (cm)")
    axes[0].set_xlabel(r"Sink (cm$^3$/day)")
    axes[0].legend()

plt.show()