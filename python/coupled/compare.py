import numpy as np
import matplotlib.pyplot as plt


def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


SMALL_SIZE = 22
MEDIUM_SIZE = 22
BIGGER_SIZE = 22

plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

cols = ["g", "r", "b", "m", "r:"] * 10
file_names = ["classical", "reference.csv", "jan", "jan_lookup2", "rhizo", ]   
labels = ["classic sink", "3d", "steady rate", "steady lookup", "cylindric", ]  # "classic sink"

plt.figure(figsize=(20, 20), dpi=80)

t_ = np.linspace(0, 7, 1000)
plt.plot(t_, 6.4 * sinusoidal(t_), "k:", label="potential transpiration")

data = []

for i, n in enumerate(file_names): 
    d = np.loadtxt(n, delimiter=';')
    plt.plot(d[0], d[1], cols[i], linewidth=2, label=labels[i])

plt.ylabel("transpiration (cm3/day)")
plt.xlabel("time (day)")    

plt.legend(loc='upper right')
plt.show()
    
