"""
total nitrate in soil over time
"""
import sys;
sys.path.append("../modules/");
sys.path.append("../../../CPlantBox/");
sys.path.append("../../../CPlantBox/src/python_modules")
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import scenario_setup as scenario
import math

name = "maize"
str_ = "_cyl3"
fname = "solute_" + name + str_
path = "results/"

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

""" load data """
data = np.load(path + fname + ".npz", allow_pickle=True)
dist = data['arr_0']
conc = data['arr_1']
l = data['arr_2']
simtime = np.shape(dist)[0]

for i in range(0, int(simtime)):
    
    dist_ = np.array(dist[i]).flatten()
    l_ = np.array(l[i]).flatten()
    conc_ = np.array(conc[i]).flatten()
    vol = np.zeros((len(conc_)))
    
    for j in range(0,len(dist_)):
        if dist_[j] == 0: 
            x_prev = 0
        else:
            x_prev = dist_[j-1]
        
        vol[j] = l_[j]*2*math.pi*(x_prev+(dist_[j]-x_prev)/2)*(dist_[j]-x_prev)

    x_ = vol
    y_ = conc_
    ind = np.argsort(-1*y_)
    x = np.cumsum(x_[ind])
    y = y_[ind]
    plt.plot(x,y,label = 'day'+str(i))

plt.legend()
plt.xlabel("Soil volume  [cm^3]")
plt.ylabel("C concentration [g / cm^3]")
plt.show()

