import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import pandas as pd
import sys


font = {'size'   : 18}
plt.rc('font', **font)
path2file2 = "../../scripts/results/Exudate"

starttime = 10 
dt = 20/60/24 #d
params = ['Exud_tot', 'Exud_liq', 'Exud_ads', 'Exud_decay']
linst = ['-', '--', ':', '-.']
scenarios = ['withdecay_withsorption', 'withdecay_nosorption', 'nodecay_nosorption']

fig, ax = plt.subplots(2,len(scenarios))
for j in range(0, len(scenarios)): 
    scenario = 'loam_4_exu_'+ scenarios[j]
    df = pd.read_csv(path2file2+"/"+scenario+"/exud.csv")

    exud_in = []
    file_in = open(path2file2+"/"+scenario+"/Q_Exud_tot.txt", 'r')
    for y in file_in.read().split('\n'):
        exud_in.append(y)
    exud_in = np.array(exud_in[:-1]).flatten()
    exud_in = exud_in.astype(float)
    x_sim = np.arange(0, dt*(len(exud_in)),dt)  


    #exudate balance input - output
    exud_out = df['Exud_tot'].loc[:].values+df['Exud_decay'].loc[:].values  #mol
    # exud_simout1 = df['Exud_liq'].loc[:].values+df['Exud_decay'].loc[:].values +df['Exud_ads_soil'].loc[:].values# +df['Exud_ads'].loc[:].values #mol
    ax[0,j].plot(x_sim, exud_in, label = 'prescribed incoming\n cumulative plant exudation') 
    ax[0,j].plot(x_sim, exud_out[:len(x_sim)], label = 'Exud_tot + Exud_decay')
    # ax[0].plot(x_sim, exud_simout1, label = 'sims output')
    ax[0,j].legend()
    ax[0,j].set_title(scenarios[j])
    
    #plot balance simulation output
    #Exud_liq +Exud_ads = Exud_tot
    for i in range(0, len(params)): 
        y = df[params[i]].loc[:].values
        x = np.arange(0, len(y))*dt
        ax[1,j].plot(x,y,linestyle = linst[i], label = params[i])
        # ax[0].set_ylim([0,1e-5])
    ax[1,j].legend()

fig.supxlabel('Time (d)')
fig.supylabel('Exudate content in soil domain (mol)', x = 0.05)

plt.show()