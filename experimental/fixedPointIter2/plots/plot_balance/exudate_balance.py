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



path2file2 = "../../scripts/results/Exudate"
scenario = "loam_4_exu_nodecay"
df = pd.read_csv(path2file2+"/"+scenario+"/exud.csv")
starttime = 10 
dt = 20/60/24 #d
params = ['Exud_tot', 'Exud_liq', 'Exud_ads', 'Exud_decay']
linst = ['-', '--', ':', '-.']

#balance with computed exudation rates and Q_exud_i.txt
path2file = "../../../../../CPlantBox/tutorial/parameterization_exudatepaper/Optimization_RS_realistic/"
data = np.load(path2file+"plant_exud_rates.npz", allow_pickle=True)
g_per_mol_C = 12.0107 #g/mol
exud_calc_ = np.mean(data['arr_0'][(starttime-1):(starttime+1),2]/g_per_mol_C) #g/d/plant --> mol/d/plant
# exud_input_ = (data['arr_0'][starttime-1,2]/g_per_mol_C) #g/d/plant --> mol/d/plant


exud_in = []
file_in = open(path2file2+"/"+scenario+"/Q_Exud_tot.txt", 'r')
for y in file_in.read().split('\n'):
    exud_in.append(y)
exud_in = np.array(exud_in[:-1]).flatten()
exud_in = exud_in.astype(float)
x_sim = np.arange(0, dt*(len(exud_in)),dt)  
exud_calc = np.repeat(exud_calc_*dt, len(exud_in))

fig, ax = plt.subplots(3,1)
#plot simulated balance
#Exud_liq +Exud_ads = Exud_tot
for i in range(0, len(params)): 
    y = df[params[i]].loc[:].values
    x = np.arange(0, len(y))*dt
    ax[0].plot(x,y,linestyle = linst[i], label = params[i])
    # ax[0].set_ylim([0,1e-5])
ax[0].legend()

#balance input - partition
#Exud_tot + Exud_decay = exud 
exud_out = df['Exud_tot'].loc[:].values+df['Exud_decay'].loc[:].values  #mol
# exud_simout1 = df['Exud_liq'].loc[:].values+df['Exud_decay'].loc[:].values +df['Exud_ads_soil'].loc[:].values# +df['Exud_ads'].loc[:].values #mol
ax[1].plot(x_sim, exud_in, label = 'Exud in') 
ax[1].plot(x_sim, exud_out[:-1], label = 'Exud out')
# ax[1].plot(x_sim, exud_simout1, label = 'sims output')
ax[1].legend()

#balance with computed exudation rates and Q_exud_tot.txt
ax[2].plot(x_sim, exud_in, label = 'sims') 
ax[2].plot(x_sim, np.cumsum(exud_calc), label = 'input')
ax[2].legend()
plt.show()