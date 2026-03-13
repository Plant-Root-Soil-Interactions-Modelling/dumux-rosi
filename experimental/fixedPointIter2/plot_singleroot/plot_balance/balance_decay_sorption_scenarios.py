import sys
sys.path.append("../../../../../CPlantBox")
sys.path.append("../../../../../CPlantBox/src")
import plantbox as pb
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from scipy import interpolate
import pandas as pd
import sys


font = {'size'   : 16}
plt.rc('font', **font)
path2file2 = "../../scripts/results_singleroot/Exudate"


params = ['Exud_tot', 'Exud_liq', 'Exud_ads', 'Exud_decay']
linst = ['-', '--', ':', '-.']
scenarios = [['loam', 'sand'],['ini100_trans0', 'ini100_trans0.6'],['ini100_trans0.6', 'ini1000_trans0.6']]
cols= ['r', 'b', 'g', 'y']
linst = ['-', '--']

#Scenarios 
soil1 =['loam', 'loam','loam', 'loam']
soil2 =['sand', 'loam','loam', 'loam']
sorption1 = ['low', 'low', 'low', 'low']
sorption2 = ['low', 'medium', 'low', 'low']
swp1 =['100','100','100', '100'] 
swp2 =['100','100', '100', '1000'] 
trans1 = ['0', '0', '0', '0.6']
trans2 = ['0', '0', '0.6', '0.6']
titles = ['Soil type', 'Sorption', 'Transpiration', 'Initial soil water potential']
add_labels = [['Loam', 'Sand'],['Low sorption', 'Medium sorption'], ['Transpiration = 0 cm/d', 'Transpiration = 0.6 cm/d'], ['Initial SWP = -100 cm', 'Initial SWP = -1000 cm']]
dt = 20/60/24 #day

fig, ax = plt.subplots(int(len(soil1)/2),int(len(soil1)/2))
i = 0
for m in range(0, int(len(soil1)/2)): 
    for k in range(0, int(len(soil1)/2)): 
        scenario1 = soil1[i]+'_res1_sorption'+sorption1[i]+'_SWP_ini'+swp1[i]+'_trans'+trans1[i]
        scenario2 = soil2[i]+'_res1_sorption'+sorption2[i]+'_SWP_ini'+swp2[i]+'_trans'+trans2[i]
        df1 = pd.read_csv(path2file2+"/"+scenario1+"/exud.csv")
        df2 = pd.read_csv(path2file2+"/"+scenario2+"/exud.csv")
        df = [df1, df2]
        
        #exudate balance input - output
        x_sim = np.arange(0, dt*len(df1['Exud_tot'].loc[:].values),dt)  
        ax[m,k].plot(x_sim, df1['Exud_tot'].loc[:].values,color = cols[0])
        ax[m,k].plot(x_sim, df1['Exud_ads'].loc[:].values,  color = cols[1])
        ax[m,k].plot(x_sim,  df1['Exud_liq'].loc[:].values, color = cols[2])
        ax[m,k].plot(x_sim, df1['Exud_decay'].loc[:].values, color = cols[3])
        
        ax[m,k].plot(x_sim, df2['Exud_tot'].loc[:].values,linestyle = '--',color = cols[0])
        ax[m,k].plot(x_sim, df2['Exud_ads'].loc[:].values,linestyle = '--',color = cols[1])
        ax[m,k].plot(x_sim, df2['Exud_liq'].loc[:].values, linestyle = '--',color = cols[2])
        ax[m,k].plot(x_sim, df2['Exud_decay'].loc[:].values,linestyle = '--',color = cols[3])
        
        
        if i == 0: 
            lines = ax[m,k].get_lines()
            legend1 = ax[m,k].legend([lines[i] for i in np.arange(0,4)], ['Total amount of exudates', 'Sorbed exudates','Dissolved exudates', 'Decomposed exudates'], loc='upper right')
            ax[m,k].add_artist(legend1)
        
        for j in range(0, 2): 
            ax[m,k].plot([], [], color = 'k', linestyle = linst[j], label = add_labels[i][j])
        ax[m,k].legend(loc = 'upper left')

        
        # ax[m,k].set_xlabel('Time (d)')
        # ax[m,k].set_ylabel('Exudate content in soil domain (mol)', x = 0.05)
        # ax[m,k].text(0.05*np.max(x_sim), 0.6*np.max(df1['Exud_tot'].loc[:].values), titles[i])
        i+=1
fig.supxlabel('Time (d)')
fig.supylabel('Exudate content in soil domain (mol)', x = 0.05)
plt.show()

