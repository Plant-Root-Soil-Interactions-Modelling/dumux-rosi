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

left  = 0.08  # the left side of the subplots of the figure
right = 0.9    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.95      # the top of the subplots of the figure
wspace = 0.3   # the amount of width reserved for blank space between subplots
hspace = 0.4   # the amount of height reserved for white space between subplots


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

fig, ax = plt.subplots(len(soil1),2)
for i in range(0, len(soil1)): 
    scenario1 = soil1[i]+'_res1_sorption'+sorption1[i]+'_SWP_ini'+swp1[i]+'_trans'+trans1[i]
    scenario2 = soil2[i]+'_res1_sorption'+sorption2[i]+'_SWP_ini'+swp2[i]+'_trans'+trans2[i]
    df1 = pd.read_csv(path2file2+"/"+scenario1+"/exud.csv")
    df2 = pd.read_csv(path2file2+"/"+scenario2+"/exud.csv")
    df = [df1, df2]
    
    #exudate balance input - output
    for j in range(0,2): 
        a = df[j]['Exud_liq'].loc[:].values/df[j]['Exud_tot'].loc[:].values*100
        b = df[j]['Exud_decay'].loc[:].values/df[j]['Exud_tot'].loc[:].values*100
        c = df[j]['Exud_ads'].loc[:].values/df[j]['Exud_tot'].loc[:].values*100
        a[np.isnan(a)] = 100
        b[np.isnan(b)] = 0
        c[np.isnan(c)] = 0
        
        data = np.vstack((a,b,c)).T
        index_ = np.arange(0, dt*len(df1['Exud_tot'].loc[:].values),dt) 
        index = ["%.2f" % x for x in index_]
        df_fin = pd.DataFrame(data,index=index, columns=['dissolved exudates', 'decomposed exudates', 'sorbed exudates'])

        df_fin.plot.area(ax=ax[i,j])
        if i>0 or j == 1: 
            ax[i,j].get_legend().remove()
        ax[i,j].set_title(add_labels[i][j])   

fig.supxlabel('Time (d)', y = 0.01)
fig.supylabel('Share of total root exudats (%)', x = 0.03)
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
plt.show()

