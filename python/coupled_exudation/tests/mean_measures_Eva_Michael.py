import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import pandas as pd
import numpy as np

font = {'size'   : 18}
matplotlib.rc('font', **font)


df = pd.read_csv("../data_magda/mean_measures_Eva_Michael.csv")
print(df.columns.values.tolist()) 
DAS = df["DAS"].loc[:].values


fig, ax = plt.subplots(3,2, figsize = (18, 10))

#optimized = np.array([11147.26, 20865.86, 39905.78, 56972.49])
#optimized1 = np.array([10002.89, 19584.2,  38741.21, 58314.23])
#optimized2 = np.array([13594.6,  23614.25, 42636.57, 53824.92])
#optimized3 = np.array([13348.97, 22165.34, 38868.31,38868.31 ])
#optimized3 = np.array([15032.78, 26361.75, 43329.03, 43329.03])
#optimized3 = np.array([19063.93, 28156.36, 40717.52, 40717.52])
optimized3 = np.array([ 3076.67707917, 13560.32198841, 59895.71241296, 59895.71241296])
optimized2 = np.array([1223.83740793, 13125.69800546, 55855.62879929,55855.62879929 ])

#optimized_diam = np.array([0.41786875, 0.33700683, 0.28101047,0.28101047 ])
optimized_diam = np.array([0.41310032, 0.32377058, 0.2718736, 0.2718736 ])
optimized_diam1 = np.array([0.41582958, 0.32561655, 0.26814852, 0.26814852])
optimized_diam2 = np.array([0.41008175, 0.3496782,  0.29133455, 0.29133455])




ax[0,0].plot(df['DAS'], df['Length']/100, color = 'r', label = 'mean of measured values')
ax[0,0].fill_between(df['DAS'], (df['Length']-df['Length std'])/100, (df['Length']+df['Length std'])/100, color = 'r', alpha = 0.2, label = 'standard deviation of measured values')
ax[0,0].plot(df['DAS'], optimized2/100, color = 'r', linestyle = '--', label = 'simulation')
#ax[0,0].plot(df['DAS'], optimized2/100, color = 'r', linestyle = ':')
#ax[0,0].plot(df['DAS'], optimized3/100, color = 'r', linestyle = '-.')
ax[0,0].set_ylabel('Total root length\n (m / plant)')
ax[0,0].text(50, 600, '(a)', fontsize=16, va='top', ha='right')

ax[0,1].plot(df['DAS'], df['RootDiameter'], color = 'g')
ax[0,1].fill_between(df['DAS'], df['RootDiameter']-df['RootDiameter std'], df['RootDiameter']+df['RootDiameter std'], color = 'g', alpha = 0.2)
#ax[1].plot(df['DAS'], optimized_diam, color = 'g', linestyle = '--')
#ax[1].plot(df['DAS'], optimized_diam1, color = 'g', linestyle = ':')
ax[0,1].plot(df['DAS'], optimized_diam2, color = 'g', linestyle = '--')
ax[0,1].set_ylabel('Root diameter\n (mm)')
ax[0,1].text(50, 0.45, '(b)', fontsize=16, va='top', ha='right')

ax[1,0].plot(df['DAS'], df['Phenolics']*10**3, color = 'b')
ax[1,0].fill_between(df['DAS'], df['Phenolics']*10**3-df['Phenolics std']*10**3, df['Phenolics']*10**3+df['Phenolics std']*10**3, color = 'b', alpha = 0.2)
ax[1,0].set_ylabel('Phenolics\n ($\mu$mol / plant / h)')
ax[1,0].text(50, 33, '(c)', fontsize=16,  va='top', ha='right')

ax[1,1].plot(df['DAS'], df['Sugars']*10**3, color = 'c')
ax[1,1].fill_between(df['DAS'], df['Sugars']*10**3-df['Sugars std']*10**3, df['Sugars']*10**3+df['Sugars std']*10**3, color = 'c', alpha = 0.2)
ax[1,1].set_ylabel('Sugars\n ($\mu$mol / plant / h)')
ax[1,1].text(50, 280, '(d)', fontsize=16, va='top', ha='right')

ax[2,0].plot(df['DAS'], df['AminoAcids']*10**3, color = 'm')
ax[2,0].fill_between(df['DAS'], df['AminoAcids']*10**3-df['AminoAcids std']*10**3, df['AminoAcids']*10**3+df['AminoAcids std']*10**3, color = 'm', alpha = 0.2)
ax[2,0].set_ylabel('Amino acids\n ($\mu$mol / plant / h)')
ax[2,0].text(50, 33, '(e)', fontsize=16,  va='top', ha='right')

# zip joins x and y coordinates in pairs
#label = label_exud[i+2*j]
#kk = 0
#for x,y in zip(y['DAS'], y['exud']):
#    ax[2].annotate(label[kk],(x,y),textcoords="offset points",xytext=(0,10),ha='center')
#    kk = kk+1
        
for i in range(0,3):
    for j in range(0,2): 
        ax[i,j].set_xlabel('Time (d)')


handles, labels = ax[0,0].get_legend_handles_labels()
order = [0,2,1]
ax[2,1].legend([handles[idx] for idx in order],[labels[idx] for idx in order])
ax[2,1].axis('off')
    
#ax[0].legend()
plt.tight_layout()
plt.show()



