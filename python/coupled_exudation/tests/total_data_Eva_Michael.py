import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
#import statsmodels.stats.multicomp as mct_ind
import pandas as pd
import numpy as np
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import sys

df = pd.read_csv("../data_magda/alldata_Eva_Michael.csv")
#print(df.columns.values.tolist()) 
DAS = df["DAS"].loc[:].values

soils = ['L', 'S']
gt = ['WT', 'rth3']
times = np.unique(DAS)
print(times) 
params = ['RootDiameter', 'Length', 'Total exudation', 'Phenolics', 'Sugars','AminoAcids']
p = 0.05

for i in range(5,6): #len(params)):
    print(params[i]) 
    for j in range(0, len(times)):

        L_WT = df[(df['DAS']==times[j]) & (df['Substrate']=='L') & (df['Genotype']=='WT')][params[i]]
        L_rth3 = df[(df['DAS']==times[j]) & (df['Substrate']=='L') & (df['Genotype']=='rth3')][params[i]]
        S_WT = df[(df['DAS']==times[j]) & (df['Substrate']=='S') & (df['Genotype']=='WT')][params[i]]
        S_rth3 = df[(df['DAS']==times[j]) & (df['Substrate']=='S') & (df['Genotype']=='rth3')][params[i]]
        
        result = f_oneway(L_WT, L_rth3, S_WT, S_rth3)
        if result[1]<=p:
            print('DAS= ', times[j])
            values = np.hstack((L_WT.values,L_rth3.values,S_WT.values,S_rth3.values))
            df1 = pd.DataFrame({'score': values,
                   'group': np.repeat(['L_WT', 'L_rth3', 'S_WT','S_rth3'], repeats=len(L_WT.values))})
            
            tukey = pairwise_tukeyhsd(endog=df1['score'],groups=df1['group'], alpha=p)
            print(tukey)

        




