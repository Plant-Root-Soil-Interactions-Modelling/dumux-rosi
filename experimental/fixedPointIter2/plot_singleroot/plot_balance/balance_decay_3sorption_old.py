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

def comp_exu(rs,f,sim_end):
    exu_tot = np.zeros((sim_end))
    for j in range(0,sim_end):
        
        rs.simulate(1, True);
        
        polylengths = rs.getParameter("length")
        radii = rs.getParameter("radius")
        types = rs.getParameter("type")
        polylines = rs.getPolylines()

        kex_ = f(j)
        sf = []
        for i in range(0, len(polylines)):
            a = radii[i]
            roottype = int(types[i])
            l_ = 0
            for k in range(0, len(polylines[i])-1):
                m = polylines[i][-1-k]
                n = polylines[i][-2-k]
                p0 = np.array([m.x, m.y, m.z])
                p1 = np.array([n.x, n.y, n.z])
                l  = np.linalg.norm(p0 - p1)
                l_ = l_+l
                #tip exudation rate 
                if l_<3.5:
                    kexu = kex_
                    c = 2
                #base exudation rate 
                else:
                    kexu = kex_/2
                    c = 1
                #if growth has already stopped (95% of total length reached) 
                if lmax[roottype]*0.99>= polylengths[i] or j>times[2]:
                    #print('REACHED')
                    kexu = kex_/2
                    c = 1
                #if artificial shoot 
                if roottype == 0:
                    kexu = 0
                    c = 0

                sf.append(2 * np.pi * a * l * 1.e-4 * kexu * 1.e3 / g_per_mol) # mol/day/root
        exu_tot[j] = np.sum(sf)
    return exu_tot

font = {'size'   : 18}
plt.rc('font', **font)
path2file2 = "../../scripts/results_singleroot/Exudate"

starttime = 1
endtime = 10
g_per_mol = 12.01
dt = 20/60/24 #d
params = ['Exud_tot', 'Exud_liq', 'Exud_ads', 'Exud_decay']
linst = ['-', '--', ':', '-.']
scenarios = ['low', 'medium', 'high']
cols= ['r', 'b', 'g']
linst = ['-', '--', ':']

#Root system direct 
times = [0, 40, 63, 98, 154]
exu_prop = np.array([0.001,0.001,0.00055,0.00039,0.00045])#[kg C/(m2 root surface  day)] 
f = interpolate.interp1d(times, exu_prop)  
rs = pb.MappedRootSystem()
rs.readParameters("../../../../../CPlantBox/modelparameter/structural/rootsystem/single_root.xml")
lmax = []
for pp in rs.getRootRandomParameter():
    lmax.append(pp.lmax)
rs.setSeed(1)
rs.initializeLB(5,4)
exu_tot = comp_exu(rs,f,endtime)
x_direct = np.linspace(starttime, endtime,endtime-starttime+1)-starttime
y_direct = list(np.cumsum(exu_tot[starttime-1:endtime-1]))
num = [0]
y_direct = np.array([0]+y_direct)


for i in range(0, len(scenarios)): 
    scenario = 'day1_loam_1_exu_withdecay_withsorption'+scenarios[i]
    df = pd.read_csv(path2file2+"/"+scenario+"/exud.csv")
    
    exud_in = [0.]
    file_in = open(path2file2+"/"+scenario+"/Q_Exud_tot.txt", 'r')
    for y in file_in.read().split('\n'):
        exud_in.append(y)
    exud_in = np.array(exud_in[:-2]).flatten()
    exud_in = exud_in.astype(float)
    x_sim = np.arange(0, dt*(len(exud_in)),dt)  
    
    #exudate balance input - output
    exud_tot = df['Exud_tot'].loc[:].values 

    if i == 0: 
        plt.plot(x_sim-dt+starttime, exud_in, '.', label = 'prescribed incoming\n cumulative plant exudation') 
        plt.plot(x_sim-dt+starttime, exud_tot[:len(x_sim)],label = 'Exud_tot')
    plt.plot(x_sim-dt+starttime, df['Exud_ads'].loc[:].values[:len(x_sim)],  color = cols[0], linestyle = linst[i], label = 'Exud_sorbed')
    plt.plot(x_sim-dt+starttime,  df['Exud_liq'].loc[:].values[:len(x_sim)], color = cols[1], linestyle = linst[i], label = 'Exud_liq')
    plt.plot(x_sim-dt+starttime, df['Exud_decay'].loc[:].values[:len(x_sim)], color = cols[2], linestyle = linst[i], label = 'Exud_decay')
    if i == 0: 
        plt.plot(np.nan, np.nan, color = 'k', linestyle = linst[0], label = 'Carboxylates')
        plt.plot(np.nan, np.nan, color = 'k', linestyle = linst[1], label = 'Sugars')
        plt.plot(np.nan, np.nan, color = 'k', linestyle = linst[2], label = 'Amino acids')
        plt.legend()
        
plt.xlabel('Time (d)')
plt.xlim(starttime,dt*(len(exud_in))*1.1+starttime)
plt.ylim(-0.1*max( df['Exud_tot'].loc[:].values[:len(x_sim)]),max( df['Exud_tot'].loc[:].values[:len(x_sim)])*1.1)
plt.ylabel('Exudate content in soil domain (mol)', x = 0.05)
plt.show()

