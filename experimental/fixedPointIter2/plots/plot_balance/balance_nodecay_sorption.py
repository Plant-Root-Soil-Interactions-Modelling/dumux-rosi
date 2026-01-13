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
path2file2 = "../../scripts/results/Exudate"

starttime = 3
endtime = 10
g_per_mol = 12.01
dt = 20/60/24 #d
params = ['Exud_tot', 'Exud_liq', 'Exud_ads', 'Exud_decay']
linst = ['-', '--', ':', '-.']

scenario = 'loam_2_exu_nodecay_withsorption'
df = pd.read_csv(path2file2+"/"+scenario+"/exud.csv")

exud_in = []
file_in = open(path2file2+"/"+scenario+"/Q_Exud_tot.txt", 'r')
for y in file_in.read().split('\n'):
    exud_in.append(y)
exud_in = np.array(exud_in[:-1]).flatten()
exud_in = exud_in.astype(float)
x_sim = np.arange(0, dt*(len(exud_in)),dt)  

#Root system direct 
times = [0, 40, 63, 98, 154]
exu_prop = np.array([0.001,0.001,0.00055,0.00039,0.00045])#[kg C/(m2 root surface  day)] 
f = interpolate.interp1d(times, exu_prop)  
rs = pb.MappedRootSystem()
rs.readParameters("../../../../../CPlantBox/modelparameter/structural/rootsystem/RS_optimized_field_loam.xml")
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

#exudate balance input - output
exud_out = df['Exud_tot'].loc[:].values+df['Exud_decay'].loc[:].values  #mol
print('exu_tot')
print('exud_in',exud_in)
print('exud_out[:len(x_sim)]',exud_out[:len(x_sim)])
print('Exud_tot', df['Exud_tot'].loc[:].values,'Exud_decay',df['Exud_decay'].loc[:].values)
print('Exud_ads', df['Exud_ads'].loc[:].values,'Exud_liq',df['Exud_liq'].loc[:].values)
print('diff', df['Exud_tot'].loc[:].values - df['Exud_ads'].loc[:].values - df['Exud_liq'].loc[:].values)
print('y_direct',y_direct*dt)

plt.plot(x_sim, exud_in, '.', label = 'prescribed incoming\n cumulative plant exudation') 
plt.plot(x_sim, exud_out[:len(x_sim)], '.', label = 'Exud_tot + Exud_decay')
plt.plot(x_sim, df['Exud_ads'].loc[:].values[:len(x_sim)], '.', label = 'Exud_sorbed')
plt.plot(x_sim,  df['Exud_liq'].loc[:].values[:len(x_sim)], '.', label = 'Exud_liq')
plt.plot(x_direct,y_direct*dt*1e6, label = 'direct calculation')
plt.legend()
plt.title('no decay, sorption')
plt.xlabel('Time (d)')
plt.xlim(0,0.5)
plt.ylim(0,max( df['Exud_tot'].loc[:].values)*1.1)
plt.ylabel('Exudate content in soil domain (mol)', x = 0.05)
plt.show()
