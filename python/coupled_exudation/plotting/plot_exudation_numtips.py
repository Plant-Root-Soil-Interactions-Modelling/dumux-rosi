"""
exudation  plot
"""
import sys; sys.path.append("../../../CPlantBox/");
sys.path.append("../../../CPlantBox/src/python_modules")
import numpy as np
import matplotlib.pyplot as plt


fname = "carbon_maize_exudate_2019"
path = "../results/"

""" exudation plot """

# load data
print(fname) 
data = np.load(path + fname + ".npy", allow_pickle=True) 

fig, ax = plt.subplots(1,2, figsize = (18, 8))

print(np.shape(data))
t = data[0]
c = data[1]
#nt = data[2]

lns1 = ax[0].plot(t, c, 'r', label = "Carbon exudation [g/day]")
ax[0].set_ylabel("Carbon exudation [g/day]")
dt = np.diff(t)
#print(c[:-1])
#sys.exit()
cuc = np.cumsum(np.multiply(c[:-1], dt))
ax2 = ax[0].twinx()
ax2.set_ylabel("Cumulative carbon exudation [g]")
lns2 = ax2.plot(t[1:],  cuc, 'r--', label = "Cumulative carbon exudation [g]")  
lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax[0].legend(lns, labs, loc=0)
ax[0].set_xlabel("Time [day]")

lns1 = ax[1].plot(t, c, 'r', label = "Carbon exudation [g/day]")
ax[1].set_ylabel("Carbon exudation [g/day]")
#ax2 = ax[1].twinx()
#ax2.set_ylabel("Number of root tips [-]")
#lns2 = ax2.plot(t,  nt, 'r--', label = "Number of root tips [-]")  
#lns = lns1+lns2
#labs = [l.get_label() for l in lns]
#ax[1].legend(lns, labs, loc=0)
ax[1].set_xlabel("Time [day]")

plt.tight_layout(pad = 5)
plt.show()
