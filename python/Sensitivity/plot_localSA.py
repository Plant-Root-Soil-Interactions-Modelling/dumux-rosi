"""
plots results of local sensitivity analysis
"""
import sys; sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

import numpy as np
import matplotlib.pyplot as plt

import run_SA as sa

def start_index(ind, ranges):
    s = 0
    for i in range(0,ind):
        if len(ranges)>1:        
            s += len(ranges[i])
    return s
    

""" def SA """
file_name = "local_soybean"
#file_name = "local_maize"
path = "results/"
not_xlog = ["theta1", "src"]

analysis_time = 100 # days

names, ranges = sa.read_ranges(path+file_name)
trans_ = np.load(path + "transpiration_" + file_name + "1" + ".npy")
times = trans_[0,:]
dt_ = np.diff(times)
dt_ = np.array([0., *dt_])
print("Simulation time from", min(times), "to ", max(times), "days")

ind_ = np.argwhere(times>analysis_time)
if len(ind_)>0:
    ind = ind_[0][0]
    print("plots for day", times[ind])
    ind += 1 
    ind10 = ind//10 # TODO check 
    print(ind)
else: 
    ind = -1
    ind10 = -1

""" font sizes """
SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title

""" make plots """
fig, ax = plt.subplots(3, 3, figsize = (16, 16))

ac = 0
for lind in range(0,len(names)):
    
    file_start_ind = 2 + start_index(lind, ranges) # 1 is initial simulation
    sa_len = len(ranges[lind])
    print("\nproducing sub-plot", names[lind], "from", file_start_ind, "to", file_start_ind+sa_len)    
    
    if sa_len>1:
        trans = np.zeros((sa_len,))
        vol = np.zeros((sa_len,))
        krs = np.zeros((sa_len,))
        for k in range(0, sa_len):
            try:
                trans_ = np.load(path + "transpiration_" + file_name + str(file_start_ind + k) + ".npy")
                trans[k] = np.sum(np.multiply(trans_[1,:ind], dt_[:ind]))
                print(trans[k])
            except:
                trans[k] = np.nan
                print("skipping file", file_name + str(file_start_ind + k))
            try:
                vol_ = np.load(path + "vol_" + file_name + str(file_start_ind + k) + ".npy")
                vol_ = vol_[:, ind10]
                vol[k] = np.sum(vol_) # sum over sub-types                
            except:
                vol[k] = np.nan
            try:
                krs_ = np.load(path + "krs_" + file_name + str(file_start_ind + k) + ".npy")
                krs[k] = krs_[ind10]
            except:
                krs[k] = np.nan
        trans = trans / trans[sa_len // 2]  # nondimensionalize
        vol = vol / vol[sa_len // 2]  # nondimensionalize
        krs = krs / krs[sa_len // 2]  # nondimensionalize
        x = ranges[lind][sa_len // 2]
        ax.flat[ac].plot(ranges[lind], trans, label = "uptake")
        #ax.flat[ac].plot(ranges[lind], vol, '-.', label = "volume")
        ax.flat[ac].plot(ranges[lind], krs, ':', label = "krs")        
        ax.flat[ac].plot([x], [1.], 'r*')
        ax.flat[ac].legend()
        ax.flat[ac].set_title(names[lind])
        #ax.flat[ac].set_ylim(0.5, 2)
        #ax.flat[ac].set_yscale('log', base = 2)
        if not names[lind] in not_xlog:
            ax.flat[ac].set_xscale('log', base = 2)
        ac += 1

plt.tight_layout(pad = 4.)
plt.show()
