"""
    Dynamic:

    plots results of local sensitivity analysis
    (see run_SA)
    
    Daniel Leitner, 2025       
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import run_SA as sa


def start_index(ind, ranges):
    s = 0
    for i in range(0, ind):
        if len(ranges) > 1:
            s += len(ranges[i])
    return s


def plot_local_SA(file_name, env_str, analysis_time, not_xlog, show = True):

    exp_name = file_name + env_str + "_time_{:g}".format(analysis_time)

    path = file_name+ env_str+"/"
    
    names, ranges = sa.read_ranges(path)
    print(len(names), "variables:")
    print(names)
    for i, r in enumerate(ranges):
        print(names[i], r)
    print()
    
    
    alldata = np.load(path + file_name + env_str +"_{:03}".format(1)+ ".npz")
    times = alldata["times"]
    times_lr = alldata["times_lr"]
    print("Simulation time from", min(times), "to ", max(times), "days")
    
    ind_ = np.argwhere(times > analysis_time)
    if len(ind_) > 0:
        ind = ind_[0][0]
        ind += 1
        ind10 = ind // 10  
    else:
        ind = -1
        ind10 = -1
    
    print("Index", ind)
    print("Plotting for day", times[ind])
    dt_ = np.diff(times[:ind])
    
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
    final_nd, final, final_trans, final_collar, final_names = [], [], [], [], []
    
    n = len(names)
    m = int(np.ceil(np.sqrt(n)))
    fig, ax = plt.subplots(int(np.ceil(n/m)), m, figsize = (16, 16))
    
    ac = 0
    for lind in range(0, len(names)):
    
        file_start_ind = 2 + start_index(lind, ranges)  # 1 is initial simulation
        sa_len = len(ranges[lind])
        # print("\nproducing sub-plot", names[lind], "from", file_start_ind, "to", file_start_ind + sa_len - 1)  # to make sure it works
    
        if sa_len > 1:
            
            trans = np.zeros((sa_len,))
            vol = np.zeros((sa_len,))
            krs = np.zeros((sa_len,))
            carbon = np.zeros((sa_len,))
            collar = np.zeros((sa_len,))
            
            for k in range(0, sa_len):
                try:
                    file_ = path + file_name + env_str +"_{:03}".format(file_start_ind + k)
                    alldata = np.load(file_+ ".npz")
        
                    trans_ = alldata["act_trans"]
                    trans[k] = -np.sum(np.multiply(trans_[:ind - 1], dt_))  # cumulative
                    
                    vol_ = alldata["vol"]
                    vol_ = vol_[:, ind10]
                    vol[k] = np.sum(vol_)  # sum over sub-types
        
                    krs_ = alldata["krs"]
                    krs[k] = krs_[ind10]
        
                    carbon_ = alldata["carbon"]
                    carbon[k] = carbon_[ind10][0] # only proportional to volume (since no underlying anatomy model)
                    
                    collar_ = alldata["collar_pot"]
                    collar[k] = np.mean(collar_) #[ind]

                except:
                    trans[k] = np.nan
                    krs[k] = np.nan
                    carbon[k] = np.nan
                    collar[k]= np.nan
                    print("*** skipping file", file_)
                       
            trans_nd = trans / trans[sa_len // 2]  # nondimensionalize
            vol_nd = vol / vol[sa_len // 2]  # nondimensionalize
            krs_nd = krs / krs[sa_len // 2]  # nondimensionalize
            carbon_nd = carbon / carbon[sa_len//2] # nondimensionalize
            collar_nd = collar / collar[sa_len//2] # nondimensionalize
        
            ax.flat[ac].plot(ranges[lind], trans_nd, '*-', label = "Water")  #  / ranges[lind][sa_len // 2]
            ax.flat[ac].plot(ranges[lind], collar_nd, '-.', label = "Leaf potential")        
            print(names[lind], collar)
            # ax.flat[ac].plot(ranges[lind], vol_nd, '-.', label = "volume")
            ax.flat[ac].plot(ranges[lind], carbon_nd, '-.', label = "Carbon")
            #ax.flat[ac].plot(ranges[lind], krs_nd, ':', label = "$K_{rs}$")
            
            x = ranges[lind][sa_len // 2]
            ax.flat[ac].plot([x], [1.], 'r*')# center red dot         
            ax.flat[ac].set_title(names[lind])
            ax.flat[ac].set_ylim(0.8, 1.2)
            # ax.flat[ac].set_yscale('log', base = 2)
            if not lind in not_xlog:
                ax.flat[ac].set_xscale('log', base = 2)
            if ac==0:
                ax.flat[ac].legend()
    
            final_names.append(names[lind]+"+")
            final_names.append(names[lind]+"-")
            
            final_nd.append(trans_nd[-1] / carbon_nd[-1])
            final_nd.append(trans_nd[0] / carbon_nd[0])
            
            final.append(trans[-1] / carbon[-1])
            final.append(trans[0] / carbon[0])
            
            final_trans.append(trans[-1])
            final_trans.append(trans[0])
            
            final_collar.append(collar_nd[-1])
            final_collar.append(collar_nd[0])
                        
            ac += 1
    
    plt.tight_layout(pad = 4.)  
    fig.suptitle(exp_name) 
    plt.savefig(exp_name+".png", dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    
    df = pd.DataFrame([final_nd, final, final_trans, final_collar], columns = final_names, 
                      index = ["trans/carbon","trans/carbon", "trans", "leaf pot"])
    df = df.round(3)    
    print()
    print(df.to_string(index=True, justify='center'))    
    
    
    I0 = np.argsort(final_nd)[::-1]
    I1 = np.argsort(final)[::-1]
    I2 = np.argsort(final_trans)[::-1]
    I3 = np.argsort(final_collar)[::-1]

    print("\nSorted")
    final_names = np.array(final_names)
    df = pd.DataFrame([final_names[I0], final_names[I1], final_names[I2], final_names[I3]], index = ["trans/carbon","trans/carbon", "trans", "leaf pot"])
    print(df.to_string(index=True, justify='left'))      
    
    
    print("Collar potential ", sa_len, collar)
    
if __name__ == "__main__":

    
    # plot_local_SA("local_soybean", env_str = "_0", analysis_time = 40., not_xlog = [5,7])
    # plot_local_SA("local_soybean", env_str = "_0", analysis_time = 87., not_xlog = [5,7])
    #plot_local_SA("local_soybean", env_str = "_1", analysis_time = 40, not_xlog = [5,7])        
    plot_local_SA("local_soybean", env_str = "_1", analysis_time = 87, not_xlog = [5,7])        

    
