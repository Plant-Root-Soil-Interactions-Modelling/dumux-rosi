import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import *
from xylem_flux import sinusoidal2
from scipy import interpolate

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title


def get_transpiration(year, sim_time, area, f, Kc):
    """ calculates transpiration with Beer's law from start_date from a pickle file"""

    sim_time += 1  # root system starts to grow one day earlier (+max dat)

    """ 0. load ET, ET0 -> ETc """
    df2 = pd.read_csv("data_magda/ET0.csv")  # evapotranspiration  in cm / day 
    yd2 = df2["ET0_"+str(year)].loc[0:sim_time].values # cm / d
    et0 = yd2
    etc = et0 * Kc[0:sim_time+1]

    """ 1. ETc -> Tpot, Evap """
    k = 0.45
    t_ = np.linspace(0, len(etc) - 1, len(etc))
    tpot = np.multiply(etc, [(1. - np.exp(-k * f(t_[i]))) for i in range(0, len(etc))])
    evap = etc - tpot
    trans = lambda t, dt:-tpot[int((t + dt / 2))] * area * sinusoidal2(t, dt)

    #t2_ = np.linspace(0, sim_time + 0.4, len(etc) * 24)  # more points more fun...
    #tpot2 = [-trans(t, 0.) / area for t in t2_]
    #plt.bar(t_, etc, label = "Crop evapotranspiration")
    #plt.bar(t_, tpot, align = 'edge', label = "Potential transpiration")
    #plt.plot(t2_, tpot2, 'r', label = "Potential transpiration")
    #plt.bar(t_, evap,  label = "Evaporation")
    #plt.xlabel("time")
    #plt.ylabel("cm / day")
    #plt.legend()
    #plt.show()

    return trans

def net_infiltration(year, soil_type, genotype, sim_time, Kc):
    """ calculates net infiltration with Beer's law from start_date from csv file"""

    sim_time += 1

    """ 0. load infiltration (precipitation + irrigation) """
    df = pd.read_csv("data_magda/Inf.csv")  # precipitation data
    yd1  = df["Inf_"+str(year)].loc[0: sim_time].values*0.1 #nmm/d - cm/d
    t_ = np.linspace(0, sim_time, yd1.shape[0] * 24)  # relative time in hours
    precip = np.array([ yd1[int(t)] * sinusoidal2(t, 0.) for t in t_ ])
    
    """ 1. load ET, ET0 -> ETc """
    df2 = pd.read_csv("data_magda/ET0.csv")  # evapotranspiration  in cm / day 
    yd2 = df2["ET0_"+str(year)].loc[0:sim_time].values # cm / d
    et0 = np.array([ yd2[int(t)] * sinusoidal2(t, 0.) for t in t_ ])

    """2. interpolate Kc """
    f_Kc = interpolate.interp1d(np.linspace(0,sim_time, sim_time+1), Kc[0:sim_time+1])
    Kc_hours = f_Kc(t_)
    etc = et0 * Kc_hours

    """3. load LAI """
    if soil_type == 'loam':
        st = 'L'
    else:
        st = 'S'
        
    df3 = pd.read_csv("data_magda/LAI_"+str(year)+".csv")  # LAI
    LAI =  df3[st+"_"+genotype].loc[0:sim_time].values
    f = interpolate.interp1d(np.linspace(0,sim_time, sim_time+1), LAI)

    """ 4. ETc -> Tpot, Evap """
    k = 0.45
    tpot = np.multiply(etc, [1. - np.exp(-k *f(t)) for t in t_])
    evap = etc - tpot
    evap = -evap
    net_inf = precip + evap
 
    y_ = []
    x_ = []
    for i in range(0, t_.shape[0] - 1):
        x_.extend([t_[i], t_[i + 1]])  # hour -> day
        y_.extend([net_inf[i], net_inf[i]])


    #fig, ax = plt.subplots(3, figsize = (18, 10))
    #ax[0].bar(t_, precip * 10, width = 1. / 24., label = "precipiation")
    #ax[0].legend()
    #ax[0].set_ylabel("[mm/day]")
    # # ax[0].set_xlabel("time")
    #ax[1].plot([0., t_[-1]], [0., 0.], 'k')
    # # ax[1].plot(t_, et0, 'k', label = "evapotranspiration")
    #ax[1].plot(t_, etc * 10, 'k:', label = "crop evapotranspiration")
    #ax[1].plot(t_, tpot*10 , 'r', label = "potential transpiration")
    #ax[1].plot(t_, -evap * 10, 'g', label = "evaporation")
    #ax[1].legend()
    # # ax[1].set_xlabel("time")
    #ax[1].set_ylabel("[mm/day]")
    #ax[2].bar(t_, net_inf * 10, width = 1. / 24., label = 'net infiltration')
    #ax[2].set_xlabel("time [day]")
    #ax[2].set_ylabel("[mm/day]")
    #ax[2].legend()
    #plt.show()

    return np.array(x_), np.array(y_), f


def decay(soil_sol_fluxes, dt, decay = 0.):
    """ removes exudates at the given decay rate  [1/day] """

    for k in soil_sol_fluxes.keys():
        #print('before', soil_sol_fluxes[k], dt, decay) 
        soil_sol_fluxes[k] = (1-decay*dt)*soil_sol_fluxes[k]
        #print('after', soil_sol_fluxes[k]) 

    return soil_sol_fluxes
