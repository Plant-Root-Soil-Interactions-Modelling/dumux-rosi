import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import *
import sys
sys.path.append("../../../../CPlantBox/");
sys.path.append("../../../../CPlantBox/src/")
sys.path.append("../../../../CPlantBox/src/functional/");
sys.path.append("../inputDataExudate/data/");
sys.path.append("data/");
from xylem_flux import sinusoidal2
from scipy import interpolate
import os



def get_transpiration(sim_time, area, Kc, soil_type):
    """ calculates transpiration with Beer's law from start_date from a pickle file"""

    year = 2019
    sim_time += 1  # root system starts to grow one day earlier (+max dat)

    """ 0. load ET, ET0 -> ETc """
    df2 = pd.read_csv("../inputDataExudate/data//ET0.csv")  # evapotranspiration  in cm / day 
    yd2 = df2["ET0_"+str(year)].loc[0:sim_time].values # cm / d
    et0 = yd2
    etc = et0 * Kc[0:sim_time+1]
    
    """3. load LAI """
    if soil_type == 'loam':
        st = 'L'
    else:
        st = 'S'
    
    df3 = pd.read_csv("../inputDataExudate/data/LAI_"+str(year)+".csv")  # LAI
    LAI =  df3[st+"_WT"].loc[0:sim_time].values
    f = interpolate.interp1d(np.linspace(0,sim_time, int(sim_time+1)), LAI)

    """ 1. ETc -> Tpot, Evap """
    k = 0.45
    t_ = np.linspace(0, len(etc) - 1, len(etc))
    tpot = np.multiply(etc*10, [(1. - np.exp(-k * f(t_[i]))) for i in range(0, len(etc))])
    evap = (etc*10 - tpot)
    trans = lambda t, dt:-tpot[int((t + dt / 2))] * area * sinusoidal2(t, dt)

    return trans

def net_infiltration(soil_type, sim_time, Kc):
    """ calculates net infiltration with Beer's law from start_date from csv file"""

    year = 2019
    genotype = 'WT'
    sim_time += 1

    """ 0. load infiltration (precipitation + irrigation) """
    print(os.getcwd()) 
    df = pd.read_csv("../inputDataExudate/data/Inf.csv")  # precipitation data
    yd1  = df["Inf_"+str(year)].loc[0: sim_time].values*0.1 #nmm/d - cm/d
    t_ = np.linspace(0, sim_time, yd1.shape[0] * 24)  # relative time in hours
    precip = np.array([ yd1[int(t)] * sinusoidal2(t, 0.) for t in t_ ])
    
    """ 1. load ET, ET0 -> ETc """
    df2 = pd.read_csv("../inputDataExudate/data//ET0.csv")  # evapotranspiration  in cm / day 
    yd2 = df2["ET0_"+str(year)].loc[0:sim_time].values # cm / d
    et0 = np.array([ yd2[int(t)] * sinusoidal2(t, 0.) for t in t_ ])

    """2. interpolate Kc """
    f_Kc = interpolate.interp1d(np.linspace(0,sim_time, int(sim_time+1)), Kc[0:int(sim_time+1)])
    Kc_hours = f_Kc(t_)
    etc = et0 * Kc_hours

    """3. load LAI """
    if soil_type == 'loam':
        st = 'L'
    else:
        st = 'S'
        
    df3 = pd.read_csv("../inputDataExudate/data/LAI_"+str(year)+".csv")  # LAI
    LAI =  df3[st+"_"+genotype].loc[0:sim_time].values
    #print('LAI', np.max(LAI),np.min(LAI))
    #sys.exit()
    f = interpolate.interp1d(np.linspace(0,sim_time, int(sim_time+1)), LAI)

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

    return np.array(x_), np.array(y_)
