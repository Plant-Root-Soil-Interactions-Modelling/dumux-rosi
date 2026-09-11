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



def get_transpiration(simtime, area, Kc, soil_type):
    """ calculates transpiration with Beer's law"""
    year = 2019
    hours_per_day = 24
    k = 0.45

    lai_data = pd.read_csv("../inputDataExudate/data/LAI.csv") 
    et0_hourly = pd.read_csv("../inputDataExudate/data/ET0.csv")[f"ET0_{year}"].values[:simtime * hours_per_day]  # cm/d
    etc_hourly = et0_hourly * Kc[:len(et0_hourly)]
    
    """3. load LAI """
    if soil_type == 'loam':
        st = 'L'
    else:
        st = 'S'
    
    lai_daily = lai_data[st+"_WT"].values[:simtime]
    lai_hourly = np.repeat(lai_daily, hours_per_day)
    
    tpot_hourly = etc_hourly * (1 - np.exp(-k * lai_hourly))
    evap_hourly = etc_hourly - tpot_hourly
    
    # Convert hourly values back to daily mean rates
    tpot_daily = tpot_hourly.reshape(-1, hours_per_day).mean(axis=1)
    evap_daily = evap_hourly.reshape(-1, hours_per_day).mean(axis=1)

    #transpiration function 
    trans = lambda t, dt:-tpot_daily[int((t + dt / 2))] * area * sinusoidal2(t, dt)
   
    return trans

def net_infiltration(soil_type, simtime, Kc):
    """ calculates net infiltration with Beer's law"""
    year = 2019
    hours_per_day = 24
    k = 0.45

    days = np.arange(simtime)
    hours = np.arange(simtime * hours_per_day) / hours_per_day

    precip_daily = pd.read_csv("../inputDataExudate/data/Inf.csv")["Inf"].values[:simtime] * 0.1  # mm/d --> cm/d
    t_ = np.linspace(0, simtime-1, precip_daily.shape[0] * 24)  # relative time in hours
    precip_hourly = np.array([ precip_daily[int(t)] * sinusoidal2(t, 0.) for t in t_ ])
    lai_data = pd.read_csv("../inputDataExudate/data/LAI.csv") 
    et0_hourly = pd.read_csv("../inputDataExudate/data/ET0.csv")[f"ET0_{year}"].values[:simtime * hours_per_day]  # cm/d
    etc_hourly = et0_hourly * Kc[:len(et0_hourly)]
    
    if soil_type == 'loam':
        st = 'L'
    else:
        st = 'S'
    
    lai_daily = lai_data[st+"_WT"].values[:simtime]
    lai_hourly = np.repeat(lai_daily, hours_per_day)

    tpot_hourly = etc_hourly * (1 - np.exp(-k * lai_hourly))
    evap_hourly = etc_hourly - tpot_hourly
    net_inf_hourly = precip_hourly - evap_hourly
    
    # Convert hourly values back to daily mean rates
    tpot_daily = tpot_hourly.reshape(-1, hours_per_day).mean(axis=1)
    evap_daily = evap_hourly.reshape(-1, hours_per_day).mean(axis=1)
    net_inf_daily = net_inf_hourly.reshape(-1, hours_per_day).mean(axis=1)
    
    return np.array(hours), np.array(net_inf_hourly)
