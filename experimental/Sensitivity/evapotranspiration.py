import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import *
from functional.xylem_flux import sinusoidal2
from lxml.html.builder import DD

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

# Soybean
# sigmoidal
# [ 7.44990146 55.4283463   0.15200146  0.34325857]
# neg exp^2
# [8.83381831e+00 8.29555707e+01 1.15834259e-03 1.75018634e-01]
#
# Maize
# sigmoidal
# [ 3.32719912 28.89532431  0.22100482  0.08645908]
# neg exp^2
# [ 5.14321003e+00  6.65910800e+01  4.09254414e-04 -1.08617435e+00]


def sigmoid(x, L , x0, k, b):
    y = L / (1 + np.exp(-k * (x - x0))) + b
    return max(y, 0)


def exp2(x, L , x0, k, b):
    y = L * np.exp(-k * (x - x0) * (x - x0)) + b
    return max(y, 0)


def lai_soybean(x):
    return sigmoid(x, 7.44990146, 55.4283463, 0.15200146, 0.34325857)  # see test_lai.py


def lai_soybean2(x):
    return exp2(x, 8.83381831e+00, 8.29555707e+01, 1.15834259e-03, 1.75018634e-01)  # see test_lai.py


def lai_maize(x):
    return sigmoid(x, 3.32719912, 28.89532431, 0.22100482, 0.08645908)  # see test_lai.py


def lai_maize2(x):
    return exp2(x, 5.14321003e+00, 6.65910800e+01, 4.09254414e-04, -1.08617435e+00)  # see test_lai.py


def lai_noroots(x):
    return 0.


def get_transpiration_pickle(filename, range_, area):
    """ takes the transpiration of filename in time @range_ from a pickle file"""
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    y_ = data["Tpot"]
    y = y_.loc[range_[0]: range_[1]].values / 10. * 24.  # mm -> cm, /hour -> /day
    trans = lambda t, dt:-y[int((t + dt / 2) * 24)] * area  # day -> hour
    return trans


def net_infiltration_table_pickle(filename, range_):
    """ takes the net infiltration of filename in time @range_ from a pickle file"""
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    yd = data["Net_infilteration"]
    y = yd.loc[range_[0]: range_[1]].values / 10. * 24.  # mm -> cm, /hour -> /day
    y_ = []
    x_ = []
    for i in range(0, y.shape[0]):
        x_.extend([float(i) / 24, float(i + 1) / 24])  # hour -> day
        y_.extend([y[i], y[i]])
    return x_, y_


def get_transpiration_beers_pickle(filename, start_date, sim_time, area, lai_f, Kc):
    """ calculates transpiration with Beer's law from start_date from a pickle file"""
    start_date = datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
    end_date = start_date + timedelta(sim_time + 1, 3600)  # root system starts to grow one day earlier (+max dat)
    range_ = [str(start_date), str(end_date)]

    with open(filename, 'rb') as f:
        data = pickle.load(f)

    yd = data["ET0"]
    et0 = yd.loc[range_[0]: range_[1]].values / 10. * 24.  # mm -> cm, /hour -> /day

    """ 1. ET0 -> ETc """
    etc = et0 * Kc  # ETc = crop evapotranspiration

    """ 2. ETc -> Tpot, Evap """
    k = 0.6
    t_ = np.linspace(0, sim_time, len(etc))
    tpot = np.multiply(etc, [(1. - np.exp(-k * lai_f(t_[i]))) for i in range(0, len(etc))])
    evap = etc - tpot

    # plt.plot(t_, etc, label = "Crop evapotranspiration")
    # plt.plot(t_, tpot, "r", label = "Potential transpiration")
    # plt.plot(t_, evap, "g", label = "Evaporation")
    # y_ = data["Tpot"]
    # # y = y_.loc[range_[0]: range_[1]].values / 10. * 24.  # mm -> cm, /hour -> /day
    # # plt.plot(t_, y, "b:", label = "Potential transpiration Ullah")
    # plt.xlabel("time")
    # plt.ylabel("cm / day")
    # plt.legend()
    # plt.show()
    # print("hourly potential transpiration len", len(tpot))
    trans = lambda t, dt:-tpot[int((t + dt / 2) * 24)] * area  # day -> hour

    return trans


def get_transpiration_beers_csv(start_date, sim_time, area, lai_f, Kc):
    """ calculates transpiration with Beer's law from start_date from a pickle file"""

    start_date = datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
    end_date = start_date + timedelta(sim_time + 1, 3600)  # root system starts to grow one day earlier (+max dat)
    range_ = [str(start_date), str(end_date)]

    """ 1. load ET, ET0 -> ETc """
    df2 = pd.read_csv("data/K0_ET0_2020-22.csv")  # evapotranspiration kg / m2 ; kg = dm3 = 0.1 cm m2 # 0.1 cm / day (?)
    df2['Datetime'] = pd.to_datetime(df2['Date'], format = '%Y-%m-%d')
    df2 = df2.set_index(pd.DatetimeIndex(df2['Datetime']))
    df2 = df2.drop(['Date'], axis = 1)
    yd2 = df2["ET (kg/m2)"].loc[range_[0]: range_[1]].values * 0.1  # dm->cm
    et0 = yd2
    etc = et0 * Kc  # ETc = crop evapotranspiration

    """ 2. ETc -> Tpot, Evap """
    k = 0.6
    t_ = np.linspace(0, sim_time, len(etc))
    tpot = np.multiply(etc, [(1. - np.exp(-k * lai_f(t_[i]))) for i in range(0, len(etc))])
    evap = etc - tpot

    # plt.plot(t_, etc, label = "Crop evapotranspiration")
    # plt.plot(t_, tpot, "r", label = "Potential transpiration")
    # plt.plot(t_, evap, "g", label = "Evaporation")
    # plt.xlabel("time")
    # plt.ylabel("cm / day")
    # plt.legend()
    # plt.show()

    trans = lambda t, dt:-tpot[int((t + dt / 2))] * area
    return trans


def get_transpiration_beers_csvS(start_date, sim_time, area, lai_f, Kc):
    """ calculates transpiration with Beer's law from start_date from a pickle file"""

    sim_time += 1  # root system starts to grow one day earlier (+max dat)
    start_date = datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
    end_date = start_date + timedelta(sim_time + 1, 0)
    range_ = [str(start_date), str(end_date)]

    """ 1. load ET, ET0 -> ETc """
    df2 = pd.read_csv("data/K0_ET0_2020-22.csv")  # evapotranspiration kg / m2 ; kg = dm3 = 0.1 cm m2 # 0.1 cm / day (?)
    df2['Datetime'] = pd.to_datetime(df2['Date'], format = '%Y-%m-%d')
    df2 = df2.set_index(pd.DatetimeIndex(df2['Datetime']))
    df2 = df2.drop(['Date'], axis = 1)
    yd2 = df2["ET (kg/m2)"].loc[range_[0]: range_[1]].values * 0.1  # dm->cm
    et0 = yd2
    etc = et0 * Kc  # ETc = crop evapotranspiration

    """ 2. ETc -> Tpot, Evap """
    k = 0.6
    print(sim_time, len(etc))
    t_ = np.linspace(0, len(etc) - 1, len(etc))
    tpot = np.multiply(etc, [(1. - np.exp(-k * lai_f(t_[i]))) for i in range(0, len(etc))])
    evap = etc - tpot

    trans = lambda t, dt:-tpot[int((t + dt / 2))] * area * sinusoidal2(t, dt)

    # t2_ = np.linspace(0, sim_time + 0.4, len(etc) * 24)  # more points more fun...
    # tpot2 = [-trans(t, 0.) / area for t in t2_]
    # # plt.bar(t_, etc, label = "Crop evapotranspiration")
    # plt.bar(t_, tpot, align = 'edge', label = "Potential transpiration")
    # plt.plot(t2_, tpot2, 'r', label = "Potential transpiration")
    # # plt.bar(t_, evap,  label = "Evaporation")
    # plt.xlabel("time")
    # plt.ylabel("cm / day")
    # plt.legend()
    # plt.show()

    return trans


def net_infiltration_table_beers_pickle(filename, start_date, sim_time, lai_f, Kc):
    """ calculates net infiltration with Beer's law from start_date from a pickle file"""

    start_date = datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
    end_date = start_date + timedelta(sim_time + 1, 3600)  # root system starts to grow one day earlier (+max dat)
    range_ = [str(start_date), str(end_date)]

    with open(filename, 'rb') as f:
        data = pickle.load(f)

    """ 0. precipitation, evaporation """
    yd = data["Net_infilteration"]
    y = yd.loc[range_[0]: range_[1]].values / 10. * 24.  # mm -> cm, /hour -> /day
    evapd = data["Epot"]
    evap = -evapd.loc[range_[0]: range_[1]].values / 10. * 24.  # mm -> cm, /hour -> /day
    precip = y - evap

    """ 1. ET0 -> ETc """
    yd = data["ET0"]
    et0 = yd.loc[range_[0]: range_[1]].values / 10. * 24.  # mm -> cm, /hour -> /day
    etc = et0 * Kc
    """ 2. ETc -> Tpot, Evap """
    k = 0.6
    t_ = np.linspace(0, sim_time, len(etc))
    tpot = np.multiply(etc, [(1. - np.exp(-k * lai_f(t_[i]))) for i in range(0, len(etc))])
    evap = etc - tpot
    evap = -evap
    net_inf = precip + evap

    # fig, ax = plt.subplots(3)
    # ax[0].bar(t_, precip, width = 0.8 / 24., label = "precipiation")
    # ax[0].legend()
    # ax[0].set_ylabel("cm / day")
    # ax[0].set_xlabel("time")
    # ax[1].plot([0., t_[-1]], [0., 0.], 'k')
    # ax[1].plot(t_, tpot, 'r', label = "potential transpiration")
    # ax[1].plot(t_, etc, 'k:', label = "crop evapotranspiration")
    # ax[1].plot(t_, -evap, 'g', label = "evaporation")
    # ax[1].legend()
    # ax[1].set_xlabel("time")
    # ax[1].set_ylabel("cm / day")
    # ax[2].bar(t_, net_inf, width = 0.8 / 24., label = 'net infiltration')
    # ax[2].set_xlabel("time")
    # ax[2].set_ylabel("cm / day")
    # ax[2].legend()
    # plt.show()

    y_ = []
    x_ = []
    for i in range(0, precip.shape[0]):
        x_.extend([float(i) / 24, float(i + 1) / 24])  # hour -> day
        y_.extend([net_inf[i], net_inf[i]])

    return x_, y_


def net_infiltration_table_beers_csv(start_date, sim_time, lai_f, Kc):
    """ calculates net infiltration with Beer's law from start_date from INARI csv file"""

    sim_time += 1
    start_date = datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
    end_date = start_date + timedelta(sim_time, 3600)  # root system starts to grow one day earlier (+max dat)
    range_ = [str(start_date), str(end_date)]

    """ 0. load precipitation """
    df = pd.read_csv("data/WeatherSummary_K0_20-22.csv")  # precipitation data
    df['Datetime'] = pd.to_datetime(df['date'], format = '%Y-%m-%d')
    df = df.set_index(pd.DatetimeIndex(df['Datetime']))
    df = df.drop(['date'], axis = 1)
    yd = df["prcp"].loc[range_[0]: range_[1]].values * 2.54  # inches -> cm
    precip = yd
    t_ = np.linspace(0, len(precip), len(precip))  # relative time

    """ 1. load ET, ET0 -> ETc """
    df2 = pd.read_csv("data/K0_ET0_2020-22.csv")  # evapotranspiration kg / m2 ; kg = dm3 = 0.1 cm m2 # 0.1 cm / day (?)
    df2['Datetime'] = pd.to_datetime(df2['Date'], format = '%Y-%m-%d')
    df2 = df2.set_index(pd.DatetimeIndex(df2['Datetime']))
    df2 = df2.drop(['Date'], axis = 1)
    yd2 = df2["ET (kg/m2)"].loc[range_[0]: range_[1]].values * 0.1  # dm->cm
    et0 = yd2
    etc = et0 * Kc
    """ 2. ETc -> Tpot, Evap """
    k = 0.6
    t_ = np.linspace(0, sim_time, len(etc))
    tpot = np.multiply(etc, [(1. - np.exp(-k * lai_f(t_[i]))) for i in range(0, len(etc))])
    evap = etc - tpot
    evap = -evap
    net_inf = precip + evap
    print()
    print("total precip", np.sum(precip) * 10, "mm")
    print("total evaporation", np.sum(evap) * 10, "mm")
    print("total net_inf", np.sum(net_inf) * 10, "mm")
    print()

    # print(net_inf.shape)
    # np.savetxt("net_inf.csv", net_inf * 10, delimiter = ",")
    # dd

    # fig, ax = plt.subplots(3)
    # ax[0].bar(t_, precip, width = 0.8, label = "precipiation")
    # ax[0].legend()
    # ax[0].set_ylabel("cm / day")
    # ax[0].set_xlabel("time")
    # ax[1].plot([0., t_[-1]], [0., 0.], 'k')
    # # ax[1].plot(t_, et0, 'k:', label = "evapotranspiration")
    # ax[1].plot(t_, etc, 'k:', label = "crop evapotranspiration")
    # ax[1].plot(t_, tpot, 'r', label = "potential transpiration")
    # ax[1].plot(t_, -evap, 'g', label = "evaporation")
    # ax[1].legend()
    # ax[1].set_xlabel("time")
    # ax[1].set_ylabel("cm / day")
    # ax[2].bar(t_, net_inf, width = 0.8, label = 'net infiltration')
    # ax[2].set_xlabel("time")
    # ax[2].set_ylabel("cm / day")
    # ax[2].legend()
    # print("net change", np.sum(net_inf) * 10, "mm")
    # plt.show()

    y_ = []
    x_ = []
    for i in range(0, precip.shape[0]):
        x_.extend([float(i), float(i + 1)])
        y_.extend([net_inf[i], net_inf[i]])

    return np.array(x_), np.array(y_)


def net_infiltration_table_beers_csvS(start_date, sim_time, lai_f, Kc):
    """ calculates net infiltration with Beer's law from start_date from INARI csv file"""

    sim_time += 1
    start_date = datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
    end_date = start_date + timedelta(sim_time, 3600)  # root system starts to grow one day earlier (+max dat)
    range_ = [str(start_date), str(end_date)]

    """ 0. load precipitation """
    df = pd.read_csv("data/WeatherSummary_K0_20-22.csv")  # precipitation data
    df['Datetime'] = pd.to_datetime(df['date'], format = '%Y-%m-%d')
    df = df.set_index(pd.DatetimeIndex(df['Datetime']))
    df = df.drop(['date'], axis = 1)
    yd = df["prcp"].loc[range_[0]: range_[1]].values * 2.54  # inches -> cm
    t_ = np.linspace(0, sim_time, yd.shape[0] * 24)  # relative time in hours
    precip = np.array([ yd[int(t)] for t in t_ ])  # TODO precip / day

    """ 1. load ET, ET0 -> ETc """
    df2 = pd.read_csv("data/K0_ET0_2020-22.csv")  # evapotranspiration kg / m2 ; kg = dm3 = 0.1 cm m2 # 0.1 cm / day (?)
    df2['Datetime'] = pd.to_datetime(df2['Date'], format = '%Y-%m-%d')
    df2 = df2.set_index(pd.DatetimeIndex(df2['Datetime']))
    df2 = df2.drop(['Date'], axis = 1)
    yd2 = df2["ET (kg/m2)"].loc[range_[0]: range_[1]].values * 0.1  # dm->cm
    et0 = np.array([ yd2[int(t)] * sinusoidal2(t, 0.) for t in t_ ])
    etc = et0 * Kc

    """ 2. ETc -> Tpot, Evap """
    k = 0.6
    tpot = np.multiply(etc, [(1. - np.exp(-k * lai_f(t))) for t in t_])
    evap = etc - tpot
    evap = -evap
    net_inf = precip + evap

    # fig, ax = plt.subplots(3, figsize = (18, 10))
    # ax[0].bar(t_, precip * 10, width = 1. / 24., label = "precipiation")
    # ax[0].legend()
    # ax[0].set_ylabel("[mm/day]")
    # # ax[0].set_xlabel("time")
    # ax[1].plot([0., t_[-1]], [0., 0.], 'k')
    # # ax[1].plot(t_, et0, 'k', label = "evapotranspiration")
    # ax[1].plot(t_, etc * 10, 'k:', label = "crop evapotranspiration")
    # ax[1].plot(t_, tpot * 10, 'r', label = "potential transpiration")
    # ax[1].plot(t_, -evap * 10, 'g', label = "evaporation")
    # ax[1].legend()
    # # ax[1].set_xlabel("time")
    # ax[1].set_ylabel("[mm/day]")
    # ax[2].bar(t_, net_inf * 10, width = 1. / 24., label = 'net infiltration')
    # ax[2].set_xlabel("time [day]")
    # ax[2].set_ylabel("[mm/day]")
    # ax[2].legend()
    # plt.show()
    # dd

    y_ = []
    x_ = []
    for i in range(0, t_.shape[0] - 1):
        x_.extend([t_[i], t_[i + 1]])  # hour -> day
        y_.extend([net_inf[i], net_inf[i]])

    return np.array(x_), np.array(y_)


def add_nitrificatin_source(s, soil_sol_fluxes, nit_flux = 0.):
    """ adds a consant nitrate source @param nit_flux due to nitrification [g/day] """
    z_ = np.linspace(-0.5, -29.5, 30)  # top 30 cm layers
    for z in z_:
        i = s.pick([0, 0, z])  # cell index
        if i in soil_sol_fluxes:
            soil_sol_fluxes[i] += nit_flux
        else:
            soil_sol_fluxes[i] = nit_flux
