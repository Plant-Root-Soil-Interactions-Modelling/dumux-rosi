import pickle
import matplotlib.pyplot as plt
import numpy as np


def sigmoid(x, L , x0, k, b):
    y = L / (1 + np.exp(-k * (x - x0))) + b
    return (y)


def lai_soybean(x):
    return sigmoid(x, 7.74397040e+00, 5.64296005e+01, 1.10652299e-01, 4.11185215e-02)  # see test_lai.py


def lai_maize(x):
    return sigmoid(x, 3.32217629e+00, 2.78937950e+01, 2.22487303e-01, 6.42233521e-03)  # see test_lai.py


def get_transpiration(filename, range_, area):
    """ takes the transpiration of filename in time @range_ """

    with open(filename, 'rb') as f:
        data = pickle.load(f)

    y_ = data["Tpot"]
    y = y_.loc[range_[0]: range_[1]].values / 10. * 24.  # mm -> cm, /hour -> /day
    # print(type(y))
    # print(len(y))
    # print(y)

    trans = lambda t, dt:-y[int((t + dt / 2) * 24)] * area  # day -> hour

    # plt.plot(y)
    # plt.xlabel("time")
    # plt.ylabel("cm / day")
    # plt.show()

    return trans


def net_infiltration_table(filename, range_):

    with open(filename, 'rb') as f:
        data = pickle.load(f)

    yd = data["Net_infilteration"]
    y = yd.loc[range_[0]: range_[1]].values / 10. * 24.  # mm -> cm, /hour -> /day
    # plt.plot(y)
    # plt.xlabel("time")
    # plt.ylabel("cm / day")
    # plt.show()

    y_ = []
    x_ = []
    for i in range(0, y.shape[0]):
        x_.extend([float(i) / 24, float(i + 1) / 24])  # hour -> day
        y_.extend([y[i], y[i]])

    return x_, y_


def get_transpiraton_beers(filename, range_, time, lai_f, Kc):

    with open(filename, 'rb') as f:
        data = pickle.load(f)

    yd = data["ET0"]
    et0 = yd.loc[range_[0]: range_[1]].values / 10. * 24.  # mm -> cm, /hour -> /day

    """ 1. ET0 -> ETc """
    etc = et0 * Kc

    """ 2. ETc -> Tpot, Evap """
    k = 0.6  # ?????????????????
    t_ = np.linspace(0, time, len(etc))
    tpot = np.multiply(etc, [(1. - np.exp(-k * lai_f(t_[i]))) for i in range(0, len(etc))])
    evap = etc - tpot

    # plt.plot(t_, etc, label = "Crop evapotranspiration")
    # plt.plot(t_, tpot, "r", label = "Potential transpiration")
    # plt.plot(t_, evap, "g", label = "Evaporation")
    # y_ = data["Tpot"]
    # y = y_.loc[range_[0]: range_[1]].values / 10. * 24.  # mm -> cm, /hour -> /day
    # plt.plot(t_, y, "b:", label = "Potential transpiration Ullah")
    # plt.xlabel("time")
    # plt.ylabel("cm / day")
    # plt.legend()
    # plt.show()

    trans = lambda t, dt:-tpot[int((t + dt / 2) * 24)] * area  # day -> hour

    return tpot


def net_infiltration_table_beers(filename, range_, time, lai_f, Kc):

    with open(filename, 'rb') as f:
        data = pickle.load(f)

    yd = data["Net_infilteration"]
    y = yd.loc[range_[0]: range_[1]].values / 10. * 24.  # mm -> cm, /hour -> /day
    evapd = data["Epot"]
    evap = -evapd.loc[range_[0]: range_[1]].values / 10. * 24.  # mm -> cm, /hour -> /day
    precip = y - evap

    # fig, ax = plt.subplots(2)
    # t_ = np.linspace(0, len(precip) / 24, len(precip))  # hourly data
    # ax[0].plot(t_, precip, 'r', label = "precipitation")  # precipitation positive sign
    # ax[0].plot(t_, evap, 'g', label = "evaporation")  # evaporation negative sign
    # ax[0].plot(t_, y, 'k:', label = 'net infiltration')
    # ax[0].plot(t_, precip + evap, 'k:', label = 'net infiltration')
    # ax[0].legend()
    # ax[0].set_ylabel("cm / day")

    """ 1. ET0 -> ETc """
    yd = data["ET0"]
    et0 = yd.loc[range_[0]: range_[1]].values / 10. * 24.  # mm -> cm, /hour -> /day
    etc = et0 * Kc
    """ 2. ETc -> Tpot, Evap """
    k = 0.6  # ?????????????????
    t_ = np.linspace(0, time, len(etc))
    tpot = np.multiply(etc, [(1. - np.exp(-k * lai_f(t_[i]))) for i in range(0, len(etc))])
    evap = etc - tpot
    evap = -evap

    # net_inf = precip + evap
    # # ax[1].plot(t_, tpot, 'r', label = "potential transpiration")
    # # ax[1].plot(t_, etc, 'k:', label = "crop evapotranspiration")
    # ax[1].plot(t_, precip, 'r', label = "precipitation")
    # ax[1].plot(t_, evap, 'g', label = "evaporation")
    # ax[1].plot(t_, net_inf, 'k:', label = 'net infiltration')
    # ax[1].legend()
    # ax[1].set_xlabel("time")
    # ax[1].set_ylabel("cm / day")
    # plt.show()

    y_ = []
    x_ = []
    for i in range(0, y.shape[0]):
        x_.extend([float(i) / 24, float(i + 1) / 24])  # hour -> day
        y_.extend([net_inf[i], net_inf[i]])

    return x_, y_

