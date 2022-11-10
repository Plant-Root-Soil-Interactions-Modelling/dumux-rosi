import matplotlib.pyplot as plt
import pickle
import pandas as pd
import numpy as np
import datetime
from evapotranspiration import *

with open('data/95.pkl', 'rb') as f:
    data = pickle.load(f)

area_maize = 75 * 15  # cm2
area_soybean = 38 * 5  # cm2

# print(type(data))
# print(data)

# wie krieg ich das aus dem panda table????

# print(data.info())
# print()
# print(data["Tpot"].info())
#
# ndata = data.to_numpy()
# print(ndata.shape)
# print()
# print(ndata[:, -1])

print(data.info())

""" Tpot """
# x = np.array(data.index.to_pydatetime(), dtype = np.datetime64)
# y = data["Tpot"]
# ndata = data.to_numpy()
#
# plt.plot(y.loc['1995-03-15 00:00:00': '1995-06-10 12:00:00'])  # 15 + 31 + 30 + 11.5 = 87.5
# plt.show()
#
# # print(y.loc['1994-11-06 00:00:00': '1995-01-06 00:00:00'].values)
# y = y.loc['1995-03-15 00:00:00': '1995-06-10 11:00:00'].values
# print(y.shape)

# y = data["Net_infilteration"]
# y = y.loc['1995-03-15 00:00:00': '1995-06-10 11:00:00'].values
# plt.plot(y)
# plt.show()
# print(y.shape)
Kc_maize = 1.2  # book "crop evapotranspiration" Allen, et al 1998
Kc_soybean = 1.15  # book "crop evapotranspiration" Allen, et al 1998

range_soybean = ['1995-03-15 00:00:00', '1995-06-10 11:00:00']
range_maize = ['1995-03-15 00:00:00', '1995-06-17 23:00:00']
get_transpiration_beers('data/95.pkl', range_soybean, area_soybean, 87.5, lai_f = lai_soybean, Kc = Kc_soybean)
get_transpiration_beers('data/95.pkl', range_maize, area_maize, 95, lai_f = lai_maize, Kc = Kc_maize)

# net_infiltration_table_beers('data/95.pkl', range_, 87.5, lai_f = lai_soybean, Kc = Kc_soybean)
# net_infiltration_table_beers('data/95.pkl', range_, 87.5, lai_f = lai_maize, Kc = Kc_maize)
