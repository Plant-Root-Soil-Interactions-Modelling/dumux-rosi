import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver

import matplotlib.pyplot as plt
import pickle
import pandas as pd
import numpy as np
import datetime
from evapotranspiration import *

area_maize = 76 * 16  # cm2
area_soybean = 76 * 3  # cm2

maize_per_m2 = (100 * 100) / area_maize

tot_use = np.array([0.8, 1.8, 2.9, 1.8, 3.8, 3.8, 1.9, 3.8, 3.8, 1.4]) * 2.54  # cm
print(np.sum(tot_use))
print(maize_per_m2)
