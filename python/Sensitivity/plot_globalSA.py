"""
plots results of local sensitivity analysis
"""
import sys; sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import run_SA as sa

""" def SA """
file_name = "global_maize_inc2"
path = "results/"

analysis_time = 30  # days

names, ranges = sa.read_ranges(path + file_name)

kr_ = ranges[0]   # 0 is y-axis
kx_ = ranges[1]   # 1 is x-axis (in imshow)
n = len(kr_) * len(kx_)
print(n)

trans_ = np.load(path + "transpiration_" + file_name + "2" + ".npy")
times = trans_[0,:]
print("Simulation time from", min(times), "to ", max(times), "days")

ind_ = np.argwhere(times > analysis_time)
if len(ind_) > 0:
    ind = ind_[0][0]
    print("plots for day", times[ind])
    ind += 1
    ind10 = ind // 10  # TODO check
    print(ind)
else:
    ind = -1
    ind10 = -1

dt_ = np.diff(times[:ind])

print("index is", ind, "from", len(times))

sum_trans = 0
print("reading global analysis")
data = np.zeros((len(kr_), len(kx_)))
for i in range(0, len(kr_)):
    for j in range(0, len(kx_)):
        lind = i * len(kx_) + j
        print(lind)
        try:
            trans_ = np.load(path + "transpiration_" + file_name + str(2 + lind) + ".npy")  # TODO try catch
            sum_trans = np.sum(np.multiply(trans_[1,:ind - 1], dt_))
            # print(lind, sum_trans)
        except:
            print(lind, i, j)
            sum_trans = np.nan
        data[i, j] = sum_trans

print(data.shape)

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
ext_ = [np.log10(kx_[0]) , np.log10(kx_[-1]), np.log10(kr_[-1]), np.log10(kr_[0])]

fig, ax = plt.subplots(1, 1, figsize = (18, 10))
ax = [ax]
cmap = matplotlib.cm.get_cmap('nipy_spectral')
im = ax[0].imshow(data, cmap = cmap, aspect = 'auto', interpolation = 'nearest', extent = ext_)  # vmin = 0., vmax = 1.e-3, extent = [kx_[0] , kx_[-1], kr_[-1], kr_[0]],
ax[0].set_ylabel("kr scale (10$^x$)")
ax[0].set_xlabel("kx scale (10$^x$)")
cb = fig.colorbar(im, orientation = 'vertical', label = "water uptake")
plt.show()

