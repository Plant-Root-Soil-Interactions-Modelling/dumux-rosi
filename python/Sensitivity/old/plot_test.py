import sys; sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from xylem_flux import sinusoidal2

kr0 = 1.e-3
kx0 = 5.e-2
#
# data = np.array([[ 388.53798698, 964.49242371, 1148.18488776, 1223.59927337, 1242.32979242,
#   1249.5923259, 1250.18067815, 1245.15006869],
#  [ 471.76008981, 1269.60478696, 1522.74642793, 1573.59446976, 1580.43109716,
#   1573.42027068, 1544.22337964, 1550.06910862],
#  [ 504.29002821, 1314.8764898, 1580.94318415, 1611.85825646, 1611.18179056,
#   1615.18818453, 1614.67587617, 1631.04912447],
#  [ 518.05613455, 1310.32294389, 1596.15965027, 1637.14239387, 1651.3538036,
#   1663.45244807, 1672.93605471, 1667.90279049],
#  [ 532.06703766, 1337.72439665, 1621.88494405, 1639.31597018, 1654.64631893,
#   1660.13009889, 1664.89651896, 1683.91677896],
#  [ 572.61764156, 1353.91147264, 1613.63515272, 1654.95301739, 1673.607607,
#   1671.09041348, 1667.52457152, 1664.60675902],
#  [ 566.97352486, 1361.70182999, 1617.8765438, 1674.95107156, 1674.46455763,
#   1678.73114635, 1687.97536225, 1681.72359464],
#  [ 536.4547878, 1363.97145456, 1606.50137654, 1672.10487028, 1696.3166485,
#   1699.2357844, 1695.19914428, 1689.68371378],
#  [ 562.63700193, 1380.71069084, 1616.13777153, 1683.34343175, 1682.41305744,
#   1683.62605925, 1689.10204033, 1687.001129],
#  [ 549.91731697, 1379.36182427, 1635.491459, 1665.46264371, 1652.76777671,
#   1705.29593605, 1685.72931712, 1706.75044642]])  # data set with stochasticity
#
data = np.array([[ 364.33242136, 933.3508926, 1120.39411515, 1191.10346069, 1214.74748129,
  1222.71676038, 1225.44879247, 1226.37093632],
 [ 449.17698496, 1214.14637338, 1475.85472958, 1525.20977716, 1531.90868173,
  1532.27258386, 1531.34671626, 1529.95571687],
 [ 475.53934743, 1260.36632515, 1532.13352081, 1578.15825898, 1587.85940803,
  1590.40849469, 1590.17539825, 1588.67600428],
 [ 490.28013414, 1282.68333505, 1554.89528704, 1600.07373791, 1611.46036076,
  1615.64887435, 1616.9631541, 1616.96909297],
 [ 500.34594309, 1297.13251987, 1568.63920654, 1613.52591202, 1626.04613643,
  1631.40670272, 1633.51382694, 1633.96781326],
 [ 508.37520141, 1308.03145121, 1578.317851, 1622.73134656, 1635.7937544,
  1641.66742903, 1644.59620424, 1645.91725839],
 [ 514.81571514, 1316.48554181, 1585.91199509, 1629.6428138, 1643.17236563,
  1649.42792856, 1652.50344792, 1654.09707984],
 [ 520.28551366, 1323.35894604, 1592.13608087, 1635.15761718, 1648.9351202,
  1655.58635324, 1658.97346796, 1660.62097346],
 [ 525.36414743, 1329.26654562, 1597.32448341, 1639.73752061, 1653.58877389,
  1660.51707632, 1664.23580476, 1666.15030207],
 [ 530.10157721, 1334.44201739, 1601.76687995, 1643.58275816, 1657.46421874,
  1664.52939766, 1668.52200349, 1670.70523161]])  # data set without stochasticity

fig, ax = plt.subplots(1, 1, figsize = (18, 10))
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size = '5%', pad = 0.05)

cmap_reversed = matplotlib.cm.get_cmap('jet_r')
# data[0, 0] = 0 # cell is left top
im = ax.imshow(data, cmap = cmap_reversed, aspect = 'auto', extent = [0.1 * kx0, 10 * kx0, 10 * kr0, 0.1 * kr0])  # interpolation = 'bicubic'

cb = fig.colorbar(im, cax = cax, orientation = 'vertical')
cb.ax.get_yaxis().labelpad = 30
cb.set_label('Cumulative uptake [cm3]', rotation = 270)

ax.set_xlabel("kx")
ax.set_ylabel("kr")
plt.show()
print(np.argmax(data))
print(data.flat[np.argmax(data)])

name = "results/transpiration_soybean_sra_sa_{:g}{:g}".format(0, 0)
i = 4
data = [np.load("results/transpiration_soybean_sra_sa_{:g}{:g}.npy".format(i, j)) for j in range(0, 8)]

trans = 10 * 0.6 * (38 * 5)  # cm3/day  for soybean
sim_time = 0.25 * 87.5
potential_trans = lambda t, dt: trans * sinusoidal2(t, dt) * t / sim_time  # soybean

dt_ = 360 / (24 * 3600)

n = len(data)
fig, ax = plt.subplots(1, 1, figsize = (18, 8))
ax = [ax]

for i in range(0, n):
    t = data[i][0]
    y = data[i][1]
    if trans > 0:
        ax[0].plot(t, potential_trans(t, dt_ * np.ones(t.shape)), 'k', label = "potential transpiration")  # potential transpiration
    ax[0].plot(t, y, 'g', label = "actual transpiration")  # actual transpiration  according to soil model
    ax[0].set_ylabel("transpiration [cm$^3$ day$^{-1}$]")
    ax[0].legend(loc = 'upper left')
    if i == 0:
        ax2 = ax[0].twinx()
    dt = np.diff(t)
    so = np.array(y)
    cup = np.cumsum(np.multiply(so[:-1], dt))
    ax2.plot(t[1:], cup, 'c--', label = "cumulative transpiration")  # cumulative transpiration (neumann)
    ax2.set_ylabel("cumulative [cm$^3$]")
    ax2.legend(loc = 'center right')
    print("cumulative uptake", cup[-1], "cm3")

ax[0].set_xlabel("Time [d]")

plt.tight_layout()

plt.show()
