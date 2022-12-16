import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def sigmoid(x, L , x0, k, b):
    y = L / (1 + np.exp(-k * (x - x0))) + b
    return (y)


def exp2(x, L , x0, k, b):
    y = L * np.exp(-k * (x - x0) * (x - x0)) + b
    return (y)


times_soybean = [20., 30, 50, 60, 70, 80, 90 ]
lai_soybean = [0.2, 0.3, 3., 4.1, 6.5, 8., 7. ]  # Priscila et al. (2013)

times_maize = [0., 20, 40, 60, 80, 100]
lai_maize = [0.01, 0.5, 3.1, 4., 3.5, 2.5 ]  # Boedhram et al. (2001)

p0_soy = [max(lai_soybean), np.median(times_soybean), 1, min(lai_soybean)]  # this is an mandatory initial guess
p0_maize = [max(lai_maize), np.mean(times_maize), 1. / np.std(times_maize), min(lai_maize)]

popt1, pcov = curve_fit(sigmoid, times_soybean, lai_soybean, p0_soy, method = 'dogbox')
popt2, pcov = curve_fit(sigmoid, times_maize, lai_maize, p0_maize, method = 'dogbox')
popt3, pcov = curve_fit(exp2, times_maize, lai_maize, p0_maize, method = 'dogbox')

# popt2 = p0_maize

x1_ = np.linspace(np.min(times_soybean), np.max(times_soybean), 100)
y1_ = np.array([sigmoid(x, popt1[0], popt1[1], popt1[2], popt1[3]) for x in x1_])

x2_ = np.linspace(np.min(times_maize), np.max(times_maize), 100)
y2_ = np.array([sigmoid(x, popt2[0], popt2[1], popt2[2], popt2[3]) for x in x2_])
y3_ = np.array([exp2(x, popt3[0], popt3[1], popt3[2], popt3[3]) for x in x2_])

fig, ax = plt.subplots(2, 1, figsize = (8, 16))

ax[0].plot(times_soybean, lai_soybean, "*")
ax[0].plot(x1_, y1_)
ax[0].set_title("Soybean")
# ax[0].xlabel("Time")
ax[0].set_ylabel("LAI")

ax[1].plot(times_maize, lai_maize, "*")
ax[1].set_title("Maize")
ax[1].plot(x2_, y2_)
ax[1].plot(x2_, y3_)
ax[1].set_xlabel("Time")
ax[1].set_ylabel("LAI")
plt.show()

print("soybean")
print(popt1)

print("maize")
print("sigmoidal")
print(popt2)
print("neg exp^2")
print(popt3)
