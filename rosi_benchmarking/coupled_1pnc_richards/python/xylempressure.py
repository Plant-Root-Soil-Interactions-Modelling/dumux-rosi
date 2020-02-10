''' xylem pressure '''

import matplotlib.pyplot as plt
import math
import numpy as np


cL = 1
psi_x = np.linspace(-15000, 0, 100)
if psi_x < -5500:
	alpha = math.exp(-50,000*cL)*math.exp(-1.02e-6*(psi_x + 5500))
else:
	alpha = math.exp(-50,000*cL)

plt.plot(psi_x, alpha) 
plt.set_xlabel("Xylem pressure (cm)")
plt.set_ylabel("$alpha$")
plt.show()
