import matplotlib.pyplot as plt
from numpy import genfromtxt
import numpy as np

data_a_act = np.genfromtxt("a_actual.out")
data_p = np.genfromtxt("a_pw_min_root.out")

collar_p = (data_p[:,1] - 1.e5) * 100. / 1.e3 / 9.81

# prepare plot
fig, ax = plt.subplots()
ax.plot(data_a_act[:,0], data_a_act[:,1]/0.02 + collar_p, "k--", linewidth=1, label="reference $\psi_{s, eq}$")
ax.set_xlabel('Time [days]')
ax.set_ylabel('Water potential [cm]')
ax.legend(loc = 'upper right')
fig.savefig("eswp_reference.pdf", dpi=300, bbox_inches='tight')
plt.show()
