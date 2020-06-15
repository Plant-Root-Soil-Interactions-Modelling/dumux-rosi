import numpy as np

a = np.loadtxt("dumux_c12_100", delimiter=',')
np.savetxt("dumux_c12_10000", a, delimiter=';')
