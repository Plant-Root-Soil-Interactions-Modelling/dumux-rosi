""" 
    Macroscopic soil model (TODO)
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import matplotlib.pyplot as plt
import numpy as np

import scenario_setup as scenario

""" maize domain  """
min_b = [-37.5, -7.5, -200.]  # Domain Mais: 60 cm Reihe, 10 cm Pflanzen
max_b = [37.5, 7.5, 0.]
# cell_number = [75, 15, 200]  # 3D 1 cm3
cell_number = [1, 1, 200]  # 1D: 75x15x1 cm3
simtime = 95  # between 90-100 days

""" soy domain  """
# min_b = [-19, -2.5, -200.]  # Domain [38 cm Reihe, 6 cm Pflanzen]
# max_b = [19, 2.5, 0.]
# cell_number = [38, 5, 200]  # 1 cm3
# simtime = 87.5  # between 75-100 days

""" soil type """
theta_r = 0.025  # sandy loam
theta_s = 0.403
alpha = 0.0383  # (cm-1) soil
n = 1.3774
k_sat = 60.  # (cm d-1)
soil_ = [theta_r, theta_s, alpha, n, k_sat]

""" set up simulator """
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = -330, p_bot = -180)

""" simulation loop """
c = []  # resulting solute concentration
dt = 1
N = int(np.ceil(simtime / dt))
for i in range(0, N):
    t = i * dt  # current simulation time
    print(t)

    if t < 1:
        s.setSoluteTopBC([1], [1.e-4])  # put something in for the first day ...
    else:
        s.setSoluteTopBC([2], [0.])

    s.solve(dt)
    c.append(s.getSolution(1))

z = s.getDofCoordinates()

if rank == 0:
    for i in range(0, len(c)):
        plt.plot(c[i], z[:, 2])
    plt.show()
