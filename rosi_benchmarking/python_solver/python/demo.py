from dumux_rosi import RichardsSP  # C++ part (Dumux binding)
from solver.richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np

""" Benchmark M2.2 """
loam = vg.Parameters([0.08, 0.43, 0.04, 1.6, 50])
simtime = 10  # days
evap = -0.1  #  [cm/day]
Nz = 399  # resolution in z - direction

cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()

s.createGrid([-5., -5., -100.], [5., 5., 0.], [1, 1, Nz])  # [cm]
s.setHomogeneousIC(-200)  # cm pressure head
s.setTopBC("atmospheric", 0.5, [[0., 1.e10], [evap, evap]])  #  [cm/day] atmospheric is with surface run-off
s.setBotBC("freeDrainage")  # BC
s.setVGParameters([loam])
s.initializeProblem()

idx_top = s.pickCell([0.0, 0.0, 0.0])  # index to watch flux

dt = simtime / N  # time step
maxDt = 1  # maximal Dumux time step [days]

x_, y_ = [], []

for i in range(0, N):
    s.solve(dt, maxDt)
    f = s.getNeumann(idx_top)  #
    x_.append(s.simTime)
    y_.append(f)

plt(x_, y_)
plt.show()
