import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

from richardssp import RichardsSP  # C++ part (Dumux binding)
from solver.richards import RichardsWrapper  # Python part
from solver.xylem_flux import XylemFlux

import numpy as np
import matplotlib.pyplot as plt

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()

# Set up the soil problem
s.createGrid([-10, -10, -30], [10, 10, 0.], [19, 19, 29], False)
loam = [0.08, 0.43, 0.04, 1.6, 50.]
s.setVGParameters([loam])
s.setHomogeneousIC(-200, True)  # cm pressure head, hydraulic equilibrium
s.setTopBC("constantFlux", 0.)
s.setBotBC("freeDrainage")
s.initializeProblem()

# Set up the root problem
r = XylemFlux(s)
path = "../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
r.readParameters(path + name + ".xml")  # Open plant and root parameter from a file
r.setKr(1.73e-4 * np.ones((1, 1)))
r.setKx(4.32e-2 * np.ones((1, 1)))
r.simulate(10)

# q_r = 0.5
# q_r = q_r * 1  # * 1 g/cm3 = g/cm2/day
# q_r = q_r * 2 * r_root * np.pi * 1.  # g/day TODO = ?
# print("Qr as sink", q_r)
# trans = -2e-8  # kg / s
# trans = -trans * 24 * 3600 * 1000  # g /day
# print(trans)

sx_old = s.toHead(s.getSolution())
rx_old = r.solve(-15000, sx_old, False)

dt = 1. / 3600. / 24.  # 1s = days

for i in range(0, 60):

    fluxes = r.eval_fluxes(sx_old, rx_old)

    s.setSource(fluxes)
    s.solve(dt)
    sx = s.toHead(s.getSolution())

    r.simulate(dt)  # simulates root growth
    rx = r.solve(-15000, sx_old, False)  # solves water flow

    sx_old = sx
    rx_old = rx

nodes = r.convert_(r.rs.getNodes())
plt.plot(rx_old, nodes[:, 2])
plt.show()

print("fin.")

