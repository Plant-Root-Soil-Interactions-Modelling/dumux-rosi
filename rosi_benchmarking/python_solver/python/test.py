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

# Set up a soil problem
s.createGrid([-10, -10, -30], [10, 10, 0.], [19, 19, 29], periodic = True)
loam = [0.08, 0.43, 0.04, 1.6, 50.]
s.setVGParameters([loam])
s.setHomogeneousIC(-500, True)  # cm pressure head, hydraulic equilibrium
s.setTopBC("constantFlux", 0.)
s.setBotBC("freeDrainage")
s.initializeProblem()

r = XylemFlux(s)
# Open plant and root parameter from a file
path = "../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
r.readParameters(path + name + ".xml")

r.setKr(1.73e-4 * np.ones((1, 1)))
r.setKx(4.32e-2 * np.ones((1, 1)))

r.simulate(10)

trans = -2e-8  # kg / s
trans = -trans * 24 * 3600 * 1000  # g /day

sx_old = s.toHead(s.getSolution())
rx_old = r.solve(trans, sx_old, False)

dt = 1. / 3600. / 24.  # 1s = days

for i in range(0, 60):

    fluxes = r.eval_fluxes(sx_old, rx_old)
    s.setSource(fluxes)
    s.simulate(dt)
    sx = s.toHead(s.getSolution())

    r.simulate(dt)  # simulates root growth
    rx = r.solve(trans, sx_old, False)  # solves water flow

    sx_old = sx
    rx_old = rx

print("fin.")

