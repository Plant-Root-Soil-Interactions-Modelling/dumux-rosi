import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

from dumux_rosi import RichardsSP  # C++ part (Dumux binding)
from solver.richards import RichardsWrapper  # Python part

from solver.xylem_flux import XylemFluxPython  # Python hybrid solver
import solver.plantbox as pb

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
root_system = pb.MappedRootSystem()

picker = lambda x, y, z : cpp_base.pick(x, y, z)  # I had problem with a common interface over different modules (plantbox, dumux-rosi)
root_system.setSoilGrid(picker)
path = "../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
root_system.readParameters(path + name + ".xml")  # Open plant and root parameter from a file
root_system.initialize()
root_system.simulate(1)

r = XylemFluxPython(root_system)
r.setKr(1.73e-4 * np.ones((1, 1)))
r.setKx(4.32e-2 * np.ones((1, 1)))

# q_r = 0.5
# q_r = q_r * 1  # * 1 g/cm3 = g/cm2/day
# q_r = q_r * 2 * r_root * np.pi * 1.  # g/day TODO = ?
# print("Qr as sink", q_r)
# trans = -2e-8  # kg / s
# trans = -trans * 24 * 3600 * 1000  # g /day
# print(trans)

sx_old = s.toHead(s.getSolution())  # inital condition
rx_old = r.solve(-15000, False)

dt = 1. / 3600. / 24.  # 1s = days

for i in range(0, 60):

    root_system.simulate(dt, True)  # simulates root growth

    r_hom = r.solve(-15000, False)  # solves water flow
    rx = r.getSolution(r_hom, sx_old)  # make inhomogeneous solution from homogeneous one

    # new fluxes are based on sx_old and solution rx

    fluxes = r.soil_fluxes(rx, sx_old)  # soil source/sink fluxes equal root system radial fluxes
    s.setSource(fluxes)
    s.solve(dt)
    sx = s.toHead(s.getSolution())

    sx_old = sx
    rx_old = rx

nodes = r.convert_(r.rs.getNodes())
plt.plot(rx_old, nodes[:, 2])
plt.show()

print("fin.")

