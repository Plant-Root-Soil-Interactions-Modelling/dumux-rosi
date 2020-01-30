import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

from richardssp import RichardsSP  # C++ part (Dumux binding)
from solver.richards import RichardsWrapper  # Python part

from solver.xylem_flux import XylemFlux

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()

# Set up a soil problem
s.createGrid([-10, -10, -30], [10, 10, 0.], [19, 19, 29], True)
loam = [0.08, 0.43, 0.04, 1.6, 50.]
s.setVGParameters([loam])
s.setHomogeneousIC(-100, True)  # cm pressure head, hydraulic equilibrium
s.setTopBC("constantFlux", 0.)
s.setBotBC("freeDrainage")
s.initializeProblem()

print()
r = XylemFlux(s)
# Open plant and root parameter from a file
path = "../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
r.readParameters(path + name + ".xml")

r.simulate(10)

print("fin.")
