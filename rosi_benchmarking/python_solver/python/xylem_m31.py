import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

from solver.xylem_flux import XylemFluxPython  # Python hybrid solver
import solver.plantbox as pb

import numpy as np

import matplotlib.pyplot as plt

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

""" 
Steady state vertical root with the Python/cpp Hybrid solver
does not work in parallel
"""

rs = pb.MappedRootSystem()

# create a straigth root by hand (todo: friendlier setter...)
N = 100
z_ = np.linspace(0., -50., N)
nodes, node_cts, segs, radii, types = [], [], [], [], []
for z in z_:
    nodes.append(pb.Vector3d(0, 0, z))
    node_cts.append(0.)
for s in range(0, N - 1):
    segs.append(pb.Vector2i(s, s + 1))
    radii.append(0.1)
    types.append(0)
rs.nodes = nodes
rs.nodeCTs = node_cts
rs.segments = segs
rs.radii = radii
rs.types = types

r = XylemFluxPython(rs)
r.setKr(1.73e-4 * np.ones((1, 1)))
r.setKx(4.32e-2 * np.ones((1, 1)))

psi_s = -200

rx_hom = r.solve(0., True)
# rx = r.getSolution(rx_hom, [psi_s])

plt.plot(rx_hom, z_)
plt.show()

