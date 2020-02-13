import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/python_solver/")

from solver.xylem_flux import XylemFluxPython  # Python hybrid solver
import solver.plantbox as pb

from math import *
import numpy as np
from numpy.linalg.linalg import norm
from scipy import sparse
import scipy.sparse.linalg as LA
# from rsml_reader import *  # located in the same directory

import matplotlib.pyplot as plt

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

""" 
Steady state vertical root with the Python/cpp Hybrid solver
does not work in parallel
"""
# Parameters
g = 9.8065 * 100.*24.*3600.*24.*3600.  # gravitational acceleration [cm day-2]
rho = 1.  # density of water, [g/cm^3]

L = 50  # length of single straight root [cm]
a = 0.2  # radius [cm] <--------------------------------------------------------- ???
kz0 = 4.32e-2  # [cm^3/day]
kz = kz0 / (rho * g)  # axial conductivity [cm^5 s / g]
kr0 = 1.728e-4  # [1/day]
kr = kr0 / (rho * g)  # radial conductivity per root type [cm^2 s / g]

p_s = -200  # static soil pressure [cm]
p0 = -500  # dircichlet bc at top

""" Analytical solution """

c = 2 * a * pi * kr / kz
p_r = lambda z: p_s + d[0] * exp(sqrt(c) * z) + d[1] * exp(-sqrt(c) * z)  #
# Boundary conditions
AA = np.array([[1, 1], [sqrt(c) * exp(-sqrt(c) * L), -sqrt(c) * exp(sqrt(c) * L)] ])  # dirichlet top, neumann bot
bb = np.array([p0 - p_s, -1])  # -rho * g
d = np.linalg.solve(AA, bb)  # compute constants d_1 and d_2 from bc

za_ = np.linspace(0, -L, 100)  # Evaluate function
pr = list(map(p_r, za_))

plt.plot(pr, za_)

""" Numeric solution """

rs = pb.MappedSegments()

# create a straigth root by hand (todo: friendlier setter...)
N = 100  # resolution
z_ = np.linspace(0., -L, N)
nodes, node_cts, segs, radii, types = [], [], [], [], []
for z in z_:
    nodes.append(pb.Vector3d(0, 0, z))
    node_cts.append(0.)
for s in range(0, N - 1):
    segs.append(pb.Vector2i(s, s + 1))
    radii.append(a)
    types.append(0)
rs.nodes = nodes
rs.nodeCTs = node_cts
rs.segments = segs
rs.radii = radii
rs.types = types

r = XylemFluxPython(rs)
r.setKr(kr0 * np.ones((1, 1)))
r.setKx(kz0 * np.ones((1, 1)))

rx_hom = r.solve(-500 - p_s, False)
rx = r.getSolution(rx_hom, [p_s])
#
plt.plot(rx, z_, "r*")
plt.xlabel("Xylem pressure (cm)")
plt.ylabel("Depth (m)")
plt.show()
