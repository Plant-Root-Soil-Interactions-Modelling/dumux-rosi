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
Benchmark M3.1 Single root: steady state vertical root solved with the Python/cpp Hybrid solver
(does not work in parallel)
"""

""" Parameters """
g = 9.8065 * 100.*24.*3600.*24.*3600.  # gravitational acceleration [cm day-2]
rho = 1.  # density of water, [g/cm^3]
L = 50  # length of single straight root [cm]
a = 0.2  # radius [cm] <--------------------------------------------------------- ???
kz0 = 4.32e-2  # [cm^3/day]
kz = kz0 / (rho * g)  # axial conductivity [cm^5 s / g]
kr0 = 1.728e-4  # [1/day]
kr = kr0 / (rho * g)  # radial conductivity per root type [cm^2 s / g]
p_s = -200  # static soil pressure [cm]
p0 = -1000  # dircichlet bc at top

""" Analytical solution """

c = 2 * a * pi * kr / kz
p_r = lambda z: p_s + d[0] * exp(sqrt(c) * z) + d[1] * exp(-sqrt(c) * z)  #

AA = np.array([[1, 1], [sqrt(c) * exp(-sqrt(c) * L), -sqrt(c) * exp(sqrt(c) * L)] ])  # # Boundary conditions dirichlet top, neumann bot
bb = np.array([p0 - p_s, -1])  # -rho * g
d = np.linalg.solve(AA, bb)  # compute constants d_1 and d_2 from bc

za_ = np.linspace(0, -L, 100)  # Evaluate function
pr = list(map(p_r, za_))

plt.plot(pr, za_)

""" Numeric solution """

N = 100  # resolution
z_ = np.linspace(0., -L, N)

nodes, segs, radii = [], [], []
for z in z_:
    nodes.append(pb.Vector3d(0, 0, z))
for s in range(0, N - 1):
    segs.append(pb.Vector2i(s, s + 1))
    radii.append(a)

rs = pb.MappedSegments(nodes, segs, radii)

r = XylemFluxPython(rs)
r.setKr([kr0])
r.setKx([kz0])

rx_hom = r.solve(0., p0 - p_s, False)
rx = r.getSolution(rx_hom, [p_s])

plt.plot(rx, z_, "r*")
plt.xlabel("Xylem pressure (cm)")
plt.ylabel("Depth (m)")
plt.show()
