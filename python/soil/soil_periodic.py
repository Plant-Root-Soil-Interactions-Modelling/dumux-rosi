import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import os
import time
import numpy as np
import plotly.graph_objects as go
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

""" 
Peridic example with water movement in soil using a constant sink,  
and plotly for vizualisation (very nice) (pip3 install plotly)
"""

cpp_base = RichardsSP()
s = RichardsWrapper(cpp_base)
s.initialize()

s.createGrid([-25, -25, -50], [25, 25, 0.], [19, 19, 19], True)
loam = [0.08, 0.43, 0.04, 1.6, 50.]
s.setVGParameters([loam])
s.setHomogeneousIC(-100, True)  # cm pressure head, hydraulic equilibrium
s.setTopBC("constantFlux", 0.)
s.setBotBC("freeDrainage")
s.initializeProblem()

if rank == 0:
    print("\nGrid bounds", s.getGridBounds(), "\n")

points = s.getDofCoordinates()  # gathered in rank = 0
cells = s.getCellCenters()  # gathered in rank = 0

dof = 0;
if rank == 0:
    dof = points.shape[0]
dof = comm.bcast(dof, root = 0)

# Test picking
p = [0, 0, 0]
id = s.pickCell(p)
print("Total dof of rank", rank, "=", dof, "picked id", id)
comm.barrier()
if rank == 0:
    print("Picked cell ", cells[id])
    print("Distance to element center", np.linalg.norm(cells[id] - p), "cm")

print()
p = [0, 25, -40]
id2 = s.pickCell(p)
print("Total dof of rank", rank, "=", dof, "picked id", id2)
comm.barrier()
if rank == 0:
    print("Picked cell ", cells[id2])
    print("Distance to element center", np.linalg.norm(cells[id2] - p), "cm")

sources = { id: 50., id2: 50.}  # gIdx: value [g/day]
s.setSource(sources)

# # Show inital condition
# x = np.array(s.getSolution())
# if rank == 0:
#     print("Water volume", s.getWaterVolume(), "cm3")
#     plt.plot(s.toHead(x), points[:, 2], "r*")
#     plt.show()

dt = 10  # ten days
s.ddt = 1.e-5

for i in range(0, 12):
    if rank == 0:
        print(i, "*** External time step ", dt, "***", "simulation time ", s.simTime)
    s.solve(dt)
    print("Water volume", s.getWaterVolume(), "cm3")

# s.writeDumuxVTK("test_periodic")

x = np.array(s.getWaterContent())

if rank == 0:

    nmin = np.min(x)
    nmax = np.max(x)
    N = s.numberOfCells
    X, Y, Z = np.mgrid[-25:25:(N[0] * 1j), -25:25:(N[1] * 1j), 0:-50:(N[2] * 1j)]  # Domain def in Python
    V = np.zeros(X.shape)
    for i in range(0, 19):
        for j in range(0, 19):
            for k in range(0, 19):
                gIdx = s.pickCell([X[i, j, k], Y[i, j, k], Z[i, j, k]])
                V[i, j, k] = x[gIdx]

    fig = go.Figure(data = go.Isosurface(
        x = X.flatten(),
        y = Y.flatten(),
        z = Z.flatten(),
        value = V.flatten(),
        isomin = nmin,
        isomax = nmax,
        surface_count = 7,  # number of isosurfaces, 2 by default: only min and max
        colorbar_nticks = 7,  # colorbar ticks correspond to isosurface values
        caps = dict(x_show = False, y_show = False)
        ))
    fig.show()
