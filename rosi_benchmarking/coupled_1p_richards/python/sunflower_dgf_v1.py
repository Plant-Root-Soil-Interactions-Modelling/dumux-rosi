"""analysis of results using signed distance functions"""
import sys
from cmath import pi
sys.path.append("../../..")
import plantbox as pb
import numpy as np
import matplotlib.pyplot as plt

path = "../modelparameter/rootsystem/"
name = "Heliantus_Pagès_2013"  # ""

rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")

# Create and set geometry
rs.setMinDx(1.e-3)
x0 = pb.Vector3d(0., 0., -1.)
nx = pb.Vector3d(1., 0., -1.)
ny = pb.Vector3d(0., 1., -1.)
soil_layer = pb.SDF_HalfPlane(x0, nx, ny)  # there was bug, with updated CPlantBox
rs.setGeometry(soil_layer)

rs.setSeed(6)
rs.initialize()

simtime = 154.  # days
dt = 1.
N = round(simtime / dt)  # steps

fig, ax1 = plt.subplots()
fig, ax2 = plt.subplots()
fig, ax3 = plt.subplots()
# Plot length and number of segments over time
stype = "length"
v_, v1_, v2_, v3_ = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
segment_ = np.zeros(N)
for i in range(0, N):
    rs.simulate(dt)
    t = np.array(rs.getParameter("type"))
    v = np.array(rs.getParameter(stype))
    segment = np.array(rs.getNumberOfSegments())
    v_[i] = np.sum(v)
    v1_[i] = np.sum(v[t == 1])
    v2_[i] = np.sum(v[t == 2])
    v3_[i] = np.sum(v[t == 3])
    segment_[i] = np.sum(segment)

t_ = np.linspace(dt, N * dt, N)
ax1.plot(t_, v_, t_, v1_, t_, v2_, t_, v3_)
ax1.set_xlabel("time (days)")
ax1.set_ylabel(stype + " (cm)")
ax1.legend(["total", "tap root", "lateral", "2. order lateral"])
ax2.plot(t_, segment_)
ax2.set_xlabel("time (days)")
ax2.set_ylabel("number of segments")

# frequency plot of segment length
ana = pb.SegmentAnalyser(rs)
length = np.array(ana.getParameter("length"))
print(length)
print(min(length))
bins=[0, 1e-3, 1, 2]
hist, binEdges = np.histogram(length, bins)
ax3.bar(range(len(hist)),hist,width=1,align='center',tick_label=
        ['{} - {}'.format(bins[i],bins[i+1]) for i,j in enumerate(hist)])

plt.show()
