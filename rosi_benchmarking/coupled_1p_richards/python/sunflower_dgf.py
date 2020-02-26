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

rs.setSeed(3)
rs.initialize()

rs.simulate(7, True)
rs.write("results/sunflower_7days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    # print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
ana.write("results/sunflower_7days.dgf")

l = np.array(ana.getParameter("length"))
print("Min ", np.min(l))

rs.simulate(7, True)
rs.write("results/sunflower_14days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    # print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
ana.write("results/sunflower_14days.dgf")

l = np.array(ana.getParameter("length"))
print("Min ", np.min(l))

rs.simulate(7, True)
rs.write("results/sunflower_21days_2.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    # print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
ana.write("results/sunflower_21days.dgf")

l = np.array(ana.getParameter("length"))
print("Min length", np.min(l))
a = np.array(ana.getParameter("radius"))
print("Min radius", np.min(a))

rs.simulate(9)
rs.write("results/sunflower_30days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
ana.write("results/sunflower_30days.dgf")

l = np.array(ana.getParameter("length"))
print("Min length", np.min(l))
a = np.array(ana.getParameter("radius"))
print("Min radius", np.min(a))

rs.simulate(15)
rs.write("results/sunflower_45days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
ana.write("results/sunflower_45days.dgf")

l = np.array(ana.getParameter("length"))
print("Min length", np.min(l))
a = np.array(ana.getParameter("radius"))
print("Min radius", np.min(a))

rs.simulate(15)
rs.write("results/sunflower_60days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
ana.write("results/sunflower_60days.dgf")

l = np.array(ana.getParameter("length"))
print("Min length", np.min(l))
a = np.array(ana.getParameter("radius"))
print("Min radius", np.min(a))

rs.simulate(30)
rs.write("results/sunflower_90days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
ana.write("results/sunflower_90days.dgf")

l = np.array(ana.getParameter("length"))
print("Min length", np.min(l))
a = np.array(ana.getParameter("radius"))
print("Min radius", np.min(a))

rs.simulate(64)
rs.write("results/sunflower_154days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
ana.write("results/sunflower_154days.dgf")

l = np.array(ana.getParameter("length"))
print("Min length", np.min(l))
a = np.array(ana.getParameter("radius"))
print("Min radius", np.min(a))
