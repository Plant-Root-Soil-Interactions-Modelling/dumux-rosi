"""analysis of results using signed distance functions"""
import sys
from cmath import pi
sys.path.append("../../..")
import plantbox as pb
import numpy as np
import matplotlib.pyplot as plt

path = "../modelparameter/rootsystem/"
name = "Glycine_max"  # ""

rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")

# Create and set geometry
rs.setMinDx(1.e-3)
x0 = pb.Vector3d(0., 0., -1.)
nx = pb.Vector3d(1., 0., -1.)
ny = pb.Vector3d(0., 1., -1.)
soil_layer = pb.SDF_HalfPlane(x0, nx, ny)  # there was bug, with updated CPlantBox
rs.setGeometry(soil_layer)

rs.setSeed(2)
rs.initialize()

rs.simulate(1, True)
rs.write("results/Glycine_max_1days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    # print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
    vt = rs.getSummed("volume")
print("Volume (cm3)", vt)
ana.write("results/Glycine_max_1days.dgf")

l = np.array(ana.getParameter("length"))
print("Min ", np.min(l))

rs.simulate(6, True)
rs.write("results/Glycine_max_7days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    # print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
    vt = rs.getSummed("volume")
print("Volume (cm3)", vt)
ana.write("results/Glycine_max_7days.dgf")

l = np.array(ana.getParameter("length"))
print("Min ", np.min(l))

rs.simulate(7, True)
rs.write("results/Glycine_max_14days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    # print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
    vt = rs.getSummed("volume")
print("Volume (cm3)", vt)
ana.write("results/Glycine_max_14days.dgf")

l = np.array(ana.getParameter("length"))
print("Min ", np.min(l))

rs.simulate(7, True)
rs.write("results/Glycine_max_21days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    # print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
    vt = rs.getSummed("volume")
print("Volume (cm3)", vt)
ana.write("results/Glycine_max_21days.dgf")

l = np.array(ana.getParameter("length"))
print("Min length", np.min(l))
a = np.array(ana.getParameter("radius"))
print("Min radius", np.min(a))

rs.simulate(9)
rs.write("results/Glycine_max_30days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
    vt = rs.getSummed("volume")
print("Volume (cm3)", vt)
ana.write("results/Glycine_max_30days.dgf")

l = np.array(ana.getParameter("length"))
print("Min length", np.min(l))
a = np.array(ana.getParameter("radius"))
print("Min radius", np.min(a))

rs.simulate(15)
rs.write("results/Glycine_max_45days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
    vt = rs.getSummed("volume")
print("Volume (cm3)", vt)
ana.write("results/Glycine_max_45days.dgf")

l = np.array(ana.getParameter("length"))
print("Min length", np.min(l))
a = np.array(ana.getParameter("radius"))
print("Min radius", np.min(a))

rs.simulate(15)
rs.write("results/Glycine_max_60days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
    vt = rs.getSummed("volume")
print("Volume (cm3)", vt)
ana.write("results/Glycine_max_60days.dgf")

l = np.array(ana.getParameter("length"))
print("Min length", np.min(l))
a = np.array(ana.getParameter("radius"))
print("Min radius", np.min(a))

rs.simulate(30)
rs.write("results/Glycine_max_90days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
    vt = rs.getSummed("volume")
print("Volume (cm3)", vt)
ana.write("results/Glycine_max_90days.dgf")

l = np.array(ana.getParameter("length"))
print("Min length", np.min(l))
a = np.array(ana.getParameter("radius"))
print("Min radius", np.min(a))

rs.simulate(64)
rs.write("results/Glycine_max_154days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
    vt = rs.getSummed("volume")
print("Volume (cm3)", vt)
ana.write("results/Glycine_max_154days.dgf")

l = np.array(ana.getParameter("length"))
print("Min length", np.min(l))
a = np.array(ana.getParameter("radius"))
print("Min radius", np.min(a))
