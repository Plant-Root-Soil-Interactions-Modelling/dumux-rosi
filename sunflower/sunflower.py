"""analysis of results using signed distance functions"""
import sys
from cmath import pi
sys.path.append("../../..")
import plantbox as pb

import numpy as np
import matplotlib.pyplot as plt

path = "../../../modelparameter/rootsystem/"
name = "Heliantus_Pagès_2013"

rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")

# Soil core analysis
depth, layers = 150., 100
z_ = np.linspace(0, -1 * depth, layers)

# Create and set geometry

# 1. creates a square 50*33 cm containter with height 150 cm
rhizotron = pb.SDF_PlantBox(50, 33, 150)
rs.setGeometry(rhizotron)

rs.setSeed(0)
rs.initialize()
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first

rs.simulate(7)
rs.write("results/sunflower_7days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
ana.write("results/sunflower_7days.dgf")

rs.simulate(7)
rs.write("results/sunflower_14days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
ana.write("results/sunflower_14days.dgf")

rs.simulate(7)
rs.write("results/sunflower_21days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
ana.write("results/sunflower_21days.dgf")

rs.simulate(9)
rs.write("results/sunflower_30days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
ana.write("results/sunflower_30days.dgf")

rs.simulate(60)
rs.write("results/sunflower_90days.vtp")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first
ana.write("results/sunflower_90days.dgf")

#ana = pb.SegmentAnalyser(rs)
#aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
#for s in aseg:
#    print("Shoot segment", s)
#    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first

#ana.write("results/sunflower_90days.dgf")

# Make a root length distribution
#ana = pb.SegmentAnalyser(rs)
#ana.cropDomain(50, 33, depth)
#layerVolume = depth / layers * 50 * 33
#rl0_ = ana.distribution("length", 0., -depth, layers, True)
#plt.plot(np.array(rl0_) / layerVolume, z_, label = '90 days')
#plt.xlabel('RLD (cm/cm^3)')  # layer size is 1 cm
#plt.ylabel('Depth (cm)')
#plt.legend(["90 days"])


# 2. creates a square 50*33 cm containter with height 150 cm
#rhizotron = pb.SDF_PlantBox(50, 33, 150)
#rs.setGeometry(rhizotron)

#rs.setSeed(0)
#rs.initialize()
#rs.simulate(30)
#rs.write("results/sunflower_30days.vtp")

#ana = pb.SegmentAnalyser(rs)
#aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
#for s in aseg:
#    print("Shoot segment", s)
#    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first

#ana.write("results/sunflower_30days.dgf")

# Make a root length distribution
#ana = pb.SegmentAnalyser(rs)
#ana.cropDomain(50, 33, depth)
#layerVolume = depth / layers * 50 * 33
#rl0_ = ana.distribution("length", 0., -depth, layers, True)
#plt.plot(np.array(rl0_) / layerVolume, z_, label = '30 days')
#plt.xlabel('RLD (cm/cm^3)')  # layer size is 1 cm
#plt.ylabel('Depth (cm)')
#plt.legend(["30 days"])

# 3. creates a square 50*33 cm containter with height 150 cm
#rhizotron = pb.SDF_PlantBox(50, 33, 150)
#rs.setGeometry(rhizotron)

#rs.setSeed(0)
#rs.initialize()
#rs.simulate(7)
#rs.write("results/sunflower_7days.vtp")

#ana = pb.SegmentAnalyser(rs)
#aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
#for s in aseg:
#    print("Shoot segment", s)
#    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first

#ana.write("results/sunflower_7days.dgf")

# Make a root length distribution
#ana = pb.SegmentAnalyser(rs)
#ana.cropDomain(50, 33, depth)
#layerVolume = depth / layers * 50 * 33
#rl0_ = ana.distribution("length", 0., -depth, layers, True)
#plt.plot(np.array(rl0_) / layerVolume, z_, label = '7 days')
#plt.xlabel('RLD (cm/cm^3)')  # layer size is 1 cm
#plt.ylabel('Depth (cm)')
#plt.legend(["7 days"])




#fig.subplots_adjust()
#plt.legend()
#plt.savefig("results/sunflower_RLD_comparsion.pdf")
#plt.show()
