#!/usr/bin/env python
import sys
print sys.path
import os
import read_hydrus

read_hydrus.SavconToVtkMaster(filename = "Geometry_4.in", folder="./")
#read_hydrus.SavconToPointVtkMaster(filename = "Geometry_4.in")   #, folder='./')

