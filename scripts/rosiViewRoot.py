#### import the simple module from the paraview
from paraview.simple import *

# find source
dic = GetSources();
keys = dic.keys();
rootsystemtestccpvd = FindSource(keys[0][0])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1416, 860]

# create a new 'Cell Data to Point Data' and 'Tube'
cellDatatoPointData1 = CellDatatoPointData(Input=rootsystemtestccpvd)
tube1 = Tube(Input=cellDatatoPointData1)
tube1.Scalars = ['POINTS', 'rootRadius']
tube1.Vectors = [None, '']
tube1.NumberofSides = 12
tube1.VaryRadius = 'By Absolute Scalar'

# show data in view
tube1Display = Show(tube1, renderView1)

# hide the original data in view
Hide(rootsystemtestccpvd, renderView1)

# show color bar/color legend
tube1Display.SetScalarBarVisibility(renderView1, True)

# change interaction mode for render view
renderView1.InteractionMode = '3D'

# reset view to fit data
renderView1.ResetCamera()
