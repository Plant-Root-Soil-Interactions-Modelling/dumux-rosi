#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
dic = GetSources();
keys = dic.keys();
rootsystemtestccpvd = FindSource(keys[0][0])

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=rootsystemtestccpvd)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1401, 860]

# get color transfer function/color map for 'p'
pLUT = GetColorTransferFunction('p')

# show data in view
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1)
# trace defaults for the display properties.
cellDatatoPointData1Display.ColorArrayName = ['POINTS', 'p']
cellDatatoPointData1Display.LookupTable = pLUT
cellDatatoPointData1Display.OSPRayScaleArray = 'p'
cellDatatoPointData1Display.OSPRayScaleFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.GlyphType = 'Arrow'

# reset view to fit data
renderView1.ResetCamera()

# hide data in view
Hide(myGridrootpvd, renderView1)

# show color bar/color legend
cellDatatoPointData1Display.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')

# create a new 'Tube'
tube1 = Tube(Input=cellDatatoPointData1)
tube1.Scalars = ['POINTS', 'p']
tube1.Vectors = [None, '1']
tube1.Radius = 0.0026372399926185607

# Properties modified on tube1
tube1.Scalars = ['POINTS', 'rootRadius']
tube1.Vectors = [None, '']
tube1.VaryRadius = 'By Absolute Scalar'

# show data in view
tube1Display = Show(tube1, renderView1)
# trace defaults for the display properties.
tube1Display.ColorArrayName = ['POINTS', 'p']
tube1Display.LookupTable = pLUT
tube1Display.OSPRayScaleArray = 'p'
tube1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tube1Display.GlyphType = 'Arrow'

# hide data in view
Hide(cellDatatoPointData1, renderView1)

# show color bar/color legend
tube1Display.SetScalarBarVisibility(renderView1, True)

# create a new 'Threshold'
threshold1 = Threshold(Input=tube1)
threshold1.Scalars = ['POINTS', 'p']
threshold1.ThresholdRange = [95000.0, 95000.0]

# Properties modified on threshold1
threshold1.Scalars = ['POINTS', 'rootAge']
threshold1.ThresholdRange = [1.0, 10000000000.0]

# show data in view
threshold1Display = Show(threshold1, renderView1)
# trace defaults for the display properties.
threshold1Display.ColorArrayName = [None, '']
threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold1Display.GlyphType = 'Arrow'

# hide data in view
Hide(tube1, renderView1)

# Properties modified on threshold1Display
threshold1Display.OSPRayScaleArray = 'TubeNormals'

# get animation scene
animationScene1 = GetAnimationScene()

animationScene1.GoToLast()

# reset view to fit data
renderView1.ResetCamera()

