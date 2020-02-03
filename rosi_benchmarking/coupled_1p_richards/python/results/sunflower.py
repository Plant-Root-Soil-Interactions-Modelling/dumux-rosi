from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
renderView1 = GetActiveViewOrCreate('RenderView')

obj1 = Plane()
obj1.Origin = [0, 0, -2.9]
obj1.Point1 = [1, 0, -2.9]
obj1.Point2 = [0, 1, -2.9]

obj1Display = Show(obj1,renderView1)
obj1Display.Opacity = 0.2
