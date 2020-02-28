import solver.plantbox as pb
from solver.vtk_tools import *

import numpy as np
import vtk

""" 
VTK Tools, by Daniel Leitner (refurbished 12/2019) 

for vtk to numpy, and numpy to vtk conversions
reading: vtp, writing: msh, dgf, vtp, rsml
"""


def segs_to_polydata(rs, zoom_factor = 10., param_names = ["radius", "type", "creationTime"]):
    """ Creates vtkPolydata from a RootSystem or Plant using segments 
    @param rs             A RootSystem, Plant, or SegmentAnalyser
    @param zoom_factor    The radial zoom factor, since root are sometimes too thin for vizualisation
    @param param_names    Parameter names of scalar fields, that are copied to the polydata
    """
    if isinstance(rs, pb.Organism):
        ana = pb.SegmentAnalyser(rs)  # for Organism like Plant or RootSystem
    else:
        ana = rs
    nodes = np_convert(ana.nodes)
    segs = np_convert(ana.segments)
    points = vtk_points(nodes)
    cells = vtk_cells(segs)
    pd = vtk.vtkPolyData()
    pd.SetPoints(points)
    pd.SetLines(cells)  # check SetPolys
    for n in param_names:
        param = np.array(ana.getParameter(n))
        if param.shape[0] == segs.shape[0]:
            if n == "radius":
                param *= zoom_factor
            data = vtk_data(param)
            data.SetName(n)
            pd.GetCellData().AddArray(data)
        else:
            print("segs_to_polydata: Warning parameter " + n + " is sikpped because of wrong size", param.shape[0], "instead of", segs.shape[0])

    c2p = vtk.vtkCellDataToPointData()
    c2p.SetPassCellData(True)
    c2p.SetInputData(pd)
    c2p.Update()
    return c2p.GetPolyDataOutput()


def render_window_(actor, windowName = "", scalarBar = None):
    """ puts a vtk actor on the stage (an interactive window) @name is the window titel 
    @param actor         the (single) actor
    @param windowName    optional
    @param scalarBar     an optional vtkScalarBarActor
    """
    colors = vtk.vtkNamedColors()  # Set the background color
    bkg = map(lambda x: x / 255.0, [26, 51, 102, 255])
    colors.SetColor("BkgColor", *bkg)
    ren = vtk.vtkRenderer()  # Set up window with interaction
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    # iren.SetInteractorStyle(vtk.vtkInteractorStyleUnicam())  # <- better than default, but maybe we find a better one
    iren.SetRenderWindow(renWin)
    ren.AddActor(actor)  # Add the actors to the renderer, set the background and size
    if scalarBar is not None:
        ren.AddActor2D(scalarBar)
    ren.SetBackground(colors.GetColor3d("BkgColor"))
    renWin.SetSize(1000, 1000)
    renWin.SetWindowName(windowName)
    iren.Initialize()  # This allows the interactor to initalize itself. It has to be called before an event loop.
    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(1.5)
    renWin.Render()
    iren.Start()  # Start the event loop.


def plot_3surf(pd):
    pass


def plot_roots(pd, pname = "type"):
    """ renders the root system in an interactive window """

    pd.GetPointData().SetActiveScalars("radius")
    tubeFilter = vtk.vtkTubeFilter()
    tubeFilter.SetInputData(pd)
    tubeFilter.SetNumberOfSides(9)
    tubeFilter.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    tubeFilter.Update()

#     pointData = pd.GetPointData()
#     pointData2 = vtk.vtkPointData()
#     pointData2.DeepCopy(pd.GetPointData())

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tubeFilter.GetOutputPort())
    mapper.ScalarVisibilityOn();
    mapper.SetScalarModeToUseCellFieldData()  # Cell is not working
    mapper.Update()

#     mapper2 = vtk.vtkPolyDataMapper()
#     mapper2.SetInputConnection(mapper.GetOutputPort())
#     mapper2.Update()

    plantActor = vtk.vtkActor()
    mapper.SetArrayName(pname)
    mapper.SelectColorArray(pname)
    # scalarRange = pd.GetScalarRange()
    # pd.GetPointData().SetActiveScalars(pname)
    mapper.UseLookupTableScalarRangeOn()

    plantActor.SetMapper(mapper)
    plantActor.RotateX(-90.0)

#     # Make the lookup table.
#     lut.SetTableRange(scalarRange) = vtk.vtkColorSeries()
#     # Select a color scheme.
#     # colorSeriesEnum = colorSeries.BREWER_DIVERGING_BROWN_BLUE_GREEN_9
#     # colorSeriesEnum = colorSeries.BREWER_DIVERGING_SPECTRAL_10
#     # colorSeriesEnum = colorSeries.BREWER_DIVERGING_SPECTRAL_3
#     # colorSeriesEnum = colorSeries.BREWER_DIVERGING_PURPLE_ORANGE_9
#     # colorSeriesEnum = colorSeries.BREWER_SEQUENTIAL_BLUE_PURPLE_9
#     # colorSeriesEnum = colorSeries.BREWER_SEQUENTIAL_BLUE_GREEN_9
#     colorSeriesEnum = colorSeries.BREWER_QUALITATIVE_SET3
#     # colorSeriesEnum = colorSeries.CITRUS
#     colorSeries.SetColorScheme(colorSeriesEnum)
#     lut = vtk.vtkLookupTable()
#     lut.SetNumberOfTableValues(16)
#     colorSeries.BuildLookupTable(lut)
#     # lut.SetNanColor(1, 0, 0, 1)
# #     lut.SetTableRange([0, 1])

    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(16)
    lut.SetHueRange(0.0, 1.0)
    lut.Build()

    mapper.SetLookupTable(lut)
    scalarBar = vtk.vtkScalarBarActor()
    scalarBar.SetLookupTable(lut)
    scalarBar.SetTitle(pname)
    lut.SetTableRange(pd.GetPointData().GetScalars(pname).GetRange())  # DataSet. GetScalarRange
#    textProperty = vtk.vtkTextProperty()
#    scalarBar.SetTitleTextProperty(textProperty)
#    scalarBar.SetLabelTextProperty(textProperty)
#    scalarBar.SetAnnotationTextProperty(textProperty)
#    scalarBar = None

    render_window_(plantActor, pname, scalarBar)

