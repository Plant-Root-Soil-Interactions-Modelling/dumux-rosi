import vtk
import numpy as np
import matplotlib.pyplot as plt

#
# vtk_tools (more tools in converstions module roco)
#
# D. Leitner, 2018
#


#
# opens a vtp and returns the vtk polydata class
#
def read_polydata(name):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(name)
    reader.Update()
    polydata = reader.GetOutput()
    return polydata


#
# opens a vtu and returns the vtk polydata class
#
def read_vtu(name):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(name)
    reader.Update()
    polydata = reader.GetOutput()
    return polydata


#
# returns the cell or vertex data (index 0 and 2 hard coded) of vtp file
#
def read1D_vtp_data(name, cell = True):
    pd = read_polydata(name)
    if cell:
        data = pd.GetCellData()
    else:
        data = pd.GetPointData()

    nocd = data.GetNumberOfArrays()
#     print("Number of arrays", nocd)
#     for i in range(0, nocd):
#         print(data.GetArrayName(i))

    sw = data.GetArray(0)  # saturation (S_lig)
    pw = data.GetArray(2)  # pressure (p_lig)

    noa = sw.GetNumberOfTuples()
    sw_ = np.ones(noa,)
    pw_ = np.ones(noa,)
    for i in range(0, noa):
        d = sw.GetTuple(i)
        sw_[i] = d[0]
        d = pw.GetTuple(i)
        pw_[i] = d[0]

    # vertex (makes actually only sense for vertex data)
    Np = pd.GetNumberOfPoints()
    z_ = np.ones(Np,)
    points = pd.GetPoints()
    for i in range(0, Np):
        p = np.zeros(3,)
        points.GetPoint(i, p)
        z_[i] = p[0]

    return sw_, pw_, z_


#
# returns the cell or vertex data (index 0 and 2 hard coded)) of vtp file
#
def read3D_vtp_data(name, cell = True):
    pd = read_vtu(name)
    if cell:
        data = pd.GetCellData()  # vertex (makes actually only sense for vertex data)
    else:
        data = pd.GetPointData()

    nocd = data.GetNumberOfArrays()
    print("Number of arrays", nocd)
    for i in range(0, nocd):
        print(data.GetArrayName(i))
    sw = data.GetArray(10)  # saturation
    pw = data.GetArray(9)  # pressure

    noa = sw.GetNumberOfTuples()
    # print("number of data points", noa)
    sw_ = np.ones(noa,)
    pw_ = np.ones(noa,)
    for i in range(0, noa):
        d = sw.GetTuple(i)
        sw_[i] = d[0]
        d = pw.GetTuple(i)
        pw_[i] = d[0]

    # vertex (makes actually only sense for vertex data)
    Np = pd.GetNumberOfPoints()
    z_ = np.ones(Np,)
    points = pd.GetPoints()
    for i in range(0, Np):
        p = np.zeros(3,)
        points.GetPoint(i, p)
        z_[i] = p[2]

    return sw_, pw_, z_

