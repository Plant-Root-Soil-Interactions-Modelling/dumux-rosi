import vtk
import numpy as np
import matplotlib.pyplot as plt

#
# vtk_tools (more tools in conversions module roco)
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
def read1D_vtp_data(name, pwidx = 2):
    pd = read_polydata(name)

    try:  # cell
        data = pd.GetCellData()
        a = data.GetArray(2).GetTuple(0)
        cell = True
    except:
        data = pd.GetPointData()
        cell = False

    nocd = data.GetNumberOfArrays()
    sw = data.GetArray(0)  # saturation (S_lig)
    pw = data.GetArray(pwidx)  # pressure (p_lig)
    noa = sw.GetNumberOfTuples()

    sw_ = np.ones(noa,)
    pw_ = np.ones(noa,)
    for i in range(0, noa):
        d = sw.GetTuple(i)
        sw_[i] = d[0]
        d = pw.GetTuple(i)
        pw_[i] = d[0]

    points = pd.GetPoints()

    if cell:
        Nc = pd.GetNumberOfCells()
        z_ = np.zeros((Nc,))
        for i in range(0, Nc):
            c = pd.GetCell(i)
            ids = c.GetPointIds()
            p1 = np.zeros(3,)
            points.GetPoint(ids.GetId(0), p1)
            p2 = np.zeros(3,)
            points.GetPoint(ids.GetId(1), p2)
            z_[i] = 0.5 * (p1[0] + p2[0])
        return sw_, pw_, z_
    else:  # not cell
        Np = pd.GetNumberOfPoints()
        z_ = np.zeros((Np,))
        for i in range(0, Np):
            p = np.zeros(3,)
            points.GetPoint(i, p)
            z_[i] = p[0]
        return sw_, pw_, z_


#
# returns the cell or vertex data (index 0 and 2 hard coded)) of vtp file
#
def read3D_vtp_data(name, axis = 2):
    pd = read_vtu(name)

    try:  # cell
        data = pd.GetCellData()
        a = data.GetArray(2).GetTuple(0)
        cell = True
    except:
        data = pd.GetPointData()
        cell = False

    nocd = data.GetNumberOfArrays()
#     print("Number of arrays", nocd)
#     for i in range(0, nocd):
#         print(data.GetArrayName(i))
    sw = data.GetArray(0)  # saturation
    pw = data.GetArray(2)  # pressure

    noa = sw.GetNumberOfTuples()
    # print("number of data points", noa)
    sw_ = np.ones(noa,)
    pw_ = np.ones(noa,)
    for i in range(0, noa):
        d = sw.GetTuple(i)
        sw_[i] = d[0]
        d = pw.GetTuple(i)
        pw_[i] = d[0]

    points = pd.GetPoints()
    if cell:
        Nc = pd.GetNumberOfCells()
        z_ = np.zeros((Nc,))
        for i in range(0, Nc):
            c = pd.GetCell(i)
            ids = c.GetPointIds()
            n = ids.GetNumberOfIds ()
            midz = 0.
            p1 = np.zeros(3,)
            for j in range(0, n):
                points.GetPoint(ids.GetId(j), p1)
                midz += p1[2] / n
            z_[i] = midz  # sould actually be the mean of the 8 cube points (but i was too lazy)
        return sw_, pw_, z_
    else:  # not cell
        Np = pd.GetNumberOfPoints()
        z_ = np.zeros((Np,))
        for i in range(0, Np):
            p = np.zeros(3,)
            points.GetPoint(i, p)
            z_[i] = p[2]
        return sw_, pw_, z_

    return sw_, pw_, z_


#
# returns the cell or vertex data (index 0 and 2 hard coded) of parallel vtp files
#
def read3Dp_vtp_data(prename, postname, n):
    z_ = np.ones(0,)
    sw_ = np.ones(0,)
    pw_ = np.ones(0,)
    for i in range(0, n):
        n_ = prename + "{:04d}-".format(i) + postname + ".vtu"
        print("opening: ", n_)
        sw, pw, z = read3D_vtp_data(n_)
        z_ = np.hstack((z_, z))
        sw_ = np.hstack((sw_, sw))
        pw_ = np.hstack((pw_, pw))

    return sw_, pw_, z_


#
# reads a dumux output style vtu
#
def read3D_vtp(name, np = 1, axis = 2):
    if np == 1:
        return read3D_vtp_data(name + ".vtu", axis)
    else:
        return read3Dp_vtp_data("s{:04d}-p".format(np), name, np)

