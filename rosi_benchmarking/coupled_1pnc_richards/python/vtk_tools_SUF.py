import vtk
import numpy as np
import matplotlib.pyplot as plt

# opens a vtp and returns the vtk polydata class
#
def read_polydata(name):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(name)
    reader.Update()
    polydata = reader.GetOutput()
    return polydata

# returns the cell or vertex data (index 0) of vtp file
#
def read3D_vtp_soilp(name):
    polydata = read_polydata(name)

    try:  # Box
        data = polydata.GetPointData()
        nocd = data.GetNumberOfArrays()
        p = data.GetArray(8)  # pressure of soil voxel
        noa = p.GetNumberOfTuples()
        cc = False
    except:  # CCTpfa
        data = polydata.GetCellData()
        nocd = data.GetNumberOfArrays()
        p = data.GetArray(8)  # pressure of soil voxel
        noa = p.GetNumberOfTuples()
        cc = True

    p_ = np.ones(noa,)
    for i in range(0, noa):
        d = p.GetTuple(i)
        p_[i] = d[0]

    return p_
    
def read3D_vtp_xylemp(name):
    polydata = read_polydata(name)

    try:  # Box
        data = polydata.GetPointData()
        nocd = data.GetNumberOfArrays()
        p = data.GetArray(0)  # xylem pressure
        noa = p.GetNumberOfTuples()
        cc = False
    except:  # CCTpfa
        data = polydata.GetCellData()
        nocd = data.GetNumberOfArrays()
        p = data.GetArray(0)  # xylem pressure
        noa = p.GetNumberOfTuples()
        cc = True
    p_ = np.ones(noa,)
    for i in range(0, noa):
        d = p.GetTuple(i)
        p_[i] = d[0]

    return p_
        
def read3D_vtp_kx(name):
    polydata = read_polydata(name)

    try:  # Box
        data = polydata.GetPointData()
        nocd = data.GetNumberOfArrays()
        kx = data.GetArray(18)  # axial conductivity of segment
        noa = kx.GetNumberOfTuples()
        cc = False
    except:  # CCTpfa
        data = polydata.GetCellData()
        nocd = data.GetNumberOfArrays()
        kx = data.GetArray(18)  # axial conductivity of segment
        noa = kx.GetNumberOfTuples()
        cc = True

    kx_ = np.ones(noa,)
    for i in range(0, noa):
        d = kx.GetTuple(i)
        kx_[i] = d[0]

    return kx_
        
def read3D_vtp_kr(name):
    polydata = read_polydata(name)

    try:  # Box
        data = polydata.GetPointData()
        nocd = data.GetNumberOfArrays()
        kr = data.GetArray(17)  # radial conductivity of segment
        noa = kr.GetNumberOfTuples()
        cc = False
    except:  # CCTpfa
        data = polydata.GetCellData()
        nocd = data.GetNumberOfArrays()
        kr = data.GetArray(17)  # radial conductivity of segment
        noa = kr.GetNumberOfTuples()
        cc = True

    kr_ = np.ones(noa,)
    for i in range(0, noa):
        d = kr.GetTuple(i)
        kr_[i] = d[0]

    return kr_

def read3D_vtp_nodes(name):
    polydata = read_polydata(name)
    """ The cells of vtkPolyData as numpy array  """
    Nc = polydata.GetNumberOfCells()
    d = polydata.GetCell(0).GetPointIds().GetNumberOfIds()
    z_ = np.zeros((Nc, d))
    for i in range(0, Nc):
        p = np.zeros(d,)
        ids = polydata.GetCell(i).GetPointIds()
        for j in range(0, d):
            p[j] = ids.GetId(j)
        z_[i, :] = p
    return z_
        
