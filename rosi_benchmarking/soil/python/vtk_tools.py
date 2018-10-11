import vtk
import numpy as np
import matplotlib.pyplot as plt

#
# vtk_tools
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
# converts a vtp to dgf
#
def vtp2dgf(name):
    pd = read_polydata(name + ".vtp")  # read vtp

    file = open(name + ".dgf", "w")  # write dgf
    file.write("DGF\n")
    # vertex
    file.write('Vertex\n')
    Np = pd.GetNumberOfPoints()
    points = pd.GetPoints()
    pdata = pd.GetPointData()
    Npd = pdata.GetNumberOfArrays()
    file.write('parameters {:g}\n'.format(Npd))
    for i in range(0, Np):
        p = np.zeros(3,)
        points.GetPoint(i, p)
        file.write('{:g} {:g} {:g} '.format(p[0], p[1], p[2]))
        for j in range(0, Npd):  # write point data - todo lets pick ids
            pdataj = pdata.GetArray(j)
            d = pdataj.GetTuple(i)
            file.write('{:g} '.format(d[0]))
        file.write('\n')

    file.write('#\n');
    file.write('Simplex\n');
    # cells
    Nc = pd.GetNumberOfCells()
    cdata = pd.GetCellData()
    Ncd = cdata.GetNumberOfArrays()
    file.write('parameters 2\n'.format(Ncd))
    for i in range(0, Nc - 1):
        cpi = vtk.vtkIdList()
        pd.GetCellPoints(i, cpi)
        for j in range(0, cpi.GetNumberOfIds()):  # write cell ids
            file.write('{:g} '.format(cpi.GetId(j)))
        for j in range(0, Ncd):  # write cell data - todo lets pick ids
            cdataj = cdata.GetArray(j)
            d = cdataj.GetTuple(i)
            file.write('{:g} '.format(d[0]))
        file.write('\n')
    # i dont know how to use these in dumux
    file.write('#\n')
    file.write('BOUNDARYSEGMENTS\n')  # how do i get the boundary segments into DUMUX ?
    file.write('2 0\n')
    file.write('3 {:g}\n'.format(Np - 1))  # vertex id, but index starts with 0
    file.write('#\n');
    file.write('BOUNDARYDOMAIN\n')
    file.write('default 1\n');
    file.write('#\n')

    print("finished writing " + name + ".dgf")
    file.close()


#
# returns the cell or vertex data (index 0 and 2 hard coded)) of vtp file
#
def read1D_vtp_data(name, cell = True):
    polydata = read_polydata(name)

    if cell:
        data = polydata.GetCellData()
    else:
        data = polydata.GetPointData()

    nocd = data.GetNumberOfArrays()
#     print("Number of arrays", nocd)
#     for i in range(0,nocd):
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

    return sw_, pw_


#
# returns the cell or vertex data (index 0 and 2 hard coded)) of vtp file
#
def read3D_vtp_data(name, cell = True):
    polydata = read_polydata(name)

    if cell:
        data = polydata.GetCellData()
    else:
        data = polydata.GetPointData()

    nocd = data.GetNumberOfArrays()
#     print("Number of arrays", nocd)
#     for i in range(0,nocd):
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

    # vertex (makes actually only sense for vertex data)
    Np = pd.GetNumberOfPoints()
    z_ = np.ones(Np,)
    points = pd.GetPoints()
    for i in range(0, Np):
        p = np.zeros(3,)
        points.GetPoint(i, p)
        z_[i] = p[2]

    return sw_, pw_, z_

