import vtk
import numpy as np
import matplotlib.pyplot as plt



def read1D_vtp_data(name, cell = True):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(name)
    reader.Update()    
    polydata = reader.GetOutput() 
    
    if cell: 
        data = polydata.GetCellData()
    else:
        data = polydata.GetPointData()
        
    nocd = data.GetNumberOfArrays()
#     print("Number of arrays", nocd)
#     for i in range(0,nocd):
#         print(data.GetArrayName(i))
 
    sw = data.GetArray(0) # saturation   
    pw = data.GetArray(2) # pressure

    noa = sw.GetNumberOfTuples()
    # print("number of data points", noa)
    
    sw_ = np.ones(noa,)
    pw_ = np.ones(noa,)
    for i in range(0,noa):    
        d = sw.GetTuple(i)
        sw_[i] = d[0]
        d = pw.GetTuple(i)
        pw_[i] = d[0]
            
    return sw_, pw_
