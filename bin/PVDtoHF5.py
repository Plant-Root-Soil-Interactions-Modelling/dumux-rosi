#!/home/h.mai/share/tools/anaconda2/bin/python2
""" Convert simulation results from PVD to H5
"""
import xml.etree.ElementTree as ET
import numpy as np
import h5py
import sys

def PVDtoHF5(PVD_filename):
    tree = ET.parse(PVD_filename)
    root = tree.getroot()

    TimeStep=[]
    filename=[]
    for DataSet in root.iter('DataSet'):
        TimeStep.append(DataSet.get('timestep'))
        fname=DataSet.get('file')
        fname=fname[fname.rindex("/")+1:len(fname)]
        filename.append(fname)

    ###############################################################################
    treevtp = ET.parse(filename[0])
    rootvtp = treevtp.getroot()
    ### Read Cell info
    CellVar=[]
    CellVarComponents=[]
    for Cells in rootvtp.iter('CellData'):
        for DataArray in Cells.iter('DataArray'):
            CellVar.append(DataArray.get('Name'))
            nComp=int(DataArray.get('NumberOfComponents'))
            CellVarComponents.append(nComp)
    ### Read Points info
    PointVar=[]
    PointVarComponents=[]
    for Points in rootvtp.iter('Points'):
        for DataArray in Points.iter('DataArray'):
            PointVar.append(DataArray.get('Name'))
            nComp=int(DataArray.get('NumberOfComponents'))
            PointVarComponents.append(nComp)
    VarType=[1]*len(CellVar)+[0]*len(PointVar)
    ###############################################################################
    with h5py.File(PVD_filename.replace(".pvd",".h5"), 'w') as hf:
        hf.create_dataset('TimeStep', data=TimeStep)
        hf.create_dataset('VarName', data=CellVar+PointVar)
        hf.create_dataset('VarType', data=VarType)
        hf.create_dataset('nComp', data=CellVarComponents+PointVarComponents)
        for i in range(len(filename)):
            sys.stdout.write ("\rReading file "+filename[i]+' .....')
            g = hf.create_group(TimeStep[i])
            treevtp = ET.parse(filename[i])
            rootvtp = treevtp.getroot()
            for DataArray in rootvtp.iter('DataArray'):
                value = np.asarray([float(s) for s in filter(None, [w.replace('\n', '') for w in DataArray.text.split (" ")])])
                g.create_dataset(DataArray.get('Name'), data=value)
if __name__=='__main__':
    sys.exit(PVDtoHF5(sys.argv[1]))
