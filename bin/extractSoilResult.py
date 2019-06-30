#!/home/h.mai/share/tools/anaconda2/bin/python2
""" edit PVD
"""
import xml.etree.ElementTree as ET
import numpy as np
import h5py
import sys

def editPVD(PVD_filename):
    tree = ET.parse(PVD_filename)
    PVDfileOut = "ex_"+PVD_filename
    root = tree.getroot()
    TimeStep=[]
    filename=[]
    for DataSet in root.iter('DataSet'):
        TimeStep.append(DataSet.get('timestep'))
        fname=DataSet.get('file')
        fname=fname[fname.rindex("/")+1:len(fname)]
        filename.append(fname)

    In_file = open(PVD_filename,"r")
    contents = In_file.read()
    In_file.close()
    replaced_contents = contents.replace(PVD_filename[0:PVD_filename.rindex(".")],\
                                            PVDfileOut[0:PVDfileOut.rindex(".")])
    Out_file=open(PVDfileOut,"w")
    Out_file.write(replaced_contents)
    Out_file.close()
    ###############################################################################
    tree = ET.parse(filename[0])
    root = tree.getroot()

    ### Read Celss info
    CellVar=[]
    for Cells in root.iter('CellData'):
        for DataArray in Cells.iter('DataArray'):
            CellVar.append(DataArray.get('Name'))
    ### Read Points info
    PointVar=[]
    for Points in root.iter('Points'):
        for DataArray in Points.iter('DataArray'):
            PointVar.append(DataArray.get('Name'))

    VarName=CellVar+PointVar
    #print VarName
    ###############################################################################

    selectedVarIdx = [1,3,19,20]
    print (u"These Variables have been chosen:")
    for i in range(len(selectedVarIdx)):
        print ("                                ", VarName[selectedVarIdx[i]],"    (",selectedVarIdx[i],")")
    print ("")
    ###############################################################################
    for i in xrange(len(filename)):
        tree = ET.parse(filename[i])
        root = tree.getroot()
        VarObject=[]

        ### Read Varible Obj
        for Cells in root.iter('CellData'):
            for DataArray in Cells.iter('DataArray'):
                VarObject.append(DataArray)
        for Points in root.iter('Points'):
            for DataArray in Points.iter('DataArray'):
                VarObject.append(DataArray)

        ### make a list of removed Var
        removeList=[]
        for j in xrange(len(VarName)):
            if j not in selectedVarIdx:
                removeList.append(VarObject[j])
        parent_map = {c:p for p in root.iter() for c in p}

        ### Removing vars
        for Variable in removeList:
            parent_map[Variable].remove(Variable)

        ### creat ouput filename
        Outfilename=filename[i].replace(PVD_filename[0:PVD_filename.rindex(".")]\
                                            ,PVDfileOut[0:PVDfileOut.rindex(".")])
        sys.stdout.write ("\rCreating output file "+Outfilename+' .....')
        #print ("Creating output file "+Outfilename+' .....', end='\r', flush=True)
        tree.write(Outfilename)
    print "DONE !!!!"
if __name__=='__main__':
    sys.exit(editPVD(sys.argv[1]))
