#!/home/h.mai/share/tools/anaconda2/bin/python2
# -*- coding: utf-8 -*-

# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
from struct import unpack,pack
import sys
import numpy as np
from os import path
import h5py

class HF5Data(object):
    def __init__(self,fileInPath,fileOutPath,*endian):
        self.fileIn = fileInPath
        self.fileOut = fileOutPath

    def scanFileForInfosAndPositions(self):
        with h5py.File(self.fileIn, 'r') as hf:
            self.timeSteps=list(hf.get('TimeStep'))
            self.varnames=list(hf.get('VarName'))
            self.vartype=np.array(hf.get('VarType'))
            self.varcomp=np.array(hf.get('nComp'))
            self.geovar=[]
            t0=hf.get(self.timeSteps[0])
            self.totalVar= np.array(t0.keys())
            for name in self.totalVar:
                if name not in self.varnames:
                    self.geovar.append(name)
            i = 0
            while (self.vartype[i]>1):
                i=i+1
            print(self.varnames[i],self.varcomp[i])
            print(t0.get(self.varnames[i]),self.varcomp[i])
            name_= self.varnames[i]
            name_=name_[0:name_.rfind("[")]
            #print (t0.get(name_))
            nPoints=int(len(t0.get(self.varnames[i]))/self.varcomp[i])
            #nPoints=int(len(t0.get(name_))/self.varcomp[i])
            self.npoin=nPoints
            i = 0
            while (self.vartype[i]<1):
                i=i+1
            nCells=int(len(t0.get(self.varnames[i]))/self.varcomp[i])
            #nCells=int(len(t0.get(name_))/self.varcomp[i])
            self.elem=nCells

    def openFileOut(self):
        self.fileOutOpen=h5py.File(self.fileOut, 'w')

    def readVarAtTime(self,selectTimeStepNum,selectedVarNum):
        with h5py.File(self.fileIn, 'r') as hf:
            t0=hf.get(self.timeSteps[selectTimeStepNum])
            arr=np.array(t0.get(self.varnames[selectedVarNum]))
        return arr
