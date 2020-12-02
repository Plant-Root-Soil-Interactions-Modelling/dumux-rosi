''' Script to compute equivalent soil water potential using SUF from Meunier et al., 2017 '''
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

import matplotlib.pyplot as plt
from vtk_tools_SUF import *
import math
import numpy as np
import glob
import scipy.sparse as sp
from scipy.sparse.linalg import inv
from scipy.sparse.linalg import spsolve
import warnings
from scipy.sparse import (spdiags, SparseEfficiencyWarning, csc_matrix,
    csr_matrix, isspmatrix, dok_matrix, lil_matrix, bsr_matrix)
warnings.simplefilter('ignore',SparseEfficiencyWarning)

# go to the right place
name = 'soybean_Honly_2003'							# problem name
Hseq_t = []
Pseq_t = []

filelist = glob.iglob(r'../../../build-cmake/cpp/coupled_1p_richards/results_' + name + '/*.vtp')

for filepath in sorted(filelist):
    kr,_ = read3D_vtp_data(filepath, index_kr)                        # get radial conductivity from root vtp
    kx, nodes = read3D_vtp_data(filepath, index_kx)                        # get axial conductivity from root vtp
    p = read3D_vtp_data(filepath, index_p)                    # get soil voxel pressure from root vtp
    xylemP = read3D_vtp_data(filepath, index_p)	# get xylem pressure from root vtp	
    xylemP = (xylemP- 1.e5) * 100. / 1.e3 / 9.81				# convert from Pa-->cm
    
    diag_kx = sp.diags(-kx)						# create a sparse diagonal matrix of kx
    diag_kr = sp.diags(-kr)						# create a sparse diagonal matrix of kr
    diag_k = sp.block_diag((diag_kx, diag_kr))				# create a sparse diagonal matrix of kx and kr
    
    IM = sp.eye(2*kr.size)		# create IM matrix with all diagonal "1"
    IM = sp.lil_matrix(IM)
	
	for i in range(diag_kx.size):
		IM[i, i+1] = -1						# fill IM[i][j+1] with "-1" 
	
	n, m = IM.shape 							# for generality
	T = np.ones((n,1))							# first column with collar connectivity (fill all 1 since it doesn't matter)
	IM_T = sp.hstack((T,IM))						# add it to IM to create IMt
	
	C_ = sp.csr_matrix(IM)*(sp.csr_matrix(diag_k))			# multiply IM with diag(K)
	C = sp.csr_matrix(C_)*(sp.csr_matrix(IM_T))				# multiply with IMt --> connectivity matrix
	C2 = C[0:kx.size, 1:kx.size+1]					# create C2 matrix
	C3 = C[0:kx.size, kx.size+1:2*kx.size+1]				# create C3 matrix
	
	I = sp.identity(kx.size)						# create identity matrix of size of C3
	x = spsolve(C2, C3)							# compute C2^-1C3
	
	C4 = sp.csr_matrix(diag_kr)*(sp.csr_matrix(I+x))			# get C4 matrix
	SUF = [] 								# define an empty array to store SUF
	Krs = C4.sum()								# sum of all elements of C4
	Krs_t.append(Krs)							# Krs (global root conductivity)
	
	for i in range(kx.size):						# sum of all elements in a row 'i' of C4
		x = C4[i, :].sum()
		x = x/Krs							# compute SUF of each element
		SUF.append(x)							# store SUF in a array
		
	Hseq = np.multiply(SUF, p)						# compute equivalent soil water potential of each segment
	Pseq = np.multiply(SUF, xylemP)					# compute equivalent xylem potential of each segment
	Hseq_t.append(np.sum(Hseq))
	Pseq_t.append(np.sum(Pseq)) 

#    Figure
fig, ax1 = plt.subplots()
time = np.linspace(0,154,3697)						# (initial day, final day, number of vtps)
ax1.plot(time, Pseq_t, 'r-')							# time vs equivalent xylem pressure using SUF as a weighing factor
ax1.plot(time, Hseq_t, 'b-')    						# time vs equivalent soil pressure using SUF as a weighing factor 
ax1.set_xlabel("Time [days]")
ax1.set_ylabel("Water potential [cm]")
ax1.legend(["xylem", "soil"], loc = 'upper left')
fig.savefig("../../../build-cmake/rosi_benchmarking/coupled_1p_richards/results_" + name + "/" + name + "_equivalent_p.pdf", bbox_inches='tight')
plt.show()
