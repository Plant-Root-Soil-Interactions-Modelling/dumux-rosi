''' Script to compute equivalent soil water potential using SUF from Meunier et al., 2017 '''

import matplotlib.pyplot as plt
from vtk_tools_SUF import *
import math
import numpy as np
#from scipy.linalg import block_diag
import glob
from scipy.sparse import block_diag

# go to the right place
name = 'soybean_Honly_2003'					# problem name
Hseq_t = []
Pseq_t = []

filelist = glob.iglob(r'../../../build-cmake/rosi_benchmarking/coupled_1p_richards/results_' + name + '/*.vtp')

for filepath in sorted(filelist):
    kr = read3D_vtp_kr(filepath)					# get radial conductivity from root vtp
    kx = read3D_vtp_kx(filepath)					# get axial conductivity from root vtp
    nodes = read3D_vtp_nodes(filepath)				# get RA connectivity from root vtp
    p = read3D_vtp_soilp(filepath)					# get soil voxel pressure from root vtp
    xylemP = read3D_vtp_xylemp(filepath)				# get xylem pressure from root vtp
    xylemP = (xylemP- 1.e5) * 100. / 1.e3 / 9.81			# convert from Pa-->cm
    
    diag_kx = np.diag(kx)						# create a diagonal matrix of kx
    diag_kr = np.diag(kr)						# create a diagonal matrix of kr
    diag_k = block_diag(diag_kx, diag_kr)				# create a big diagonal matrix of kx and kr
    
    IM = scipy.sparse(shape=(np.shape(diag_k)))			# create IM matrix with all zeros
    np.fill_diagonal(IM, 1)						# fill all diagonal elements with "1"
    for i in range(len(diag_kx)):
        IM[i][i+1] = -1						# fill IM[i][j+1] with "-1" 
    
    n,m = IM.shape 							# for generality
    T = np.ones((n,1))							# first column with collar connectivity (fill all 1 since it doesn't matter)
    IM_T = np.hstack((T,IM))						# add it to IM to create IMt
    C = np.matmul(IM, diag_k)						# multiply IM with diag(K)
    C = np.matmul(C, IM_T)						# multiply with IMt --> connectivity matrix

    C2 = C[0:len(diag_kx), 1:len(diag_kx)+1]				# create C2 matrix
    C3 = C[0:len(diag_kx), len(diag_kx) + 1:2*len(diag_kx)+1]	# create C3 matrix
    inv_C2 = np.linalg.inv(C2)					# get inverse of C2
    C2_C3 = np.matmul(inv_C2, C3)					# get [C2^-1 * C3]
    I = np.identity(len(C3))						# create identity matrix of size of C3
    I = I + C2_C3							# I + [C2^-1C3]
    C4 = np.matmul(diag_kr, I)					# get C4 matrix

    SUF = [] 								# define an empty array to store SUF
    total = np.sum(C4)							# sum of all elements of C4
    for i in range(len(C4)):						# sum of all elements in a row 'i' of C4
        x = sum(C4[i])
        x = x/total							# compute SUF of each element
        SUF.append(x)							# store SUF in a array
#    print(np.sum(SUF))
    for j in range(len(SUF)):
        Hseq = np.multiply(SUF, p)					# compute equivalent soil water potential of each segment
        Pseq = np.multiply(SUF, xylemP)					# compute equivalent xylem potential of each segment
    
    Hseq_ = np.sum(Hseq)						# total equivalent soil water potential
    Pseq_ = np.sum(Pseq)						# total equivalent xylem potential
    
    Hseq_t.append(Hseq_)
    Pseq_t.append(Pseq_)

#print(Hseq_t)
#print("\n")
#print(Pseq_t)

#    Figure
fig, ax1 = plt.subplots()
time = np.linspace(0,154,3697)						# (initial day, final day, number of vtps)
ax1.plot(time, Pseq_t, 'r-')						# time vs equivalent xylem pressure using SUF as a weighing factor
ax1.plot(time, Hseq_t, 'b-')    					# time vs equivalent soil pressure using SUF as a weighing factor 
ax1.set_xlabel("Time [days]")
ax1.set_ylabel("Water potential [cm]")
ax1.legend(["xylem", "soil"], loc = 'upper left')
fig.savefig("../../../build-cmake/rosi_benchmarking/coupled_1p_richards/results_" + name + "/" + name + "_equivalent_p.pdf", bbox_inches='tight')
plt.show()
