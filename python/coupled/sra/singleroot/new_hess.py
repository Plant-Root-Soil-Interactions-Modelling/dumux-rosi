import numpy as np
import scipy
from scipy import sparse

""" 
    Andreas example 
    
    (pre-study for building doussan matrix)
"""
segs = [ [0, 1], [1, 2], [1, 3], [3, 4] ]
vn = len(segs)  # number of edges
nn = 5  # number of nodes

ii_, jj_, vv_ = [], [], []
for i, s in enumerate(segs):  # build incidence matrix from edges
    ii_.append(i)
    ii_.append(i)
    jj_.append(segs[i][0])
    jj_.append(segs[i][1])
    vv_.append(-1.)
    vv_.append(1.)
IM = sparse.coo_matrix((np.array(vv_), (np.array(ii_), np.array(jj_))), shape = (vn, nn))
IMt = IM.transpose()
Kx = sparse.diags(0.1 * np.ones((vn,)))  # e.g. const kx = 0.1
L = IMt @ Kx @ IM

# alternatively use sparse.laplace (!!!) segs -> adjacency matrix -> L

print("IM")
print(IM.todense())

print("\nIM tranposed")
print(IMt.todense())

print("\nKx")
print(Kx.todense())

print("\nL")
print(L.todense())

""" Dirichlet """
diag = sparse.eye(L.shape[0]).tolil()
diag[0, 0] = 0
Lbc = diag.dot(L)
Lbc[0, 0] = 1  # found no readable cheap way to do it
print("\nLbc")
print(Lbc.todense())

Lm = L[1:, 1:]
print("\nLm")
print(Lm.todense())

