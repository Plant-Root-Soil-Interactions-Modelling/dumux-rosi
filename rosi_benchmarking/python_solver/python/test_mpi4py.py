import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/soil_richards/")
from richardsyaspsolver import *

import numpy as np
import os
import time

from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# # SCATTER EXAMPLE
# if rank == 0:
#     data = [(i+1)**2 for i in range(size)]
#     print(data)
#     print("data ", len(data))
# else:
#     data = None
#     data = comm.scatter(data, root=0)
# print("now ", data)
# assert data == (rank+1)**2

# GATHERING EXAMPLE
data = [rank]*rank
print("rank", rank, "data", data)
data = comm.gather(data, root=0)
if rank == 0:
    print("process 0 now ", len(data))
    for i in range(0,4):
        print("data ", data[i])
else:
    assert data is None


