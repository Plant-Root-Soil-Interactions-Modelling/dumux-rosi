from solverbase import SolverWrapper
from richards import RichardsWrapper
from mpi4py import MPI; comm = MPI.COMM_WORLD; size = comm.Get_size(); rank = comm.Get_rank()

import numpy as np


class RichardsFlatWrapper(RichardsWrapper):
    """ 
    get the outputs as flattened arrays for 1D models, without re-writing RichardsNoMPIWrapper
    """

    def __init__(self, base):
        super().__init__(base)
