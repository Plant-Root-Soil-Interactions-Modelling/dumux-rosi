from rosi.solverbase import SolverWrapper
from rosi.richards import RichardsWrapper
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
max_rank = comm.Get_size()


class RichardsFlatWrapper(RichardsWrapper):
    """ 
    get the outputs as flattened arrays for 1D models, without re-writing RichardsNoMPIWrapper
    """

    def __init__(self, base):
        super().__init__(base)
