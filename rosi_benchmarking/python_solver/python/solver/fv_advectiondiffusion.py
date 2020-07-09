import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
from scipy.sparse import identity as I
import scipy.linalg as la

from solver.fv_grid import *
from solver.fv_solver import *
import solver.van_genuchten as vg


class FVAdvectionDiffusion(FVSolver):
    """ 
    AdvectionDiffusion solver (first draft)
    """

    def __init__(self, grid :FVGrid):
        super().__init__(grid)
        self.mid = self.grid.centers()  # precompute cell centers
        #
        n = self.n
        self.b = np.ones((n,)) * 1.  # [1] buffer power
        self.D = np.ones((n,)) * 0.864  #  [cm2/day] 1.e-5 cm2/s = 0.864 cm2 / day
        self.u = np.zeros((n, self.grid.dim))  # [cm/day] velocity field
        # TODO make local, save only A
        self.beta = np.zeros((n,))  # const part of diagonal entries [1/day]
        self.alpha = np.zeros(grid.neighbours.shape)  # secondary diagonals [1/day]
        i_ = np.array(range(0, n), dtype = np.int64)
        cols = np.ones((1, self.grid.number_of_neighbours), dtype = np.int64)
        self.alpha_i = np.outer(i_, cols)
        self.alpha_j = self.grid.neighbours.copy()
        for (i, j) in self.grid.boundary_faces:
            self.alpha_j[i, j] = 0  # corresponding alpha = 0 at these boundary cells
            self.alpha_j[i, j] = 0
        self.A = None

    def solver_initialize(self):
        """ call back function for initialization"""
        self.sim_time = 0.  # current simulation time [day]
        self.prepare_linear_system()

    def solver_proceed(self, dt):
        """ call back function for each time step """
        pass

    def solve_single_step(self, dt, verbose):
        """ solve for time step dt, called by solve with step size control 
        @param dt     time step [day]
        @param verbose            tell me more        
        """
        self.bc_to_source(dt)
        b = self.x0 + self.sources
        b[0] = max(b[0], 0.)  # TODO
        c2 = self.solve_crank_nicolson(b, dt)
        return (c2, True, 0)

    def prepare_linear_system(self):
        """ experimental ... TODO """
        self.alpha = np.zeros(self.grid.neighbours.shape)
        self.beta = np.zeros((self.n,))

        for i in range(0, self.n):  # over cells
            for  j in range(0, self.grid.number_of_neighbours):
                neighbour = self.grid.neighbours[i, j]
                if neighbour >= 0:  # -1 means there is no neighoour
                    n = self.mid[neighbour] - self.mid[i]
                    n = n / np.linalg.norm(n)
                    self.beta[i] += (1. / self.b[i]) * self.grid.area_per_volume[i, j] * (-self.D[i] / self.grid.dx[i, j] - np.inner(self.u[i], n))
                    self.alpha[i, j] += (1. / self.b[i]) * self.grid.area_per_volume[i, j] * (self.D[i] / self.grid.dx[i, j])

        A = sparse.coo_matrix((self.alpha.flat, (self.alpha_i.flat, self.alpha_j.flat)))
        B = sparse.coo_matrix((self.beta, (np.array(range(0, self.n)), np.array(range(0, self.n)))))
        self.A = A + B

    def solve_forward_euler(self, c, dt):
        """ explicit euler (matrix vector multiplication) """
        return (I(self.n) + self.A * dt) * c

    def solve_backward_euler(self, c, dt):
        """ implicit euler (solve a linear system) """
        return LA.spsolve(I(self.n) - self.A * dt, c, use_umfpack = True)

    def solve_crank_nicolson(self, c, dt):
        """ arithmetic mean of forward and backward eulers """
        return 0.5 * (self.solve_forward_euler(c, dt) + self.solve_backward_euler(c, dt))

    def bc_concentration(self, i, c, dx):
        """ flux boundary condition [g / cm2 / day] """
        q = (self.D[i]) * (c - self.x0[i]) / dx
        # print(q, "g/cm2/day")
        return q

    def bc_to_source(self, dt):
        """
        Computes a cell source term from a boundary conditions given in self.bc,
        self.bc is dictionary with keys (cell, face_id) and values is a boundary function returning a value [g / cm2 / day].

        Only one face per cell is allowed, overwrites any other sources

        dictionary with lambda functions would be nicer, but causes trouble with parallelisation
        """
        for (i, j), (type, v) in self.bc.items():
            if type == "concentration":
                bc = self.bc_concentration(i, *v[0:2])
            else:
                raise("Unkown boundary condition")
            self.sources[i] = dt * bc * self.grid.area_per_volume[i, j]  # [g / cm3]


class FVAdvectionDiffusion_richards(FVAdvectionDiffusion):
    """ 
    Links FVAdvectionDiffusion to a FVRichards model, 
    for effective Diffusivity and Darcy velocity 
    """

    def __init__(self, grid :FVGrid, richards):
        super().__init__(grid)
        self.richards = richards
        self.D0 = self.D

    def solver_initialize(self):
        """ call back function for initialization"""
        super().solver_initialize()
        self.solver_proceed(0.)

    def solver_proceed(self, dt):
        """ retrieve effective diffusion and velocity field from richards """
        self.u = self.richards.darcy_velocity()
        theta = vg.water_content(self.richards.x0, self.richards.soil)
        self.D = self.D0 * theta * 0.5

