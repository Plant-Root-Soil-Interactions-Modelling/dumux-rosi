import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
from scipy.sparse import identity as I
import scipy.linalg as la

from fv_grid import *
from fv_solver import *
import van_genuchten as vg


class FVAdvectionDiffusion(FVSolver):
    """ 
    AdvectionDiffusion solver (first draft)
    """

    def __init__(self, grid :FVGrid):
        super().__init__(grid)
        mid = self.grid.centers()  # precompute cell face normals
        cols = np.ones((1, self.grid.number_of_neighbours), dtype = np.int64)
        self.fn = mid[self.grid.neighbours] - mid[np.outer(np.arange(0, self.n), cols)]
        for i in range(0, self.grid.dim):
            self.fn[:, :, i] = self.fn[:, :, i] / np.linalg.norm(self.fn, axis = 2)  # normalize
        for (i, j) in self.grid.boundary_faces:
            self.fn[i, j, :] = 0  # set face boundary normals to zero
        self.b = np.ones((self.n,)) * 1.  # [1] buffer power
        self.D = np.ones((self.n,)) * 0.864  #  [cm2/day] effective diffusion, 1.e-5 cm2/s = 0.864 cm2 / day
        self.u = np.zeros((self.n, self.grid.dim))  # [cm/day] velocity field
        self.alpha = None
        self.beta = None
        self.A = None

    def solver_initialize(self):
        """ call back function for initialization"""
        self.sim_time = 0.  # current simulation time [day]
        self.prepare_linear_system()

    def solver_proceed(self, x, dt):
        """ call back function for each time step """
        pass

    def solve_single_step(self, dt, verbose):
        """ solve for time step dt, called by solve with step size control 
        @param dt                 time step [day]
        @param verbose            tell me more        
        """
        self.bc_to_source(dt)
        b = self.x0 + self.sources
        b[0] = max(b[0], 0.)  # TODO
        c2 = self.solve_backward_euler(b, dt)
        # c2 = self.solve_crank_nicolson(b, dt)  # <--- choose temporal discretisation
        return (c2, True, 0)

    def create_alpha_beta(self):
        """ creates diagonal (beta) and seconary diagonals (alpha) """
        cols = np.ones((1, self.grid.number_of_neighbours), dtype = np.int64)  # [1,1]
        beta_inv = np.outer(np.divide(np.ones(self.b.shape), self.b), cols)
        gamma = np.multiply(beta_inv, self.grid.area_per_volume)

#         D = 0.5 * (self.D[np.outer(np.arange(0, self.n), cols)] + self.D[self.grid.neighbours])  # boundary normals are zero TODO change u0 to surface fluxes?
#         for (i, j) in self.grid.boundary_faces:
#             D[i, j] = self.D[i]
        D = self.D[np.outer(np.arange(0, self.n), cols)]
        ddx = np.divide(D, self.grid.dx)

        beta_ = np.multiply(gamma, -ddx)
        for (i, j) in self.grid.boundary_faces:
            beta_[i, j] = 0  # zero flux over boundaries

        u = 0.5 * (self.u[np.outer(np.arange(0, self.n), cols)] + self.u[self.grid.neighbours])  # boundary normals are zero TODO change u0 to surface fluxes?
        un = np.multiply(u, self.fn)
        self.beta = np.sum(beta_ - un[:, :, 0] , axis = 1)  # TODO for nD

        self.alpha = -beta_
        for (i, j) in self.grid.boundary_faces:
            self.alpha[i, j] = 0  # zero flux over boundaries

    def prepare_linear_system(self):
        """ builds the linear system """
        self.create_alpha_beta()
        alpha_i = np.outer(np.array(range(0, self.n), dtype = np.int64), np.ones((1, self.grid.number_of_neighbours), dtype = np.int64))  # [0..n] x [1 1]
        alpha_j = self.grid.neighbours.copy()
        for (i, j) in self.grid.boundary_faces:
            alpha_j[i, j] = 0  # corresponding alpha = 0 at these boundary cells
        self.A = sparse.coo_matrix((np.hstack((self.alpha.flat, self.beta)), (np.hstack((alpha_i.flat, alpha_i[:, 0].flat)),
                                                                              np.hstack((alpha_j.flat, alpha_i[:, 0].flat)))))

    def solve_forward_euler(self, c, dt):
        """ explicit euler (matrix vector multiplication) """
        return (dt * self.A + I(self.n)) * c

    def solve_backward_euler(self, c, dt):
        """ implicit euler (solve a linear system) """
        return LA.spsolve(I(self.n) - dt * self.A, c, use_umfpack = True)

    def solve_crank_nicolson(self, c, dt):
        """ arithmetic mean of forward and backward eulers """
        return 0.5 * (self.solve_forward_euler(c, dt) + self.solve_backward_euler(c, dt))

    def bc_concentration(self, i, c, dx, n):
        """ flux boundary condition [g / cm2 / day] """
        q = (self.D[i]) * (c - self.x0[i]) / dx - np.inner(self.u[i], n)
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
                bc = self.bc_concentration(i, *v[0:3])
            else:
                raise("Unkown boundary condition")
            self.sources[i] = dt * bc * self.grid.area_per_volume[i, j]  # [g / cm3]

    def cfl(self, dt):
        """ returns the Courant-Friedrichs-Lewy number for time step @param dt """
        min_dx = np.min(self.grid.dx[:])
        max_u = np.max(np.linalg.norm(self.u, axis = 1))
        max_D = np.max(self.D[:])  # speed of diffusion is approximated with 1
        return max(max_u, max_D) * dt / min_dx

    def cfl_dt(self, c):
        """ returns a time step to achieve the Courant-Friedrichs-Lewy @param c """
        min_dx = np.min(self.grid.dx[:])
        max_u = np.max(np.linalg.norm(self.u, axis = 1))
        max_D = np.max(self.D[:])  # speed of diffusion is approximated with 1
        return c * min_dx / max(max_u, max_D)


class FVAdvectionDiffusion1D(FVAdvectionDiffusion):
    """ 
        Specialization for 1D:
        uses a direct banded solver, and builds sparse accordingly 
    """

    def __init__(self, grid :FVGrid):
        super().__init__(grid)
        self.ab = np.zeros((3, self.n))

    def prepare_linear_system(self):
        """ builds the linear system """
        self.create_alpha_beta()
        self.A = sparse.diags([self.alpha[1:, 0], self.beta, self.alpha[:-1, 1]], [-1, 0, 1])  # only needed for explicit Euler

    def solve_backward_euler(self, c, dt):
        """ implicit euler (solve a linear system) """
        self.ab[0, 0] = 0
        self.ab[2, -1] = 0
        self.ab[0, 1:] = -dt * self.alpha[:-1, 1]
        self.ab[1, :] = np.ones(self.beta.shape) - dt * self.beta
        self.ab[2, :-1] = -dt * self.alpha[1:, 0]
        return la.solve_banded ((1, 1), self.ab, c, overwrite_ab = True)


class FVAdvectionDiffusion_richards(FVAdvectionDiffusion):
    """ 
    Links FVAdvectionDiffusion to a FVRichards model, 
    for effective Diffusivity (Millington and Quirk, 1961) and Darcy velocity 
    """

    def __init__(self, grid :FVGrid, richards):
        super().__init__(grid)
        self.richards = richards

    def solver_proceed(self, x, dt):
        """ retrieve effective diffusion and velocity field from richards """
        r = self.richards

        self.u = r.darcy_velocity()
        if len(self.u.shape) == 1:  # (n,) -> (n,1)
            self.u = np.expand_dims(self.u, axis = 1)

        theta = vg.water_content(r.x0, r.soil)
        phi = r.soil.theta_S
        sw = theta / phi
        self.D = phi * (sw ** 3) * np.cbrt(theta) * self.D0
        self.b = self.b0 + theta
        self.prepare_linear_system()


class FVAdvectionDiffusion1D_richards(FVAdvectionDiffusion1D):
    """ 
    Links FVAdvectionDiffusion1D to a FVRichards model, 
    for effective Diffusivity (Millington and Quirk, 1961) and Darcy velocity 
    """

    def __init__(self, grid :FVGrid, richards):
        super().__init__(grid)
        self.richards = richards

    def solver_proceed(self, x, dt):
        """ retrieve effective diffusion and velocity field from richards """
        r = self.richards

        self.u = r.darcy_velocity()
        if len(self.u.shape) == 1:  # (n,) -> (n,1)
            self.u = np.expand_dims(self.u, axis = 1)

        theta = vg.water_content(r.x0, r.soil)
        phi = r.soil.theta_S
        sw = theta / phi
        self.D = phi * (sw ** 3) * np.cbrt(theta) * self.D0
        self.b = self.b0 + theta
        self.prepare_linear_system()

