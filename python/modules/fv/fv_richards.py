import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
import scipy.linalg as la
import timeit

from fv.fv_grid import *
from fv.fv_solver import *
import functional.van_genuchten as vg


class FVRichards(FVSolver):
    """ 
    Solves richards equation using finite volumes following van Dam and Feddes (2000) 
    
    
    van Dam, J.C., Feddes, R.A., 2000. Numerical simulation of infiltration,
    evaporation and shallow groundwater levels with the Richards equation.
    J. Hydrol. 233, 72-85.    
    
    FV_Richards is the base class of    
    FV_Richards1D(FV_Richards)         adds a special bounary conditions for 1d, uses a driect banded solver
    """

    def __init__(self, grid:FVGrid, soil):
        """ Initializes Richards solver """
        super().__init__(grid)
        self.soil = vg.Parameters(soil)  #  currently, homogeneous soil (todo)
        #
        n = self.n
        self.k = np.zeros(grid.neighbours.shape)  # contains precomputed  hydraulic conductivities [cm/day]
        self.alpha = np.zeros((n * grid.number_of_neighbours,))
        i_ = np.array(range(0, n), dtype = np.int64)
        cols = np.ones((1, self.grid.number_of_neighbours), dtype = np.int64)
        self.alpha_i = np.outer(i_, cols)
        self.alpha_j = self.grid.neighbours.copy()
        for (i, j) in self.grid.boundary_faces:
            self.alpha_j[i, j] = 0  # corresponding alpha = 0 at these boundary cells
            self.alpha_j[i, j] = 0
        self.beta_const = np.zeros((n,))  # const part of diagonal entries [1]
        self.f_const = np.zeros((n,))  # const part of load vector [1]
        #
        self.gravitation = 0  # assure dimension equals grid.dim (todo) currently gravitation is ignored
        # Solver settings
        self.max_iter = 100
        self.stop_tol = 1e-9  # stopping tolerance
        # self.NB = 10  # boundary iterations (for flux_out, rootsystem, ...TODO)s

    def solver_initialize(self):
        """ call back function for initialization"""
        self.sim_time = 0.  # current simulation time [day]

    def solver_proceed(self, dt):
        """ call back function for each successfully performed time step """
        pass

    def solve_single_step(self, dt, verbose):
        """ picard iteration (a fix point iteration) """
        # t0 = timeit.default_timer()
        # create constant parts dependent on self.x0
        self.create_k()
        # t1 = timeit.default_timer()
        self.create_f_const()
        # t2 = timeit.default_timer()
        self.create_alpha_beta(dt)
        # t3 = timeit.default_timer()
        self.prepare_linear_system()
        # t4 = timeit.default_timer()

        h_pm1 = self.x0.copy()
        for i in range(0, self.max_iter):

            c_pm1 = vg.specific_moisture_storage(h_pm1, self.soil)
            theta_pm1 = vg.water_content(h_pm1, self.soil)
            beta = c_pm1 + self.beta_const
            self.bc_to_source(dt)
            f = np.multiply(c_pm1, h_pm1) - theta_pm1 + self.f_const + self.sources
            # t5 = timeit.default_timer()

            try:
                h = self.solve_linear_system(beta, f)
#                 t6 = timeit.default_timer()
#                 if i == 0:
#                     tdt = timeit.default_timer() - t0
#                     tdt1 = t1 - t0
#                     tdt2 = t2 - t1
#                     tdt3 = t3 - t2
#                     tdt4 = t4 - t3
#                     tdt5 = t5 - t4
#                     tdt6 = t6 - t5
#                     print("Step {:g} s : {:g}, {:g}, {:g}, {:g}, {:g}, {:g}".format(tdt, tdt1 / tdt, tdt2 / tdt, tdt3 / tdt, tdt4 / tdt, tdt5 / tdt, tdt6 / tdt))
            except:
               print("Linear solver did raise something")
               return (h_pm1, False, i)

            delta_h = np.linalg.norm(h - h_pm1, np.inf)  # check if it is converging
            if delta_h < self.stop_tol:
                return (h, True, i)
            else:
                h_pm1 = h
        print("Maximal number of iterations reached")
        return (h_pm1, False, self.max_iter)

    def prepare_linear_system(self):
        """ the part of the linear system that can be constructed outside the fix-point iteration """
        self.A = sparse.coo_matrix((self.alpha.flat, (self.alpha_i.flat, self.alpha_j.flat)))

    def solve_linear_system(self, beta, f):
        """ constructs and solves linear system """
        B = sparse.coo_matrix((beta, (np.array(range(0, self.n)), np.array(range(0, self.n)))))
        # plt.spy(A + B)
        # plt.show()
        return LA.spsolve(self.A + B, f, use_umfpack = True)

    def create_k(self):
        """ sets up hydraulic conductivities from last time step solution self.x0, for each neighbour"""
        cols = np.ones((1, self.grid.number_of_neighbours))
        a = vg.hydraulic_conductivity(self.x0, self.soil)
        a_ = np.outer(a, cols)
        b = a[self.grid.neighbours]
        self.k = 2 * np.divide(np.multiply(a_, b), a_ + b)

    def create_f_const(self):
        """ sets up constant part of the load vector """
        self.f_const = vg.water_content(self.x0, self.soil)

    def create_alpha_beta(self, dt):
        """ calculate net fluxes over the cells """
        v = dt * np.divide(np.multiply(self.grid.area_per_volume, self.k), self.grid.dx)
        for (i, j) in self.grid.boundary_faces:
            v[i, j] = 0.
        self.alpha = -v
        self.beta_const = np.sum(v, axis = 1)

    def bc_flux(self, q):
        """ flux boundary condition [cm3 / cm2 / day] """
        return q

    def bc_flux_out(self, i, q, crit_p, dx):
        """ outflow (out of domain) boundary condition limits to critical pressure [cm3 / cm2 / day]"""
        k = vg.hydraulic_conductivity(self.x0[i], self.soil)
        max_q = k * (crit_p - self.x0[i]) / dx  # maximal possible outflux
        return min(max(q, max_q), 0.)

    def bc_flux_in(self, i, q, crit_p, dx):
        """ inflow (into domain) boundary condition limits to maximal possible infiltration [cm3 / cm2 / day]"""
        k = vg.hydraulic_conductivity(self.x0[i], self.soil)  # [cm / day]
        max_q = k * (0. - self.x0[i]) / dx  # maximal possible influx [cm3 / cm2 /day]
        return max(min(q, max_q), 0.)

    def bc_flux_in_out(self, i, q, crit_p, dx):
        """ inflow and outflow boundary condition limited to critical pressures [cm3 / cm2 / day]"""
        if q >= 0:
            return self.bc_flux_in(i, q, 0., dx)
        else:
            return self.bc_flux_out(i, q, crit_p, dx)

    def bc_to_source(self, dt):
        """ 
        Computes a cell source term from a boundary conditions given in self.bc, 
        self.bc is dictionary with keys (cell, face_id) and values is a boundary function returning a value [cm3 / cm2 / day].          
        
        Only one face per cell is allowed, overwrites any other sources
        
        dictionary with lambda functions would be nicer, but causes trouble with parallelisation    
        """
        for (i, j), (type, v) in self.bc.items():
            if type == "flux_in_out":
                bc = self.bc_flux_in_out(i, *v[0:2], self.grid.dx[i, j])
            elif type == "flux_in":
                bc = self.bc_flux_in(i, *v[0:2], self.grid.dx[i, j])
            elif type == "flux_out":
                bc = self.bc_flux_out(i, *v[0:2], self.grid.dx[i, j])
            elif type == "flux":
                bc = self.bc_flux(v[0])
            else:
                raise Exception("Unkown boundary condition")
            self.sources[i] = dt * bc * self.grid.area_per_volume[i, j]  # [cm3 / cm3]

    def get_flux(self, cell_id, face_id):
        """ flux [cm3/cm2/day] over the inner face given by @param face_id in cell @param cell_id """
        i = cell_id
        return self.k[i, face_id] * (self.x0[i] - self.x0[i + 1]) / self.grid.dx[i, face_id]  # [cm3/cm2/day]

    def get_water_content(self):
        """ current water content for each cell """
        return vg.water_content(self.x0, self.soil)

    def darcy_velocity(self):
        """ Darcy velocity """
        self.create_k()  # todo remove, but make sure to take current k
        dx = np.zeros(self.grid.dx.shape)
        for i in range(0, self.x0.shape[0]):  # TODO
            for j in range(0, self.grid.number_of_neighbours):
                ni = self.grid.neighbours[i, j]
                if ni >= 0:
                    dx[i, j] = self.x0[i] - self.x0[ni]
        q = np.divide(np.multiply(self.k, dx) , self.grid.dx)  # [cm3/cm2/day]
        # multiply by normal (1d) TODO (nd)

        # include boundary fluxes
        self.bc_to_source(1.)  # cm3/cm3/day # todo remove, but make sure to take current sources
        for (i, j), (type, v) in self.bc.items():
                    q[i, j] -= self.sources[i] / self.grid.area_per_volume[i, j]

        q[:, 0] *= -1.

        return np.mean(q, axis = 1)


class FVRichards1D(FVRichards):
    """ 1d sepezialisation of the general richards solverFV_Richards
        
        picks a direct banded linear solver (instead of default umfpack)
        adds the realized inner flux (left interval boundary) for cylindrical models
        adds special boundary conditions to couple to a root system
    """

    def __init__(self, grid:FVGrid, soil):
        super().__init__(grid, soil)
        delattr(self, 'alpha_i')  # indices are not longer needed
        delattr(self, 'alpha_j')
        self.ab = np.zeros((3, self.n))  # banded matrix
        self.innerFlux = 0.  # cm3/cm2/day
        # for get_inner_head
        self.dx0 = self.grid.center(0) - self.grid.nodes[0]
        self.dx1 = self.grid.center(1) - self.grid.center(0)

    def prepare_linear_system(self):
        """ the part of the linear system that can be constructed outside the fix-point iteration """
        self.ab[0, 0] = 0
        self.ab[2, -1] = 0
        self.ab[0, 1:] = self.alpha[:-1, 1]
        self.ab[2,:-1] = self.alpha[1:, 0]

    def solve_linear_system(self, beta, f):
        """ banded solver to solve tridiagonal matrix """
        self.ab[1,:] = beta  # diagonal is changed each iteration
        return la.solve_banded ((1, 1), self.ab, f, overwrite_b = True)

    def solver_initialize(self):
        """ call back function """
        self.innerFlux = 0.

    def solver_proceed(self, x, dt):
        """ call back function for each time step """
        self.innerFlux += self.sources[0] / (self.grid.area_per_volume[0, 0])  # [cm3/cm3] -> [cm3 /cm2] exact inner flux

    def get_inner_flux(self):
        """ exact flux at cell_id = 0 from last time step, todo generalize """
        return self.innerFlux

    def get_inner_head(self):
        """ extrapolates last result to root surface (specalized 1d) """
        x0, x1 = self.x0[0], self.x0[1]
        return x0 - ((x1 - x0) / self.dx1) * self.dx0

    def bc_to_source(self, dt):
        """ 
        Computes a cell source term from a boundary conditions given in self.bc, 
        self.bc is dictionary with keys (cell, face_id) and values is a boundary function returning a value [cm3 / cm2 / day].          
        
        Only one face per cell is allowed, overwrites any other sources
        
        dictionary with lambda functions would be nicer, but causes trouble with parallelisation    
        """
        for (i, j), (type, v) in self.bc.items():
            if type == "rootsystem_exact":
                bc = self.bc_rootsystem_exact(*v[0:6], self.grid.dx[i, j])
            elif type == "flux_in_out_mid":
                x = self.x0[i]
                bc = self.bc_flux_in_out(i, *v[0:2], self.grid.dx[i, j])
                self.x0[i] += (dt / 2) * bc * self.grid.area_per_volume[i, j]  # [cm3 / cm3]
                bc = self.bc_flux_in_out(i, *v[0:2], self.grid.dx[i, j])
                self.x0[i] = x
            elif type == "flux_in_out":
                bc = self.bc_flux_in_out(i, *v[0:2], self.grid.dx[i, j])
            elif type == "flux_in":
                bc = self.bc_flux_in(i, *v[0:2], self.grid.dx[i, j])
            elif type == "flux_out":
                bc = self.bc_flux_out(i, *v[0:2], self.grid.dx[i, j])
            elif type == "flux":
                bc = self.bc_flux(v[0])
            elif type == "rootsystem":
                bc = self.bc_rootsystem(*v[0:2], self.grid.dx[i, j])
            else:
                raise Exception("Unkown boundary condition")
            self.sources[i] = dt * bc * self.grid.area_per_volume[i, j]  # [cm3 / cm3]

    def bc_rootsystem(self, rx, kr, dx):
        """ flux is given by radial conductivity times difference in matric potentials 
        @param rx      root xylem matric potential [cm]
        @param kr      root radial condcutivitiy [1 / day] (intrinsic)
        """
        h = self.get_inner_head()
        k = vg.hydraulic_conductivity(h, self.soil)  # [cm / day]
        # dx = self.grid.nodes[0]  # = radius [cm]
        q = min(kr, k / dx) * (rx - h)
        return q  # [cm3 / cm2 / day]

    def bc_rootsystem_exact(self, rx0, rx1, kr, kz, a, l, dx):
        """ flux is given following Meunier et al., using exact solution over single segment  
        @param kz      root xylem matric potential [cm]
        @param kr      root radial condcutivitiy [1 / day] (intrinsic)
        """
        h = self.get_inner_head()  # value at self.grid.nodes[0]
        k = vg.hydraulic_conductivity(h, self.soil)  # [cm / day]
        dx = self.grid.nodes[0]  # = radius [cm]
        kr = min(kr, k / dx)
        f = 2 * a * np.pi * kr
        tau = np.sqrt(f / kz)  # sqrt(c) [cm-1]
        d = np.exp(-tau * l) - np.exp(tau * l)  #  det
        fExact = f * (1. / (tau * d)) * (rx0 - h + rx1 - h) * (2. - np.exp(-tau * l) - np.exp(tau * l))
        # print(10. / (24.*3600.) * fExact / (2 * a * np.pi * l))
        return fExact / (2 * a * np.pi * l)  # [cm3 / cm2 / day]

