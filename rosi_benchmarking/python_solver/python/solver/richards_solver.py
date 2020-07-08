import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
import scipy.linalg as la
import matplotlib.pyplot as plt

import solver.van_genuchten as vg
from solver.fv_grid import *


class FV_Richards:
    """ 
    Solves richards equation using finite volumes following van Dam and Feddes (2000) 
    
    
    van Dam, J.C., Feddes, R.A., 2000. Numerical simulation of infiltration,
    evaporation and shallow groundwater levels with the Richards equation.
    J. Hydrol. 233, 72-85.    
    
    FV_Richards is the base class of    
    FV_Richards1D(FV_Richards)         adds a special bounary conditions for 1d, uses a driect banded solver
    """

    def __init__(self, grid :FV_Grid, soil):
        self.grid = grid  # simplistic fv grid
        self.n = self.grid.n_cells
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
        self.sim_time = 0  # current simulation time [day]
        self.h0 = np.zeros((n,))  # solution of last time step [cm]
        self.bc = { }  # boundary conditions, map with key (cell_id, face_id) containing list of values
        self.sources = np.zeros((n,))  # [cm3 / cm3]
        #
        self.gravitation = 0  # assure dimension equals grid.dim (todo) currently gravitation is ignored

    def solve(self, output_times :list, max_dt = 0.5, verbose = True):
        """ solves richards equation for several output times 
        @param output_times       final simulation times [day]
        @param max_dt             maximal internal time step
        """
        self.solver_initialize()
        dt = max_dt  #  initially proposed time step
        k = 0
        h_out = []
        while k < len(output_times):

            dt_ = min(output_times[k] - self.sim_time, dt)  #  actual time step
            dt_ = min(dt_, max_dt)
            h2, ok, i = self.picard_iteration(dt_)

            if ok:

                self.sim_time = self.sim_time + dt_  # increase current time

                self.h0 = h2
                self.solver_success()

                if output_times[k] <= self.sim_time:  # store result
                    h_out.append(h2.copy())
                    k = k + 1

                if verbose:
                    print('Time {:g} days, iterations {:g}, last time step {:g}'.format(self.sim_time, i, dt_))

                if dt_ == dt:
                    if i < 5:
                        dt = dt * 1.25
                    else:
                        dt = dt / 1.25
            else:
                dt = dt / 10.
                print("retry with max {:g} = {:g} day".format(dt, min(output_times[k] - self.sim_time, dt)))
                if dt < 1.e-10:
                    raise Exception("I did not find a solution")

        return np.array(h_out)

    def solve_single(self, dt = 0.5, verbose = True):
        """ solves richards equation by using a single picard iteration. 
            for less overhead in sequential couplings for very small dt
        """
        self.solver_initialize()
        self.h0, ok, i = self.picard_iteration(dt)
        if verbose:
            print('Picard iterations {:g}, time step {:g}'.format(i, dt))
        if ok:
            self.solver_success()
        else:
            raise Exception("solve_single did not find a solution, decrease time step, or use solve")

    def solver_initialize(self):
        """ call back function for initialization"""
        pass

    def solver_success(self):
        """ call back function for each time step """
        pass

    def picard_iteration(self, dt):
        """ fix point iteration """
        max_iter = 100
        stop_tol = 1e-9  # stopping tolerance

        # create constant parts dependent on self.h0
        self.create_k()
        self.create_f_const()
        self.create_alpha_beta(dt)
        self.prepare_linear_system()

        h_pm1 = self.h0.copy()
        for i in range(0, max_iter):

            c_pm1 = vg.specific_moisture_storage(h_pm1, self.soil)
            theta_pm1 = vg.water_content(h_pm1, self.soil)
            beta = c_pm1 + self.beta_const
            self.bc_to_source(dt)
            f = np.multiply(c_pm1, h_pm1) - theta_pm1 + self.f_const + self.sources

            try:
                h = self.solve_linear_system(beta, f)
            except:
               print("Linear solver did raise something")
               return (h_pm1, False, i)

            delta_h = np.linalg.norm(h - h_pm1, np.inf)  # check if it is converging
            if delta_h < stop_tol:
                return (h, True, i)
            else:
                h_pm1 = h

        print("Maximal number of iterations reached")
        return (h_pm1, False, max_iter)

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
        """ sets up hydraulic conductivities from last time step solution self.h0, for each neighbour"""
        cols = np.ones((1, self.grid.number_of_neighbours))
        a = vg.hydraulic_conductivity(self.h0, self.soil)
        a_ = np.outer(a, cols)
        b = vg.hydraulic_conductivity(self.h0[self.grid.neighbours], self.soil)
        self.k = 2 * np.divide(np.multiply(a_, b), a_ + b)

    def create_f_const(self):
        """ sets up constant part of the load vector """
        self.f_const = vg.water_content(self.h0, self.soil)

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
        """ outflow boundary condition limits to out or zero flux, and to critical pressure [cm3 / cm2 / day]"""
        k = vg.hydraulic_conductivity(self.h0[i], self.soil)
        max_q = k * (crit_p - self.h0[i]) / dx  # maximal possible outflux
        # print("bc_flux_out", q, max_q, i , crit_p, self.h0[i], k, dx)
        return min(max(q, max_q), 0.)

    def bc_flux_in(self, i, q, crit_p, dx):
        """ inflow boundary condition limits to out or zero flux, and to critical pressure [cm3 / cm2 / day]"""
        k = vg.hydraulic_conductivity(self.h0[i], self.soil)  # [cm / day]
        max_q = k * (crit_p - self.h0[i]) / dx  # maximal possible influx [cm3 / cm2 /day]
        # print("bc_flux_in", q, max_q, i, crit_p, self.h0[i], k, dx)
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
                bc = self.bc_flux_in_out(i, *v[0:3])
            elif type == "flux_in":
                bc = self.bc_flux_in(i, *v[0:3])
            elif type == "flux_out":
                bc = self.bc_flux_out(i, *v[0:3])
            elif type == "flux":
                bc = self.bc_flux(v[0])
            else:
                raise("Unkown boundary condition")
            self.sources[i] = dt * bc * self.grid.area_per_volume[i, j]  # [cm3 / cm3]

    def getFlux(self, cell_id, face_id):
        """ flux [cm3/cm2/day] over the inner face given by @param face_id in cell @param cell_id """
        i = cell_id
        return self.k[i, face_id] * (self.h0[i] - self.h0[i + 1]) / self.grid.dx[i, face_id]  # [cm3/cm2/day]

    def getWaterContent(self):
        """ current water content for each cell """
        return vg.water_content(self.h0, self.soil)


class FV_Richards1D(FV_Richards):
    """ Sepecialisation of the general richards solverFV_Richards
        
        picks a direct banded linear solver (instead of default umfpack)
        adds the realized inner flux (left interval boundary) for cylindrical models
        add special boundary conditions to couple to a root system
    """

    def __init__(self, grid :FV_Grid, soil):
        super().__init__(grid, soil)
        delattr(self, 'alpha_i')  # indices are not longer needed
        delattr(self, 'alpha_j')
        self.ab = np.zeros((3, self.n))  # banded matrix
        self.innerFlux = 0.  # cm3/cm2/day
        # for getInnerHead()
        self.dx0 = self.grid.center(0)  #  (self.grid.mid[0] - self.grid.nodes[0])
        self.dx1 = self.grid.center(1) - self.grid.center(0)  # todo
        # for bc_rootsystem, bc_rootsystem_exact
        self.dx = self.grid.nodes[0]

    def prepare_linear_system(self):
        """ the part of the linear system that can be constructed outside the fix-point iteration """
        self.ab[0, 0] = 0
        self.ab[2, -1] = 0
        self.ab[0, 1:] = self.alpha[:-1, 1]
        self.ab[2, :-1] = self.alpha[1:, 0]

    def solve_linear_system(self, beta, f):
        self.ab[1, :] = beta  # diagonal is changed each iteration
        return la.solve_banded ((1, 1), self.ab, f, overwrite_b = True)

    def solver_initialize(self):
        """ call back function """
        self.innerFlux = 0.

    def solver_success(self):
        """ call back function for each time step """
        self.innerFlux += self.sources[0] / (self.grid.area_per_volume[0, 0])  # [cm3/cm3] -> [cm3 /cm2] exact inner flux

    def getInnerFlux(self):
        """ exact flux at cell_id = 0, todo generalize """
        return self.innerFlux

    def getInnerHead(self):
        """ extrapolates result to root surface (specalized 1d) """
        right_neighbour = self.grid.neighbours[0, 1]
        h0, h1 = self.h0[0], self.h0[1]
        h = h0 - ((h1 - h0) / self.dx1) * self.dx0
        # print("linear head extrapolation", h0, h1, h, dx0, dx1, 1)
        return h

    def bc_to_source(self, dt):
        """ 
        Computes a cell source term from a boundary conditions given in self.bc, 
        self.bc is dictionary with keys (cell, face_id) and values is a boundary function returning a value [cm3 / cm2 / day].          
        
        Only one face per cell is allowed, overwrites any other sources
        
        dictionary with lambda functions would be nicer, but causes trouble with parallelisation    
        """
        for (i, j), (type, v) in self.bc.items():
            if type == "rootsystem":
                bc = self.bc_rootsystem(*v[0:2])
                self.sources[i] = dt * bc * self.grid.area_per_volume[i, j]  # [cm3 / cm3]
            elif type == "rootsystem_exact":
                bc = self.bc_rootsystem_exact(*v[0:6])
                self.sources[i] = dt * bc * self.grid.area_per_volume[i, j]  # [cm3 / cm3]
            else:
                super().bc_to_source(dt)

    def bc_rootsystem(self, rx, kr):
        """ flux is given by radial conductivity times difference in matric potentials 
        @param rx      root xylem matric potential [cm]
        @param kr      root radial condcutivitiy [1 / day] (intrinsic)
        """
        h = self.getInnerHead()
        k = vg.hydraulic_conductivity(h, self.soil)  # [cm / day]
        q = min(kr, k / self.dx) * (rx - h)
        return q  # [cm3 / cm2 / day]

    def bc_rootsystem_exact(self, rx0, rx1, kr, kz, a, l):
        """ flux is given following Meunier et al., using exact solution over single segment  
        @param kz      root xylem matric potential [cm]
        @param kr      root radial condcutivitiy [1 / day] (intrinsic)
        """
        h = self.getInnerHead()
        k = vg.hydraulic_conductivity(h, self.soil)  # [cm / day]
        kr = min(kr, k / self.dx)
        f = -2 * a * np.pi * kr
        tau = np.sqrt(2 * a * np.pi * kr / kz)  # sqrt(c) [cm-1]
        d = np.exp(-tau * l) - np.exp(tau * l)  #  det
        fExact = -f * (1. / (tau * d)) * (rx0 - h + rx1 - h) * (2. - np.exp(-tau * l) - np.exp(tau * l))
        return fExact / (2 * a * np.pi * l)  # [cm3 / cm2 / day]

