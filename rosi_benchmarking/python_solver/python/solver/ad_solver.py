import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
from scipy.sparse import identity as I
import scipy.linalg as la
import matplotlib.pyplot as plt

import solver.van_genuchten as vg
from solver.fv_grid import *


class FV_AdvectionDiffusion:
    """ 
    AdvectionDiffusion solver (first draft)
    """

    def __init__(self, grid :FV_Grid):
        self.grid = grid  # simplistic fv grid
        self.n = self.grid.n_cells
        self.mid = self.grid.centers()  # precompute cell centers
        #
        n = self.n
        self.b = np.ones((n,)) * 1.  # [1] buffer power
        self.D = np.ones((n,)) * 0.864  #  [cm2/day] 1.e-5 cm2/s = 0.864 cm2 / day
        self.u = np.zeros((n, self.grid.dim))  # [cm/day] velocity field
        self.c0 = np.zeros((n,))  # [g/cm3] concentration of last time step
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
        #
        self.sim_time = 0  # current simulation time
        self.bc = { }  # boundary conditions, map with key (cell_id, face_id) containing list of values
        self.sources = np.zeros((n,))  # [cm3 / cm3]

    def solve(self, output_times :list, max_dt = 0.5, verbose = True):
        """ solves richards equation for several output times 
        @param output_times       final simulation times [day]
        @param max_dt             maximal internal time step
        """
        self.solver_initialize()
        self.prepare_linear_system()
        dt = max_dt  #  initially proposed time step
        k = 0
        h_out = []

        while k < len(output_times):

            dt_ = min(output_times[k] - self.sim_time, dt)  #  actual time step
            dt_ = min(dt_, max_dt)

            self.bc_to_source(dt)
            # print("sources", self.sources[0])
            b = self.c0 + self.sources
            b[0] = max(b[0], 0.)
            c2 = self.solve_crank_nicolson(b, dt)
            ok = True
#             try:
#                 c2 = self.solve_backward_euler(self.c, dt)
#                 ok = True
#             except:
#                 ok = False

            if ok:

                self.sim_time = self.sim_time + dt_  # increase current time

                self.c0 = c2
                self.solver_success()

                if output_times[k] <= self.sim_time:  # store result
                    h_out.append(c2.copy())
                    k = k + 1

                if verbose:
                    print('Time {:g} days, last time step {:g}'.format(self.sim_time, dt_))

                if dt_ == dt:
                    dt = dt * 1.25

            else:
                dt = dt / 10.
                print("retry with max {:g} = {:g} day".format(dt, min(output_times[k] - self.sim_time, dt)))
                if dt < 1.e-10:
                    raise Exception("I did not find a solution")

        return np.array(h_out)

    def solver_initialize(self):
        """ call back function for initialization"""
        pass

    def solver_success(self):
        """ call back function for each time step """
        pass

    def prepare_linear_system(self):
        """ experimental """
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
        return 0.5 * (self.solve_forward_euler(c, dt) + self.solve_backward_euler(c, dt))

    def bc_concentration(self, i, c, dx):
        """ flux boundary condition [g / cm2 / day] """
        q = (1. / self.b[i]) * (self.D[i] / dx) * (c - self.c0[i])
        # print(q, "g/cm2/day")
        return q

#     def bc_flux(self, q):
#         """ flux boundary condition [cm3 / cm2 / day] """
#         return q
#
#     def bc_flux_out(self, i, q, crit_p, dx):
#         """ outflow boundary condition limits to out or zero flux, and to critical pressure [cm3 / cm2 / day]"""
#         k = vg.hydraulic_conductivity(self.h0[i], self.soil)
#         max_q = k * (crit_p - self.h0[i]) / dx  # maximal possible outflux
#         # print("bc_flux_out", q, max_q, i , crit_p, self.h0[i], k, dx)
#         return min(max(q, max_q), 0.)
#
#     def bc_flux_in(self, i, q, crit_p, dx):
#         """ inflow boundary condition limits to out or zero flux, and to critical pressure [cm3 / cm2 / day]"""
#         k = vg.hydraulic_conductivity(self.h0[i], self.soil)  # [cm / day]
#         max_q = k * (crit_p - self.h0[i]) / dx  # maximal possible influx [cm3 / cm2 /day]
#         # print("bc_flux_in", q, max_q, i, crit_p, self.h0[i], k, dx)
#         return max(min(q, max_q), 0.)
#
#     def bc_flux_in_out(self, i, q, crit_p, dx):
#         """ inflow and outflow boundary condition limited to critical pressures [cm3 / cm2 / day]"""
#         if q >= 0:
#             return self.bc_flux_in(i, q, 0., dx)
#         else:
#             return self.bc_flux_out(i, q, crit_p, dx)
#
    def bc_to_source(self, dt):
        """
        Computes a cell source term from a boundary conditions given in self.bc,
        self.bc is dictionary with keys (cell, face_id) and values is a boundary function returning a value [cm3 / cm2 / day].

        Only one face per cell is allowed, overwrites any other sources

        dictionary with lambda functions would be nicer, but causes trouble with parallelisation
        """
        for (i, j), (type, v) in self.bc.items():
            if type == "concentration":
                bc = self.bc_concentration(i, *v[0:2])
            else:
                raise("Unkown boundary condition")
            self.sources[i] = dt * bc * self.grid.area_per_volume[i, j]  # [g / cm3]

