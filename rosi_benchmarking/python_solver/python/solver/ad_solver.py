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

    def __init__(self, grid :FV_Grid, u):
        self.grid = grid  # simplistic fv grid
        self.n = self.grid.n_cells
        self.mid = self.grid.centers()  # precompute cell centers
        #
        self.b = 1  # buffer power
        self.D = 1e-4  # cm2 / day
        self.u = u  # velocity field
        n = self.n
        #
        self.beta = np.zeros((n,))  # const part of diagonal entries [1/day]
        self.alpha = np.zeros(grid.neighbours.shape)  # secondary diagonals [1/day]
        self.alpha_i = np.outer(i_, cols)  # construtor
        self.alpha_j = self.grid.neighbours
        for (i, j) in self.grid.boundary_faces:
            self.alpha_j[i, j] = 0  # corresponding alpha = 0 at these boundary cells
            self.alpha_j[i, j] = 0
        #
        self.c0 = np.zeros((n,))  # concentration of last time step [g/cm3]
        self.bc = { }  # boundary conditions, map with key (cell_id, face_id) containing list of values
        self.sources = np.zeros((n,))  # [cm3 / cm3]

    def solve(self, output_times :list, max_dt = 0.5, verbose = True):
        """ solves richards equation for several output times 
        @param output_times       final simulation times [day]
        @param max_dt             maximal internal time step
        """
        self.solver_initialize()
        dt = max_dt  #  initially proposed time step
        sim_time = 0  # current simulation time
        k = 0
        h_out = []
        while k < len(output_times):

            dt_ = min(output_times[k] - sim_time, dt)  #  actual time step
            dt_ = min(dt_, max_dt)

            self.prepare_linear_system()

            try:
                c2 = self.solve_backward_euler(self.c)
                ok = True
            except:
                ok = False

            if ok:

                sim_time = sim_time + dt_  # increase current time

                self.c0 = c2
                self.solver_success()

                if output_times[k] <= sim_time:  # store result
                    h_out.append(c2.copy())
                    k = k + 1

                if verbose:
                    print('Time {:g} days, last time step {:g}'.format(sim_time, dt_))

                if dt_ == dt:
                    dt = dt * 1.25

            else:
                dt = dt / 10.
                print("retry with max {:g} = {:g} day".format(dt, min(output_times[k] - sim_time, dt)))
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
        self.alpha = np.zeros(grid.neighbours.shape)
        for i in range(0, self.n):  # over cells
            for  j in range(0, self.grid.number_of_neighbours):
                if j >= 0:  # -1 means there is no neighoour
                    n = mid[neighbors[j]] - mid[i]
                    n = n / np.linalg.norm(n)
                    self.beta[i] += (1. / beta) * self.grid.area_per_volume[i, j] * (self.D / self.dx[i, j] - np.dot(self.u, n)) * self.c0[i]
                    self.alpha[i, j] += -(1. / beta) * self.grid.area_per_volume[i, j] * (self.D / self.dx[i, j]) * self.c0[neighbours[j]]
        A = sparse.coo_matrix((self.alpha.flat, (self.alpha_i.flat, self.alpha_j.flat)))
        B = sparse.coo_matrix((beta, (np.array(range(0, self.n)), np.array(range(0, self.n)))))
        self.A = dt * (A + B)

    def solve_forward_euler(self, c):
        """ todo """
        return (I(self.n) + self.A) * c

    def solve_backward_euler(self, c):
        """ solves a linear system """
        return LA.spsolve(I(self.n) - self.A, c, use_umfpack = True)

#     def solve_crank_nicolson(self, c):
#         return 0.5 * (self.solve_forward_euler(c) + self.solve_backward_euler(c))
#
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
#     def bc_to_source(self, dt):
#         """
#         Computes a cell source term from a boundary conditions given in self.bc,
#         self.bc is dictionary with keys (cell, face_id) and values is a boundary function returning a value [cm3 / cm2 / day].
#
#         Only one face per cell is allowed, overwrites any other sources
#
#         dictionary with lambda functions would be nicer, but causes trouble with parallelisation
#         """
#         for (i, j), (type, v) in self.bc.items():
#             if type == "rootsystem":
#                 bc = self.bc_rootsystem(*v[0:2])
#             elif type == "flux_in_out":
#                 bc = self.bc_flux_in_out(i, *v[0:3])
#             elif type == "rootsystem_exact":
#                 bc = self.bc_rootsystem_exact(*v[0:6])
#             elif type == "flux_in":
#                 bc = self.bc_flux_in(i, *v[0:3])
#             elif type == "flux_out":
#                 bc = self.bc_flux_out(i, *v[0:3])
#             elif type == "flux":
#                 bc = self.bc_flux(v[0])
#             else:
#                 raise("Unkown boundary condition")
#             self.sources[i] = dt * bc * self.grid.area_per_volume[i, j]  # [cm3 / cm3]

