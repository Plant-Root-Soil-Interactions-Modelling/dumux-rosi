import numpy as np

from solver.fv_grid import *
import solver.van_genuchten as vg


class FVSolver:
    """ 
    Base class for FVRichards, FVAdvectionDiffusion    
    """

    def __init__(self, grid :FVGrid):
        """ Initializes the solver """
        self.grid = grid
        self.n = self.grid.n_cells
        self.x0 = np.zeros((self.n,))  # solution of last time step [cm]
        self.sim_time = 0.
        self.bc = { }  # boundary conditions, map with key (cell_id, face_id) containing list of values
        self.sources = np.zeros((self.n,))  # [cm3 / cm3]

    def solve(self, output_times :list, max_dt = 0.5, verbose = True):
        """ solves richards equation for several output times 
        @param output_times       final simulation times [day]
        @param max_dt             maximal time step [day]
        @param verbose            tell me more
        """
        self.solver_initialize()
        self.sim_time = 0.  # current simulation time [day]
        dt = max_dt  #  initially proposed time step
        k = 0
        h_out = []
        while k < len(output_times):

            dt_ = min(output_times[k] - self.sim_time, dt)  #  actual time step
            dt_ = min(dt_, max_dt)
            x, ok, i = self.solve_single_step(dt_, verbose)

            if ok:

                self.sim_time = self.sim_time + dt_  # increase current time

                self.x0 = x
                self.solver_proceed(dt_)

                if output_times[k] <= self.sim_time:  # store result
                    h_out.append(x.copy())
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

    def solver_initialize(self):
        """ call back function for initialization"""
        pass

    def solver_proceed(self, dt):
        """ call back function after each successfully performed time step, e.g. to update coefficientsf from coupled equations
        @param dt     time step [day]
        """
        pass

    def solve_single_step(self, dt, verbose):
        """ solve for time step dt, called by solve with step size control 
        @param dt     time step [day]
        @param verbose            tell me more        
        """
        pass
