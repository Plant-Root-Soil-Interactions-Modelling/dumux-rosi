import numpy as np
import copy

from fv_grid import *
from fv_solver import *
import van_genuchten as vg


class FVSystem(FVSolver):
    """ 
    A Simulation loop for several equations (FVRichards, or FVAdvectionDiffusion)
    """

    def __init__(self, grid :FVGrid):
        """ Initializes the solver """
        super().__init__(grid)
        self.solvers = []

    def add_eqn(self, solver):
        """ add a solver to the system of equations, one solver per equation """
        self.solvers.append(solver)

    def solve(self, output_times :list, max_dt = 0.5, verbose = True):
        """ solves richards equation for several output times 
        @param output_times       final simulation times [day]
        @param max_dt             maximal time step [day]
        @param verbose            tell me more
        """
        for s in self.solvers:
            s.solver_initialize()
        self.solver_initialize()

        dt = max_dt  #  initially proposed time step
        k = 0
        h_out = []
        while k < len(output_times):

            dt_ = min(output_times[k] - self.sim_time, dt)  #  actual time step
            dt_ = min(dt_, max_dt)

            x, ok, i = self.solve_single_step(dt_, verbose)

            if ok:

                self.sim_time = self.sim_time + dt_  # increase current time
                for s in self.solvers:
                    s.sim_time = self.sim_time + dt_

                self.solver_proceed(x, dt_)

                if output_times[k] <= self.sim_time:  # store result
                    h_out.append(copy.deepcopy(x))
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
        self.sim_time = 0.

    def solver_proceed(self, x, dt):
        """ call back function after each successfully performed time step, e.g. to update coefficientsf from coupled equations
        @param dt     time step [day]
        """
        for j, s in enumerate(self.solvers):
            s.solver_proceed(x[j], dt)
            s.x0 = x[j]

    def solve_single_step(self, dt, verbose):
        """ solve for time step dt, called by solve with step size control 
        @param dt     time step [day]
        @param verbose            tell me more        
        """
        x = [None] * len(self.solvers)
        ok = True
        i = 0
        for j, s in enumerate(self.solvers):
            x[j], ok_, i_ = s.solve_single_step(dt, verbose)
            i = max(i, i_)
            ok = ok and ok_
            if not ok:
                break
        return (x, ok, i)

