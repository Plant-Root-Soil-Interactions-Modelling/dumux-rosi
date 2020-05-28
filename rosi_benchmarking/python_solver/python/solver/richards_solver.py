import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
import matplotlib.pyplot as plt

import solver.van_genuchten as vg
from solver.fv_grid import *


class FV_Richards:
    """ 
    Solves richards equation using finite volumes following van Dam and Feddes (2000) 
    
    
    van Dam, J.C., Feddes, R.A., 2000. Numerical simulation of infiltration,
    evaporation and shallow groundwater levels with the Richards equation.
    J. Hydrol. 233, 72-85.    
    """
    
    def __init__(self, grid :FV_Grid, soil):
        self.grid = grid  # simplistic fv grid           
        self.soil = vg.Parameters(soil)  #  currently, homogeneous soil (todo)                
        # matrix
        self.k = np.zeros(grid.neighbours.shape)  # contains precomputed  hydraulic conductivities [cm/day]
        self.alpha = np.zeros((grid.n_cells * grid.number_of_neighbours(),))
        self.alpha_i = np.zeros((grid.n_cells * grid.number_of_neighbours(),), dtype=np.int64)
        self.alpha_j = np.zeros((grid.n_cells * grid.number_of_neighbours(),), dtype=np.int64)      
        self.beta_const = np.zeros((grid.n_cells,))  # const part of diagonal entries [1] 
        self.f_const = np.zeros((grid.n_cells,))  # const part of load vector [1]           
        # 
        self.bc = { }  # boundary conditions, map with key (cell_id, face_id) containing list of values
        self.sources = np.zeros((grid.n_cells,))  # [cm3 / cm3]
        self.h0 = np.zeros((grid.n_cells,))  # [cm]   
        self.innerFlux = 0.  # cm3/cm2/day 
        #
        self.gravitation = 0  # assure dimension equals grid.dim (todo)
        
    def create_k(self):
        """ sets up hydraulic conductivities from last time step solution @param h , for each neighbour"""
        for i in range(0, self.grid.n_cells):                        
            for j, ni in enumerate(self.grid.neighbours[i, :]):                        
                    if ni > -1:
                        a = vg.hydraulic_conductivity(self.h0[i], self.soil)
                        b = vg.hydraulic_conductivity(self.h0[ni], self.soil)
                        self.k[i, j] = 2 * a * b / (a + b)  # harmonic 
#                        self.k[i, j] = 0.5 * (a + b)  # arithmetic
#                        self.k[i, j] = vg.hydraulic_conductivity(0.5 * (h[i] + h[ni]), self.soil)  # arithmetic mean2
    
    def create_f_const(self):
        """ sets up constant part of the load vector """
        self.f_const = vg.water_content(self.h0, self.soil)   
                        
    def create_alpha_beta(self, dt):  
        """ calculate net fluxes over the cells """
        c = 0
        for i in range(0, self.grid.n_cells):                        
            self.beta_const[i] = 0.
            for j, ni in enumerate(self.grid.neighbours[i, :]):                        
                if ni > -1:
                    dx = np.linalg.norm(self.grid.mid[i] - self.grid.mid[ni])  # we could precompute these
                    v = dt * self.k[i, j] * (self.grid.area[i, j] / self.grid.volume[i]) / dx
                    self.alpha[c] = -v 
                    self.alpha_i[c] = i
                    self.alpha_j[c] = ni
                    c += 1
                    self.beta_const[i] += v                    
    
    def bc_flux(self, i, j, v):
        """ flux boundary condition """
        return v[0] * (self.grid.area[i, j] / self.grid.volume[i])  # [1/day]   
    
    def bc_flux_out(self, i, j, v):
        """ outflow boundary condition limits to out or zero flux, and to critical pressure"""
        q, crit_p, dx = v[0], v[1], v[2]  
        q = min(q, 0.)  # limit to negative (out) flux  
        k = vg.hydraulic_conductivity(self.h0[i], self.soil)
        max_q = k * (crit_p - self.h0[i]) / dx  # maximal possible outflux   
        return min(max(q, max_q), 0.) * (self.grid.area[i, j] / self.grid.volume[i])  # [1/day]    

    def bc_flux_in(self, i, j, v):
        """ outflow boundary condition limits to out or zero flux, and to critical pressure"""
        q, crit_p, dx = v[0], v[1], v[2]  
        q = max(q, 0.)  # limit to positive (in) flux  
        k = vg.hydraulic_conductivity(self.h0[i], self.soil)
        max_q = k * (crit_p - self.h0[i]) / dx  # maximal possible outflux   
        # print("bc_flux_in", q, max_q, i , j, crit_p, self.h0[i], k, dx)
        return max(min(q, max_q), 0.) * (self.grid.area[i, j] / self.grid.volume[i])  # [1/day]    

    def bc_rootsystem(self, i, j, v):
        """ flux is given by radial conductivity times difference in matric potentials """
        rx = v[0]  # root xylem matric potential [cm]
        kr = v[1]  # root condcutivitiy cm3 / day          
        k = vg.hydraulic_conductivity(self.h0[i], self.soil)            
        q = min(kr, k) * (rx - self.h0[i]) * (self.grid.area[i, j] / self.grid.volume[i])  #
        # print("root system flux", rx, self.h0[i], kr, k, q)            
        return min(q, 0.)  # limit to outflux  [1/day]

    def bc_to_source(self, dt):
        """ 
        Computes a cell source term from a boundary conditions given in self.bc, 
        self.bc is dictionary with keys (cell, face_id) and 
        values is a list, where the first argument is a function, i.e. v[0](cell_id, face_id, v[1:]), 
        and the other list values are passed to this boundary funcntion    
        """          
        for k, v in self.bc.items(): 
            i, j = k[0], k[1]  # cell_id, face_id        
            bc = v[0]  # fetch lambda bc function            
            self.sources[i] = dt * bc(i, j, v[1:])
    
    def getFlux(self, cell_id, face_id):
        """ flux [cm3/cm2/day] over the inner face given by @param face_id in cell @param cell_id """
        i = cell_id
        return self.k[i, face_id] * (self.h0[i] - self.h0[i + 1]) / dx  # [cm3/cm2/day]

    def getWaterContent(self):
        """ current water content for each cell """
        return vg.water_content(self.h0, self.soil)
    
    def getInnerFlux(self):
        return self.innerFlux
#         q = self.bc[(0, 0)][0]
#         crit_p = self.bc[(0, 0)][1]  
#         max(q, 0.)
#         dx = self.grid.mid[0] - self.grid.nodes[0]  # TODO only holds for 1D, missing concept    
#         a = vg.hydraulic_conductivity(h0, self.soil)
#         b = vg.hydraulic_conductivity(crit_p, self.soil)
#         # k = 2 * a * b / (a + b)  
#         # k = (a + b) / 2             
#         k = a 
#         max_q = k * (crit_p - h0) / dx  # SIGN        
#         # print("InnerFlux", q, max_q, h0)
#         return min(max(q, max_q), 0.)    

    def picard_iteration(self, dt):        
        """ fix point iteration """
        n = self.grid.n_cells    
        max_iter = 100
        stop_tol = 1e-9  # stopping tolerance
        
        # create constant parts dependent on self.h0
        self.create_k()
        self.create_f_const() 
        self.create_alpha_beta(dt) 

        h_pm1 = self.h0.copy()        
        for i in range(0, max_iter):                        
            c_pm1 = vg.specific_moisture_storage(h_pm1, self.soil)
            theta_pm1 = vg.water_content(h_pm1, self.soil)            
            beta = c_pm1 + self.beta_const
            self.bc_to_source(dt)
            f = np.multiply(c_pm1, h_pm1) - theta_pm1 + self.f_const + self.sources           
            A = sparse.coo_matrix((self.alpha, (self.alpha_i, self.alpha_j)))
            B = sparse.coo_matrix((beta, (np.array(range(0, n)), np.array(range(0, n)))))  
#             plt.spy(A + B)
#             plt.show()
            try:
                h = LA.spsolve(A + B, f, use_umfpack=True)
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
    
    def solve(self, output_times :list, max_dt=0.5, verbose=True):
        """ solves richards equation for several output times 
        @param output_times       final simulation times [day]
        @param max_dt             maximal internal time step
        """ 
        dt = 0.1  #  time step
        
        self.innerFlux = 0.  # for this solving step
        sim_time = 0  # current simulation time
        k = 0
        h_out = []
        while k < len(output_times):
            
            dt_ = min(output_times[k] - sim_time, dt)  #  actual time step            
            dt_ = min(dt_, max_dt)
            h2, ok, i = self. picard_iteration(dt_)

            if ok:          
                self.innerFlux += self.sources[0] / (self.grid.area[0, 0] / self.grid.volume[0])  # exact inner flux    
                sim_time = sim_time + dt_  # increase current time
                
                if output_times[k] <= sim_time:
                    h_out.append(h2.copy())
                    k = k + 1
                if verbose:
                    print('Time {:g} days, iterations {:g}, last time step {:g}'.format(sim_time, i, dt_))
                
                self.h0 = h2  
                
                if dt_ == dt:
                    if i < 5:
                        dt = dt * 1.25
                    else:
                        dt = dt / 1.25                
            else:         
                dt = dt / 10
                print("retry with max {:g} = {:g} day".format(dt, min(output_times[k] - sim_time, dt)))
                if dt < 1.e-10:
                    print("\nI failed you\n")
                    return
    
        return np.array(h_out)

