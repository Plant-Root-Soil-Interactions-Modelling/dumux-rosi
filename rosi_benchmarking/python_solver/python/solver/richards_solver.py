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
        i_ = np.array(range(0, self.grid.n_cells), dtype=np.int64)
        cols = np.ones((1, self.grid.number_of_neighbours()), dtype=np.int64)  
        self.alpha_i = np.outer(i_, cols)  # construtor 
        self.alpha_j = self.grid.neighbours  # rename  
        # matrix
        self.k = np.zeros(grid.neighbours.shape)  # contains precomputed  hydraulic conductivities [cm/day]
        self.alpha = np.zeros((grid.n_cells * grid.number_of_neighbours(),))
        self.beta_const = np.zeros((grid.n_cells,))  # const part of diagonal entries [1] 
        self.f_const = np.zeros((grid.n_cells,))  # const part of load vector [1]           
        # 
        self.bc = { }  # boundary conditions, map with key (cell_id, face_id) containing list of values
        self.sources = np.zeros((grid.n_cells,))  # [cm3 / cm3]
        self.h0 = np.zeros((grid.n_cells,))  # [cm]   
        self.innerFlux = 0.  # cm3/cm2/day 
        #
        self.gravitation = 0  # assure dimension equals grid.dim (todo)
#                 self.bc_f = { "flux": lambda: self.bc_flux, "flux_in": lambda: self.bc_flux_in, "flux_out":lambda: self.bc_flux_ot,
#                      "rootsystem": lambda: self.bc_rootsystem }        
        
    def create_k(self):
        """ sets up hydraulic conductivities from last time step solution @param h , for each neighbour"""
        cols = np.ones((1, self.grid.number_of_neighbours()))
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
        self.alpha = -v
        self.beta_const = np.sum(v, axis=1)
    
    def bc_flux(self, q):
        """ flux boundary condition """
        return q  # [cm3 / cm2 / day]   
    
    def bc_flux_out(self, i, q, crit_p, dx):
        """ outflow boundary condition limits to out or zero flux, and to critical pressure"""
        k = vg.hydraulic_conductivity(self.h0[i], self.soil)
        max_q = k * (crit_p - self.h0[i]) / dx  # maximal possible outflux 
        # print("bc_flux_out", q, max_q, i , crit_p, self.h0[i], k, dx)  
        return min(max(q, max_q), 0.)  # [cm3 / cm2 / day]   

    def bc_flux_in(self, i, q, crit_p, dx):
        """ inflow boundary condition limits to out or zero flux, and to critical pressure"""
        k = vg.hydraulic_conductivity(self.h0[i], self.soil)  # [cm / day]
        max_q = k * (crit_p - self.h0[i]) / dx  # maximal possible outflux [cm3 / cm2 /day]   
        # print("bc_flux_in", q, max_q, i, crit_p, self.h0[i], k, dx)
        return max(min(q, max_q), 0.)  # [cm3 / cm2 / day]    

    def bc_rootsystem(self, rx, kr):
        """ flux is given by radial conductivity times difference in matric potentials 
        @param rx      root xylem matric potential [cm]
        @param kr      root radial condcutivitiy [1 / day] (intrinsic)
        """
        h = self.getInnerHead()
        k = vg.hydraulic_conductivity(h, self.soil)  # [cm / day]              
        dx = self.grid.nodes[0]  # self.grid.mid[0]
        q = min(kr, k / dx) * (rx - h)  
        return min(q, 0.)  # limit to outflux  # [cm3 / cm2 / day] 

    def bc_to_source(self, dt):
        """ 
        Computes a cell source term from a boundary conditions given in self.bc, 
        self.bc is dictionary with keys (cell, face_id) and values is a boundary function returning a value [cm3 / cm2 / day].          
        
        Only one face per cell is allowed, overwrites any other sources
        
        dictionary with lambda functions would be nicer, but causes trouble with parallelisation    
        """                                     
        for (i, j), (type, p1, p2, p3) in self.bc.items(): 
            if type == "rootsystem":
                bc = self.bc_rootsystem(p1, p2)
            elif type == "flux_in":
                bc = self.bc_flux_in(i, p1, p2, p3)
            elif type == "flux_out":
                bc = self.bc_flux_out(i, p1, p2, p3)
            elif type == "flux":
                bc = self.bc_flux_out()
            else:
                raise("Unkown boundary condition")
            self.sources[i] = dt * bc * self.grid.area_per_volume[i, j]
    
    def getFlux(self, cell_id, face_id):
        """ flux [cm3/cm2/day] over the inner face given by @param face_id in cell @param cell_id """
        i = cell_id
        return self.k[i, face_id] * (self.h0[i] - self.h0[i + 1]) / dx  # [cm3/cm2/day]

    def getWaterContent(self):
        """ current water content for each cell """
        return vg.water_content(self.h0, self.soil)
    
    def getInnerHead(self):
        """ extrapolates result to root surface (specalized 1d) """        
        right_neighbour = self.grid.neighbours[0, 1]  
        dx0 = self.grid.mid[0]  #  (self.grid.mid[0] - self.grid.nodes[0])         
        dx1 = self.grid.mid[1] - self.grid.mid[0]  
        h0, h1 = self.h0[0], self.h0[1]
        h = h0 - ((h1 - h0) / dx1) * dx0    
        # print("linear head extrapolation", h0, h1, h, dx0, dx1, 1) 
        return h             
    
    def getInnerFlux(self):
        """ exact flux at cell_id = 0, todo generalize """
        return self.innerFlux

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
            A = sparse.coo_matrix((self.alpha.flat, (self.alpha_i.flat, self.alpha_j.flat)))
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
                sim_time = sim_time + dt_  # increase current time
                
                if output_times[k] <= sim_time:
                    h_out.append(h2.copy())
                    k = k + 1
                if verbose:
                    print('Time {:g} days, iterations {:g}, last time step {:g}'.format(sim_time, i, dt_))
                
                self.h0 = h2  

                self.bc_to_source(dt_)
                self.innerFlux += self.sources[0] / (self.grid.area_per_volume[0, 0])  # exact inner flux    
                
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

