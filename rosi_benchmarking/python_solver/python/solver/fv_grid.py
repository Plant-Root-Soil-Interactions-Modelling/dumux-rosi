import numpy as np


class FV_Grid: 
    """ abstract class """
    
    def __init__(self):
            
        self.dim = 0  # coordiante dimensions

        self.nodes = None  # [cm] coordinates, np.array 
        self.cells = None  # node indices, np.array(dtype=np.int64)
        self.n_nodes = 0  # number of nodes 
        self.n_cells = 0  # number of cells
    
        self.neighbours = None  # cell indices, np.array(dtype=np.int64)        
        self.area_per_volume = None  # [cm2/cm3] cell face area to neigbour cell divided by volume. shape of neigbours 
        self.mid = None  # cell mid points [cm]
        self.dx = None  # distances to neighbours. shape of neigbours

    def number_of_neighbours(self):
        """ returns the maximal number of neighbouring cells sharing a commonface (normally const) """
        return 0


class FV_Grid1D(FV_Grid):
    """ 1d grid """
    
    def __init__(self, nodes :np.array):
        n = nodes.shape[0] - 1  # number of cells
        self.n_nodes = n + 1
        self.n_cells = n
        self.dim = 1
                
        self.nodes = nodes        
        self.cells = np.transpose(np.vstack((np.array(range(0, n), dtype=np.int64), np.array(range(1, n + 1), dtype=np.int64))))
        self.neighbours = np.transpose(np.vstack((np.array(range(-1, n - 1), dtype=np.int64), np.array(range(1, n + 1), dtype=np.int64))))
        self.neighbours[0, 0] = 0              
        self.neighbours[n - 1, 1] = n - 1      
           
        cols = np.ones((1, self.number_of_neighbours()))           
        volume = self.nodes[self.cells[:, 1]] - self.nodes[self.cells[:, 0]]                
        volumes = np.outer(volume, cols)        
        area = np.ones(self.neighbours.shape)
        self.area_per_volume = np.divide(area, volumes)
                         
        self.mid = 0.5 * (self.nodes[self.cells[:, 0]] + self.nodes[self.cells[:, 1]])
        
        self.dx = np.ones(self.neighbours.shape)
        for i in range(0, self.n_cells):                        
            for j, ni in enumerate(self.neighbours[i, :]):                        
                if i != ni:
                    self.dx[i, j] = np.linalg.norm(self.mid[i] - self.mid[ni]) 
                else:  
                    self.dx[i, j] = np.linalg.norm(self.mid[i] - self.nodes[i])  # half cell width    
                    
    def number_of_neighbours(self):
        """ returns the maximal number of neighbouring cells sharing a common face (normally const) """
        return 2


class FV_Grid1Dcyl(FV_Grid1D):
    """ 1d cylindrical grid """
    
    def __init__(self, nodes :np.array):
        super().__init__(nodes)
        cols = np.ones((1, self.number_of_neighbours()))
        volume = np.pi * (self.nodes[self.cells[:, 1]] ** 2 - self.nodes[self.cells[:, 0]] ** 2)  # correct volume  
        volumes = np.outer(volume, cols) 
        area = 2 * np.pi * self.nodes[self.cells]  # correct area
        self.area_per_volume = np.divide(area, volumes)        

