import numpy as np


class FV_Grid: 
    """ abstract class """
    
    def __init__(self):
            
        self.nodes = None  # coordinates, np.array
        self.cells = None  # node indices, np.array(dtype=np.int64)
        self.neighbours = None  # cell indices, np.array(dtype=np.int64)
        self.dim = 0  # coordiante dimensions
        self.n_nodes = 0  # number of nodes
        self.n_cells = 0  # number of cells
        
        # precompute
        self.area = None  # cell face area to neigbour cell. np.array, shape of neigbours 
        self.volume = None  # cell volume    
        self.mid = None  # cell mid points

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
        self.neighbours[n - 1, 1] = -1  #  -1 means no neighbour
        self.area = np.ones(self.neighbours.shape)
        self.volume = self.nodes[self.cells[:, 1]] - self.nodes[self.cells[:, 0]]        
        self.mid = 0.5 * (self.nodes[self.cells[:, 0]] + self.nodes[self.cells[:, 1]])

    def number_of_neighbours(self):
        """ returns the maximal number of neighbouring cells sharing a common face (normally const) """
        return 2


class FV_Grid1Dcyl(FV_Grid1D):
    """ 1d cylindrical grid """
    
    def __init__(self, nodes :np.array):
        super().__init__(nodes)
        self.area = 2 * np.pi * np.multiply(self.area, self.nodes[self.cells])  # correct area
        self.volume = np.pi * (self.nodes[self.cells[:, 1]] ** 2 - self.nodes[self.cells[:, 0]] ** 2)  # correct volume  

