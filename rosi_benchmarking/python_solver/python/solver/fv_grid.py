import numpy as np


class FVGrid:
    """ simplistic FV base class """

    def __init__(self):
        self.dim = 0  # coordiante dimensions
        self.number_of_neighbours = 0  # the maximal number of neighbouring cells

        self.nodes = None  # [cm] coordinates, np.array
        self.cells = None  # node indices, np.array(dtype=np.int64)
        self.n_nodes = 0  # number of nodes
        self.n_cells = 0  # number of cells

        self.neighbours = None  # cell indices, np.array((n_cells, number_of_neighbours), dtype=np.int64), -1 indicates no neigbour
        self.boundary_faces = []  # tuples defining the boundary faces of the grid (cell_id, face_id)
        self.area_per_volume = None  # [cm2/cm3] cell face area to neigbour cell divided by volume. shape of neigbours
        self.dx = None  # distances to neighbours. shape of neigbours

    def centers(self):
        """ cell centers as numpy array"""
        c = np.zeros((self.n_cells, self.dim))
        for i in range(0, self.n_cells):
            c[i, :] = self.center(i)
        return c

    def center(self, i):
        """ the cell center of the cell with index i """
        m = np.zeros((self.dim,))
        for j in self.cells[i, :]:
            m += self.nodes[j, :] / self.number_of_neighbours  # TODO works only for dim 1 ?!, how to (n,) -> (n,1)
        return m


class FVGrid1D(FVGrid):
    """ 1d grid """

    def __init__(self, nodes :np.array):
        super().__init__()

        self.dim = 1
        self.number_of_neighbours = 2

        if len(nodes.shape) == 1:  # (n,) -> (n,1)
            self.nodes = np.expand_dims(nodes, axis = 1)
        else:
            self.nodes = nodes
        n = nodes.shape[0] - 1  # number of cells
        self.cells = np.transpose(np.vstack((np.array(range(0, n), dtype = np.int64), np.array(range(1, n + 1), dtype = np.int64))))
        self.n_nodes = n + 1
        self.n_cells = n

        self.neighbours = np.transpose(np.vstack((np.array(range(-1, n - 1), dtype = np.int64), np.array(range(1, n + 1), dtype = np.int64))))
        self.neighbours[0, 0] = -1  # indicates there is no neighbour
        self.neighbours[n - 1, 1] = -1  # indicates there is no neighbou
        self.boundary_faces.append((0, 0))
        self.boundary_faces.append((n - 1, 1))

        cols = np.ones((1, self.number_of_neighbours))
        volume = self.nodes[self.cells[:, 1]] - self.nodes[self.cells[:, 0]]
        volumes = np.outer(volume, cols)
        area = np.ones(self.neighbours.shape)
        self.area_per_volume = np.divide(area, volumes)

        center = self.centers()
        self.dx = np.ones(self.neighbours.shape)
        for i in range(0, self.n_cells):  # TODO get rid of loops, use boundary_faces...
            for j, ni in enumerate(self.neighbours[i, :]):
                if ni >= 0:
                    self.dx[i, j] = np.linalg.norm(center[i] - center[ni])
                else:
                    self.dx[i, j] = np.linalg.norm(center[i] - self.nodes[i])  # half cell width


class FVGrid1Dcyl(FVGrid1D):
    """ 1d cylindrical grid """

    def __init__(self, nodes :np.array):
        super().__init__(nodes)
        volume = np.pi * (self.nodes[self.cells[:, 1]] ** 2 - self.nodes[self.cells[:, 0]] ** 2)  # correct volume
        volumes = np.outer(volume, np.ones((1, self.number_of_neighbours)))  # volume x [1,1]
        area = 2 * np.pi * self.nodes[self.cells]  # correct area
        self.area_per_volume = np.divide(area[:, :, 0], volumes)
