import timeit
import math
import numpy as np
import matplotlib.pylab as plt
from scipy import sparse
import scipy.sparse.linalg as LA

import solver.rsml_reader as rsml
import solver.plantbox as pb
from solver.richards import RichardsWrapper


class XylemFlux:
    """  Hybrid flux solver (following Meunier et al)"""

    def __init__(self, soil :RichardsWrapper):  # remove type later for richardsnc stuff
        self.soil = soil
        self.rs = None
        self.nodes = None
        self.seg = None
        self.node_cts = None
        self.radii = None
        self.types = None
        self.node2cell = None
        self.cell2node = None
        self.kx_ = None
        self.kx_t_ = None
        self.kr_ = None
        self.kr_t_ = None
        self.kx = None
        self.kr = None

        # units:
        self.rho = 1.  # [cm3/g]
        self.g = 9.8065 * 100. / (24. *3600.*24. *3600.)  # [cm/day^2]

#     def openRSML(self, file_name):
#         polylines, props, funcs = rsml.read_rsml(file_name)
#         nodes, segs = rsml.get_segments(polylines, props)
#         self.nodes = np.array(nodes) * 1.e-3  # [mm] -> [m], and convert_ from list to numpy array
#         self.segs = np.array(segs, dtype = np.int64)  # convert_ from list to numpy array
#         radii, cts, types = rsml.get_parameter(polylines, funcs, props)
#         self.radiir = np.array(radii) * 1.e-3  # [mm]->[m]
#         self.cts = np.array(cts)  # [days]

    def readParameters(self, file_name):
        """ opens a .xml parameter file, and initializes the root system """
        self.rs = pb.RootSystem()
        self.rs.readParameters(file_name)
        self.rs.initialize()
        self.nodes = self.convert_(self.rs.getNodes())  # contains the initial nodes of tap, basal and shootborne roots
        self.seg = self.convert_(self.rs.getShootSegments(), np.int64)
        self.node_cts = np.zeros(self.nodes.shape[0])
        self.radii = 0.2 * np.ones(self.seg.shape[0])
        self.types = np.zeros(self.seg.shape[0], np.int64)
        self.node2cell = {}  # reset the mappers
        self.cell2node = {}
        self.update_maps(self.seg)

    def simulate(self, dt):
        """ simulates the root system growth and maps the new segments to the soil cells"""

        assert self.rs != None, "simulate: Initialize a root system first with readParameters "

        # simulate
        self.rs.simulate(dt, True)

        # move nodes
        uni = np.array(self.rs.getUpdatedNodeIndices(), dtype = np.int64)
        unodes = self.convert_(self.rs.getUpdatedNodes())
        if len(uni) != 0:
            self.nodes[uni] = unodes  # do the update
        # add nodes
        newnodes = self.convert_(self.rs.getNewNodes())
        newnode_cts = self.convert_(self.rs.getNewNodeCTs())
        if len(newnodes) != 0:
            self.nodes = np.vstack((self.nodes, newnodes))
            self.node_cts = np.hstack((self.node_cts, newnode_cts))
        # add segments
        newsegs = self.convert_(self.rs.getNewSegments(), np.int64)
        if len(newsegs) != 0:
            self.seg = np.vstack((self.seg, newsegs))
        # add radii and type
        newseg_origins = self.rs.getNewSegmentOrigins()  # todo improve (no append for new_...)
        new_radii = np.array([])
        new_types = np.array([], np.int64)
        for o in newseg_origins:
            new_radii = np.append(new_radii, o.getParam().a)
            new_types = np.append(new_types, o.getParam().subType)
        if len(new_radii) != 0:
            self.radii = np.append(self.radii, new_radii)
            self.types = np.append(self.types, new_types)

        self.update_maps(newsegs)

    def update_maps(self, seg):
        """ maps node indices to soil cell indices, and vice versa """
        for s in seg:
            nodeIdx = s[1]
            node = self.nodes[s[1], :]
            cellIdx = self.soil.pickCell(list(node))
            # print("node idx", nodeIdx, " mapped to cell idx", cellIdx)
            assert cellIdx > -1, "updateMaps: root node " + str(node) + " i out of the soil domain"
            self.node2cell[nodeIdx] = cellIdx
            self.cell2node[cellIdx] = nodeIdx

    def eval_fluxes(self, soil, root):
        """ evaluate fluxes per segment index (= nodeIdx-1) [cm3 day-1], 
            where @param soil and @param root are the solutions for the last time step, 
            (todo -> c++) """
        simTime = self.rs.getSimTime()  # CARE that we use the right time (todo)
        fluxes = { }
        for i, s in enumerate(self.seg):
            nodeIdx = s[1]
            cellIdx = self.node2cell[s[1]]
            age = simTime - self.node_cts[nodeIdx]
            type = self.types[i]
            kr = float(self.kr(age, type))
            l = np.linalg.norm(self.nodes[s[1], :] - self.nodes[s[0], :])
            surf = float(2 * math.pi * self.radii[i] * l)
            # print(soil[cellIdx], root[nodeIdx], kr, surf)
            qr = -surf * kr * float(soil[cellIdx] - root[nodeIdx])
            if cellIdx in fluxes:
                 fluxes[cellIdx] += qr
            else:
                fluxes[cellIdx] = qr
        return fluxes

    def solve(self, trans, soilP, neumann = True) :
        """ solves the flux equations, with neumann or dirichlet boundary condtion,
            @param trans    [cm day-1]
            @param soilP    [cm]
            @parm neumann   Neumann or Dirichlet
         """
        start = timeit.default_timer()
        Q, b = self.linear_system_(soilP)

        if neumann:
            Q, b = self.bc_neumann(Q, b, np.array([0]), np.array([trans]))
        else:
            Q, b = self.bc_dirichlet(Q, b, np.array([0]), np.array([trans]))

        # print ("assembled in", timeit.default_timer() - start, " s")
        # start = timeit.default_timer()
        x = LA.spsolve(Q, b, use_umfpack = True)  # direct
        print ("linear system solved in", timeit.default_timer() - start, " s")

        return x

    def solve_wp(self, trans, soilP, wiltingPoint = -15000):
        """ solves the flux equations using neumann and switching to dirichlet in case, 
            (todo assembling once should be enough) 
            @param trans        [cm day-1]
            @param soilP        [cm]
            @parm wiltingPoint  wiltingPoint            
        """
        try:
            x = solve(trans, soilP)
        except:
            x = [-15001]

        if x[0] < -15000:
            x = solve (wiltingPoint, soilP, False)

        return x

    def linear_system_(self, soilP):
        """ assembles the linear system,
         (todo -> c++) """
        simTime = self.rs.getSimTime()  # to calculate age from ct
        Ns = self.seg.shape[0]
        N = self.nodes.shape[0]

        I = np.zeros(4 * Ns, dtype = np.int64)
        J = np.zeros(4 * Ns, dtype = np.int64)
        V = np.zeros(4 * Ns)
        b = np.zeros(N)
        k = 0
        for s in range(0, Ns):

            i, j = self.seg[s, 0], self.seg[s, 1]

            a = self.radii[s]
            age = -self.node_cts[j]
            type = self.types[s]
            p_s = soilP[self.node2cell[s + 1]]  # node index is segment index +1

            kx = self.kx(age, type)
            kr = self.kr(age, type)

            n1, n2 = self.nodes[i, :], self.nodes[j, :]
            v = n2 - n1
            l = np.linalg.norm(v)
            vz = v[2] / l  # normed direction

            c = 2.*a * math.pi * kr / kx  # Eqn (2)
            d = math.exp(-math.sqrt(c) * l) - math.exp(math.sqrt(c) * l)  # Eqn (5)
            di = 1. / d

            cii = -kx * di * math.sqrt(c) * (math.exp(-math.sqrt(c) * l) + math.exp(math.sqrt(c) * l))  # Eqn 16
            cij = 2 * kx * di * math.sqrt(c)  # Eqn 17
            bi = kx * self.rho * self.g * vz  # Eqn 18

            b[i] += bi
            I[k], J[k], V[k] = i, i, cii
            k += 1
            I[k], J[k], V[k] = i, j, cij
            k += 1

            i, j = j, i  # edge ji
            b[i] -= bi  # Eqn 14 with changed sign
            I[k], J[k], V[k] = i, i, cii
            k += 1
            I[k], J[k], V[k] = i, j, cij
            k += 1

        Q = sparse.coo_matrix((V, (I, J)))
        Q = sparse.csr_matrix(Q)
        return (Q, b)

    def setKr(self, values :np.array, age :np.array = np.array([])):
        """ Sets the radial conductivity in [1 day-1], 
            converts to [cm2 day g-1] by dividing by rho*g """
        self.kr_ = values / (self.rho * self.g)
        if age.size == 0:
            if values.shape[0] == 1:  # constant
                self.kr = lambda age, type: self.kr_
            else:  # constant per type
                self.kr = lambda age, type: self.kr_[type]
        else:
            self.kx_t_ = age
            if values.shape[0] == 1:  # age dependent
                self.kr = lambda age, type: np.interp(age, self.kr_t_, self.kr_)
            else:  # table per type
                self.kr = lambda age, type: np.interp(age, self.kr_t_[type], self.kr_[type])

    def setKx(self, values :np.array, age :np.array = np.array([])):
        """ Sets the axial conductivity [cm3 day-1], 
            converts to [cm5 day g-1] by dividing by rho*g """
        self.kx_ = values / (self.rho * self.g)
        if age.size == 0:
            if values.shape[0] == 1:  # constant
                self.kx = lambda age, type: self.kx_
            else:  # constant per type
                self.kx = lambda age, type: self.kx_[type]
        else:
            self.kx_t_ = age
            if values.shape[0] == 1:  # age dependent
                self.kx = lambda age, type: np.interp(age, self.kx_t_, self.kx_)
            else:  # table per type
                self.kx = lambda age, type: np.interp(age, self.kx_t_[type], self.kx_[type])

    @staticmethod
    def bc_dirichlet(Q, b, n0, d):
        c = 0
        for c in range(0, len(n0)):
            i = n0[c]
            e0 = np.zeros((1, Q.shape[1]))  # build zero vector
            Q[i, :] = sparse.csr_matrix(e0)  # replace row i with ei
            Q[i, i] = 1
            b[i] = d[c]
        return Q, b

    @staticmethod
    def bc_neumann(Q, b, n0, f):
        c = 0
        for c in range(0, len(n0)):
            i = n0[c]  # print("Neumann BC at node "+str(i))
            b[i] += f[c]
        return Q, b

    @staticmethod
    def convert_(x, dtype = np.float64):
        return np.array(list(map(lambda x: np.array(x, dtype), x)), dtype)  # is there a better way?
