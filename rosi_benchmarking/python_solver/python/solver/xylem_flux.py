from math import *
import numpy as np
from numpy.linalg.linalg import norm
import matplotlib.pylab as plt
from scipy import sparse
import scipy.sparse.linalg as LA

import solver.rsml_reader as rsml
import solver.plantbox as pb
from solver.richards import RichardsWrapper


class XylemFlux:
    """  Hybrid flux solver (following Meunier et al)"""

    def __init__(self, soil :RichardsWrapper):
        self.rs = None
        self.nodes = None
        self.seg = None
        self.radii = None
        self.node_cts = None
        self.types = None
        self.node2cell = {}
        self.cell2node = {}
        self.soil = soil
        self.soilP = None

#     def openRSML(self, file_name):
#         polylines, props, funcs = rsml.read_rsml(file_name)
#         nodes, segs = rsml.get_segments(polylines, props)
#         self.nodes = np.array(nodes) * 1.e-3  # [mm] -> [m], and convert from list to numpy array
#         self.segs = np.array(segs, dtype = np.int64)  # convert from list to numpy array
#         radii, cts, types = rsml.get_parameter(polylines, funcs, props)
#         self.radiir = np.array(radii) * 1.e-3  # [mm]->[m]
#         self.cts = np.array(cts)  # [days]

    def readParameters(self, file_name):
        self.rs = pb.RootSystem()
        self.rs.readParameters(file_name)
        self.rs.initialize()
        self.nodes = self.convert(self.rs.getNodes())  # contains the initial nodes of tap, basal and shootborne roots
        self.seg = self.convert(self.rs.getShootSegments(), np.int64)
        self.node_cts = np.zeros(self.nodes.shape[0])
        self.radii = 0.2 * np.ones(self.seg.shape[0])
        self.types = np.zeros(self.seg.shape[0])
        self.update_maps(self.seg)

    def simulate(self, dt):
        """ simulates the root system growth and registers new segments """

        assert self.rs != None, "simulate: Initialize a root system first with readParameters "

        # simulate
        self.rs.simulate(dt, True)

        # move nodes
        uni = np.array(self.rs.getUpdatedNodeIndices(), dtype = np.int64)
        unodes = self.convert(self.rs.getUpdatedNodes())
        if len(uni) != 0:
            nodes[uni] = unodes  # do the update
        # add nodes
        newnodes = self.convert(self.rs.getNewNodes())
        newnode_cts = self.convert(self.rs.getNewNodeCTs())
        if len(newnodes) != 0:
            self.nodes = np.vstack((self.nodes, newnodes))
            self.node_cts = np.hstack((self.node_cts, newnode_cts))
        # add segments
        newsegs = self.convert(self.rs.getNewSegments(), np.int64)
        if len(newsegs) != 0:
            self.seg = np.vstack((self.seg, newsegs))
        # add radii and type
        new_radii = []
        new_types = []
        newseg_origins = self.rs.getNewSegmentOrigins()
        for o in newseg_origins:
            new_radii.append(o.getParam().a)
            new_types.append(o.getParam().subType)
        if len(new_radii) != 0:
            self.radii = np.hstack((self.radii, new_radii))
            self.types = np.hstack((self.types, new_types))

        self.update_maps(newsegs)  # ALL (do incremential)

    def update_maps(self, seg):
        """ maps node indices to soil cell indices, and vice versa """
        for s in seg:
            nodeIdx = s[1]
            node = self.nodes[s[1], :]
            cellIdx = self.soil.pickCell(list(node))
            assert cellIdx > -1, "updateMaps: root node " + str(node) + " i out of the soil domain"
            self.node2cell[nodeIdx] = cellIdx
            self.cell2node[cellIdx] = nodeIdx

    def eval_fluxes(self, soil, root):
        """ evaluate fluxes per cellId, where @param soil and @param root are the solutions for the last time step """
        fluxes = {}
        for i, s in enumerate(self.seg):
            nodeIdx = s[1]
            cellIdx = self.node2cell(s[1])
            qr = -2 * pi * self.radii[i] * kr[i] * (soil[cellId] - root[nodeIdx])
            fluxes[cellIdx] = qr
        return fluxes

    def solve(self):
        """ solves the flux equations """
        Q, b = linar_system_()

        pass

    def linear_system_(self, soil):

        Ns = self.seg.shape[0]
        N = self.nodes.shape[0]

        # TODO move I,J,V loop to C++
        I = np.zeros(4 * Ns, dtype = np.int64)
        J = np.zeros(4 * Ns, dtype = np.int64)
        V = np.zeros(4 * Ns)
        b = np.zeros(N)
        k = 0

        for s in range(0, Ns):

            i, j = seg[s, 0], seg[s, 1]
            n1, n2 = nodes[i, :], nodes[j, :]
            p_s = soil(self.seg2cell[s])
            v = n2 - n1
            l = norm(v)
            vz = v[2] / l  # normed direction
            a = radius[s]

            c = 2.*a * pi * kr[s] / kz[s]  # Eqn (2)
            d = exp(-sqrt(c) * l) - exp(sqrt(c) * l)  # Eqn (5)
            di = 1. / d

            cii = -kz[s] * di * sqrt(c) * (exp(-sqrt(c) * l) + exp(sqrt(c) * l))  # Eqn 16
            cij = 2 * kz[s] * di * sqrt(c)  # Eqn 17
            bi = kz[s] * rho * g * vz  # Eqn 18

            b[i] += bi
            I[k], J[k], V[k] = i, i, cii
            k += 1
            I[k], J[k], V[k] = i, j, cij
            k += 1

            # edge ji
            i, j = j, i
            b[i] -= bi  # Eqn 14 with changed sign
            I[k], J[k], V[k] = i, i, cii
            k += 1
            I[k], J[k], V[k] = i, j, cij
            k += 1

        Q = sparse.coo_matrix((V, (I, J)))
        Q = sparse.csr_matrix(Q)

        return (Q, b)

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
    def convert(x, dtype = np.float64):
        return np.array(list(map(lambda x: np.array(x, dtype), x)), dtype)  # is there a better way?
