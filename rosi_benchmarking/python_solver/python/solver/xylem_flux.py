import timeit
import math

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA

from solver.plantbox import MappedRootSystem
from solver.plantbox import XylemFlux
import solver.rsml_reader as rsml  # todo


class XylemFluxPython(XylemFlux):
    """  Hybrid flux solver (following Meunier et al)"""

    def __init__(self, rs):
        super().__init__(rs)
        self.kx_ = None
        self.kx_t_ = None
        self.kr_ = None
        self.kr_t_ = None
        self.kx = None
        self.kr = None

    def solve(self, value, neumann = True) :
        """ solves the flux equations, with neumann or dirichlet boundary condtion,
            @param value    [cm day-1] or [cm] pressure head
            @parm neumann   Neumann or Dirichlet
         """
        start = timeit.default_timer()

        # self.linearSystem()
        I, J, V, b = self.linear_system()
        Q = sparse.coo_matrix((V, (I, J)))
        # Q = sparse.coo_matrix((self.aV, (self.aI, self.aJ)))
        Q = sparse.csr_matrix(Q)

        if neumann:
            self.aB = b
            print("neumann")
            Q, b = self.bc_neumann(Q, self.aB, [0], [value])
        else:
            self.aB = b
            print("dirichlet")
            Q, b = self.bc_dirichlet(Q, self.aB, [0], [value])

        x = LA.spsolve(Q, b, use_umfpack = True)  # direct
        print(value)
        print(x[0])
        print(Q)

        print ("linear system solved in", timeit.default_timer() - start, " s")
        return x

    def solve_wp(self, trans, wiltingPoint = -15000):
        """ solves the flux equations using neumann and switching to dirichlet in case, 
            (todo assembling once should be enough) 
            @param trans        [cm day-1]
            @parm wiltingPoint  wiltingPoint            
        """
        try:
            x = solve(trans)
        except:
            x = [-15001]

        if x[0] < -15000:
            x = solve (wiltingPoint, False)

        return x

    def getSolution_(self, rx_, sx_):
        """ creates the inhomogenneous solution from the homogeneous one"""
        cIdx = self.rs.seg2cell[0]
        rx_[0] += sx_[cIdx]
        for i in range(1, len(rx_)):
            cIdx = self.rs.seg2cell[i - 1]
            rx_[i] += sx_[cIdx]
        return rx_

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

    def setKr(self, values :np.array, age :np.array = np.array([])):
        """ Sets the radial conductivity in [1 day-1], 
            converts to [cm2 day g-1] by dividing by rho*g """
        self.kr_ = values / (self.rho * self.g)
        print("Kr", values, "day-1 -> ", self.kr_, "cm2 day g-1")
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
        self.setKrF(self.kr)  # to the cpp

    def setKx(self, values :np.array, age :np.array = np.array([])):
        """ Sets the axial conductivity [cm3 day-1], 
            converts to [cm5 day g-1] by dividing by rho*g """
        self.kx_ = values / (self.rho * self.g)
        print("Kx", values, "day-1 -> ", self.kx_, "cm5 day g-1")
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
        self.setKxF(self.kx)  # to the cpp

    def linear_system(self):
        """ assembles the linear system,
         (todo -> c++) """
        simTime = self.rs.getSimTime()  # to calculate age from ct
        Ns = len(self.rs.segments)
        N = len(self.rs.nodes)

        I = np.zeros(4 * Ns, dtype = np.int64)
        J = np.zeros(4 * Ns, dtype = np.int64)
        V = np.zeros(4 * Ns)
        b = np.zeros(N)
        k = 0
        for s in range(0, Ns):

            i, j = self.rs.segments[s].x, self.rs.segments[s].y

            age = simTime - self.rs.nodeCTs[j]
            type = self.rs.types[j - 1]
            a = self.rs.radii[j - 1]

            kx = self.kx(age, type)  # j-1 = segment index
            kr = self.kr(age, type)

            n1 = np.array([self.rs.nodes[i].x, self.rs.nodes[i].y, self.rs.nodes[i].z])
            n2 = np.array([self.rs.nodes[j].x, self.rs.nodes[j].y, self.rs.nodes[j].z])
            v = n2 - n1
            l = np.linalg.norm(v)
            vz = v[2] / l  # normed direction

            c = 2.*a * math.pi * kr / kx  # Eqn (2)
            d = math.exp(-math.sqrt(c) * l) - math.exp(math.sqrt(c) * l)  # Eqn (5)
            # print(i, j, n1, n2, l, d)
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

        return I, J, V, b

