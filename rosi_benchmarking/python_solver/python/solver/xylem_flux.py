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

    def __init__(self, rs :pb.MappedRootSystem):  # remove type later for richardsnc stuff
        self.rs = rs  #
        self.kx_ = None
        self.kx_t_ = None
        self.kr_ = None
        self.kr_t_ = None
        self.kx = None
        self.kr = None

        # units:
        self.rho = 1.  # [cm3/g]
        self.g = 9.8065 * 100. / (24. *3600.*24. *3600.)  # [cm/day^2]

    def solve(self, trans, neumann = True) :
        """ solves the flux equations, with neumann or dirichlet boundary condtion,
            @param trans    [cm day-1]
            @parm neumann   Neumann or Dirichlet
         """
        start = timeit.default_timer()
        Q, b = self.linear_system_()

        if neumann:
            Q, b = self.bc_neumann(Q, b, np.array([0]), np.array([trans]))
        else:
            Q, b = self.bc_dirichlet(Q, b, np.array([0]), np.array([trans]))

        # print ("assembled in", timeit.default_timer() - start, " s")
        # start = timeit.default_timer()
        x = LA.spsolve(Q, b, use_umfpack = True)  # direct
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

    def linear_system_(self):
        """ assembles the linear system,
         (todo -> c++) """
        simTime = self.rs.getSimTime()  # to calculate age from ct
        Ns = self.rs.segments.shape[0]
        N = self.rs.nodes.shape[0]

        I = np.zeros(4 * Ns, dtype = np.int64)
        J = np.zeros(4 * Ns, dtype = np.int64)
        V = np.zeros(4 * Ns)
        b = np.zeros(N)
        k = 0
        for s in range(0, Ns):

            i, j = self.rs.segments[s, 0], self.rs.segments[s, 1]

            kx = rs.kx(j - 1)  # j-1 = segment index
            kr = self.kr(j - 1)
            a = self.rs.radii[j - 1]

            n1, n2 = self.rs.nodes[i, :], self.rs.nodes[j, :]
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

    def getSolution(self, x_):
        """ creates the inhomogenneous solution from the homogeneous one"""
        x_[0] += rs.soilPressure(0)
        for i in range(1, len(x_)):
            x_[i] += rs.soilPressure(cIdx - 1)
        return x_

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

