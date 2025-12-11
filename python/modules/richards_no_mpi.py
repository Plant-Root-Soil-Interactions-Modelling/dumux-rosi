from rosi.solverbase import SolverWrapper
from rosi.richards import RichardsWrapper

import numpy as np


class RichardsNoMPIWrapper(RichardsWrapper):
    """ 
    rewrites all methods using MPI to single process ones
    """

    def __init__(self, base):
        super().__init__(base)

    def initialize(self, args_ = [""], verbose = True, doMPI_ = False):
        """ Writes the Dumux welcome message, and creates the global Dumux parameter tree """
        self.base.initialize(args_, verbose, doMPI = doMPI_)

    def solve(self, dt:float, maxDt = -1.):
        """ Simulates the problem, the internal Dumux time step ddt is taken from the last time step 
        @param dt      time span [days] 
        @param mxDt    maximal time step [days] 
        """
        self.base.solveNoMPI(dt * 24.*3600., maxDt * 24.*3600.)  # days -> s

    def getAllNeumann(self, eqIdx = 0):
        """ Gathers the neuman fluxes into rank 0 as a map with global index as key [cm / day]"""
        dics = self.base.getAllNeumann(eqIdx)
        # for 1D
        for key, value in dics.items():
            dics[key] = value / 1000 * 24 * 3600 * 100.  # [kg m-2 s-1] / rho = [m s-1] -> cm / day
        return dics

        # flat_dic = {}
        # for d in dics:
        #     flat_dic.update(d)
        # for key, value in flat_dic:
        #     flat_dic[key] = value / 1000 * 24 * 3600 * 100.  # [kg m-2 s-1] / rho = [m s-1] -> cm / day
        return flat_dic

    def getNetFlux(self, eqIdx = 0):
        """ Gathers the net fluxes fir each cell into rank 0 as a map with global index as key [cm3 / day]"""
        self.checkGridInitialized()
        return self._map((self.base.getNetFlux(eqIdx)), 0) * 1000. *24 * 3600  # kg/s -> cm3/day

    def _map(self, x, type, dtype = np.float64):
        """Converts rows of x to numpy array and maps it to the right indices         
        @param type 0 dof indices, 1 point (vertex) indices, 2 cell (element) indices   
        """
        if type == 0:  # auto (dof)
            indices = (self.base.getDofIndices())
        elif type == 1:  # points
            indices = (self.base.getPointIndices())
        elif type == 2:  # cells
            indices = (self.base.getCellIndices())
        else:
            raise Exception('PySolverBase._map: type must be 0, 1, or 2.')
        if indices:  # only for rank 0 not empty
            assert len(indices) == len(x), "_map: indices and values have different length"
            ndof = max(indices) + 1
            if isinstance(x[0], list):
                m = len(x[0])
            else:
                m = 1
            p = np.zeros((ndof, m), dtype = dtype)
            for i in range(0, len(indices)):  #
                p[indices[i],:] = np.array(x[i], dtype = dtype)
            return p
        else:
            return 0

    def gather(self, data2gather, root = 0):
        """ dummy gather 
            to make the overloading of functions from richards and solverbase easier
            @see solverbase::gather()
        """
        return data2gather

    def _flat0(self, xx):
        """flattens the gathered list in rank 0, empty list for other ranks """
        if isinstance(xx[0], (list, type(np.array([])))):
            return np.array([item for sublist in xx for item in sublist])
        else:
            return xx

    def distributeSources(self, dt, source, inner_fluxs, eqIdx, plantM):
        """ when we have a source in for a 1d model (exchange with bulk soil),
            divid the source between the cell of the 1d model according to the water
            or solute content
            @param sources: sources to divide between the cells for water [list, cm3/day] and/or solutes [list, mol/day]
            @param inner_fluxs: negative or positive flux at the root surface for water [list, cm3] and/or solutes [list, mol]
            @param eqIdx: index of the components [list of int]
        """
        splitVals = list()
        # iterate through each component
        assert len(source) == len(inner_fluxs)
        for i, src in enumerate(source):  # [cm3/day] or [mol/day]
            splitVals.append(self.distributeSource(dt, src, inner_fluxs[i],
                                    eqIdx[i], plantM = plantM))
        return np.array(splitVals, dtype = object)

    def distributeSource(self, dt, source: float, inner_flux:float, eqIdx: int, plantM):
        """ when we have a source in for a 1d model (exchange with bulk soil),
            divid the source between the cell of the 1d model according to the water
            or solute content
            @param sources: sources to divide between the cells for water [cm3/day] or solute [mol/day]
            @param inner_flux: negative or positive flux at the root surface for water [cm3] or solute [mol]
            @param eqIdx: index of the components [int]
        """
        assert self.dimWorld != 3
        assert eqIdx < self.numFluidComp

        # distribute the value between the cells
        splitVals = self.distributeVals(dt, source, inner_flux, eqIdx, plantM)
        # send the source data to dumux
        values = list(splitVals.copy())
        res = {i: values[i] for i in range(len(values))}
        self.setSource(res.copy(), eq_idx = eqIdx)  # [mol/day], in modules/richards.py

        return splitVals

    def distributeVals(self, dt, source: float, inner_flux: float, eqIdx: int, plantM):
        """ when we have a source in for a 1d model (exchange with bulk soil),
            divid the source between the cell of the 1d model according to the water
            or solute content
            @param sources: sources to divide between the cells for water [cm3/day] or solute [mol/day]
            @param inner_flux: negative or positive flux at the root surface for water [cm3/day] or solute [mol/day]
            @param eqIdx: index of the components [int]
        """
        # will this still work even if theta < wilting point during solve() calls?

        splitVals = np.array([0. for i in range(self.numberOfCellsTot)])

        if source != 0.:  # [cm3/day] or [mol/day]

            if eqIdx == 0:  # compute the amount of water potentially available in each cell
                splitVals = self.distributeValWater(dt, source, inner_flux, plantM)
            else:
                splitVals = self.distributeValSolute(eqIdx, dt, source, inner_flux, plantM)
            source_ = sum(splitVals)
            try:
                assert (((splitVals >= 0).all()) or ((splitVals <= 0).all()))
                assert abs(sum(abs(splitVals)) - abs(source_)) < 1e-13
                assert len(splitVals) == self.numberOfCellsTot
            except:
                print('splitVals', splitVals)
                print('troubleshoot data', 'source=', source, 'source_=', source_,
                      ';dt=', dt, ';inner_flux=', inner_flux,
                ';theta=', repr(self.getWaterContent()), ';cylVol=', repr(self.getCellVolumes()),
                [self.vg_soils[0].theta_R, self.vg_soils[0].theta_S,
                self.vg_soils[0].alpha, self.vg_soils[0].n, self.vg_soils[0].Ksat], self.theta_wilting_point)
                raise Exception

        return splitVals

    def distributeValSolute(self, eqIdx, dt, source: float, inner_flux: float, plantM):
        seg_values_content = np.maximum(self.getContent(eqIdx), 0.)  # during the solve() loop, we might still get seg_values <0

        seg_values_content[0] += inner_flux  # add root solute release or uptake
        cylVol = self.getCellVolumes()

        return np.array(plantM.distributeValSolute_(seg_values_content.copy(), cylVol.copy(), source, dt))  # mol/day

    def distributeValWater(self, dt, source: float, inner_flux: float, plantM):
        cylVol = self.getCellVolumes()
        seg_values_perVol_ = np.maximum(self.getWaterContent(), 0.)  # cm3/cm3

        seg_values_perVol_[0] += inner_flux / cylVol[0]  # add root water release or uptake

        return np.array(plantM.distributeValWater_(seg_values_perVol_.copy(), cylVol.copy(), source, dt,
                                            self.vg_soils[0].theta_S, self.theta_wilting_point))
