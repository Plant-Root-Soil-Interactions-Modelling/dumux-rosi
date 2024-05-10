from solverbase import SolverWrapper
from richards import RichardsWrapper

from mpi4py import MPI; 
comm = MPI.COMM_WORLD; size = comm.Get_size(); rank = comm.Get_rank()
max_rank = comm.Get_size()

import numpy as np


class RichardsNoMPIWrapper(RichardsWrapper):
    """ 
    rewrites all methods using MPI to single process ones
    """

    def __init__(self, base, usemoles):
        super().__init__(base, usemoles)
        self.useMoles = usemoles

    def initialize(self, args_ = [""], verbose = False, doMPI = False):
        """ Writes the Dumux welcome message, and creates the global Dumux parameter tree """
        self.base.initialize(args_, verbose, doMPI)
        
        
    def initializeProblem(self, maxDt = -1.):
        """ After the grid is created, the problem can be initialized """
        self.base.initializeProblem(maxDt * 24.*3600.)
        
        self.dofIndices_   = self.base.getDofIndices()
        self.pointIndices_ = self.base.getPointIndices()
        self.cellIndices_  = self.base.getCellIndices()
        
        self.dofIndices   = self.dofIndices_
        self.pointIndices = self.pointIndices_
        self.cellIndices  = self.cellIndices_
        
        self.CellVolumes_ = np.array( self.base.getCellVolumes()).reshape(-1) * 1.e6  # m2 -> cm2
        self.CellVolumes = self._map(self.CellVolumes_ , 0)  # m2 -> cm2
        
    def solve(self, dt:float, saveInnerDumuxValues_ = True):
        """ Simulates the problem, the internal Dumux time step ddt is taken from the last time step 
        @param dt      time span [days] 
        @param mxDt    maximal time step [days] 
        """                            
        self.base.solveNoMPI(dt * 24.*3600., saveInnerDumuxValues = saveInnerDumuxValues_)  # days -> s

    def getSolutionHead(self, eqIdx=0):
        """Gathers the current solution into rank 0, and converts it into a numpy array (Ndof, neq), 
        model dependent units, [Pa, ...]"""
        self.checkGridInitialized()
        return self._map((self.base.getSolutionHead(eqIdx)), 0)        
    
    def getKrw(self):
        """Gathers the current solution's saturation into rank 0, and converts it into a numpy array (Nc, 1) [1]"""
        self.checkGridInitialized()
        return self._map((self.base.getKrw()), 0)

    def getContentCyl(self,idComp, isDissolved,gId = None ):
        return self.getContent(idComp, isDissolved)

    def getWaterContent(self):
        """Gathers the current solution's saturation into rank 0, and converts it into a numpy array (Nc, 1) [1]"""
        self.checkGridInitialized()
        return self._map((self.base.getWaterContent()), 2)
    
    def getWaterContent_(self):
        return self.getWaterContent()
    
    def getPoints(self):
        """Gathers vertices into rank 0, and converts it into numpy array (Np, 3) [cm]"""
        self.checkGridInitialized()
        return self._map((self.base.getPoints()), 1) * 100.  # m -> cm

    def getCellCenters(self):
        """Gathers cell centers into rank 0, and converts it into numpy array (Nc, 3) [cm]"""
        self.checkGridInitialized()
        return self._map((self.base.getCellCenters()), 2) * 100.  # m -> cm

    def getDofCoordinates(self):
        """Gathers dof coorinates into rank 0, and converts it into numpy array (Ndof, 3) [cm]"""
        self.checkGridInitialized()
        return self._map((self.base.getDofCoordinates()), 0) * 100.  # m -> cm

    def getCells(self):
        """ Gathers dune elements (vtk cells) as list of list of vertex indices (vtk points) (Nc, Number of corners per cell) [1]"""
        return self._map((self.base.getCells()), 2, np.int64)

    def getCellVolumesCyl(self):
        """ Gathers element volumes (Nc, 1) [cm3] """
        return self.CellVolumes

    def getCellSurfacesCyl(self):
        """ Gathers element volumes (Nc, 1) [cm3] """
        return self._map((self.base.getCellSurfacesCyl()), 2) * 1.e4  # m3 -> cm3

    def getCellCenters(self):
        """Gathers cell centers into rank 0, and converts it into numpy array (Nc, 3) [cm]"""
        self.checkGridInitialized()
        return self._map((self.base.getCellCenters()), 2) * 100.  # m -> cm

    def getSolution(self, eqIdx=0):
        """Gathers the current solution into rank 0, and converts it into a numpy array (dof, neq), 
        model dependent units [Pa, ...]"""
        self.checkGridInitialized()
        return self._map((self.base.getSolution(eqIdx)), 0)
    
    
    def getSavedBC(self,  rIn, rOut ):  
        
        """returns the total boundary flow simulated during the last dumux solumation
            @param rIn: position of the inner boundary == root radius [cm]
            @param rOut: position of the outer boundary of the perihirzal zone [cm]
        """
        assert self.dimWorld != 3
        length = self.segLength
        
        # values saved for each dumux sub-timestep
        BC_in_vals = np.array(self.base.BC_in_vals) #[ mol / (m^2 \cdot s)]
        BC_out_vals = np.array(self.base.BC_out_vals) #[ mol / (m^2 \cdot s)]
        BC_ddt = np.array(self.base.BC_ddt)# s
        # total flows
        q_out_vals = BC_out_vals* BC_ddt[:, None] #mol / m^2
        q_in_vals = BC_in_vals* BC_ddt[:, None] #mol / m^2        
        q_out = np.sum(q_out_vals,axis = 0) #mol / m^2
        q_in = np.sum(q_in_vals,axis = 0) #mol / m^2
        # mean flow rate during the dumux simulaiton
        q_out_m = q_out/ sum(BC_ddt)  #mol / m^2/s
        q_in_m = q_in/ sum(BC_ddt) #mol / m^2/s
        
        if not self.useMoles:
            raise Exception # unitConversions = 1/  10000 * 24 * 3600 * 100.  # [kg m-2 s-1] / rho = [m s-1] -> cm / day
        else:
            molarMassWat = 18. # [g/mol]
            densityWat = 1. #[g/cm3]
            # [mol/cm3] = [g/cm3] /  [g/mol] 
            molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 
            unitConversionW =- (1/1e4) * (24.* 3600.)  / molarDensityWat; # [m2/cm2]*[s/d] * [cm3/mol]
            unitConversion =- (1/1e4) * (24.* 3600.) ; # [m2/cm2]*[s/d]
            unitConversions = np.full(len(q_in_m), unitConversion)
            unitConversions[0] = unitConversionW
        # water: [mol m-2 s-1]*[m2/cm2]*[s/d] * [cm3/mol] * cm2-> [cm3/d] 
        # solute: [mol m-2 s-1]*[m2/cm2]*[s/d] * cm2 -> [mol/d] 
        Q_in_m  = q_in_m * unitConversions * (2 * np.pi * rIn  * length)
        Q_out_m = q_out_m * unitConversions * (2 * np.pi * rOut * length) 
        
        return(Q_in_m, Q_out_m)

    def getAllNeumann(self, eqIdx=0):
        """ Gathers the neuman fluxes into rank 0 as a map with global index as key [cm / day]"""
        verbose = False
        if not self.useMoles:
            assert  eqIdx == 0
            unitConversion = 1/  1000 * 24 * 3600 * 100.  # [kg m-2 s-1] / rho = [m s-1] -> cm / day
        else:
            if eqIdx == 0:
                molarMassWat = 18. # [g/mol]
                densityWat = 1. #[g/cm3]
                # [mol/cm3] = [g/cm3] /  [g/mol] 
                molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 
                unitConversion =- (1/10000) * (24.* 3600.)  / molarDensityWat; # [mol m-2 s-1]*[cm2/m2]*[s/d] * [cm3/mol]-> [cm/d] 
            else:
                unitConversion =- (1/10000) * (24.* 3600.) ; # [mol m-2 s-1] -> [mol cm-2 d-1] 
        dics = self.base.getAllNeumann(eqIdx) #[mol or kg /(m^2 * s)]
        
        if self.dimWorld != 1:#??
            flat_dic = {}
            for d in dics:
                flat_dic.update(d)
        else:
            flat_dic = dics
        try:
            for key, value in flat_dic.items():          
                posFace = 1.
                if (self.dimWorld == 1) and (key == 0):
                    posFace /= (self.getPoints()[key]/100) #axyssymmetric, was multiplied by r coord of face (in cm) in neumann
                elif (self.dimWorld == 1) and (key != 0):
                    posFace /= (self.getPoints()[-1]/100) #axyssymmetric
                flat_dic[key]  = value * unitConversion *posFace # [kg m-2 s-1] / rho /m= [s-1] -> 1 / day 
                    
        except:
            print('error unpacking dict',flat_dic,type(flat_dic))
            raise Exception
        return flat_dic

    def getNetFlux(self, eqIdx=0):
        """ Gathers the net fluxes fir each cell into rank 0 as a map with global index as key [cm3 / day]"""
        assert not self.useMoles
        self.checkGridInitialized()
        return self._map((self.base.getNetFlux(eqIdx)), 0) * 1000. *24 * 3600  # kg/s -> cm3/day

    def _map(self, x, idType, dtype=np.float64):
        """Converts rows of x to numpy array and maps it to the right indices         
        @param idType 0 dof indices, 1 point (vertex) indices, 2 cell (element) indices   
        """
        if idType == 0:  # auto (dof)
            indices = self.getDofIndices_()
        elif idType == 1:  # points
            indices = self.getPointIndices_()
        elif idType == 2:  # cells
            indices = self.getCellIndices_()
        else:
            raise Exception('PySolverBase._map: idType must be 0, 1, or 2.')
        if len(indices) >0:  # only for rank 0 not empty
            try:
                assert len(indices) == len(x), "_map: indices and values have different length"
            except:
                print(len(indices) , len(x), indices)
                raise Exception
            ndof = max(indices) + 1
            
            try:
                if isinstance(x[0], (list,type(np.array([])))) :
                    m = len(x[0])
                    p = np.zeros((ndof, m), dtype = dtype)
                    for i in range(0, len(indices)):  #
                        p[indices[i],:] = np.array(x[i], dtype = dtype)
                    if m == 1:
                        p = p.flatten()
                else:
                    p = np.zeros(ndof, dtype = dtype)
                    p[indices] = np.asarray(x, dtype = dtype)
            except:
                print('indices',indices, 'x',x)
                print('x[0]',x[0])
                raise Exception
            return p
        else:
            return np.array([])

        
        
    def distributeSources(self, source, eqIdx, numDissolvedSoluteComp: int, selectCell = None):
        """ when we have a source in for a 1d model (exchange with bulk soil),
            divid the source between the cell of the 1d model according to the water
            or solute content
            @param sources: sources to divide between the cells for water [list, cm3/day] and/or solutes [list, mol/day]
            @param eqIdx: index of the components [list of int]
            @param numDissolvedSoluteComp: number of solutes which are disolved
            @param selectCell: optional, manually set the index of the cell the source will be applied to [int]
        """
        splitVals = list()
        # iterate through each component
        for i, src in enumerate(source):# [cm3/day] or [mol/day]
            splitVals.append(self.distributeSource( src, eqIdx[i], numDissolvedSoluteComp, selectCell))
        return np.array(splitVals, dtype = object)
            
    def distributeSource(self, source: float, eqIdx: int, numDissolvedSoluteComp: int, selectCell = None):
        """ when we have a source in for a 1d model (exchange with bulk soil),
            divid the source between the cell of the 1d model according to the water
            or solute content
            @param sources: sources to divide between the cells for water [cm3/day] or solute [mol/day]
            @param eqIdx: index of the components [int]
            @param numDissolvedSoluteComp: number of solutes which are disolved
            @param selectCell: optional, manually set the index of the cell the source will be applied to [int]
        """
        assert self.dimWorld != 3
        
        # distribute the value between the cells
        splitVals = self.distributeVals(source, eqIdx, numDissolvedSoluteComp, selectCell)       
        # send the source data to dumux
        if source != 0.:# [cm3/day] or [mol/day]
            test_values = list(splitVals.copy())
            test_keys = np.array([i for i in range(len(test_values))])
            res = {}
            for key in test_keys:
                for value in test_values:
                    res[key] = value
                    test_values.remove(value)
                    break                        
            self.setSource(res.copy(), eq_idx = eqIdx)  # [mol/day], in modules/richards.py
        else:
            res = dict()
            res[0] = 0.
            self.setSource(res.copy(), eq_idx = eqIdx)  # [mol/day], in modules/richards.py
        return splitVals
    
    def distributeVals(self, source: float, eqIdx: int, numDissolvedSoluteComp: int, selectCell = None):
        """ when we have a source in for a 1d model (exchange with bulk soil),
            divid the source between the cell of the 1d model according to the water
            or solute content
            @param sources: sources to divide between the cells for water [cm3/day] or solute [mol/day]
            @param eqIdx: index of the components [int]
            @param numDissolvedSoluteComp: number of solutes which are disolved
            @param selectCell: optional, manually set the index of the cell the source will be applied to [int]
        """
        splitVals = np.array([0.])
        verbose = False
        if source != 0.:# [cm3/day] or [mol/day]
            if selectCell == None:
                if eqIdx == 0:
                    seg_values_ = self.getSolutionHead()#self.getWaterVolumes()#getWaterVolumesCyl(length)
                    seg_values = seg_values_ - min(seg_values_) +1e-14
                else:
                    isDissolved = (eqIdx <= numDissolvedSoluteComp)
                    seg_values = self.getContentCyl(eqIdx, isDissolved)
                # during the solve() loop, we might still get seg_values <0 
                # assert min(seg_values) >= 0.
                seg_values = np.maximum(seg_values,0.)

                if (sum(seg_values) == 0.):# should normally only happen with source >= 0
                    weightVals =np.full(len(seg_values), 1 /len(seg_values))
                elif source < 0:# goes away from the 1d models
                    weightVals = seg_values /sum(seg_values)
                    if verbose:
                        print("sum(seg_values[segIds])", seg_values, weightVals)
                else:# goes toward  the 1d models
                    if min(abs(seg_values)) == 0.:# at least 1 seg_values = 0 but not all
                        seg_values = np.maximum(seg_values,1.e-14)
                        assert min(abs(seg_values)) != 0.
                    weightVals = (1 / seg_values) / sum(1/seg_values)
            else:
                assert isinstance(selectCell, int) and (selectCell >= 0) and (selectCell < self.numberOfCellsTot)
                weightVals = np.zeros(self.numberOfCellsTot)
                weightVals[selectCell] = 1.
            splitVals = weightVals * source
            try:
                assert (sum(weightVals) - 1.) < 1e-13
                assert len(splitVals) == self.numberOfCellsTot
            except:
                print('(sum(weightVals) - 1.) < 1e-13',rank,weightVals, sum(weightVals),(sum(weightVals) - 1.) ,(sum(weightVals) - 1.) < 1e-13)
                print('splitVals',splitVals, self.numberOfCellsTot, len(splitVals) == self.numberOfCellsTot)
                raise Exception
            if verbose:
                print(rank,'distributeVals',eqIdx,'splitVals',splitVals)   
        return splitVals
