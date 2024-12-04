import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import griddata
from mpi4py import MPI; comm = MPI.COMM_WORLD; size = comm.Get_size(); rank = comm.Get_rank()

from helpfull import StdoutRedirector

class SolverWrapper():
    """ Additional functionality to the C++ Python binding. 
    
        The C++ class derived from SolverBase is wrapped, 
        SolverBase functions are passed to C++ 
        
        Additionally, contains mainly methods that are easier to write in Python, 
        e.g. MPI communication, writeVTK, interpolate        
    """

    def __init__(self, base, usemoles):
        """ 
            @param base is the C++ base class that is wrapped. 
            @param usemoles: do we use moles (mol) or mass (g) in dumux
            other parameters are created and will be set in
            @see scenario_setup::setSoilParam()
        """
        self.base = base
        self.useMoles = usemoles
        self.solidDensity = 0
        self.solidMolarMass=0
        self.solidMolDensity=0
        self.bulkDensity_m3 =0 # mol soil minerals / m3 bulk soil density
        self.molarMassWat = 18. # [g/mol]
        self.densityWat_m3 = 1e6 #[g/m3]
        self.m3_per_cm3 = 1e-6 #m3/cm3
        self.cm3_per_m3 = 1e6 #cm3/m3
        # [mol/m3] = [g/m3] /  [g/mol] 
        self.molarDensityWat_m3 =  self.densityWat_m3 / self.molarMassWat # [mol wat/m3 wat] 
        self.molarDensityWat_cm3 =  self.molarDensityWat_m3 /1e6 # [mol wat/cm3 wat] 
        self.results_dir = "./results/"
        self.pindx = 0
        self.periodic = False
        

    def initialize(self, args_ = [""], verbose = False,doMPI_=True):
        """ Writes the Dumux welcome message, and creates the global Dumux parameter tree 
        
            @params verbose: level of verbosity of the dumux object [bool]
            @params doMPI_: do we use multi processes in the dumux object [bool]
            If the python code is run in parallel, we can use doMPI_=False for the
            1d models only (FoamGrid). Does not work for the other (3d) grids.
        """
        try:
            with StdoutRedirector() as redirector:
                self.base.initialize(args_, verbose, doMPI=doMPI_)
        except Exception as e:
            target_filepath = self.results_dir + 'stdcout_cpp'+str(self.pindx)+"_"+str(rank)+'.txt'
            with open(target_filepath, 'w') as f:
                f.write(redirector.buffer)
            raise Exception
            
    @property
    def dimWorld(self):
        """Get the dimention of the dumux domain (1d, 2d, 3d) """
        return self.base.dimWorld
        
    @property
    def numSoluteComp(self):
        """Get the number of components evaluated (not counting water)."""
        return self.base.numComp() - 1 
    
    @property
    def numComp(self):
        """Get the number of components evaluated (including water) ."""
        return self.base.numComp()  
        
    @property
    def numFluidComp(self):
        """Get the number of components in the fluid phase (including water)."""
        return self.base.numFluidComp()
    
    @property
    def numDissolvedSoluteComp(self):
        """Get the number of disolved solutes."""
        return self.base.numFluidComp() - 1
        
    def setMaxTimeStepSize(self, maxDt):
        """
            change the maximum inner time step which can be tested by dumux
        """
        self.base.setMaxTimeStepSize(maxDt)
        
    def createLinearSolver(self):
        """
            manually (re)create nonlinear solver. 
            useful to implement new solver parameters
        """
        self.base.createLinearSolver()
    
    def createNewtonSolver(self):
        """
            manually (re)create nonlinear solver. 
            useful to implement new solver parameters
        """
        try:
            with StdoutRedirector() as redirector:
                self.base.createNewtonSolver()
        except Exception as e:
            target_filepath = self.results_dir + 'stdcout_cpp'+str(self.pindx)+"_"+str(rank)+'.txt'
            with open(target_filepath, 'w') as f:
                f.write(redirector.buffer)
            raise Exception
            
            
    def save(self):
        """ 
            saves the current value of the solution vector.
            that saves will be over-written when calling
            @see solve() by the solution vector before the soling.
            The solution vector is reset to that value when calling @see reset()
            usefull for the fixed point iteraiton
        """
        self.base.save()
        
        
    def saveManual(self):
        """ saves the current value of the solution vector. 
            Allows to save a different solution vector from @see save()
            The solution vector is reset to that value when calling @see resetManual()
            usefull for the fixed point iteraiton
        """
        self.base.saveManual()
        
    def reset(self):
        """ 
            reset solution vector to its value before solve function or to
            the value when calling @see save (whichever came last)
        """
        self.base.reset()
        
    def resetManual(self):
        """ 
            reset solution vector to its value when calling @see saveManual()
        """
        self.base.resetManual()
        
        
    def setVerbose(self, verbose:int):
        """ set verbose level """
        self.base.setVerbose(verbose)

    def gather(self, data2gather):
        return comm.gather(data2gather, root=0) 
    
    def allgatherv(self,data2share, keepShape = False, data2share_type_default = float): 
        """ own allgather() function for vectors
            use it if the vectors exchanged between the threads might be of
            different size
            @param data2share: (local) vector of the thread to share
            @param keepShape: the output vector has to be of the same shape as data2share
            @data2share_type_default: type of the data in data2share
        """
        try:
            assert isinstance(data2share, (list, type(np.array([]))))
        except:
            error = ('on rank {}: '.format(rank)+
                'solverbase::allgatherv() the data to share should be of type list or np.array.'+
                r'not {}.'.format(type(data2share))     
            ) 
            print(error)
            raise Exception
            
        data2share = np.array(data2share)
            
        # get shape and size of local array
        if len((data2share).shape) == 2:
            local_size = (data2share).shape[0] * (data2share).shape[1]
            shape0 = (data2share).shape[0]
            shape1 = (data2share).shape[1]
            data2share_type = type(data2share[0][0])
        elif len((data2share).shape) == 1:
            local_size = (data2share).shape[0] 
            shape0 = (data2share).shape[0]
            shape1 = 0
            if len(data2share) > 0:
                data2share_type = type(data2share[0])
            else:
                data2share_type = data2share_type_default
        else:
            raise Exception  
        
        # data2share needs to use floats for Allgatherv with MPI.DOUBLE to work.
        data2share = np.array(data2share, dtype = np.float64)
        
        # other data needed by comm.Allgatherv
        all_sizes = np.array(comm.allgather(local_size))
        work_size = sum(all_sizes)
        all_data2share = np.zeros(work_size)

        offsets = np.zeros(len(all_sizes), dtype=np.int64)
        offsets[1:]=np.cumsum(all_sizes)[:-1]
        all_sizes =tuple(all_sizes)
        offsets =tuple( offsets)      
        
        # share the vectors
        comm.Allgatherv( [data2share.reshape(-1), MPI.DOUBLE],[all_data2share,all_sizes,offsets,MPI.DOUBLE])
        
        all_data2share = np.array(all_data2share, dtype =data2share_type)
        
        if keepShape:
            if shape1 > 0:
                all_data2share = all_data2share.reshape( size,shape0,shape1)
            else:
                all_data2share = all_data2share.reshape(size, shape0)
        else:
            if shape1 > 0:
                all_data2share = all_data2share.reshape(-1,shape1)
            else:
                all_data2share = all_data2share.reshape(-1)
        
        return all_data2share

    def createGridFromInput(self, modelParamGroup = ""):
        """ Creates the Grid and gridGeometry from the global DuMux parameter tree """
        self.base.createGrid(modelParamGroup)

    def createGrid(self, boundsMin, boundsMax, numberOfCells, periodic = False):
        """ Creates a rectangular grid with a given resolution             
            @param boundsMin        domain corner [cm]
            @param boundsMax        domain corner [cm]        
            @param numberOfCells    resoultion [1]
            @param periodic         If true, the domain is periodic in x and y, not in z 
        """
        
        self.periodic = periodic
        try:
            with StdoutRedirector() as redirector:
                self.base.createGrid(np.array(boundsMin) / 100., np.array(boundsMax) / 100., 
                                     np.array(numberOfCells), periodic
                                    )  # cm -> m
        except Exception as e:
            target_filepath = self.results_dir + 'stdcout_cpp'+str(self.pindx)+"_"+str(rank)+'.txt'
            with open(target_filepath, 'w') as f:
                f.write(redirector.buffer)
            raise Exception
            
        self.numberOfCellsTot = np.prod(self.numberOfCells)
        if self.dimWorld == 3:
            self.numberOfFacesTot = self.numberOfCellsTot * 6
        elif self.dimWorld == 1:
            self.numberOfFacesTot = self.numberOfCellsTot * 2
        else:
            raise Exception

    def createGrid1d(self, points):
        """ todo
        """
        p = []
        for v in points:
            p.append([v / 100.])  # cm -> m
        try:
            with StdoutRedirector() as redirector:
                self.base.createGrid1d(p)
        except Exception as e:
            target_filepath = self.results_dir + 'stdcout_cpp'+str(self.pindx)+"_"+str(rank)+'.txt'
            with open(target_filepath, 'w') as f:
                f.write(redirector.buffer)
            raise Exception
            
        self.numberOfCellsTot = len(points) -1
        self.numberOfFacesTot = self.numberOfCellsTot * 2


    def readGrid(self, file:str):
        """ Creates a grid from a file (e.g. dgf or msh)"""
        self.base.readGrid(file)

    def getGridBounds(self):
        """  Returns a rectangular bounding box around the grid geometry [cm] """
        return np.array(self.base.getGridBounds()) * 100.  # m -> cm

        
    def setParameter(self, key:str, value:str):
        """ Writes a parameter into the global Dumux parameter map """
        self.base.setParameter(key, value)

    def getParameter(self, key:str):
        """ Reads a parameter from the global Dumux parameter map, returns an empty string if value is not set """
        return self.base.getParameter(key)

    def initializeProblem(self, maxDt = -1.):
        """ After the grid is created, the problem can be initialized 
            @param: maxDt : maximum inner-time step that dumux can use [d]            
        """
        try:
            with StdoutRedirector() as redirector:
                self.base.initializeProblem(maxDt * 24.*3600.)
        except Exception as e:
            target_filepath = self.results_dir + 'stdcout_cpp'+str(self.pindx)+"_"+str(rank)+'.txt'
            with open(target_filepath, 'w') as f:
                f.write(redirector.buffer)
            raise Exception
            
        ##       saves shape data of the grid  
        # local indices
        self.dofIndices_   = self.base.getDofIndices()
        self.pointIndices_ = self.base.getPointIndices()
        self.cellIndices_  = self.base.getCellIndices()
        # all indices
        self.dofIndices   = np.asarray(self._flat0(self.gather(self.dofIndices_)), np.int64)
        self.pointIndices = np.asarray(self._flat0(self.gather(self.pointIndices_)), np.int64) 
        self.cellIndices  = np.asarray(self._flat0(self.gather(self.cellIndices_)), np.int64) 
        self.cellsVertex = self._map(self._flat0(self.gather(self.base.getCells())), 2, np.int64)
            
        # volumes, surface
        self.CellVolumes_ =np.array( self.base.getCellVolumes()) * 1.e6  # m3 -> cm3
        self.CellVolumes = self._map(self._flat0(self.gather(self.CellVolumes_)), 2)   # m2 -> cm2
        # coordinates
        self.pointCoords = self._map(self._flat0(self.gather(self.base.getPoints())), 1) * 100.  # m -> cm
        self.cellCenters = self._map(self._flat0(self.gather(self.base.getCellCenters())), 2) * 100.  # m -> cm
        self.dofCoordinates = self._map(self._flat0(self.gather(self.base.getDofCoordinates())), 0) * 100.  # m -> cm

        #self.cellSurfacesCyl = self._map(self._flat0(self.gather(self.base.getCellSurfacesCyl())), 2) * 1.e4  # m2 -> cm2

    def setInitialCondition(self, ic, eqIdx = 0):
        """ Sets the initial conditions for all global elements, processes take from the shared @param ic """
        self.base.setInitialCondition(ic, eqIdx)

    def setInitialConditionHead(self, ic):
        """ Sets the initial conditions for all global elements, processes take from the shared @param ic """
        self.base.setInitialConditionHead(ic)

    def solve(self, dt:float, doMPIsolve_=True, saveInnerFluxes_ = True):
        """ Simulates the problem, the internal Dumux time step ddt is taken from the last time step 
        @param dt      time span [days] 
        @param doMPIsolve_    use MPI when solving? [bool] only usefull when using 1d models ( doMPIsolve_= False)
        with python code in parallel
        @param saveInnerFluxes_ save the inter-cell and boundary flows computed by dumux [bool]
        usefull to compute the error rate and for 1d-3d coupling
        """
        self.base.solve(dt * 24.*3600.,doMPIsolve=doMPIsolve_, saveInnerDumuxValues = saveInnerFluxes_)  # days -> s

    def solveSteadyState(self):
        """ Finds the steady state of the problem """
        self.base.solveSteadyState()

    def getPoints(self):
        """Gathers vertices into rank 0, and converts it into numpy array (Np, 3) [cm]"""
        self.checkGridInitialized()
        return self.pointCoords

    def getPoints_(self):
        """nompi version of """
        self.checkGridInitialized()
        return np.array(self.base.getPoints()) * 100.  # m -> cm

    def getCellCenters(self):
        """Gathers cell centers into rank 0, and converts it into numpy array (Nc, 3) [cm]"""
        self.checkGridInitialized()
        return self.cellCenters

    def getCellCenters_(self):
        """nompi version of getCellCenters"""
        self.checkGridInitialized()
        return np.array(self.base.getCellCenters()) * 100.  # m -> cm

    def getDofCoordinates(self):
        """Gathers dof coorinates into rank 0, and converts it into numpy array (Ndof, 3) [cm]"""
        return self.dofCoordinates

    def getDofCoordinates_(self):
        """nompi version of getDofCoordinates"""
        self.checkGridInitialized()
        return np.array(self.base.getDofCoordinates()) * 100.  # m -> cm

    def getCells(self):
        """ Gathers dune elements (vtk cells) as list of list of vertex indices (vtk points) (Nc, Number of corners per cell) [1]"""
        return self.cellsVertex

    def getCells_(self):
        """nompi version of getCells"""
        return np.array(self.base.getCells(), dtype = np.int64)

    #def getCellSurfacesCyl(self):
    #    """ Gathers element volumes (Nc, 1) [cm3] """
    #    raise Exception
    #    #return self.cellSurfacesCyl 

    #def getCellSurfacesCyl_(self):
    #    """nompi version of getCellSurfacesCyl"""
    #    raise Exception
    #    #return np.array(self.base.getCellSurfacesCyl()) * 1.e4  # m2 -> cm2
        
    def getCellVolumes(self):
        """ Gathers element volumes (Nc, 1) [cm3] """
        return self.CellVolumes # m2 -> cm2

    def getCellVolumes_(self):
        """nompi version of  """
        return self.CellVolumes_ # m3 -> cm3

    def getCellVolumesCyl(self):
        """ Gathers element volumes (Nc, 1) [cm3] """
        return self.getCellVolumes()

    def getCellVolumesCyl_(self):
        """nompi version of  getCellVolumesCyl"""
        return self.getCellVolumes_()#np.array(self.base.getCellVolumesCyl()) * 1.e6  # m3 -> cm3

    def getCellIndices(self):
        """Gathers dof indicds into rank 0, and converts it into numpy array (dof, 1)"""
        if rank > 0:
            return []
        else:
            return self.cellIndices
        
    def getPointIndices(self):
        """Gathers dof indicds into rank 0, and converts it into numpy array (dof, 1)"""
        if rank > 0:
            return []
        else:
            return self.pointIndices
        
    def getDofIndices(self):
        """Gathers dof indicds into rank 0, and converts it into numpy array (dof, 1)"""
        if rank > 0:
            return []
        else:
            return self.dofIndices

    def getCellIndices_(self):
        """nompi version of  """
        self.checkGridInitialized()
        return self.cellIndices_
    
    def getPointIndices_(self):
        """nompi version of  """
        self.checkGridInitialized()
        return self.pointIndices_
    
    def getDofIndices_(self):
        """nompi version of  """
        self.checkGridInitialized()
        return self.dofIndices_

    def getSolution(self, eqIdx = 0):
        """Gathers the current solution into rank 0, and converts it into a numpy array (dof, neq), 
        model dependent units [Pa, ...]"""
        self.checkGridInitialized()
        return self._map(self._flat0(self.gather(self.base.getSolution(eqIdx))),0)

    def getAvgDensity_(self):
        """nompi version of  """
        self.checkGridInitialized()
        return np.array(self.base.getAvgDensity())
        
    def getSolution_(self, eqIdx = 0):
        """nompi version of  """
        self.checkGridInitialized()
        return np.array(self.base.getSolution(eqIdx))

    def getSolutionAt(self, gIdx, eqIdx = 0):
        """Returns the current solution at a cell index, model dependent units [Pa, ...]"""
        return self.base.getSolutionAt(gIdx, eqIdx)

    def getNeumann(self, gIdx, eqIdx = 0):
        """ Gathers the neuman fluxes into rank 0 as a map with global index as key [cm / day]
            ATT: only gives the current flux, NOT the mean flux computed by dumux during the last solve() call
        """
        assert not self.useMoles
        return self.base.getNeumann(gIdx, eqIdx) / 1000 * 24 * 3600 * 100.  # [kg m-2 s-1] / rho = [m s-1] -> cm / day

    def getAllNeumann(self, eqIdx = 0):
        """ Gathers the neuman fluxes into rank 0 as a map with global index as key [cm / day]
            ATT: only gives the current flux, NOT the mean flux computed by dumux during the last solve() call
        """
        assert not self.useMoles # need the update the unit change before using when useMoles==True
        dics = self.gather(self.base.getAllNeumann(eqIdx))
        flat_dic = {}
        for d in dics:
            flat_dic.update(d)
        for key, value in flat_dic:
            flat_dic[key] = value / 1000 * 24 * 3600 * 100.  # [kg m-2 s-1] / rho = [m s-1] -> cm / day
        return flat_dic

    def getAllNeumann_(self, eqIdx = 0):
        """nompi version of (TODO is that working?)"""
        assert not self.useMoles
        dics = self.base.getAllNeumann(eqIdx)
        flat_dic = {}
        for d in dics:
            flat_dic.update(d)
        for key, value in flat_dic:
            flat_dic[key] = value / 1000 * 24 * 3600 * 100.  # [kg m-2 s-1] / rho = [m s-1] -> cm / day
        return flat_dic

    def getNetFlux(self, eqIdx = 0):
        """ Gathers the net fluxes fir each cell into rank 0 as a map with global index as key [cm3 / day]"""
        assert not self.useMoles
        self.checkGridInitialized()
        return self._map(self._flat0(self.gather(self.base.getNetFlux(eqIdx))), 0) * 1000. *24 * 3600  # kg/s -> cm3/day

    def getNetFlux_(self, eqIdx = 0):
        """nompi version of """
        assert not self.useMoles
        self.checkGridInitialized()
        return np.array(self.base.getNetFlux(eqIdx)) * 1000. *24 * 3600  # kg/s -> cm3/day

    def pickCell(self, pos):
        """ Picks a cell and returns its global element cell index """
        return self.base.pickCell(np.array(pos) / 100.)  # cm -> m

    def pick(self, x):
        """ Picks a cell and returns its global element cell index """
        return self.base.pick(np.array(x) / 100.)  # cm -> m
        
        
    def pick_(self,coordCell):
        """ non MPI version of pick, to use when the plant object is not computed on all threads
            coordCell: cell coordinates in cm
        """
        bounds = self.getGridBounds(); # cm
        min_b = bounds[:3]; max_b = bounds[3:]
        cell_number_ = np.array(self.numberOfCells)
        
        if (self.isPeriodic) :
            for i in [0,1]: # for x and y, not z
                minx = min_b[i];
                xx = max_b[i]-minx; # unit of demain distance
                if (not np.isinf(xx)) :# periodic in x
                    coordCell[i] -= minx # start at 0 (relative coordinates)
                    if (coordCell[i]>=0) :
                        coordCell[i] = coordCell[i] - int(coordCell[i]/xx)*xx;
                    else :
                        coordCell[i] = coordCell[i] + int((xx-coordCell[i])/xx)*xx;
                    
                    coordCell[i] += minx# back to absolute coordinates             
            
        # is in domain according to x,y,z?
        ratioDist = (coordCell - min_b)/(max_b - min_b)
        if ((ratioDist > 1.) |(ratioDist < 0.) ).any():#not in the domain
            return -1
                 
        id_rows = ratioDist*cell_number_
        onCellOuterFace = np.where((id_rows == np.round(id_rows)) & (id_rows > 0) )[0]
        id_rows[onCellOuterFace] -= 1
        id_rows= np.floor(id_rows)
        id_cell = id_rows[2]*(cell_number_[0]*cell_number_[1])+id_rows[1]*cell_number_[0]+id_rows[0]
        return int(id_cell) 

    def __str__(self):
        """ Solver representation as string """
        return str(self.base)

    def checkGridInitialized(self):
        """ Checks if the problem was initialized, i.e. initializeProblem() was called """
        return self.base.checkGridInitialized()

    def checkProblemInitialized(self, throwError = True):
        """ Checks if the problem was created, i.e. initializeProblem() was called """
        if self.base.checkProblemInitialized():
            return True
        elif throwError:
            raise Exception('problem was not created, call initializeProblem()')
        else:
            return False
    
    @property
    def simTime(self):
        """ Current simulation time (read only) [days]"""
        return self.base.simTime / 24. / 3600.  # s->days

    @property
    def rank(self):
        """ MPI rank (read only) """
        return self.base.rank

    @property
    def maxRank(self):
        """ number of MPI processes (read only) """
        return self.base.maxRank

    @property
    def ddt(self):
        """ last internal time step, i.e. Dumux time step [days]"""
        return self.base.ddt / 24. / 3600.  # s-> days

    @ddt.setter
    def ddt(self, value):
        """ sets internal time step, i.e. Dumux time step [days]"""
        self.base.ddt = value * 24.*3600.  # days -> s

    def interpolate(self, xi, eq = 0):
        """ interpolates the solution at position ix [cm],
        model dependent units
        todo: test"""
        self.checkGridInitialized()
        points = self.getDofCoordinates()
        solution = self.getSolution()
        if rank == 0:
            yi = np.zeros((xi.shape[0]))
            for i in range(0, xi.shape[0]):
                yi[i] = griddata(points, solution[:, eq], xi / 100., method = 'linear', rescale = True)  # cm -> m
            return yi
        else:
            return []

    def interpolate_(self, xi, points, solution, eq = 0):
        """ interpolates the solution at position ix [cm] (same es interpolate_ but passes points and solution),
        model dependent units """
        self.checkGridInitialized()
        if rank == 0:
            yi = np.zeros((xi.shape[0]))
            for i in range(0, xi.shape[0]):
                yi[i] = griddata(points, solution[:, eq], xi / 100., method = 'linear')  # cm -> m
            return yi
        else:
            return []

    def interpolateNN(self, xi, eq = 0):
        """ solution at the points xi (todo currently works only for CCTpfa)"""
        self.checkGridInitialized()
        solution = self.getSolution()
        if rank == 0:
            y = np.zeros((xi.shape[0]))
        else:
            y = []
        for i in range(0, xi.shape[0]):
            idx = self.pickCell(xi[i,:])  # cm -> m
            if rank == 0:
                y[i] = solution[idx, eq]
        return y

    def writeVTK(self, file:str, small:bool = False):
        """writes a vtk file (todo additional fields) 
        @param file 
        @param small Determines if data are compressed and stroed binary  TODO
        """
        points = self.getDofCoordinates()
        if self.rank == 0:
            pd = vtk.vtkPolyData()
            pd.SetPoints(_vtkPoints(points))
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName(file)
            writer.SetInputData(pd)
            writer.SetDataModeToBinary()
            writer.SetCompressorTypeToZLib()
            writer.Write()

    def _map(self, x, type_, dtype = np.float64):
        """Converts rows of x to numpy array and maps it to the right indices         
        @param type_ 0 dof indices, 1 point (vertex) indices, 2 cell (element) indices   
        """
        verbose_ = False
        if type_ == 0:  # auto (dof)
            indices = self.getDofIndices() 
        elif type_ == 1:  # points
            indices = self.getPointIndices() 
        elif type_ == 2:  # cells
            indices =  self.getCellIndices() 
        else:
            raise Exception('PySolverBase._map: type_ must be 0, 1, or 2.')
            
        if rank == 0:
            try:
                assert len(indices) == len(x), "_map: indices and values have different length"
            except:
                print('assert len(indices) == len(x)',
                      len(indices) , len(x), indices)
                raise Exception
            ndof = max(indices) + 1
            if isinstance(x[0], (list,type(np.array([])))) :
                m = len(x[0])
                p = np.zeros((ndof, m), dtype = dtype)
                p[indices,:] = np.asarray(x, dtype = dtype)
            else:
                p = np.zeros(ndof, dtype = dtype)
                p[indices] = np.asarray(x, dtype = dtype)
            return p
        
        else:
            #return array instead of float to be able to have same object type no matter what
            return np.array([])

    def _flat0(self, xx, dtype_=None):
        """flattens the gathered list in rank 0, empty list for other ranks """
        if rank == 0:
            return np.array([item for sublist in xx for item in sublist],dtype = dtype_)
        else:
            return np.array([])

    @staticmethod
    def _vtkPoints(p):
        """ Creates vtkPoints from an numpy array"""
        assert(p.shape(1) == 3)
        da = vtk.vtkDataArray.CreateDataArray(vtk.VTK_DOUBLE)
        da.SetNumberOfComponents(3)  # vtk point dimension is always 3
        da.SetNumberOfTuples(p.shape[0])
        for i in range(0, p.shape[0]):
            da.InsertTuple3(i, p[i, 0], p[i, 1], p[i, 2])
        points = vtk.vtkPoints()
        points.SetData(da)
        return points

    @staticmethod
    def _vtkCells(t):
        """ Creates vtkCells from an numpy array"""
        cellArray = vtk.vtkCellArray()
        if t.shape[1] == 2:
            Simplex = vtk.vtkLine
        elif t.shape[1] == 4:
            Simplex = vtk.vtkTetra
        elif t.shape[1] == 8:
            Simplex = vtk.Hexaedron
        else:
            raise Exception('PySolverBase._vtkCells: do not know what to do with {} vertices'.format(t.shape[1]))
        for vert in t:
            tetra = Simplex()
            for i, v in enumerate(vert):
                tetra.GetPointIds().SetId(i, int(v))
            cellArray.InsertNextCell(tetra)
        return cellArray

    @staticmethod
    def _vtkData(data):
        """ todo """
        pass

    @staticmethod
    def to_pa(ph):
        """ Converts cm pressure head to Pascal [kg/ (m s^2)] """
        return 1.e5 + ph / 100 * 1000. * 9.81;

    @staticmethod
    def to_head(p):
        """ Converts Pascal [kg/ (m s^2)] to cm pressure head """
        return (p - 1.e5) * 100. / 1000. / 9.81;

    @property
    def segLength(self):
        """ 
            return the legnth of the root segment (when using the perirhizal zone)
            as saved by dumux [cm] 
            returns None I think for dimWorld != 1
        """
        self.checkGridInitialized() # is this the correct check? we need to have called initializeProblem before
        return self.base.segLength()*100 # m to cm