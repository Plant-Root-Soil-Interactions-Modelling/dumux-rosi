from solverbase import SolverWrapper

import numpy as np
from mpi4py import MPI; 
comm = MPI.COMM_WORLD; size = comm.Get_size(); rank = comm.Get_rank()
max_rank = comm.Get_size()

class RichardsWrapper(SolverWrapper):
    """ 
    Adds functionality specifically for Richards equation             
    Wraps passing parameters to Dumux for VanGenuchtenParameter, IC, BC, 
    """

    def __init__(self, base, usemoles):
        """ 
            @param base is the C++ base class that is wrapped. 
            @param usemoles: do we use moles (mol) or mass (g) in dumux
        """
        super().__init__(base, usemoles)
        self.soils = []
        self.param_group = "Soil."
        self.useMoles = usemoles
        # self.molarMassC = 12.011 # g/mol

    def setParameterGroup(self, group:str):
        """ sets the DuMux paramter group, must end with a dot, e.g. 'Soil.' """
        self.param_group = group
        

    @staticmethod
    def dumux_str(l):
        """ to pass lists to dumux parameter tree as string """
        if type(l).__name__ == "float" or type(l).__name__ == "int":
            return str(l)
        elif type(l).__name__ == 'float64':
            return str(l)
        else:
            s = ""
            for l_ in l:
                s += " " + str(l_)
            return s[1:]

    def setVGParameters(self, soils:list):
        """ Sets a list of vg parameter lists, 5 values per soil 
            set before creating the problem with SolverBase.initializeProblem 
        
            # todo doc params
        """
        if self.checkProblemInitialized(throwError = False):
            raise Exception('use setVGParameters before calling initialize()')
            
        qr, qs, alpha, n, ks = [], [], [], [], []

        for soil in soils:
            assert len(soil) == 5, "setVGParameters, soil need to be a list of 5 parameters: qr, qs, alpha, n, ks "
            qr.append(soil[0])
            qs.append(soil[1])
            alpha.append(soil[2])
            n.append(soil[3])
            ks.append(soil[4])

        self.setParameter(self.param_group + "VanGenuchten.Qr", self.dumux_str(qr))
        self.setParameter(self.param_group + "VanGenuchten.Qs", self.dumux_str(qs))
        self.setParameter(self.param_group + "VanGenuchten.Alpha", self.dumux_str(alpha))
        self.setParameter(self.param_group + "VanGenuchten.N", self.dumux_str(n))
        self.setParameter(self.param_group + "VanGenuchten.Ks", self.dumux_str(ks))

        self.soils = soils.copy()

    def setLayersZ(self, number:list, z:list = []):
        """ sets depth dependent layers 
        
        @param number     list of layer numbers at the z-positions (if given), or per soil layer, [cm] pressure head.      
        @param z          list of z-positions [cm].  Between the sampling points linear interpolation is applied.                              
        """
        self.setParameter(self.param_group + "Layer.Number", self.dumux_str(number))
        if z:
            assert len(number) == len(z), "setLayersZ: sample point values and z coordinates have unequal length"
            self.setParameter(self.param_group + "Layer.Z", self.dumux_str(np.array(z) / 100.))  # cm -> m

    def setICZ(self, p:list, z:list = []):
        """ sets depth dependent initial condtions 
        
        @param p     list of pressures at the z-positions (if given), or per soil layer, [cm] pressure head.      
        @param z     list of z-positions [cm].  Between the sampling points linear interpolation is applied.                              
        """
        self.setParameter(self.param_group + "IC.P", self.dumux_str(p))
        if z:
            assert(len(p) == len(z))  # sample points
            self.setParameter(self.param_group + "IC.Z", self.dumux_str(np.array(z) / 100.))

    def setICZ_solute(self, c:list, z:list = []):
        """ sets depth dependent initial condtions for solutes 
        
        @param p     list of concentrations at the z-positions (if given), or per soil layer, [g/cm3].      
        @param z     list of z-positions [cm].  Between the sampling points linear interpolation is applied.                              
        """
        if isinstance(c, float):
            c = [float(c)]
        self.setParameter(self.param_group + "IC.C", self.dumux_str(c))
        if z:
            assert(len(c) == len(z))  # sample points
            self.setParameter(self.param_group + "IC.CZ", self.dumux_str(np.array(z) / 100.))

    def setHomogeneousIC(self, p:float, equilibrium = False):
        """ sets homogeneous initial condions 
        
        @param p              mean matric potential [cm] pressure head
        @param equilibrium    in hydrostatic equilibrium (True) or with a static matric potential (False)
                              for hydrostatic equilibrium the grid must be created before  
        """
        if equilibrium:
            bounds = self.getGridBounds()
            z = [bounds[2], bounds[5]]
            m = (z[1] - z[0]) / 2.
            p = [p + m, p - m]
            self.setICZ(p, z)
        else:
            self.setICZ([p])

    def setLinearIC(self, top, bot):
        """ sets linear initial conditions from @param top matric potential to @param bot matric potential """
        bounds = self.getGridBounds()
        z = [bounds[2], bounds[5]]  # min, max
        p = [bot, top]
        self.setICZ(p, z)

    def setTopBC(self, type_top:str, value_top:float = 0., climate:list = []):
        """ Top boundary conditions are set before creating the problem with SolverBase.initializeProblem 
        
        @param type_top:
        type_top ==  "constantPressure", value_top is a constant pressure [cm] pressure head
        type_top ==  "constantFlux", value_top is the constant flux [cm/day]
        type_top ==  "constantFluxCyl", value_top is the constant flux [cm/day] for cylindrical coordinates                          
        type_top ==  "atmospheric", value_top is given by climatic data describing evapotranspiration [cm/day], 
                     Data are given in @param climate, the value of value_top is ignored.  
                     Minus denotes evaporation, plus transpiraton.                                            
                     Evaporation stops at a critical pressure of -10000 cm, infiltration is with run off.     

        @param climate:  Two lists are expected, first a list of times [day], second a list of evapotranspirations [cm/day], 
                         between the values linear interpolation is applied.                                                                          
        """
        if type_top == "constantPressure" or type_top == "pressure":
            t = 1
        elif type_top == "constantFlux" or type_top == "flux":
            t = 2
        elif type_top == "constantFluxCyl" or type_top == "fluxCyl":
            t = 3
        elif type_top == "atmospheric":
            t = 4
        elif type_top == "noflux" or type_top == "noFlux" or type_top == "no-flux":
            t = 2
            assert value_top == 0., "setTopBC: value_top must be zero in case of no flux"
        else:
            raise Exception('richards.setTopBC(): Top type should be "constantPressure", "constantFlux", constantCyl", "noFlux", or "atmospheric" unknown top type {}'.format(type_top))

        self.setParameter(self.param_group + "BC.Top.Type", str(t))
        self.setParameter(self.param_group + "BC.Top.Value", str(value_top))

        if t == 4:  # atmospheric
            if climate:
                assert(len(climate[0]) == len(climate[1]))  # sample points
                self.setParameter("Climate.Time", self.dumux_str(climate[0]))
                self.setParameter("Climate.Precipitation", self.dumux_str(climate[1]))  # TODO confusing name (should be Evapotranspiration)

            else:
                raise Exception('richards.setTopBC(): Atmospheric boundary conditions where set, but no climatic data were given')

    def setTopBC_solute(self, type_top:str, value_top:float = 0., managed:list = []):
        if type_top == "constantConcentration":
            t = 1
        elif type_top == "constantFlux":
            t = 2
        elif type_top == "constantFluxCyl" or type_top == "fluxCyl":
            t = 3
        elif type_top == "outflow":
            t = 6
        elif type_top == "linear":
            t = 7
        elif type_top == "michaelisMenten":
            t = 8
        elif type_top == "managed":
            t = 9
        else:
            raise Exception('richards.setTopBC_solute(): Top solute type should be "constantConcentration", "constantFlux", "constantFluxCyl",  or "managed" unknown top type {}'.format(type_top))

        self.setParameter(self.param_group + "BC.Top.SType", str(t))
        self.setParameter(self.param_group + "BC.Top.CValue", str(value_top))

        if t == 9:  # managed (nitrogen input)
            if managed:
                assert(len(managed[0]) == len(managed[1]))  # sample points
                self.setParameter("Managed.Time", self.dumux_str(managed[0]))
                self.setParameter("Managed.Input", self.dumux_str(managed[1]))

            else:
                raise Exception('Managed boundary conditions where set, but no managment data were given')

    def setBotBC(self, type_bot:str, value_bot = 0.):
        """ Top boundary conditions are set before creating the problem with SolverBase.initializeProblem 
        
        @param type_bot:
        type_bot ==  "constantPressure", value_bot is a constant pressure [cm] pressure head
        type_bot ==  "constantFlux", value_bot is the constant flux [cm/day] 
        type_bot ==  "constantFluxCyl", value_bot is the constant flux [cm/day] for cylindrical coordinates        
        type_bot ==  "freeDrainage", free drainage, the value of value_bot is ignored                       
        """

        if type_bot == "constantPressure" or type_bot == "pressure":
            b = 1
        elif type_bot == "constantFlux" or type_bot == "flux":
            b = 2
        elif type_bot == "constantFluxCyl" or type_bot == "fluxCyl":
            b = 3
        elif type_bot == "freeDrainage" or type_bot == "fleeFlow" or type_bot == "free":
            b = 5
        elif type_bot == "rootSystem" or type_bot == "rootSystemExact":
            b = 6
        elif type_bot == "noflux"  or type_bot == "noFlux" or type_bot == "no-flux":
            b = 2
            assert value_bot == 0., "setBotBC: value_bot must be zero in case of no flux"
        else:
            raise Exception('richards.setBotBC(): Bottom type should be "constantPressure", "constantFlux", "constantFluxCyl", "noFlux", or "freeDrainage", unknown bottom type {}'.format(type_bot))

        self.setParameter(self.param_group + "BC.Bot.Type", str(b))
        self.setParameter(self.param_group + "BC.Bot.Value", str(value_bot))

    def setBotBC_solute(self, type_bot:str, value_bot:float = 0.):
        """ Top boundary conditions are set before creating the problem with SolverBase.initializeProblem                    
        """
        if type_bot == "constantConcentration":
            b = 1
        elif type_bot == "constantFlux":
            b = 2
        elif type_bot == "constantFluxCyl" or type_bot == "fluxCyl":
            b = 3
        elif type_bot == "outflow":
            b = 6
        elif type_bot == "linear":
            b = 7
        elif type_bot == "michaelisMenten":
            b = 8
        else:
            raise Exception('richards.setBotBC_solute(): Bottom type should be "constantConcentration", "constantFlux", "constantFluxCyl", "outflow", "lineaer" or "michaelisMenten", unknown bottom type {}'.format(type_bot))

        self.setParameter(self.param_group + "BC.Bot.SType", str(b))
        self.setParameter(self.param_group + "BC.Bot.CValue", str(value_bot))

    def setInnerFluxCyl(self, flux):
        """ 
        Sets the flux directly in the problem (problem must be initialized), 
        calls base.setBotBC, @see setBotBC in richards.hh
        
        @param flux      [cm/day] negative means outflow          
        """
        self.base.setBotBC(3, flux)

    def setInnerMatricPotential(self, x):
        """ 
        Sets the dirichlet BC directly in the problem (problem must be initialized), 
        calls base.setBotBC, @see setBotBC in richards.hh
        
        @param flux      [cm] xylem matric potentail          
        """
        self.base.setBotBC(1, x)

    def getInnerFlux(self, eq_idx = 0):
        """ [cm3 / cm2 / day] """
        assert not self.useMoles
        return self.base.getInnerFlux(eq_idx) * 24 * 3600 * 10.  # [kg m-2 s-1] = g / cm2 / day

    def setOuterFluxCyl(self, flux):
        """ 
        sets the flux directly in the problem (problem must be initialized), calls base.setToptBC, @see setToptBC in richards.hh
        @param flux      [cm/day] negative means outflow    
        """
        self.base.setTopBC(3, flux)

    def getOuterFlux(self, eq_idx = 0):
        """ [cm / day]"""
        assert not self.useMoles
        return self.base.getOuterFlux(eq_idx) / 1000 * 24 * 3600 * 100.  # [kg m-2 s-1] / rho = [m s-1] -> cm / day

    def setOuterPressure(self, value_top = 0.):
        """ sets the flux directly in the problem (problem must be initialized), calls base.setToptBC, @see setToptBC in richards.hh"""
        self.base.setTopBC(1, value_top)

    def setInnerBC(self, type_bot:str, value_bot = 0.):
        """ more appropriate name for cylindrical coordinates, calls setBotBC, @see setBotBC """
        self.setBotBC(type_bot, value_bot)

    def setOuterBC(self, type_top:str, value_top = 0.):
        """ more appropriate name for cylindrical coordinates, calls setToptBC, @see setToptBC """
        self.setTopBC(type_top, value_top)

    def setInnerBC_solute(self, type_bot:str, value_bot = 0.):
        """ more appropriate name for cylindrical coordinates, calls setBotBC, @see setBotBC """
        self.setBotBC_solute(type_bot, value_bot)

    def setOuterBC_solute(self, type_top:str, value_top = 0.):
        """ more appropriate name for cylindrical coordinates, calls setToptBC, @see setToptBC """
        self.setTopBC_solute(type_top, value_top)

    def setInnerBCRootSystem(self, params):
        self.base.botValues_ = []

    def setRootSystemBC(self, params):
        """Couples the inner boundary of a cylindrical model to an exact root surface flow 
        @params [x0, x1, kr, kx, l]
               x0, x1        matric potential  
               kr, kx        root system conductivities
               l             length [cm]
        """
        self.base.setRootSystemBC(params)

    def setSoluteTopBC(self, type_top, value_top):
        """  sets the flux directly in the problem (problem must be initialized), calls base.setSToptBC in richards.hh"""
        self.base.setSTopBC(type_top, value_top)
        
    def setSoluteBotBC(self, type_bot, value_bot):
        """  sets the flux directly in the problem (problem must be initialized), calls base.setSToptBC in richards.hh"""
        self.base.setSBotBC(type_bot, value_bot)

    def getInnerHead(self):
        """Gets the pressure head at the inner boundary [cm] """
        shift = 0
        return self.base.getInnerHead(shift)  # -> richards_cyl.hh

    def getInnerSolutes(self, compId = 1):
        """Gets the concentration at the inner boundary [mol/cm3] """    
        CC = np.array(self.getSolution_(compId)).flatten()[0] #mol/mol
        isDissolved = compId <= self.numDissolvedSoluteComp # is the component in the water phase (True) [bool]
        if self.useMoles:
            if isDissolved:
                CC *= self.molarDensityWat_m3/1e6 #mol/cm3 scv
            else:
                CC *= self.bulkDensity_m3/1e6 #mol/cm3 scv
        return CC

    def setSource(self, source_map, eq_idx = 0):
        """Sets the source term as map with global cell index as key, and source as value [cm3/day] """
        self.checkGridInitialized()
        # useMole fraction or mass fraction? 
        if not self.useMoles:
            unitConversion = 1/ 24. / 3600. / 1.e3;  # [cm3/day] -> [kg/s] (richards.hh)
        else:
            if eq_idx == 0:
                molarMassWat = 18. # [g/mol]
                densityWat = 1. #[g/cm3]
                # [mol/cm3] = [g/cm3] /  [g/mol] 
                molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 
                unitConversion = 1/ 24. / 3600.  * molarDensityWat; # [cm3/day] -> [mol/s] 
            else:
                unitConversion = 1/ 24. / 3600. ; # [mol/day] -> [mol/s] 

        for cellId, value in source_map.items(): 
            source_map[cellId] = value * unitConversion     
        self.base.setSource(source_map, eq_idx) # 

    # def applySource(self, dt, source_map, crit_p):
    #     """Sets the source term as map with global cell index as key, and source as value [cm3/day] """
    #     self.checkGridInitialized()
    #     assert not self.useMoles
    #     for key, value in source_map.items():
    #         source_map[key] = value / 24. / 3600. / 1.e3;  # [cm3/day] -> [kg/s]
    #     self.base.applySource(dt * 24.*3600., source_map, self.to_pa(crit_p))

    def setCriticalPressure(self, critical):
        """ Sets the critical pressure to limit flow for boundary conditions constantFlow, constantFlowCyl, and atmospheric """
        self.base.setCriticalPressure(critical)

    def getSolutionHead(self, eqIdx = 0):
        """Gathers the current solution into rank 0, and converts it into a numpy array (Ndof, neq), 
        model dependent units, [Pa, ...]"""
        self.checkGridInitialized()
        return self._map(self._flat0(self.gather(self.base.getSolutionHead(eqIdx))), 0)#.flatten()

    def getSolutionHead_(self, eqIdx = 0):
        """ no mpi version of getSolutionHead() """
        self.checkGridInitialized()
        return np.array(self.base.getSolutionHead(eqIdx))

    def getSolutionHeadAt(self, gIdx, eqIdx = 0):
        """Returns the current solution at a cell index"""
        return self.base.getSolutionHeadAt(gIdx, eqIdx)

    def getKrw(self):
        """Gathers the current solution's saturation into rank 0, and converts it into a numpy array (Nc, 1) [1]"""
        self.checkGridInitialized()
        return self._map(self._flat0(self.gather(self.base.getKrw())), 0)
        
        
    def getSaturation(self):
        """Gathers the current solution's saturation into rank 0, and converts it into a numpy array (Nc, 1) [1]"""
        self.checkGridInitialized()
        return self._map(self._flat0(self.gather(self.base.getSaturation())), 0)
        
    def getWaterContent(self):
        """Gathers the current solution's saturation into rank 0, and converts it into a numpy array (Nc, 1) [1]"""
        self.checkGridInitialized()
        theta= self._map(self._flat0(self.gather(self.base.getWaterContent())), 2)
        
        try:
            if len(theta) > 0:
                assert min(theta) >= self.vg_soil.theta_R
                assert max(theta) <= self.vg_soil.theta_S
        except:
            print('getting out of range theta')
            print(repr(theta), 'miin and max', min(theta), max(theta))
            sat = self.getSaturation_()
            print('saturation', repr(sat), min(sat), max(sat))
            phead = self.getSolutionHead_()
            print('phead', repr(self.getSolutionHead_()), min(phead), max(phead))
            write_file_array('thetaerror',theta, directory_ =self.results_dir, fileType = '.csv')
            raise Exception
        return theta

    def getWaterContent_(self):
        """no mpi version of getWaterContent() """
        self.checkGridInitialized()
        theta = np.array(self.base.getWaterContent())
        try:
            if len(theta) > 0:
                assert min(theta) >= self.vg_soil.theta_R
                assert max(theta) <= self.vg_soil.theta_S
        except:
            print('getting out of range theta', rank)
            print(repr(theta), 'min and max', min(theta), max(theta))
            sat = self.getSaturation_()
            print('saturation', repr(sat), min(sat), max(sat))
            phead = self.getSolutionHead_()
            print('phead', repr(self.getSolutionHead_()), min(phead), max(phead))
            raise Exception
        return theta

    def getWaterVolume(self):
        """Returns total water volume of the domain [cm3]"""
        self.checkGridInitialized()
        assert self.dimWorld != 1
        return self.base.getWaterVolume() * 1.e6  # m3 -> cm3
        
    def getWaterVolumesCyl(self, verbose = False):
        """Returns total water volume of the domain [cm3]"""
        self.checkGridInitialized()
        assert self.dimWorld != 3
        vols = self.getCellVolumes() #cm3 scv
        watCont = self.getWaterContent()# # cm3 wat/cm3 scv
        
        return np.multiply(vols , watCont  )  
        
    def getWaterVolumes(self):
        """Returns total water volume of the domain [cm3]"""
        self.checkGridInitialized()
        vols = self.getCellVolumes() #cm3 scv 
        watCont = self.getWaterContent()  # cm3 wat/cm3 scv
        return np.multiply(vols , watCont  )
            
    def getTotCContent_each(self):
        """ copmute the total content per solute class in each cell of the domain """
        #vols = self.getCellVolumes()#.flatten() #cm3 scv   
        totC = np.array([self.getContent(i+1) for i in range(self.numSoluteComp)])
            
        
        if rank == 0:
            try:
                assert np.array(totC).shape == (self.numSoluteComp,self.numberOfCellsTot)
            except:
                print('totC',totC,totC.shape , (self.numSoluteComp, self.numberOfCellsTot))
                raise Exception
            
        return totC
        
    def getTotCContent(self):
        """ copmute the total content of all solutes in each cell of the domain """
        return self.getTotCContent_each().sum(axis=0)
        
        
    def phaseDensity(self, isDissolved):
        """ returns the density of the phase the solute is in [mol / m3]
            i.e., water phase for disolved solutes, solide phase otherwise
            @param isDissolved: is the solute disolved in water [bool]
        """
        if isDissolved: #mol wat / m3 wat
            return self.molarDensityWat_m3
        else:   # mol scv / m3 scv
            return self.bulkDensity_m3
        
    def getFace2CellIds_(self):
        return np.array(self.base.face2CellIds).max(axis = 0)
         
    # move that to c++
    def getFlux_10cBU(self): 
        """ returns the total inter-cell flux of water [cm3] per cell
            and solutes [mol] during the last @see solve() call
            for one thread 
        """
        assert self.dimWorld == 3
        ff10c_ = self.getFlux_10c_() # get the flux of the local thread for each cell
        f2cidx_ = self.getFace2CellIds_() # see to which cell correspond each face
        ff10c = np.array([sum(ff10c_[np.where(f2cidx_ == idx_)[0]]) for idx_ in set(f2cidx_) if list(f2cidx_).count(idx_) == 6]) # sum values per cell
        # only keep the value is all 6 faces of the cell belong to this thread. Otherwise the sum value is incorrect
        f2cidx = np.array([idx_ for idx_ in set(f2cidx_) if list(f2cidx_).count(idx_) == 6])
        # get global index
        dofind = np.array(self.base.getDofIndices())
        f2cidx_g = [] # give correct inner shape 
        if len(f2cidx) > 0:
            f2cidx_g = dofind[f2cidx] 
            
            
        # gather
        f2cidx_gAll = self._flat0(self.gather(f2cidx_g)) # 
        ff10c_All = self._flat0(self.gather(ff10c)) #  
        if rank == 0:
            # get arrays 
            # ff10c_All = np.vstack([i for i in ff10c_All if len(i) >0]) # 'if len(i) >0' in case some threads have no cells
            # f2cidx_gAll = list( np.concatenate([i for i in f2cidx_gAll],dtype=object)) # empty lists ([], for threads with no cells) are automatically taken out
            # some cells are computed on several threads, so take out duplicates
            f2cidx_gAll_unique = np.array(list(set(f2cidx_gAll)),dtype=int)
            ff10c_All_unique = np.array([ff10c_All[
              max(np.where(
                np.array(f2cidx_gAll, dtype=int) == idx_)[0])
                ] for idx_ in f2cidx_gAll_unique ]) # select value from one of the threads which simulate all the faces of the dumux cell
            
            flux10cCell = np.transpose(ff10c_All_unique) # [comp][cell]
            assert flux10cCell.shape == (self.numComp ,self.numberOfCellsTot)
            import sys
            np.set_printoptions(threshold=sys.maxsize)
            print('f2cidx_g',f2cidx_g)
            print('ff10c',ff10c)
            print('f2cidx_gAll',f2cidx_gAll)
            print('flux10cCell',flux10cCell)
            # mol to cm3 for water
            molarMassWat = 18. # [g/mol]
            densityWat = 1. #[g/cm3]
            # [mol/cm3] = [g/cm3] /  [g/mol] 
            molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 
            flux10cCell[0] /=  molarDensityWat  
        
        else:
            flux10cCell = None
                
        #raise Exception
        return flux10cCell
    
    def getFlux_10c(self): 
        """ returns the total inter-cell flux of water [cm3] per cell
            and solutes [mol] during the last @see solve() call
            for one thread 
        """
        assert self.dimWorld == 3
        ff10c_ = self.getFlux_10c_() # get the flux of the local thread for each cell
        f2cidx_ = self.getFace2CellIds_() # see to which cell correspond each face
        ff10c = np.array([sum(ff10c_[np.where(f2cidx_ == idx_)[0]]) for idx_ in set(f2cidx_) if list(f2cidx_).count(idx_) == 6]) # sum values per cell
        # only keep the value is all 6 faces of the cell belong to this thread. Otherwise the sum value is incorrect
        f2cidx = np.array([idx_ for idx_ in set(f2cidx_) if list(f2cidx_).count(idx_) == 6])
        # get global index
        dofind = np.array(self.base.getDofIndices())
        f2cidx_g = [] # give correct inner shape 
        if len(f2cidx) > 0:
            f2cidx_g = dofind[f2cidx] 
        # gather
        f2cidx_gAll = self.gather(f2cidx_g) # 
        ff10c_All = self.gather(ff10c) #  
        if rank == 0:
            # get arrays 
            ff10c_All = np.vstack([i for i in ff10c_All if len(i) >0]) # 'if len(i) >0' in case some threads have no cells
            f2cidx_gAll = list( np.concatenate([i for i in f2cidx_gAll],dtype=object)) # empty lists ([], for threads with no cells) are automatically taken out
            # some cells are computed on several threads, so take out duplicates
            f2cidx_gAll_unique = np.array(list(set(f2cidx_gAll)),dtype=int)
            ff10c_All_unique = np.array([ff10c_All[
              max(np.where(
                np.array(f2cidx_gAll, dtype=int) == idx_)[0])
                ] for idx_ in f2cidx_gAll_unique ]) # select value from one of the threads which simulate all the faces of the dumux cell
            
            flux10cCell = np.transpose(ff10c_All_unique) # [comp][cell]
            assert flux10cCell.shape == (self.numComp ,self.numberOfCellsTot)
            
            # mol to cm3 for water
            molarMassWat = 18. # [g/mol]
            densityWat = 1. #[g/cm3]
            # [mol/cm3] = [g/cm3] /  [g/mol] 
            molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 
            flux10cCell[0] /=  molarDensityWat  
        
        else:
            flux10cCell = None
                
        
        return flux10cCell   
    
    def getFlux_10cBU2(self): 
        """ returns the total inter-cell flux of water [cm3] per cell
            and solutes [mol] during the last @see solve() call
            for one thread 
        """
        assert self.dimWorld == 3
        ff10c_ = self.getFlux_10c_() # get the flux computed on local thread for cells , [face][comp]
        f2cidx_ = self.getFace2CellIds_() # see to which cell correspond each face (local index) , [face]
        ff10c = np.array([sum(ff10c_[np.where(f2cidx_ == idx_)[0]]) for idx_ in set(f2cidx_) if list(f2cidx_).count(idx_) == 6]) # sum values per cell, [cell][comp]
        # only keep the value is all 6 faces of the cell belong to this thread. Otherwise the sum value is incorrect
        f2cidx = np.array([idx_ for idx_ in set(f2cidx_) if list(f2cidx_).count(idx_) == 6])#, [cell]
        # get global index
        dofind = np.array(self.base.getDofIndices())
        f2cidx_g = [] 
        if len(f2cidx) > 0:
            f2cidx_g = dofind[f2cidx] # local to global index , [cell]
            
        flux10cCell_ = np.transpose(ff10c) # [comp][cell]
        
        # gather
        f2cidx_gAll = self._flat0(self.gather(f2cidx_g) ,dtype_=object)# 
        #print('before gather',rank, flux10cCell_)
        flux10cCell = np.array([self._flat0(self.gather(flux10cCell__)
                                           ) for flux10cCell__ in flux10cCell_])#  
        if rank == 0:
            #print('f2cidx_gAll',rank, f2cidx_gAll.shape,len(f2cidx_gAll),
            #      len(np.unique(f2cidx_gAll)))#,np.array(flux10cCell).shape)
            #print('flux10cCell',flux10cCell,flux10cCell[0])
            # no need for _map: each threads has only data for cells in Dune::Partitions::interior
            assert len(f2cidx_gAll) == self.numberOfCellsTot
            #print('(np.unique(f2cidx_gAll) == f2cidx_gAll).all()',(np.unique(f2cidx_gAll) == f2cidx_gAll).all(),f2cidx_gAll[np.unique(f2cidx_gAll) != f2cidx_gAll])
            #print('np.unique(f2cidx_gAll)',np.unique(f2cidx_gAll))
            #print('f2cidx_gAll',f2cidx_gAll)
            #assert (np.unique(f2cidx_gAll) == f2cidx_gAll).all()
            #import sys
            #np.set_printoptions(threshold=sys.maxsize)
            
            f2cidx_gAll = np.array(f2cidx_gAll,dtype=int)
            try:
                flux10cCell =np.array([ flux10cCell_[f2cidx_gAll] for flux10cCell_ in flux10cCell])
            except:
                print(type(f2cidx_gAll),f2cidx_gAll)
                raise Exception
            # get arrays 
            #ff10c_All = np.vstack([i for i in ff10c_All if len(i) >0]) # 'if len(i) >0' in case some threads have no cells
            #f2cidx_gAll = list( np.concatenate([i for i in f2cidx_gAll],dtype=object)) # empty lists ([], for threads with no cells) are automatically taken out
            # some cells are computed on several threads, so take out duplicates
            #f2cidx_gAll_unique = np.array(list(set(f2cidx_gAll)),dtype=int)
            #ff10c_All_unique = np.array([ff10c_All[
            #  max(np.where(
            #    np.array(f2cidx_gAll, dtype=int) == idx_)[0])
            #    ] for idx_ in f2cidx_gAll_unique ]) # select value from one of the threads which simulate all the faces of the dumux cell
            
            #flux10cCell = np.transpose(ff10c_All_unique) # [comp][cell]
            assert flux10cCell.shape == (self.numComp ,self.numberOfCellsTot)
            
            import sys
            np.set_printoptions(threshold=sys.maxsize)
            print('f2cidx_g',f2cidx_g)
            print('ff10c',ff10c)
            print('f2cidx_gAll',f2cidx_gAll)
            print('flux10cCell',flux10cCell)
            # mol to cm3 for water
            molarMassWat = 18. # [g/mol]
            densityWat = 1. #[g/cm3]
            # [mol/cm3] = [g/cm3] /  [g/mol] 
            molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 
            flux10cCell[0] /=  molarDensityWat  
        
        else:
            flux10cCell = None
                
        #raise Exception
        return flux10cCell
        
    
    # move that to c++
    def getFlux_10cBU3(self): 
        """ returns the total inter-cell flux of water [cm3] per cell
            and solutes [mol] during the last @see solve() call
            for one thread 
        """
        assert self.dimWorld == 3
        ff10c_ = self.getFlux_10c_() # get the flux of the local thread for each cell
        f2cidx_ = self.getFace2CellIds_() # see to which cell correspond each face
        ff10c = np.array([sum(ff10c_[np.where(f2cidx_ == idx_)[0]]) for idx_ in set(f2cidx_) if list(f2cidx_).count(idx_) == 6]) # sum values per cell
        # only keep the value is all 6 faces of the cell belong to this thread. Otherwise the sum value is incorrect
        f2cidx = np.array([idx_ for idx_ in set(f2cidx_) if list(f2cidx_).count(idx_) == 6])
        # get global index
        dofind = np.array(self.base.getDofIndices())
        f2cidx_g = [] # give correct inner shape 
        if len(f2cidx) > 0:
            f2cidx_g = dofind[f2cidx] 
            
            
        # gather
        f2cidx_gAll = self._flat0(self.gather(f2cidx_g)) # 
        ff10c_All = self._flat0(self.gather(ff10c)) #  
        if rank == 0:
            # get arrays 
            # ff10c_All = np.vstack([i for i in ff10c_All if len(i) >0]) # 'if len(i) >0' in case some threads have no cells
            # f2cidx_gAll = list( np.concatenate([i for i in f2cidx_gAll],dtype=object)) # empty lists ([], for threads with no cells) are automatically taken out
            # some cells are computed on several threads, so take out duplicates
            f2cidx_gAll_unique = np.array(list(set(f2cidx_gAll)),dtype=int)
            ff10c_All_unique = np.array([ff10c_All[
              max(np.where(
                np.array(f2cidx_gAll, dtype=int) == idx_)[0])
                ] for idx_ in f2cidx_gAll_unique ]) # select value from one of the threads which simulate all the faces of the dumux cell
            
            flux10cCell = np.transpose(ff10c_All_unique) # [comp][cell]
            assert flux10cCell.shape == (self.numComp ,self.numberOfCellsTot)
            
            # mol to cm3 for water
            molarMassWat = 18. # [g/mol]
            densityWat = 1. #[g/cm3]
            # [mol/cm3] = [g/cm3] /  [g/mol] 
            molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 
            flux10cCell[0] /=  molarDensityWat  
        
        else:
            flux10cCell = None
                
        
        return flux10cCell
    
    def getFlux_10c_(self):#  mol        
        """ returns the total inter-cell flux of water and solutes [mol] per face
            during the last @see solve() call
            for one thread 
            used by @see getFlux_10c()
        """
        # get the data saved during solve()
        inFluxes = np.array(self.base.inFluxes) #[ mol / s], flux rate at each dumux sub-timestep  [time][face][comp]
        inFluxes_ddt = np.array(self.base.inFluxes_ddt)# s, duration of each sub timestep
        
        
        # go from [ mol / s] to [mol] 
        inFluxes_tot = np.array([ np.array([xxx * inFluxes_ddt[isrcs] for xxx in inFluxes[isrcs] ]) for isrcs in range(len(inFluxes_ddt))])
        
        # sum to get the total flow per element during the whole solve() simulation
        inFluxes_tot = inFluxes_tot.sum(axis = 0) # cm3 or mol, [face][comp]
        # important => need minus sign
        return - inFluxes_tot
        
        
    def getSource_10c(self):
        """ returns the total inter-cell flux of water [cm3] per cell
            and solutes [mol] during the last @see solve() call
            for one thread 
        """
        # get total source and cell volumes for each thread
        src10c =  self._map(self._flat0(self.gather(self.getSource_10c_())), 0) * 1e-6 # [ mol / m^3]   => [ mol / cm^3]   
        vols = self.getCellVolumes() # [ cm^3]   
        if rank == 0:
            # from mol/m^3 to mol
            src10c = np.array([src10_ * vols[cellidx] for cellidx, src10_ in enumerate(src10c)])
            src10c = np.transpose(src10c) # [comp][cell]

            molarMassWat = 18. # [g/mol]
            densityWat = 1. #[g/cm3]
            # [mol/cm3] = [g/cm3] /  [g/mol] 
            molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 
            src10c[0] /=  molarDensityWat # mol to cm3 for water 
        
            assert src10c.shape == (self.numComp ,self.numberOfCellsTot)#numComp + water 
        return src10c
    
    def getSource_10c_(self):# [ mol / m^3]                 
        """ returns the total source of water or solute [mol/m^3] per cell
            during the last @see solve() call
            for one thread 
            used by @see getSource_10c()
        """
        # get values saved by dumux for each sub time step
        inSources = np.array(self.base.inSources) #[ mol / (m^3 \cdot s)] 
        inFluxes_ddt = np.array(self.base.inFluxes_ddt)# s
        
        # sum to go from rates to total source
        inSources_tot = np.array([ np.array([xxx * inFluxes_ddt[isrcs] for idc, xxx in enumerate(inSources[isrcs] )]) for isrcs in range(len(inFluxes_ddt))])
        # sum over all the sub time steps
        inSources_tot = inSources_tot.sum(axis = 0) # mol, [cell][comp]
        
        return inSources_tot
            
    def getConcentration(self,idComp):      
        """ returns the concentraiton of a component (solute)
            @param idcomp: index of the component (> 0) [int]
        """
        isDissolved = idComp <= self.numDissolvedSoluteComp # is the component in the water phase (True) [bool]
        C_ = self.getSolution(idComp) # mol/mol wat or mol/mol scv
        if not isDissolved:
            if self.useMoles:
                C_ *= self.bulkDensity_m3 #mol/m3 scv
            return C_  /1e6  #mol/cm3 scv
        if self.useMoles:
            C_ *= self.molarDensityWat_m3 # mol/m3 wat            
        return C_ /1e6 #mol/cm3 scv
        
    def getContent(self,idComp):
        """ returns the content of a component (solute)
            @param idcomp: index of the component (> 0) [int]
        """
        assert idComp > 0 # do not use for water content
        isDissolved = idComp <= self.numDissolvedSoluteComp # is the component in the water phase (True) [bool]
        vols = (1/  1e6)*self.getCellVolumes_() #m3 scv            
        
        if idComp <= self.numComp:
            C_ = self.getSolution_(idComp)#.flatten() # mol/mol or g/g 
        else:
            print('wrong idComp', idComp)
            raise Exception
            
        try:
            assert (C_ >= 0.).all()
        except:
            print('getContent',idComp, isDissolved,min( vols),max( vols),min( C_),max( C_))
            raise Exception
            
        if not isDissolved:
            if self.useMoles:
                C_ *= self.bulkDensity_m3 #mol/m3 scv
            return self._map(self._flat0(self.gather(np.multiply(vols , C_  ))),0)
        
        watCont = self.getWaterContent_()#.flatten() # m3 wat/m3 scv
        if self.useMoles:
            C_ *= self.molarDensityWat_m3 # mol/m3 wat
            
        return self._map(self._flat0(self.gather(np.multiply(np.multiply(vols , watCont) , C_ ))),0)
          

    def setSolution_(self, sol, eqIdx = 1):
        """nompi version of  """
        self.checkGridInitialized()
        assert max_rank == 1
        return np.array(self.base.setSolution(sol, eqIdx))
        
    def getVeclocity1D(self):
        """Returns the Darcy velocities [cm/day] TODO not working! """
        self.checkGridInitialized()
        assert not self.useMoles
        return np.array(self.base.getVelocity1D()) * 100.*24 * 3600  # m/s -> cm/day

    def setRegularisation(self, pcEps, krEps):
        self.checkGridInitialized()
        """ Van Genuchten regularisation parameters"""
        self.base.setRegularisation(pcEps, krEps)

    def writeDumuxVTK(self, file_name):
        """Uses the Dumux VTK writer to write the current model output"""
        self.checkGridInitialized()
        return self.base.writeDumuxVTK(file_name)

    @staticmethod
    def to_pa(ph):
        return 1.e5 + ph / 100 * 1000. * 9.81;

    @staticmethod
    def to_head(p):
        return (p - 1.e5) * 100. / 1000. / 9.81;

    @property
    def numberOfCells(self):
        """ In case of a rectangular domain, the number of cells in each direction, 
            not valid if the grid was given by a file (read only) """
        return self.base.numberOfCells

