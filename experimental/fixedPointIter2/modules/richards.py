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

    def __init__(self, base):
        """ 
            @param base is the C++ base class that is wrapped. 
            @param usemoles: do we use moles (mol) or mass (g) in dumux
        """
        super().__init__(base)
        self.soils = []
        self.param_group = "Soil."
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
            
        self.qr, self.qs, self.alpha, self.n, self.ks = [], [], [], [], []

        for soil in soils:
            assert len(soil) == 5, "setVGParameters, soil need to be a list of 5 parameters: qr, qs, alpha, n, ks "
            self.qr.append(soil[0])
            self.qs.append(soil[1])
            self.alpha.append(soil[2])
            self.n.append(soil[3])
            self.ks.append(soil[4])

        self.setParameter(self.param_group + "VanGenuchten.Qr", self.dumux_str(self.qr))
        self.setParameter(self.param_group + "VanGenuchten.Qs", self.dumux_str(self.qs))
        self.setParameter(self.param_group + "VanGenuchten.Alpha", self.dumux_str(self.alpha))
        self.setParameter(self.param_group + "VanGenuchten.N", self.dumux_str(self.n))
        self.setParameter(self.param_group + "VanGenuchten.Ks", self.dumux_str(self.ks))

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
        try:
            assert len(type_top) == len(value_top)
        except:
            raise ValueError(f'richards::setSoluteTopBC: Length of the solute boundary condition arrays for types and values ({len(type_top)}, {len(value_top)}) different from number of dissolved solutes {self.numDissolvedSoluteComp}.')
        print('setSoluteTopBC',type_top, value_top)
        self.base.setSTopBC(type_top, value_top)
        
    def setSoluteBotBC(self, type_bot, value_bot):
        """  sets the flux directly in the problem (problem must be initialized), calls base.setSToptBC in richards.hh"""
        try:
            assert len(type_bot) == len(value_bot)
        except:
            raise ValueError(f'richards::setSoluteBotBC: Length of the solute boundary condition arrays for types and values ({len(type_bot)} and {len(value_bot)}) different from number of dissolved solutes {self.numDissolvedSoluteComp}.')
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
            if eq_idx == 0:
                unitConversion = 1/ 24. / 3600. / 1.e3;  # [cm3/day] -> [kg/s] (richards.hh)
            else:
                unitConversion = 1. / 24. / 3600.  
            
        else:
            if eq_idx == 0:
                molarMassWat = self.molarMassWat # [g/mol]
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
                assert min(theta) >= self.qr[0]
                assert max(theta) <= self.qs[0]
        except:
            print('getting out of range theta')
            print(repr(theta), 'miin and max', min(theta), max(theta))
            #sat = self.getSaturation_()
            #print('saturation', repr(sat), min(sat), max(sat))
            phead = self.getSolutionHead_()
            print('phead', repr(self.getSolutionHead_()), min(phead), max(phead))
            #write_file_array('thetaerror',theta, directory_ =self.results_dir, fileType = '.csv')
            raise Exception
        return theta

    def getWaterContent_(self):
        """no mpi version of getWaterContent() """
        self.checkGridInitialized()
        theta = np.array(self.base.getWaterContent())
        try:
            if len(theta) > 0:
                assert min(theta) >= self.qr[0]
                assert max(theta) <= self.qs[0]
        except:
            print('getting out of range theta', rank)
            print(repr(theta), 'min and max', min(theta), max(theta))
            #sat = self.getSaturation_()
            #print('saturation', repr(sat), min(sat), max(sat))
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
            
    
    def getInnerFlow(self, eqIdx = 0, length = None):
        """ mean flow rate at the inner boundary [cm3 / day] """
        assert self.dimWorld == 1
        pos = np.array(self.base.getCylFaceCoordinates())
        flux = self.getBoundaryFluxesPerFace_(eqIdx, length)
        return flux[pos == min(pos)]
        
    def getOuterFlow(self, eqIdx = 0, length = None):
        """ mean flow rate at the outer boundary [cm3 / day] """
        assert self.dimWorld == 1
        pos = np.array(self.base.getCylFaceCoordinates())
        flux = self.getBoundaryFluxesPerFace_(eqIdx, length)
        return flux[pos == max(pos)]
    
    
    def getCss1_(self): # mol/cm3 scv
        ''' att: only valid for the 4c problem
            todo: update
        '''
        C_S_W = self.getConcentration_(1)
        C_S_S2 = self.getSolution_(2)
        css1 = np.array([self.base.computeCSS1(self.bulkDensity_m3,C_S_W[i]*1e6, i) for i, idx in enumerate(self.getCellIndices_())])*1e-6
        css1[C_S_S2>0] = 0.
        return css1
    
    def getCss1(self): # mol/cm3 scv
        ''' att: only valid for the 4c problem
            todo: update
        '''        
        return self._map(self._flat0(self.gather(self.getCss1_())), 2)
            
    def getTotCContent_each(self):
        """ copmute the total content per solute class in each cell of the domain """
        #vols = self.getCellVolumes()#.flatten() #cm3 scv   
        totC = np.array([self.getContent(i+1) for i in range(self.numSoluteComp)])
        vols = self.getCellVolumes() 
        css1 = self.getCss1() * vols
        idTemps = css1 != 0
        totC[1][idTemps] = css1[idTemps]
        #if self.dimWorld==1:
        #    assert totC[1].sum() == 0.
        #    totC[1] = self.getCss1() * vols
        
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

    def getSource(self, eqIdx = 0):
        """ Gathers the net sources for each cell into rank 0 as a map with global index as key [kg / cm3 / day]"""
        self.checkGridInitialized()
        return self._map(self._flat0(comm.gather(self.getSource_(eqIdx), root = 0)), 2) 

    def getSource_(self, eqIdx = 0):
        """nompi version of """
        self.checkGridInitialized()
        unitChange = 24 * 3600 / 1e6 # [ kgOrmol/m3/s -> kgOrmol/cm3/day]
        if eqIdx == 0:
            if self.useMoles:
                unitChange *=  self.molarMassWat # [mol/cm3/day -> cm3/cm3/day]
            else:
                unitChange *= 1e3  # [kg/cm3/day -> cm3/cm3/day]
        return np.array(self.base.getScvSources()[eqIdx]) * unitChange     
        
    def getReactions(self, eqIdx = 0):
        """ Gathers the net sources for each cell into rank 0 as a map with global index as key [kg / cm3 / day]"""
        self.checkGridInitialized()
        return self._map(self._flat0(comm.gather(self.getReactions_(eqIdx), root = 0)), 2) 

    def getReactions_(self, eqIdx = 0):
        """nompi version of """
        self.checkGridInitialized()
        unitChange = 24 * 3600 / 1e6 # [ kgOrmol/m3/s -> kgOrmol/cm3/day]
        if eqIdx == 0:
            if self.useMoles:
                unitChange *=  self.molarMassWat # [mol/cm3/day -> cm3/cm3/day]
            else:
                unitChange *= 1e3  # [kg/cm3/day -> cm3/cm3/day]
        return np.array(self.base.getScvReactions()[eqIdx]) * unitChange
    
    def getFluxesPerCell(self, eqIdx = 0, length = None):
        """ Gathers the net sources fir each cell into rank 0 as a map with global index as key [cm3/ day or kg/ day or mol / day]"""
        self.checkGridInitialized()
        return self._flat0(comm.gather(self.getFluxesPerCell_(eqIdx, length), root = 0))#self._map(, 2)
 
    def getFluxesPerCell_(self, eqIdx = 0, length = None):
        """nompi version of """
        self.checkGridInitialized()
        scvfFluxes = np.array(self.getFluxesPerFace_(eqIdx, length)  ) # [cm3/day] 
        face2CellIdx = np.array(self.getFace2CellIdx_())
        scvFluxes = np.array([
             sum(scvfFluxes[face2CellIdx == cellindx]) for cellindx in np.unique(face2CellIdx) 
            ]) 
        return scvFluxes
            
    def getFluxesPerFace_(self, eqIdx = 0, length = None):
        """nompi version of 
            [cm3/day] 
        """
        self.checkGridInitialized() 
        scvfFluxes = self.getBoundaryFluxesPerFace_(eqIdx, length)      
        scvfFluxes += self.getCell2CellFluxesPerFace_(eqIdx, length)  
        return scvfFluxes 
        
    def getBoundaryFluxesPerFace_(self, eqIdx = 0, length = None):
        """  [cm3/day] """
        self.checkGridInitialized()
        unitChange = 24 * 3600 / 1e4 # [ kgOrmol/m2/s -> kgOrmol/cm2/day]
        if eqIdx == 0:
            if self.useMoles:
                unitChange *=  self.molarMassWat # [mol/cm2/day -> cm3/cm2/day]
            else:
                unitChange *= 1e3  # [kg/cm2/day -> cm3/cm2/day]
        surfaces = self.getFaceSurfaces_(length)  # cm2 
        notGhost = surfaces > 0
        return (np.array(self.base.getScvfBoundaryFluxes()[eqIdx]) * surfaces)[notGhost] * unitChange
        
    def getCell2CellFluxesPerFace_(self, eqIdx = 0, length = None):
        """  [cm3/day] """
        self.checkGridInitialized()
        unitChange = 24 * 3600 # [ kgOrmol/s -> kgOrmol/day]
        if eqIdx == 0:
            if self.useMoles:
                unitChange *=  self.molarMassWat # [mol/day -> cm3/day]
            else:
                unitChange *= 1e3  # [kg/day -> cm3/day]
        surfaces = self.getFaceSurfaces_(length)        
        notGhost = surfaces > 0
        return np.array(self.base.getScvfInnerFluxes()[eqIdx])[notGhost] * unitChange
        
    def getFaceSurfaces_(self, length = None):
        if self.dimWorld == 1:
            if length is None:
                raise Exception('getFaceSurfaces_: length parameter required to get flux of 1d domain')
            return np.array(self.base.getCylFaceCoordinates()) * 100 * 2 * np.pi * length
        else:
            return np.array( self.base.getFaceSurfaces()) * 1e4 # cm2
        
    def getFace2CellIdx_(self):
        face2CellIdx = np.array(self.base.getFace2CellIdx())
        surfaces = np.array( self.base.getFaceSurfaces()) 
        notGhost = surfaces > 0
        return face2CellIdx[notGhost]
            
    def getConcentration_(self,idComp):      
        """ returns the concentraiton of a component (solute)
            @param idcomp: index of the component (> 0) [int]
        """
        isDissolved = idComp <= self.numDissolvedSoluteComp # is the component in the water phase (True) [bool]
        C_ = self.getSolution_(idComp) # mol/mol wat or mol/mol scv
        if not isDissolved:
            if self.useMoles:
                C_ *= self.bulkDensity_m3 #mol/m3 scv
            return C_  /1e6  #mol/cm3 scv
        if self.useMoles:
            C_ *= self.molarDensityWat_m3 # mol/m3 wat            
        return C_ /1e6 #mol/cm3 scv
        
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

