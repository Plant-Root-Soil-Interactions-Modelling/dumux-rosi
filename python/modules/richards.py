import sys
sys.path.append("../../../CPlantBox/src")
import functional.van_genuchten as vg  
from solverbase import SolverWrapper
import numpy as np
from mpi4py import MPI; comm = MPI.COMM_WORLD; size = comm.Get_size(); rank = comm.Get_rank()


class RichardsWrapper(SolverWrapper):
    """ 
    Adds functionality specifically for Richards equation             
    Wraps passing parameters to Dumux for VanGenuchtenParameter, IC, BC, 
    """

    def __init__(self, base):
        super().__init__(base)
        self.soils = []
        self.param_group = "Soil."

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
        self.vg_soils = []
        for soil in soils:
            assert len(soil) == 5, "setVGParameters, soil need to be a list of 5 parameters: qr, qs, alpha, n, ks "
            qr.append(soil[0])
            qs.append(soil[1])
            alpha.append(soil[2])
            n.append(soil[3])
            ks.append(soil[4])
            self.vg_soils.append(vg.Parameters(soil))

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
        if self.checkProblemInitialized(throwError = False):
            raise Exception('use setLayersZ before calling initialize()')

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
            
        if self.checkProblemInitialized(throwError = False):
            raise Exception('use setICZ_solute before calling initialize()')

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

    def setTopBC(self, type_top, value_top:float = 0., climate:list = []):
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
        assert isinstance(type_top, (str, int))
        if isinstance(type_top, str):
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
            type_top = t

        if self.checkProblemInitialized(throwError = False):
            self.base.setTopBC(type_top, value_top)
        else:
            self.setParameter(self.param_group + "BC.Top.Type", str(type_top))
            self.setParameter(self.param_group + "BC.Top.Value", str(value_top))

            if t == 4:  # atmospheric
                if climate:
                    assert(len(climate[0]) == len(climate[1]))  # sample points
                    self.setParameter("Climate.Time", self.dumux_str(climate[0]))
                    self.setParameter("Climate.Precipitation", self.dumux_str(climate[1]))  # TODO confusing name (should be Evapotranspiration)

                else:
                    raise Exception('richards.setTopBC(): Atmospheric boundary conditions where set, but no climatic data were given')

    def setTopBC_solute(self, type_top:list, value_top:list = [], managed:list = []):
        if not isinstance(type_top, (list, type(np.array([])))):
            type_top = [type_top]
        if not isinstance(value_top, (list, type(np.array([])))):
            value_top = [value_top]
        if len(value_top) == 0:
            value_top = [0. for i in type_top]
            
        assert isinstance(type_top[0], (str, int))
        for id_tt, tt in enumerate(type_top):
            if isinstance(tt, str):
                if tt == "constantConcentration":
                    t = 1
                elif tt == "constantFlux":
                    t = 2                
                elif tt == "noflux"  or tt == "noFlux" or tt == "no-flux":
                    t = 2
                    if not value_top[id_tt] == 0:
                        print("RichardsWrapper.setTopBC_solute() warning: value_top must be zero in case of no flux")
                    value_top[id_tt] = 0.
                elif tt == "constantFluxCyl" or tt == "fluxCyl":
                    t = 3
                elif tt == "outflow":
                    t = 6
                elif tt == "linear":
                    t = 7
                elif tt == "michaelisMenten":
                    t = 8
                elif tt == "managed":
                    t = 9
                else:
                    raise Exception('richards.setTopBC_solute(): Top solute type should be "constantConcentration", "constantFlux", "constantFluxCyl",  or "managed" unknown top type {}'.format(type_top))
                type_top[id_tt] = t

        if self.checkProblemInitialized(throwError = False):

            self.base.setSTopBC(type_top, value_top)
        else:
            self.setParameter(self.param_group + "BC.Top.SType", self.dumux_str(type_top))
            self.setParameter(self.param_group + "BC.Top.CValue", self.dumux_str(value_top))

            if t == 9:  # managed (nitrogen input)
                if managed:
                    assert(len(managed[0]) == len(managed[1]))  # sample points
                    self.setParameter("Managed.Time", self.dumux_str(managed[0]))
                    self.setParameter("Managed.Input", self.dumux_str(managed[1]))

                else:
                    raise Exception('Managed boundary conditions where set, but no managment data were given')

    def setBotBC(self, type_bot, value_bot = 0.):
        """ Top boundary conditions are set before creating the problem with SolverBase.initializeProblem 
        
        @param type_bot:
        type_bot ==  "constantPressure", value_bot is a constant pressure [cm] pressure head
        type_bot ==  "constantFlux", value_bot is the constant flux [cm/day] 
        type_bot ==  "constantFluxCyl", value_bot is the constant flux [cm/day] for cylindrical coordinates        
        type_bot ==  "freeDrainage", free drainage, the value of value_bot is ignored                       
        """
        assert isinstance(type_bot, (str, int))
        if isinstance(type_bot, str):
            if type_bot == "constantPressure" or type_bot == "pressure" or type_bot == "matricPotential" or type_bot == "potential":
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
                if not value_bot == 0:
                    print("RichardsWrapper.setBotBC() warning: value_bot must be zero in case of no flux")
                value_bot = 0.
                # assert value_bot == 0., "setBotBC: value_bot must be zero in case of no flux"
            else:
                raise Exception('richards.setBotBC(): Bottom type should be "constantPressure", "constantFlux", "constantFluxCyl", "noFlux", or "freeDrainage", unknown bottom type {}'.format(type_bot))
            type_bot = b
        if self.checkProblemInitialized(throwError = False):
            self.base.setBotBC(type_bot, value_bot)
        else:
            self.setParameter(self.param_group + "BC.Bot.Type", str(type_bot))
            self.setParameter(self.param_group + "BC.Bot.Value", str(value_bot))

    def setBotBC_solute(self, type_bot:list, value_bot:list = []):
        """ Top boundary conditions are set before creating the problem with SolverBase.initializeProblem                    
        """
        if not isinstance(type_bot, (list, type(np.array([])))):
            type_bot = [type_bot]
        if not isinstance(value_bot, (list, type(np.array([])))):
            value_bot = [value_bot]
        if len(value_bot) == 0:
            value_bot = [0. for i in type_bot]
        assert isinstance(type_bot[0], (str, int))
        for id_tb, tb in enumerate(type_bot):
            if isinstance(tb, str):
                if tb == "constantConcentration":
                    b = 1
                elif tb == "constantFlux":
                    b = 2
                elif tb == "constantFluxCyl" or tb == "fluxCyl":
                    b = 3
                elif tb == "noflux"  or type_bot == "noFlux" or type_bot == "no-flux":
                    b = 2
                    if not value_bot[id_tb] == 0:
                        print("RichardsWrapper.setBotBC() warning: value_bot must be zero in case of no flux")
                    value_bot[id_tb] = 0.
                elif tb == "outflow":
                    b = 6
                elif tb == "linear":
                    b = 7
                elif tb == "michaelisMenten":
                    b = 8
                else:
                    raise Exception('richards.setBotBC_solute(): Bottom type should be "constantConcentration", "constantFlux", "constantFluxCyl", "outflow", "lineaer" or "michaelisMenten", unknown bottom type {}'.format(type_bot))
                type_bot[id_tb] = b
        

        if self.checkProblemInitialized(throwError = False):
            self.base.setSBotBC(type_bot, value_bot)
        else:
            self.setParameter(self.param_group + "BC.Bot.SType", self.dumux_str(type_bot))
            self.setParameter(self.param_group + "BC.Bot.CValue", self.dumux_str(value_bot))

    def setInnerFluxCyl(self, flux):
        """ 
        Sets the flux directly in the problem (problem must be initialized),  
        
        @param flux      [cm/day] negative means outflow          
        """
        self.setBotBC(3, flux)

    def setInnerMatricPotential(self, x):
        """ 
        Sets the dirichlet BC directly in the problem (problem must be initialized),  
        
        @param flux      [cm] xylem matric potentail          
        """
        self.setBotBC(1, x)

    def getBCFlows(self ):          
        """returns the total boundary flow simulated during the last dumux solumation
        """
        assert self.dimWorld == 1
        # water: [cm3/d] 
        # solute: [g/d] or [mol/d] 
        Q_in_m  = np.array([self.getInnerFlow(nc, self.length) for nc in range(self.numFluidComp)])
        Q_out_m = np.array([self.getOuterFlow(nc, self.length) for nc in range(self.numFluidComp)]) 
        
        return(Q_in_m, Q_out_m)
        
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

    def getInnerFlux(self, eq_idx = 0):
        """ instantaneous flow rate at the inner boundary [cm3 / cm2 / day] """
        return self.base.getInnerFlux(eq_idx) * 24 * 3600 * 10.  # [kg m-2 s-1] = g / cm2 / day

    def setOuterFluxCyl(self, flux):
        """ 
        sets the flux 
        @param flux      [cm/day] negative means outflow    
        """
        self.setTopBC(3, flux)

    def getOuterFlux(self, eq_idx = 0):
        """ instantaneous flow rate at the outer boundary  [cm / day]"""
        return self.base.getOuterFlux(eq_idx) / 1000 * 24 * 3600 * 100.  # [kg m-2 s-1] / rho = [m s-1] -> cm / day

    def setOuterPressure(self, value_top = 0.):
        """ sets the flux """
        self.setTopBC(1, value_top)

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
        """  sets the flux directly in the problem """
        self.setTopBC_solute(type_bot, value_bot)

    def setSoluteBotBC(self, type_bot, value_bot):
        """  sets the flux directly in the problem """
        self.setBotBC_solute(type_bot, value_bot)

    def getInnerHead(self, shift = 0):
        """Gets the pressure head at the inner boundary [cm] """
        return self.base.getInnerHead(shift)  # -> richards_cyl.hh

    def getInnerSolutes(self, shift = 0):
        """Gets the concentration at the inner boundary [cm] """
        return self.base.getInnerSolutes(shift)

    def getSource(self, eqIdx = 0):
        """ Gathers the net sources for each cell into rank 0 as a map with global index as key [kg / cm3 / day]"""
        self.checkGridInitialized()
        return self._map(self._flat0(self.gather(self.getSource_(eqIdx), root = 0)), 2)

    def getSource_(self, eqIdx = 0):
        """nompi version of """
        self.checkGridInitialized()
        unitChange = 24 * 3600 / 1e6  # [ kgOrmol/m3/s -> kgOrmol/cm3/day]
        if eqIdx == 0:
            if self.useMoles:
                unitChange *= self.molarDensityWat  # [mol/cm3/day -> cm3/cm3/day]
            else:
                unitChange *= 1e3  # [kg/cm3/day -> cm3/cm3/day]
        else:
            if not self.useMoles:
                unitChange *= 1e3 # [kg/cm3/day -> g/cm3/day]
            
        return np.array(self.base.getScvSources()[eqIdx]) * unitChange

    def setSource(self, source_map, eq_idx = 0):
        """Sets the source term as map with global cell index as key, and source as value [cm3/day] """
        self.checkGridInitialized()
        for key, value in source_map.items():
            if eq_idx == 0:
                if self.useMoles:
                    source_map[key] = value / 24. / 3600. / self.molarDensityWat;  # [cm3/day] -> [mol/s]
                else:
                    source_map[key] = value / 24. / 3600. / 1.e3;  # [cm3/day] -> [kg/s] (richards.hh)
            else:
                source_map[key] = value / 24. / 3600. / 1.e3;  # [g/day] -> [kg/s] or [mol/day] -> [mol/s]
        self.base.setSource(source_map, eq_idx)

    def applySource(self, dt, source_map, crit_p):
        """Sets the source term as map with global cell index as key, and source as value [cm3/day] """
        self.checkGridInitialized()
        for key, value in source_map.items():
            source_map[key] = value / 24. / 3600. / 1.e3;  # [cm3/day] -> [kg/s]
        self.base.applySource(dt * 24.*3600., source_map, self.to_pa(crit_p))

    def setCriticalPressure(self, critical):
        """ Sets the critical pressure to limit flow for boundary conditions constantFlow, constantFlowCyl, and atmospheric """
        self.base.setCriticalPressure(critical)

    def getSolutionHead(self, eqIdx = 0):
        """Gathers the current solution into rank 0, and converts it into a numpy array (Ndof, neq), 
        model dependent units, [Pa, ...]"""
        self.checkGridInitialized()
        return self._map(self._flat0(self.gather(self.base.getSolutionHead(eqIdx), root = 0)), 0)

    def getSolutionHead_(self, eqIdx = 0):
        """ no mpi version of getSolutionHead() """
        self.checkGridInitialized()
        return np.array(self.base.getSolutionHead(eqIdx))

    def getSolutionHeadAt(self, gIdx, eqIdx = 0):
        """Returns the current solution at a cell index"""
        return self.base.getSolutionHeadAt(gIdx, eqIdx)

    def getWaterContent(self):
        """Gathers the current solution's saturation into rank 0, and converts it into a numpy array (Nc, 1) [1]"""
        self.checkGridInitialized()
        return self._map(self._flat0(self.gather(self.base.getWaterContent(), root = 0)), 2)

    def getWaterContent_(self):
        """no mpi version of getWaterContent() """
        self.checkGridInitialized()
        return np.array(self.base.getWaterContent())

    def getWaterVolume(self):
        """Returns total water volume of the domain [cm3]"""
        self.checkGridInitialized()
        return self.base.getWaterVolume() * 1.e6  # m3 -> cm3

    def getVeclocity1D(self):
        """Returns the Darcy velocities [cm/day] TODO not working! """
        self.checkGridInitialized()
        return np.array(self.base.getVelocity1D()) * 100.*24 * 3600  # m/s -> cm/day

    def setRegularisation(self, pcEps, krEps):
        self.checkGridInitialized()
        """ Van Genuchten regularisation parameters"""
        self.base.setRegularisation(pcEps, krEps)

    def addVanGenuchtenDomain(self, min_, max_, layerIndex):
        """ add Van Genuchten domain """
        self.checkGridInitialized()
        self.base.addVanGenuchtenDomain(min_[0] / 100., min_[1] / 100., min_[2] / 100., max_[0] / 100., max_[1] / 100., max_[2] / 100., layerIndex)

    def changeVanGenuchtenSet(self, vgIndex, qr, qs, alpha, n, ks):
        """ add Van Genuchten domain """
        self.checkGridInitialized()
        self.base.changeVanGenuchtenSet(vgIndex, qr, qs, alpha, n, ks)

    def getLayerCellIndx(self, layerIndex):
        """ get indexes of cells in a specific soil layer (Nc, 1) 
            defined via @see addVanGenuchtenDomain
        """
        self.checkGridInitialized()
        return np.unique(self._flat0(self.gather(self.base.getLayerCellIndx(layerIndex), root = 0)))

    def getConcentration_(self,idComp):      
        """ non mpi version of
        """
        # only valid for solutes (idComp > 0). 
        if (idComp <= 0)|(idComp > self.numComp):
            raise Exception("richards::getConcentration_: Requested element id {idComp} does not correspond to a solute.")
        isDissolved = idComp <= self.numDissolvedSoluteComp # is the component in the water phase (True) [bool]
        C_ = self.getSolution_(idComp) # mol/mol wat or mol/mol scv
        return C_ * self.phaseDensity(isDissolved)
        
    def getConcentration(self,idComp):      
        """ returns the concentration of a component (solute) [mol/cm3 phase] or [g/cm3 phase]
            with phase either the liquid (water) or the solid (bulk soil) phase
            @param idcomp: index of the component (> 0) [int]
        """
        # only valid for solutes (idComp > 0). 
        if (idComp <= 0)|(idComp > self.numComp):
            raise Exception("richards::getConcentration: Requested element id {idComp} does not correspond to a solute.")
        isDissolved = idComp <= self.numDissolvedSoluteComp # is the component in the water phase (True) [bool]
        C_ = self.getSolution(idComp) # [mol/mol phase] or [g/g phase]   
        return C_ * self.phaseDensity(isDissolved)
        
    def getContent(self,idComp):
        """ returns the content of a component (solute) [mol] or [g]
            @param idcomp: index of the component (> 0) [int]
        """
        # only valid for solutes (idComp > 0). 
        if (idComp <= 0)|(idComp > self.numComp):
            raise Exception("richards::getContent: Requested element id {idComp} does not correspond to a solute.")
            
        isDissolved = idComp <= self.numDissolvedSoluteComp # is the component in the water phase (True) [bool]
        vols = self.getCellVolumes_() # cm3 scv            
        
        C_ = self.getSolution_(idComp) # concentration in [mol solute/mol phase] or [g solute/g phase] 
        C_ *= self.phaseDensity(isDissolved) # [mol solute/mol phase] * [mol phase/cm3 phase] or [g solute/g phase] * [g phase/cm3 phase]
            
        if not isDissolved: # is in the soil phase
            return self._map(self._flat0(self.gather(np.multiply(vols , C_  ))),0)
        
        # is in the water phase: needs to be multiplied by the water content to get the total concentration in the subcrontrol volume
        watCont = self.getWaterContent_() # cm3 wat/cm3 scv
        return self._map(self._flat0(self.gather(np.multiply(np.multiply(vols , watCont) , C_ ))),0)
          
        
    def phaseDensity(self, isDissolved):
        """ returns the density of the phase the solute is in [mol / m3]
            i.e., water phase for disolved solutes, solide phase otherwise
            @param isDissolved: is the solute disolved in water [bool]
        """
        if isDissolved: 
            if self.useMoles:
                return self.molarDensityWat #mol wat / cm3 wat
            else:
                return self.densityWat # g wat / cm3 wat
        else:   
            if self.useMoles:
                return self.molarBulkDensity # mol soil minerals / cm3 scv
            else:
                return self.bulkDensity # g soil minerals / cm3 scv
            
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

