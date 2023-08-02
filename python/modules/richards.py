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
        return self.base.getInnerFlux(eq_idx) * 24 * 3600 * 10.  # [kg m-2 s-1] = g / cm2 / day

    def setOuterFluxCyl(self, flux):
        """ 
        sets the flux directly in the problem (problem must be initialized), calls base.setToptBC, @see setToptBC in richards.hh
        @param flux      [cm/day] negative means outflow    
        """
        self.base.setTopBC(3, flux)

    def getOuterFlux(self, eq_idx = 0):
        """ [cm / day]"""
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
        #
        # print("setRootSystemBC", params)
        self.base.setRootSystemBC(params)

    def setSoluteTopBC(self, type_top, value_top):
        """  sets the flux directly in the problem (problem must be initialized), calls base.setSToptBC in richards.hh"""
        self.base.setSTopBC(type_top, value_top)

    def getInnerHead(self, shift = 0):
        """Gets the pressure head at the inner boundary [cm] """
        return self.base.getInnerHead(shift)  # -> richards_cyl.hh

    def getInnerSolutes(self, shift = 0):
        """Gets the concentration at the inner boundary [cm] """
        return self.base.getInnerSolutes(shift)

    def setSource(self, source_map, eq_idx = 0):
        """Sets the source term as map with global cell index as key, and source as value [cm3/day] """
        self.checkInitialized()
        for key, value in source_map.items():
            source_map[key] = value / 24. / 3600. / 1.e3;  # [cm3/day] -> [kg/s] (richards.hh)
        self.base.setSource(source_map, eq_idx)

    def applySource(self, dt, source_map, crit_p):
        """Sets the source term as map with global cell index as key, and source as value [cm3/day] """
        self.checkInitialized()
        for key, value in source_map.items():
            source_map[key] = value / 24. / 3600. / 1.e3;  # [cm3/day] -> [kg/s]
        self.base.applySource(dt * 24.*3600., source_map, self.to_pa(crit_p))

    def setCriticalPressure(self, critical):
        """ Sets the critical pressure to limit flow for boundary conditions constantFlow, constantFlowCyl, and atmospheric """
        self.base.setCriticalPressure(critical)

    def getSolutionHead(self, eqIdx = 0):
        """Gathers the current solution into rank 0, and converts it into a numpy array (Ndof, neq), 
        model dependent units, [Pa, ...]"""
        self.checkInitialized()
        return self._map(self._flat0(comm.gather(self.base.getSolutionHead(eqIdx), root = 0)), 0)

    def getSolutionHead_(self, eqIdx = 0):
        """ no mpi version of getSolutionHead() """
        self.checkInitialized()
        return np.array(self.base.getSolutionHead(eqIdx))

    def getSolutionHeadAt(self, gIdx, eqIdx = 0):
        """Returns the current solution at a cell index"""
        return self.base.getSolutionHeadAt(gIdx, eqIdx)

    def getWaterContent(self):
        """Gathers the current solution's saturation into rank 0, and converts it into a numpy array (Nc, 1) [1]"""
        self.checkInitialized()
        return self._map(self._flat0(comm.gather(self.base.getWaterContent(), root = 0)), 2)

    def getWaterContent_(self):
        """no mpi version of getWaterContent() """
        self.checkInitialized()
        return np.array(self.base.getWaterContent())

    def getWaterVolume(self):
        """Returns total water volume of the domain [cm3]"""
        self.checkInitialized()
        return self.base.getWaterVolume() * 1.e6  # m3 -> cm3

    def getVeclocity1D(self):
        """Returns the Darcy velocities [cm/day] TODO not working! """
        self.checkInitialized()
        return np.array(self.base.getVelocity1D()) * 100.*24 * 3600  # m/s -> cm/day

    def setRegularisation(self, pcEps, krEps):
        self.checkInitialized()
        """ Van Genuchten regularisation parameters"""
        self.base.setRegularisation(pcEps, krEps)

    def addVanGenuchtenDomain(self, min_, max_, layerIndex):
        """ add Van Genuchten domain """
        self.checkInitialized()
        self.base.addVanGenuchtenDomain(min_[0] / 100., min_[1] / 100., min_[2] / 100., max_[0] / 100., max_[1] / 100., max_[2] / 100., layerIndex)

    def changeVanGenuchtenSet(self, vgIndex, qr, qs, alpha, n, ks):
        """ add Van Genuchten domain """
        self.checkInitialized()
        self.base.changeVanGenuchtenSet(vgIndex, qr, qs, alpha, n, ks)

    def writeDumuxVTK(self, file_name):
        """Uses the Dumux VTK writer to write the current model output"""
        self.checkInitialized()
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

