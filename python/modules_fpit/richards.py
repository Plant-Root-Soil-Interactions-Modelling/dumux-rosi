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
        super().__init__(base, usemoles)
        self.soils = []
        self.param_group = "Soil."
        self.useMoles = usemoles
        self.mpiVerbose = False
        self.molarMassC = 12.011 # g/mol

    def setParameterGroup(self, group:str):
        """ sets the DuMux paramter group, must end with a dot, e.g. 'Soil.' """
        self.param_group = group
        
    def setComputeDtCSS2(self,DtCSS2):        
        return self.base.setComputeDtCSS2(DtCSS2)  
        
    def computeDtCSS2(self,CSS1, CSW, CSS2):    
        # mol/m3/s * m3/cm3 * s/d => mol/cm3/d 
        return  self.base.computeDtCSS2(CSS1, CSW, CSS2) * self.m3_per_cm3 * (24*3600)    
    
    def computeInitCSS2(self,CSS1, CSW):    
        # mol/m3 * m3/cm3 => mol/cm3 
        return  self.base.computeInitCSS2(CSS1, CSW) * self.m3_per_cm3 

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

    def setHomogeneousIC(self, p:float, equilibrium = False, doReturn = False):
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
        #
        # print("setRootSystemBC", params)
        self.base.setRootSystemBC(params)

    def setSoluteTopBC(self, type_top, value_top):
        """  sets the flux directly in the problem (problem must be initialized), calls base.setSToptBC in richards.hh"""
        self.base.setSTopBC(type_top, value_top)
        
    def setSoluteBotBC(self, type_bot, value_bot):
        """  sets the flux directly in the problem (problem must be initialized), calls base.setSToptBC in richards.hh"""
        self.base.setSBotBC(type_bot, value_bot)

    def getInnerHead(self, shift = 0):
        """Gets the pressure head at the inner boundary [cm] """
        return self.base.getInnerHead(shift)  # -> richards_cyl.hh

    def getInnerSolutes(self, shift = 0, compId = 1, isDissolved = True):
        """Gets the concentration at the inner boundary [mol/cm3] """    
        CC = np.array(self.getSolution_(compId)).flatten()[0] #mol/mol
        if self.useMoles:
            if isDissolved:
                CC *= self.molarDensityWat_m3/1e6 #mol/cm3 scv
            else:
                CC *= self.bulkDensity_m3/1e6 #mol/cm3 scv
        return CC

    def setSource(self, source_map, eq_idx = 0):
        """Sets the source term as map with global cell index as key, and source as value [cm3/day] """
        self.checkInitialized()
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

        # to go from [mol/s] to [mol/m3/s].
        # for dimWorld==3, scv.volume() computed in dumux. otherwise, do it here
        #if self.dimWorld == 1:
        #    vols = self.getCellSurfacesCyl() / 1e4 * cyl_length / 100 #m3 scv
        #else:
        #    vols = np.ones(max(list(source_map.keys())) +1) # volume computed within dumux

        # print("setSourceA", source_map, eq_idx)
        for cellId, value in source_map.items(): 
            source_map[cellId] = value * unitConversion# / vols[cellId]         
        #mol/s
        #print('setsource',source_map, eq_idx )
        self.base.setSource(source_map, eq_idx) # 

    def applySource(self, dt, source_map, crit_p):
        """Sets the source term as map with global cell index as key, and source as value [cm3/day] """
        self.checkInitialized()
        assert not self.useMoles
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
        if (self.mpiVerbose and (size > 1)):
            comm.barrier()
            print("richards::getSolutionHead", rank)
            comm.barrier()
        return (self._map(self.allgatherv(self.base.getSolutionHead(eqIdx)), 0))#.flatten()

    def getSolutionHead_(self, eqIdx = 0):
        """ no mpi version of getSolutionHead() """
        self.checkInitialized()
        #assert max_rank == 1
        return np.array(self.base.getSolutionHead(eqIdx))

    def getSolutionHeadAt(self, gIdx, eqIdx = 0):
        """Returns the current solution at a cell index"""
        return self.base.getSolutionHeadAt(gIdx, eqIdx)

    def getKrw(self):
        """Gathers the current solution's saturation into rank 0, and converts it into a numpy array (Nc, 1) [1]"""
        self.checkInitialized()
        if (self.mpiVerbose and (size > 1)):
            comm.barrier()
            print("richards::getKrw", rank)
            comm.barrier()
        return self._map(self.allgatherv(self.base.getKrw()), 0)
        
        
    def getSaturation(self):
        """Gathers the current solution's saturation into rank 0, and converts it into a numpy array (Nc, 1) [1]"""
        self.checkInitialized()
        if (self.mpiVerbose and (size > 1)):
            comm.barrier()
            print("richards::getSaturation", rank)
            comm.barrier()
        return self._map(self.allgatherv(self.base.getSaturation()), 0)
        
    def getWaterContent(self):
        """Gathers the current solution's saturation into rank 0, and converts it into a numpy array (Nc, 1) [1]"""
        self.checkInitialized()
        if (self.mpiVerbose and (size > 1)):
            comm.barrier()
            print("richards::getWaterContent", rank)
            comm.barrier()
        theta= (self._map(self.allgatherv(self.base.getWaterContent()), 2))#.flatten()
        if len(theta) > 0:
            assert min(theta) >= self.vg_soil.theta_R # better install a systematic heck no? and/or put it at the end of each solve function
            assert max(theta) <= self.vg_soil.theta_S
        return theta

    def getWaterContent_(self):
        """no mpi version of getWaterContent() """
        self.checkInitialized()
        #assert max_rank == 1
        theta = np.array(self.base.getWaterContent())
        if len(theta) > 0:
            assert min(theta) >= self.vg_soil.theta_R
            assert max(theta) <= self.vg_soil.theta_S
        return theta

    def getWaterVolume(self):
        """Returns total water volume of the domain [cm3]"""
        self.checkInitialized()
        assert self.dimWorld != 1
        return self.base.getWaterVolume() * 1.e6  # m3 -> cm3
        
    def getWaterVolumesCyl(self, verbose = False):
        """Returns total water volume of the domain [cm3]"""
        self.checkInitialized()
        assert self.dimWorld != 3
        vols = self.getCellVolumes() #cm3 scv
        watCont = self.getWaterContent()#.flatten() # cm3 wat/cm3 scv
        if(verbose):
            print("getWaterVolumesCyl")
            #print(length)
            print(vols , watCont )
            print(np.multiply(vols , watCont  ) )
        
        return np.multiply(vols , watCont  )  
        
    def getWaterVolumes(self):
        """Returns total water volume of the domain [cm3]"""
        self.checkInitialized()
        vols = self.getCellVolumes()#.flatten() #cm3 scv 
        watCont = self.getWaterContent()#.flatten()  # cm3 wat/cm3 scv
        return np.multiply(vols , watCont  )
        
    def getCSS1_th(self):
        raise Exception
        C_S_W = self.molarDensityWat_m3*np.array(self.getSolution(1))#.flatten()
        vols = self.getCellVolumes()/1e6#m3
        #mol C / cm3 scv
        css1_th = self.CSSmax * (C_S_W/(C_S_W+ self.k_sorp*1e6)) * self.f_sorp # if cell is empty, can get it directly from the solver.
        return  np.array(css1_th*vols)
    
    def getTotCContent_each(self):
        if self.dimWorld == 3:
            vols = self.getCellVolumes()#.flatten() #cm3 scv   
            totC = np.array([self.getContent(i+1, (i < 2)) for i in range(self.numComp +1)])
            
        elif self.dimWorld == 1:
            vols = self.getCellVolumes()#self.getCellSurfacesCyl()  * l  #cm3 scv
            totC = np.array([self.getContentCyl(i+1, (i < 2)) for i in range(self.numComp+1)])
        # mol/mol * (mol/m3) = mol/m3 
        
        # css1 = self.getCSS1_out()  #mol C / cm3 scv
        # print('css1',css1.shape,self.getCellVolumes_().shape)
        # print('rank',rank,totC.shape, vols.shape )
        # totC = np.vstack((totC, css1*vols))
        # raise Exception
        if rank == 0:
            try:
                assert np.array(totC).shape == (self.numComp +1,self.numberOfCellsTot)
            except:
                print('totC',totC,totC.shape , (self.numComp +1, self.numberOfCellsTot))
                raise Exception
            
        return totC
        
    def getTotCContent(self):
        return self.getTotCContent_each().sum(axis=0)
        
    def getCSS1_out_th(self):#mol C / cm3 scv
        raise Exception
        if (self.css1Function == 0) or (self.css1Function == 4) :
            C_  = self.molarDensityWat_m3*np.array(self.getSolution(1))#.flatten()
            return self.CSSmax * (C_ /(C_ + self.k_sorp*1e6)) * self.f_sorp
        elif (self.css1Function == 1) or (self.css1Function == 3):
            return 0.
        elif self.css1Function == 2:
            C_  = self.molarDensityWat_m3*np.array(self.getSolution(1))#.flatten()
            return self.CSSmax * C_ /(self.k_sorp*1e6) * self.f_sorp
        elif self.css1Function == 5:
            vols = self.getCellVolumes()/1e6#m3
            watCont = self.getWaterContent()
            C_  = self.molarDensityWat_m3*np.array(self.getSolution(1)) * vols * watCont
            return self.CSSmax * (C_ /(C_ + self.k_sorp*1e6)) * self.f_sorp
        elif self.css1Function == 6:#cssmax is content
            vols = self.getCellVolumes()/1e6#m3
            watCont = self.getWaterContent()
            C_  = self.molarDensityWat_m3*np.array(self.getSolution(1)) * vols * watCont
            return self.CSSmax *1e6 * (C_ /(C_ + self.k_sorp*1e6)) * self.f_sorp / (vols *1e6)
        elif self.css1Function == 7:
            vols = self.getCellVolumes()/1e6#m3
            watCont = self.getWaterContent()
            C_  = self.molarDensityWat_m3*np.array(self.getSolution(1)) * vols * watCont
            return self.CSSmax * C_ /( self.k_sorp*1e6) * self.f_sorp
        elif self.css1Function == 8:
            vols = self.getCellVolumes()/1e6#m3
            watCont = self.getWaterContent()# m3/m3
            # mol = 
            C_  = self.molarDensityWat_m3*np.array(self.getSolution(1)) * vols * watCont
            return self.CSSmax * C_ /( self.k_sorp)  / (vols *1e6) #* self.f_sorp
        else:
            raise Exception
            
    #def getCSS1_out_real(self):#mol C / cm3 scv
    #    return np.array(self.base.getCSS1_out())/1e6
        
    def getCSS1_out_(self):#mol C / cm3 scv
        temp = np.array(self.base.computeCSS1s())#mol C / m3 scv zone 1
        #print('getCSS1_out_', temp, self.f_sorp,self.CSSmax)
        return self.f_sorp * temp /1e6
        
    def getCSS1_out(self):#mol C / cm3 scv
        temp = self._map(self.allgatherv(self.getCSS1_out_()),0)
        return temp
        #else:
        #    return self.getCSS1_out_real()
        
    def getContentCyl_deprecated(self,idComp, isDissolved,gId = None ):
        if False:
            assert self.dimWorld != 3
            assert idComp > 0 # do not use for water content
            vols = self.getCellVolumes() /1e6#/ 1e4 * length / 100 #m3 scv
            
            if idComp <= self.numComp:
                C_ = self.getSolution(idComp)#.flatten() # mol/mol or g/g 
            elif (idComp == (self.numComp +1)):
                C_ = self.getCSS1_out() *1e6#  mol C / m3 scv
                #print('getCSS1_out_th',self.CSSmax * (C_S_W/(C_S_W+ self.k_sorp*1e6)) * self.f_sorp)
            else:
                print('wrong idComp', idComp)
                raise Exception
            
            
            try:
                assert (C_ >= 0.).all()
            except:
                print('getContentCyl',idComp, isDissolved, vols, C_)
                raise Exception
                
            if not isDissolved:
                if ((self.useMoles) and (idComp != (self.base.numComp()))):
                    C_ *= self.bulkDensity_m3 #mol/m3 scv
                return np.multiply(vols , C_  ) # mol
                
            watCont = self.getWaterContent()#.flatten() # m3 wat/m3 scv
            if self.useMoles:
                C_ *= self.molarDensityWat_m3 # mol/mol wat* mol wat/m3 wat
            #print("np.multiply(vols , watCont)", sum(np.multiply(vols , watCont)))    
            return np.multiply(np.multiply(vols , watCont) , C_ )
        
    def phaseDensity(self, isDissolved):# mol / m3
        if isDissolved: #mol wat / m3 wat
            return self.molarDensityWat_m3
        else:   # mol scv / m3 scv
            return self.bulkDensity_m3
         
    def getFlux_10c(self):
        assert self.dimWorld == 3
        ff10c_ = self.getFlux_10c_()
        f2cidx_ = self.getFace2CellIds_()
        #print('ff10c_',rank, ff10c_)
        #print('f2cidx_',rank, list(f2cidx_))
        #occurences = [list(f2cidx_).count(idx_) for idx_ in set(f2cidx_)]
        #print('occurences',rank,[list(f2cidx_).count(idx_) for idx_ in set(f2cidx_)])
        ff10c = np.array([sum(ff10c_[np.where(f2cidx_ == idx_)[0]]) for idx_ in set(f2cidx_) if list(f2cidx_).count(idx_) == 6]) # sum values per cell
        f2cidx = np.array([idx_ for idx_ in set(f2cidx_) if list(f2cidx_).count(idx_) == 6])# keep those local elem  indexes
        dofind = np.array(self.base.getDofIndices())
        f2cidx_g = dofind[f2cidx] # get global index
        #print('f2cidx',rank, f2cidx,f2cidx_)
        #raise Exception
        #print('dofind',rank,dofind)
        f2cidx_gAll = self.allgatherv(f2cidx_g)
        ff10c_All = self.allgatherv(ff10c)
        if False:
            print(f2cidx_gAll)
            print(ff10c_All)
            print(len(f2cidx_gAll), len(ff10c_All), len(f2cidx_g), self.numberOfCells)
            print(f2cidx_gAll)
        f2cidx_gAll_unique = np.array(list(set(f2cidx_gAll)))
        
        if False:
            for jjj in f2cidx_gAll_unique:
                tempff = ff10c_All[np.where(f2cidx_gAll == jjj)]
                try:
                    assert(tempff == tempff[0]).all()
                except:
                    print('jjj',jjj,tempff )
                    raise Exception
            #print(np.where(f2cidx_gAll == jjj), max(np.where(f2cidx_gAll == jjj)[0]))
        ff10c_All_unique = np.array([ff10c_All[max(np.where(f2cidx_gAll == idx_)[0])] for idx_ in f2cidx_gAll_unique])
        
        #print('ff10c',rank,ff10c)
        #print('f2cidx',rank,f2cidx)
        #setf2cidx__ = np.array(list(set(f2cidx)))
        #print('setf2cidx__',rank, setf2cidx__, setf2cidx__.shape)
        #setf2cidx_ = self.allgatherv(setf2cidx__)
        #setf2cidx = self._map( setf2cidx_, 0)
        #print('setf2cidx_',rank, setf2cidx_)
        #print('setf2cidx',rank, setf2cidx)
        #print('faceGIdxs',self.base.faceGIdxs.keys(),len(self.base.faceGIdxs.keys()))
        #print('facemap',self.base.facemap.keys(),len(self.base.facemap.keys()), self.numberOfFacesTot)
        
        # raise Exception
        # flux10c = self.allgatherv(ff10c)#, keepShape =True) #TODO: check that gathergin works  
        # face2CellIds = self.allgatherv(f2cidx)#, keepShape =True)
        

        # print('getFlux_10c()',rank,flux10c.shape, np.array(self.getCellVolumes_()).shape ,len(inFluxes_ddt), self.numberOfFacesTot,self.base.numComp())
        #flux10c_ = flux10c[0]
        #face2CellIds = face2CellIds.transpose((1,0))
        #flux10c = flux10c.transpose((2,0,1))
        # print('face2CellIds',face2CellIds.shape, flux10c.shape)
        # print([(nf, np.where(valThreads == max(valThreads))[0][0],  max(valThreads)) for nf, valThreads in enumerate(face2CellIds)])
        #for nf, valThreads in enumerate(face2CellIds):
        #    flux10c_[nf][:] = flux10c[np.where(valThreads == max(valThreads))[0]][nf][:]
        #flux10c = np.array([flux10c[np.where(valThreads == max(valThreads))[0][0]][nf][:] for nf, valThreads in enumerate(face2CellIds)])
        #face2CellIds = face2CellIds.max(axis = 1)
        # print('face2CellIds',face2CellIds)
        
        #flux10c = flux10c_
        # print(flux10c.shape,np.array([np.where(face2CellIds == nCell) for nCell in range(self.numberOfCellsTot)]),
        #     np.array([flux10c[np.where(face2CellIds == nCell)] for nCell in range(self.numberOfCellsTot)]))
        #flux10cCell = np.array([flux10c[np.where(face2CellIds == nCell)].sum(axis=0) for nCell in range(self.numberOfCellsTot)])  
        # print(flux10cCell)
        flux10cCell = np.transpose(ff10c_All_unique) # [comp][cell]
        if rank == 0:
            assert flux10cCell.shape == (self.base.numComp() ,self.numberOfCellsTot)
        
        molarMassWat = 18. # [g/mol]
        densityWat = 1. #[g/cm3]
        # [mol/cm3] = [g/cm3] /  [g/mol] 
        molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 
        flux10cCell[0] /=  molarDensityWat # mol to cm3 for water 
        if rank == 0:
            assert flux10cCell.shape == (self.base.numComp() ,self.numberOfCellsTot)#numComp + water 
            
        # print('flux shape',flux10cCell.shape)
        #print('getFlux_10c()',flux10cCell)
        #raise Exception
        return flux10cCell
        
    def getFace2CellIds_(self):
        return np.array(self.base.face2CellIds).max(axis = 0)
    
    def getFlux_10c_(self):# cm3 or mol
        verbose = False
        inFluxes = np.array(self.base.inFluxes) #[ mol / s]
        inFluxes_ddt = np.array(self.base.inFluxes_ddt)# s
        
        
        try:
            if size == 1:# with MPI, could be that axis 1 size < self.numberOfCellsTot
                assert inFluxes.shape == (len(inFluxes_ddt), self.numberOfFacesTot,self.base.numComp())
            else:
                assert inFluxes.shape[0] == len(inFluxes_ddt)
                assert inFluxes.shape[2] == (self.base.numComp())
                
        except:
            print('shape failed, inFluxes',size, inFluxes, 'inFluxes.shape',inFluxes.shape,(len(inFluxes_ddt), self.numberOfFacesTot,self.base.numComp()))
            raise Exception
        
        inFluxes_tot = np.array([ np.array([xxx * inFluxes_ddt[isrcs] for idc, xxx in enumerate(inFluxes[isrcs] )]) for isrcs in range(len(inFluxes_ddt))])
        
        #inFluxes_tot = np.array([inFx * inFluxes_ddt[idx] for idx, inFx in enumerate(inFluxes)]) # cm3 or mol at each dumux sub-time step
        inFluxes_tot = inFluxes_tot.sum(axis = 0) # cm3 or mol, [cell][comp]
        #css1_flux = np.full((self.numberOfCellsTot,1) , 0.)
        
        #inFluxes_tot = np.hstack((inFluxes_tot, css1_flux))
        #print('inFluxes_tot',inFluxes_tot[:],inFluxes_tot.shape)
        #print('getFlux_10c_()',rank,inFluxes_tot.shape, np.array(self.getCellVolumes_()).shape ,len(inFluxes_ddt), self.numberOfFacesTot,self.base.numComp())
        return inFluxes_tot
        
    def getFlux_10c_Old(self):# cm3 or mol
        verbose = False
        inFluxes = np.array(self.base.inFluxes) #[ mol / s]
        inFluxes_ddt = np.array(self.base.inFluxes_ddt)# s
        
        
        try:
            if size == 1:# with MPI, could be that axis 1 size < self.numberOfCellsTot
                assert inFluxes.shape == (len(inFluxes_ddt), self.numberOfCellsTot*6,self.base.numComp())
            else:
                assert inFluxes.shape[0] == len(inFluxes_ddt)
                assert inFluxes.shape[2] == (self.base.numComp())
                
        except:
            print('shape failed, inFluxes',size, inFluxes, inFluxes.shape,(len(inFluxes_ddt), self.numberOfCellsTot,self.base.numComp()))
            raise Exception
        
        inFluxes_tot = np.array([ np.array([xxx * inFluxes_ddt[isrcs] for idc, xxx in enumerate(inFluxes[isrcs] )]) for isrcs in range(len(inFluxes_ddt))])
        
        #inFluxes_tot = np.array([inFx * inFluxes_ddt[idx] for idx, inFx in enumerate(inFluxes)]) # cm3 or mol at each dumux sub-time step
        inFluxes_tot = inFluxes_tot.sum(axis = 0) # cm3 or mol, [cell][comp]
        #css1_flux = np.full((self.numberOfCellsTot,1) , 0.)
        
        #inFluxes_tot = np.hstack((inFluxes_tot, css1_flux))
        print('inFluxes_tot',inFluxes_tot[:],inFluxes_tot.shape)
        print('self.base.getDofIndices()',self.base.getDofIndices()[:])
        return inFluxes_tot
        
    def getSource_10c(self):
        src10c =  comm.bcast(self._map(self.allgatherv(self.getSource_10c_()), 0), root = 0)#css1_before, css1_after)), 0) #TODO: check that gathergin works  
        vols = comm.bcast(self.getCellVolumes(), root = 0)/1e6
        
        #print('src10cA',np.array([src10_ * vols[cellidx] for cellidx, src10_ in enumerate(src10c)]))
        src10c = np.array([src10_ * vols[cellidx] for cellidx, src10_ in enumerate(src10c)])
        src10c = np.transpose(src10c) # [comp][cell]
        #print('src10c',src10c)
        
        molarMassWat = 18. # [g/mol]
        densityWat = 1. #[g/cm3]
        # [mol/cm3] = [g/cm3] /  [g/mol] 
        molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 
        src10c[0] /=  molarDensityWat # mol to cm3 for water 
        if rank == 0:
            assert src10c.shape == (self.base.numComp() ,self.numberOfCellsTot)#numComp + water 
        return src10c
    
    def getSource_10c_(self):#, css1_before = None, css1_after = None):# cm3 or mol
        verbose = False
        inSources = np.array(self.base.inSources) #[ mol / (m^3 \cdot s)] 
        inFluxes_ddt = np.array(self.base.inFluxes_ddt)# s
        #if self.dimWorld == 1:
        #    vols = self.getCellVolumes()/1e6 # / 1e4 * length / 100 #m3 scv
        #elif self.dimWorld == 3:
        #    vols = (1/  1e6)*self.getCellVolumes_() #m3 scv
        #else:
        #    raise Exception
        if False:    
            try:
                assert inSources.shape == (len(inFluxes_ddt), self.numberOfCellsTot,self.base.numComp())
            except:
                print('shape failed, inSources',inSources, inSources.shape,(len(inFluxes_ddt), self.numberOfCellsTot,self.base.numComp()))
                raise Exception
          
        # inSources_tot = np.array([ np.array([xxx * vols[idc] * inFluxes_ddt[isrcs] for idc, xxx in enumerate(inSources[isrcs] )]) for isrcs in range(len(inFluxes_ddt))])
        inSources_tot = np.array([ np.array([xxx * inFluxes_ddt[isrcs] for idc, xxx in enumerate(inSources[isrcs] )]) for isrcs in range(len(inFluxes_ddt))])
        
        #inSources_tot = np.array([inFx * inFluxes_ddt[idx] * vols for idx, inFx in enumerate(inSources)]) # cm3 or mol at each dumux sub-time step
        inSources_tot = inSources_tot.sum(axis = 0) # cm3 or mol, [cell][comp]
        #d_css1 = np.full( (self.numberOfCellsTot,1),0.)
        #if css1_before is not None:
        #    d_css1 = (css1_after - css1_before).reshape(self.numberOfCellsTot,1)
            
        #inSources_tot = np.hstack((inSources_tot,d_css1))
        
        return inSources_tot
            
    def getConcentration(self,idComp, isDissolved):
        C_ = self.getSolution(idComp)#.flatten()  # mol/mol wat or mol/mol scv
        if not isDissolved:
            if self.useMoles:
                C_ *= self.bulkDensity_m3 #mol/m3 scv
            return C_  /1e6  #mol/cm3 scv
        if self.useMoles:
            C_ *= self.molarDensityWat_m3 # mol/m3 wat            
        return C_ /1e6 #mol/cm3 scv
        
    def getContent(self,idComp, isDissolved):
        #assert self.dimWorld != 1
        assert idComp > 0 # do not use for water content
        vols = (1/  1e6)*self.getCellVolumes_()#.flatten() #m3 scv            
        
        if idComp <= self.numComp:
            C_ = self.getSolution_(idComp)#.flatten() # mol/mol or g/g 
        elif (idComp == (self.numComp +1)):
            C_ = self.getCSS1_out_()*1e6 # mol/m3
        else:
            print('wrong idComp', idComp)
            raise Exception
            
        try:
            assert (C_ >= 0.).all()
        except:
            print('getContent',idComp, isDissolved,min( vols),max( vols),min( C_),max( C_))
            raise Exception
            
        if not isDissolved:
            if ((self.useMoles) and (idComp != (self.numComp +1))):
                C_ *= self.bulkDensity_m3 #mol/m3 scv
            return self._map(self.allgatherv(np.multiply(vols , C_  ) ),0)
        watCont = self.getWaterContent_()#.flatten() # m3 wat/m3 scv
        if self.useMoles:
            C_ *= self.molarDensityWat_m3 # mol/m3 wat
            
        return self._map(self.allgatherv(np.multiply(np.multiply(vols , watCont) , C_ )),0)
          

    def setSolution_(self, sol, eqIdx = 1):
        """nompi version of  """
        self.checkInitialized()
        assert max_rank == 1
        return np.array(self.base.setSolution(sol, eqIdx))
        
    def getVeclocity1D(self):
        """Returns the Darcy velocities [cm/day] TODO not working! """
        self.checkInitialized()
        assert not self.useMoles
        return np.array(self.base.getVelocity1D()) * 100.*24 * 3600  # m/s -> cm/day

    def setRegularisation(self, pcEps, krEps):
        self.checkInitialized()
        """ Van Genuchten regularisation parameters"""
        self.base.setRegularisation(pcEps, krEps)

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

