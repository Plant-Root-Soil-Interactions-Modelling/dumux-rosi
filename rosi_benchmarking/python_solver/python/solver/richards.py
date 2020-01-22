import os
import sys
sys.path.append("../../../../build-cmake/rosi_benchmarking/python_solver/")

from solver.solverbase import SolverWrapper


class RichardsWrapper(SolverWrapper):
    """ 
    Adds functionality specifically for Richards equation             
    Wraps passing parameters to Dumux for VanGenuchtenParameter, IC, BC, 
    """

    def __init__(self, base):
        super().__init__(base)
        self.param_group = "Soil."

    def  setParameterGroup(self, group :str):
        """ sets the DuMux paramter group, must end with a dot, e.g. 'Soil.' """
        self.param_group = group

    @staticmethod
    def toHead(pa):
        """ Converts Pascal (kg/ (m s^2)) to cm pressure head """
        g , rho, ref = 9.81, 1.e3, 1.e5  # (m/s^2), (kg/m^3), Pa
        return (pa - ref) * 100 / rho / g

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

    def setVGParameters(self, soils :list):
        """ Sets a list of vg parameter lists, 5 values per soil 
            set before creating the problem with SolverBase.initializeProblem 
        
            # todo doc params
        """
        qr, qs, alpha, n, ks = [], [], [], [], []

        for soil in soils:
            assert(len(soil) == 5)
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

    def setLayersZ(self, number :list, z :list = []):
        """ sets depth dependent layers 
        
        @param number     list of layer numbers at the z-positions (if given), or per soil layer, [cm] pressure head.      
        @param z          list of z-positions [cm].  Between the sampling points linear interpolation is applied.                              
        """
        self.setParameter(self.param_group + "Layer.Number", self.dumux_str(number))
        if z:
            assert(len(number) == len(z))  # sample points
            self.setParameter(self.param_group + "Layer.Z", self.dumux_str(z))

    def setICZ(self, p :list, z :list = []):
        """ sets depth dependent initial condtions 
        
        @param p     list of pressures at the z-positions (if given), or per soil layer, [cm] pressure head.      
        @param z     list of z-positions [cm].  Between the sampling points linear interpolation is applied.                              
        """
        self.setParameter(self.param_group + "IC.P", self.dumux_str(p))
        if z:
            assert(len(p) == len(z))  # sample points
            self.setParameter(self.param_group + "IC.Z", self.dumux_str(z))

    def setHomogeneousIC(self, p :float, equilibrium = False):
        """ sets homogeneous initial condions 
        
        @param p              mean matric potential [cm] pressure head
        @param equilibrium    in hydrostatic equilibrium (True) or with a static matric potential (False)
                              for hydrostatic equilibrium the grid must be created before  
        """
        if equilibrium:
            bounds = self.getGridBounds()
            z = [bounds[2], bounds[5]]
            m = 100. * (z[1] - z[0]) / 2.
            p = [p + m, p - m]
            self.setICZ(p, z)
        else:
            self.setICZ(p)

    def setTopBC(self, type_top :str, value_top :float, climate :list = []):
        """ Top boundary conditions are set before creating the problem with SolverBase.initializeProblem 
        
        @param type_top:
        type_top ==  "constantPressure", value_top is a constant pressure [cm] pressure head
        type_top ==  "constantFlux", value_top is the constant flux [cm/day]         
        type_top ==  "atmospheric", value_top is given by climatic data describing evapotranspiration [cm/day], 
                     Data are given in @param climate, the value of value_top is ignored.  
                     Minus denotes evaporation, plus transpiraton.                                            
                     Evaporation stops at a critical pressure of -10000 cm, infiltration is with run off.     

        @param climate:  Two lists are expected, first a list of times [day], second a list of evapotranspirations [cm/day], 
                         between the values linear interpolation is applied.                                                                          
        """
        if type_top == "constantPressure":
            t = 1
        elif type_top == "constantFlux":
            t = 2
        elif type_top == "atmospheric":
            t = 4
        else:
            raise Exception('Top type should be "constantPressure", "constantFlux", or "atmospheric", unknown top type {}'.format(type_top))

        self.setParameter(self.param_group + "BC.Top.Type", str(t))
        self.setParameter(self.param_group + "BC.Top.Value", str(value_top))

        if t == 4:  # atmospheric
            if climate:
                assert(len(climate[0]) == len(climate[1]))  # sample points
                self.setParameter("Climate.Time", self.dumux_str(climate[0]))
                self.setParameter("Climate.Precipitation", self.dumux_str(climate[1]))  # TODO confusing name (should be Evapotranspiration)

            else:
                raise Exception('Atmospheric boundary conditions where set, but no climatic data where given')

    def setBotBC(self, type_bot :str, value_bot = 0.,):
        """ Top boundary conditions are set before creating the problem with SolverBase.initializeProblem 
        
        @param type_bot:
        type_bot ==  "constantPressure", value_bot is a constant pressure [cm] pressure head
        type_bot ==  "constantFlux", value_bot is the constant flux [cm/day] 
        type_bot ==  "freeDrainage", free drainage, the value of value_bot is ignored                       
        """

        if type_bot == "constantPressure":
            b = 1
        elif type_bot == "constantFlux":
            b = 2
        elif type_bot == "freeDrainage":
            b = 5
        else:
            raise Exception('Bottom type should be "constantPressure", "constantFlux", or "freeDrainage", unknown bottom type {}'.format(type_bot))

        self.setParameter(self.param_group + "BC.Bot.Type", str(b))
        self.setParameter(self.param_group + "BC.Bot.Value", str(value_bot))

    def getSaturation(self):
        """Gathers the current solution's saturation into rank 0, and converts it into a numpy array (dof, 1) """
        self.checkInitialized()
        return self._map(self._flat0(MPI.COMM_WORLD.gather(super().getSaturation(), root = 0)))

