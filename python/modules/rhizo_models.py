# import sys; sys.path.append("../modules/"); sys.path.append("../modules/fv/"); sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src");
# sys.path.append("../../build-cmake/cpp/python_binding/")

import plantbox as pb
import functional.xylem_flux as xylem_flux

from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding) of cylindrcial model
from rosi_richardsnc_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial model (a single cylindrical model is not allowed to run in parallel)
from fv.fv_grid import *
import fv.fv_richards as rich  # local pure Python cylindrical models
import functional.van_genuchten as vg

import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()


class RhizoMappedSegments(pb.MappedSegments):
    """
        Adds 1-dimensional rhizosphere models to each root segment of a MappedSegments (or later MappedPlant)    
        
        modes:        
        "dumux"              inner boundary is FluxCyl,                    <- currently best
        "dumux_exact"        inner boundary is rootSystemExact (same as "dumux", but slower, for experiments...) 
        "python"             inner boundary is flux_in,
        "python_exact"       inner boundary is rootsystem_exact (slower, for experiments...)
    """

    # TODO copy mapped segments constructors (!)...

    def __init__(self, x, wilting_point, NC, logbase, mode):
        """ @param file_name is either a pb.MappedRootSystem, pb.MappedSegments, or a string containing a rsml filename"""
        if isinstance(x, str):
            ms = xylem_flux.XylemFluxPython.read_rsml(x)
        elif isinstance(x, xylem_flux.XylemFluxPython):
            ms = x.rs
            self.set_xylem_flux(x)
        elif isinstance(x, pb.MappedSegments):  # should also be true for MappedRootSystem (since MappedSegments is base class)
            ms = x
        super().__init__(ms.nodes, ms.nodeCTs, ms.segments, ms.radii, ms.subTypes)

        # TODO replace by a copy constructor or copy() at some point
        self.seg2cell = ms.seg2cell
        self.cell2seg = ms.cell2seg
        self.soil_index = ms.soil_index
        self.minBound = ms.minBound
        self.maxBound = ms.maxBound
        self.resolution = ms.resolution

        # self.cutAtGrid = ms.cutAtGrid # TODO export in pybind11

        # additional variables
        self.cyls = []
        self.wilting_point = wilting_point
        self.NC = NC
        self.logbase = logbase
        self.mode = mode  # more precise RichardsCylFoam, mode="dumux"
        self.last_dt = 0.
        self.outer_radii = None

    def initialize(self, soil, x, eidx = None, cc = None):
        """ calls the specific initializer according to @param mode (no parallelisation)  
        @param soil     van genuchten parameters as list        
        @param x        is the solution (or initial condition) of the soil model
        @þaram eidx     in case of mpi the cylinder indices per process
        """
        if eidx is None:
            eidx = np.array(range(0, len(self.radii)), np.int64)  # segment indices for process
        self.eidx = np.array(eidx, np.int64)
        self.outer_radii = np.array(self.segOuterRadii())  # not parallel yet
        # self.outer_radii = np.minimum(np.array(self.segOuterRadii()), 1.)  # not parallel yet
        self.seg_length = self.segLength()
        self.soil = soil
        self.vg_soil = vg.Parameters(soil)
        vg.create_mfp_lookup(self.vg_soil, -1.e5, 1000)
        self.cyls = []
        self.dx2 = []
        if self.mode == "dumux":
            for i in eidx:
                self.cyls.append(self.initialize_dumux_(i, x[self.seg2cell[i]], False, False))
        elif self.mode == "dumux_exact":
            for i in eidx:
                self.cyls.append(self.initialize_dumux_(i, x[self.seg2cell[i]], True, False))
        elif self.mode == "dumux_dirichlet":
            for i in eidx:
                self.cyls.append(self.initialize_dumux_(i, x[self.seg2cell[i]], False, True))
        elif self.mode == "dumux_dirichlet_nc":
            for i in eidx:
                self.cyls.append(self.initialize_dumux_nc_(i, x[self.seg2cell[i]], cc[self.seg2cell[i]]))
        elif self.mode == "dumux_nc":
            for i in eidx:
                self.cyls.append(self.initialize_dumux_(i, x[self.seg2cell[i]], False, True))
        elif self.mode == "python" or self.mode == "python_exact":
            for i in eidx:
                x0 = x[self.seg2cell[i]]
                self.cyls.append(self.initialize_python_(i, x0))
        else:
            raise Exception("RhizoMappedSegments.initialize: unknown solver {}".format(self.mode))

    def initialize_dumux_nc_(self, i, x, c):
        a_in = self.radii[i]
        a_out = self.outer_radii[i]
        if a_in < a_out:
            cyl = RichardsNoMPIWrapper(RichardsNCCylFoam())  # only works for RichardsCylFoam compiled without MPI
            cyl.initialize()
            cyl.setVGParameters([self.soil])
            lb = self.logbase
            points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), self.NC, base = lb)
            cyl.createGrid1d(points)
            self.dx2.append(0.5 * (points[1] - points[0]))
            cyl.setHomogeneousIC(x)  # cm pressure head
            cyl.setICZ_solute(c)  # [kg/m2]
            cyl.setInnerBC("pressure", 0.)  # [cm/day]
            cyl.setInnerBC_solute("constantConcentration", c[0])  # [kg/m2]
            cyl.setOuterBC("fluxCyl", 0.)
            cyl.setOuterBC_solute("fluxCyl", 0.)

            cyl.setParameter("Component.MolarMass", "6.2e-2")  # TODO no idea, where this is neeeded, i don't want to use moles ever (nitrate 62,0049 g/mol)
            cyl.setParameter("Component.LiquidDiffusionCoefficient", "1.7e-9")  # m2 s-1 # nitrate = 1700 um^2/sec
            # # cyl.setParameter("Component.BufferPower", "140")  # buffer power = \rho * Kd [1]

            cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
            cyl.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
            cyl.initializeProblem()
            cyl.setCriticalPressure(self.wilting_point)  # cm pressure head
            return cyl
        else:
            print("RhizoMappedSegments.initialize_dumux_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))
            return []

    def initialize_dumux_(self, i, x, exact, dirichlet = False):
        """ Dumux RichardsCylFoam solver"""
        a_in = self.radii[i]
        a_out = self.outer_radii[i]
        if a_in < a_out:
            cyl = RichardsNoMPIWrapper(RichardsCylFoam())  # only works for RichardsCylFoam compiled without MPI
            cyl.initialize()
            lb = self.logbase
            points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), self.NC, base = lb)
            cyl.createGrid1d(points)
            self.dx2.append(0.5 * (points[1] - points[0]))
            cyl.setHomogeneousIC(x)  # cm pressure head
            cyl.setVGParameters([self.soil])
            cyl.setOuterBC("fluxCyl", 0.)  # [cm/day]
            if exact:
                cyl.setInnerBC("rootSystemExact")  # parameters are passed with cyl.setRootSystemBC (defined in richards_cyl.hh)
            else:
                if dirichlet:
                    cyl.setInnerBC("pressure", 0.)  # [cm/day]
                else:
                    cyl.setInnerBC("fluxCyl", 0.)  # [cm/day]
            cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
            cyl.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
            cyl.initializeProblem()
            cyl.setCriticalPressure(self.wilting_point)  # cm pressure head
            return cyl
        else:
            print("RhizoMappedSegments.initialize_dumux_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))
            return []

    def initialize_python_(self, i, x):
        """ Python home grown richards fv solver"""
        a_in = self.radii[i]
        a_out = self.outer_radii[i]
        if a_in < a_out:
            lb = self.logbase
            points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), self.NC, base = lb)
            grid = FVGrid1Dcyl(points)
            cyl = rich.FVRichards1D(grid, self.soil)
            cyl.x0 = np.ones((self.NC - 1,)) * x
            return cyl
        else:
            print("RhizoMappedSegments.initialize_python_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))
            return []

    def set_xylem_flux(self, rs):
        """ we need XylemFlux for the kr and kx call back functions (for python)"""
        self.rs = rs

    def get_inner_heads(self, shift = 0):
        """ matric potential at the root surface interface [cm]"""
        rsx = np.zeros((len(self.cyls),))
        if self.mode.startswith("dumux"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                # print("inner head", cyl.getInnerHead())
                # print("innerIdx", cyl.base.innerIdx)
                rsx[i] = cyl.getInnerHead(shift)  # [cm] (in richards.py, then richards_cyl.hh)
        elif self.mode.startswith("python"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                rsx[i] = cyl.get_inner_head()  # [cm]
        else:
            print("RhizoMappedSegments.get_inner_heads: Warning, mode {:s} unknown".format(self.mode))
        return self._map(self._flat0(comm.gather(rsx, root = 0)))  # gathers and maps correctly

    def get_inner_solutes(self, shift = 0):
        """ matric potential at the root surface interface [cm]"""
        rsx = np.zeros((len(self.cyls),))
        for i, cyl in enumerate(self.cyls):  # run cylindrical models
            rsx[i] = cyl.getInnerSolutes(shift)  # [cm]
        return self._map(self._flat0(comm.gather(rsx, root = 0)))  # gathers and maps correctly

    def get_soil_k(self, rx):
        """ TODO """
        soil_k = np.zeros((len(self.cyls),))
        if self.mode.startswith("dumux"):
            for i, idx in enumerate(self.eidx):  # run cylindrical models
                rsx = self.cyls[i].getInnerHead()
                nidx = idx + 1  # segment index+1 = node index
                try:
                    soil_k[i] = ((vg.fast_mfp[self.vg_soil](rsx) - vg.fast_mfp[self.vg_soil](rx[nidx])) / (rsx - rx[nidx])) / self.dx2[i]
                except:
                    print(rsx, rx[nidx])
        else:
            print("RhizoMappedSegments.get_soil_k: Warning, mode {:s} not implemented or unknown".format(self.mode))
        return self._map(self._flat0(comm.gather(soil_k, root = 0)))

    def get_dx2(self):
        """ TODO doc me AND only for mode="dumux" yet (set in initialize)"""
        return self._map(self._flat0(comm.gather(self.dx2, root = 0)))

    def get_inner_fluxes(self):
        """ fluxes [cm3/day] at the root surface, i.e. inner boundary  """
        fluxes = np.zeros((len(self.cyls),))
        if self.mode.startswith("dumux"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                fluxes[i] = -float(cyl.getInnerFlux()) * (2 * np.pi * self.radii[j] * self.seg_length[j]) / self.radii[j]  # [cm/day] -> [cm3/day], ('/inner_radii' comes from cylindrical implementation)
        elif self.mode.startswith("python"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                fluxes[i] = cyl.get_inner_flux() * (2 * np.pi * self.radii[j] * self.seg_length[j]) / self.last_dt  # divide by dt is correct here! getInnerFlux only gives the source in cm3/cm2
        else:
            print("RhizoMappedSegments.get_inner_fluxes: Warning, mode {:s} unknown".format(self.mode))
        return self._map(self._flat0(comm.gather(fluxes, root = 0)))  # gathers and maps correctly

    def get_inner_concentrations(self):  # TODO
        """ solute concentration at the root surface interface [g / cm3]"""
        rsx = np.zeros((len(self.cyls),))
        if self.mode.startswith("dumux"):
#             for i, cyl in enumerate(self.cyls):  # run cylindrical models
#                 rsx[i] = cyl.getInnerHead()  # [cm]
            pass  # currently zero flux !!!!!!
        else:
            print("RhizoMappedSegments.get_inner_concentrations: Warning, mode {:s} unknown".format(self.mode))
        return self._map(self._flat0(comm.gather(rsx, root = 0)))  # gathers and maps correctly

    def get_inner_mass_fluxes(self):  # TODO check
        """ fluxes [cm3/day] at the root surface, i.e. inner boundary  """
        fluxes = np.zeros((len(self.cyls),))
        if self.mode.startswith("dumux"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                fluxes[i] = -float(cyl.getInnerFlux(1)) * (2 * np.pi * self.radii[j] * self.seg_length[j]) / self.radii[j]  # [cm/day] -> [cm3/day], ('/inner_radii' comes from cylindrical implementation)
        else:
            print("RhizoMappedSegments.get_inner_fluxes: Warning, mode {:s} unknown".format(self.mode))
        return self._map(self._flat0(comm.gather(fluxes, root = 0)))  # gathers and maps correctly

    def get_water_volumes(self):
        """ fluxes [cm3/day] at the root surface, i.e. inner boundary  """
        volumes = np.zeros((len(self.cyls),))
        if self.mode.startswith("dumux"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                volumes[i] = cyl.getWaterVolume()
        else:
            print("RhizoMappedSegments.get_water_volumes: Warning, mode {:s} unknown".format(self.mode))
        return self._map(self._flat0(comm.gather(volumes, root = 0)))  # gathers and maps correctly

    def solve(self, dt, *argv):
        """ set bc for cyl and solves it """
        self.last_dt = dt
        if self.mode == "dumux":
            proposed_inner_fluxes = argv[0]
            proposed_outer_fluxes = argv[1]
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                l = self.seg_length[j]
                cyl.setInnerFluxCyl(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l))  # [cm3/day] -> [cm /day]
                cyl.setOuterFluxCyl(proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l))  # [cm3/day] -> [cm /day]
                try:
                    cyl.solve(dt)
                except:
                    str = "RhizoMappedSegments.solve: dumux exception with boundaries in flow {:g} cm3/day, out flow {:g} cm3/day, segment radii [{:g}-{:g}] cm"
                    str = str.format(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l), proposed_outer_fluxes[j] / (2 * np.pi * self.radii[j] * l), self.radii[j], self.outer_radii[j])
                    print("node ", self.nodes[self.segments[j].y])
                    self.plot_cylinder(j)
                    self.plot_cylinders()
                    raise Exception(str)
        elif self.mode == "dumux_dirichlet":
            rx = argv[0]
            proposed_outer_fluxes = argv[1]
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                l = self.seg_length[j]
                cyl.setInnerMatricPotential(rx[i])
                cyl.setOuterFluxCyl(proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l))  # [cm3/day] -> [cm /day]
                try:
                    cyl.solve(dt)
                except:
                    # str = "RhizoMappedSegments.solve: dumux exception with boundaries in flow {:g} cm3/day, out flow {:g} cm3/day, segment radii [{:g}-{:g}] cm"
                    # str = str.format(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l), proposed_outer_fluxes[j] / (2 * np.pi * self.radii[j] * l), self.radii[j], self.outer_radii[j])
                    # print("node ", self.nodes[self.segments[j].y])
                    self.plot_cylinder(j)
                    self.plot_cylinders()
                    raise Exception(str)
        elif self.mode == "dumux_dirichlet_nc":
            rx = argv[0]
            proposed_outer_fluxes = argv[1]
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                l = self.seg_length[j]
                cyl.setInnerMatricPotential(rx[i])
                cyl.setOuterFluxCyl(proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l))  # [cm3/day] -> [cm /day]
                try:
                    cyl.solve(dt)
                except:
                    # str = "RhizoMappedSegments.solve: dumux exception with boundaries in flow {:g} cm3/day, out flow {:g} cm3/day, segment radii [{:g}-{:g}] cm"
                    # str = str.format(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l), proposed_outer_fluxes[j] / (2 * np.pi * self.radii[j] * l), self.radii[j], self.outer_radii[j])
                    # print("node ", self.nodes[self.segments[j].y])
                    self.plot_cylinder(j)
                    self.plot_cylinders()
                    raise Exception(str)
        elif self.mode == "dumux_exact":
            rx = argv[0]
            proposed_outer_fluxes = argv[1]
            rsx = argv[2]
            soil_k = argv[3]
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                l = self.seg_length[j]
                seg_ct = self.nodeCTs[self.segments[j].y]
                age = 1.  # self.seg_ages[j] + t  # age  = (maxCT -nodeCT[j]) + t = (maxCT+t) - nodeCT[j] = rs_sim_time - nodeCT[j]
                type_ = self.subTypes[j]
                x0 = rx[self.segments[j].x]
                x1 = rx[self.segments[j].y]
                kr = self.rs.kr_f(age, type_)
                if len(soil_k) > 0:
                    kr = min(kr, soil_k[j])
                kx = self.rs.kx_f(age, type_)
                cyl.setRootSystemBC([x0, x1, kr, kx, l, rsx[j]])
                cyl.setOuterFluxCyl(proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l))  # [cm3/day] -> [cm /day]
                # print("solve dumux_exact", x0, x1, self.rs.kr_f(age, type_), self.rs.kx_f(age, type_), l, self.radii[j])
                try:
                    cyl.solve(dt)
                except:
                    str = "RhizoMappedSegments.solve: dumux exception with boundaries in flow out flow {:g} cm3/day, segment radii [{:g}-{:g}] cm"
                    str = str.format(proposed_outer_fluxes[j] / (2 * np.pi * self.radii[j] * l), self.radii[j], self.outer_radii[j])
                    raise Exception(str)
        elif self.mode == "python":
            proposed_inner_fluxes = argv[0]
            proposed_outer_fluxes = argv[1]
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                l = self.seg_length[j]
                proposed_inner_fluxes = argv[0]
                proposed_outer_fluxes = argv[1]
                q_inner = proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l)
                cyl.bc[(0, 0)] = ("flux_in_out", [q_inner , self.wilting_point])
                q_outer = proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l)
                ndof = self.NC - 1
                cyl.bc[(ndof - 1, 1)] = ("flux_in_out", [q_outer , self.wilting_point])
                cyl.solve([dt], dt, False)
                # print("solve python", x0, x1, self.rs.kr_f(age, type_), self.rs.kx_f(age, type_), l, self.radii[j])
        elif self.mode == "python_exact":
            rx = argv[0]
            proposed_outer_fluxes = argv[1]
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                l = self.seg_length[j]
                seg_ct = self.nodeCTs[self.segments[j].y]
                age = 1.  # self.seg_ages[j] + t  # age  = (maxCT -nodeCT[j]) + t = (maxCT+t) - nodeCT[j] = rs_sim_time - nodeCT[j]
                type_ = self.subTypes[j]
                x0 = rx[self.segments[j].x]
                x1 = rx[self.segments[j].y]
                #         rx_approx = 0.5 * (x0 + x1)
                #         cyl.bc[(0, 0)] = ("rootsystem", [rx_approx, self.kr_f(age, type_)])
                cyl.bc[(0, 0)] = ("rootsystem_exact", [x0, x1, self.rs.kr_f(age, type_), self.rs.kx_f(age, type_), self.radii[j], l])
                # print("solve python", x0, x1, self.rs.kr_f(age, type_), self.rs.kx_f(age, type_), l, self.radii[j])
                q_outer = proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l)
                ndof = self.NC - 1
                cyl.bc[(ndof - 1, 1)] = ("flux_in_out", [q_outer , self.wilting_point])
                cyl.solve([dt], dt, False)
        elif self.mode == "dumux_nc":
            proposed_inner_fluxes = argv[0]
            proposed_outer_fluxes = argv[1]
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                l = self.seg_length[j]
                cyl.setInnerFluxCyl(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l))  # [cm3/day] -> [cm /day]
                cyl.setOuterFluxCyl(proposed_outer_fluxes[j] / (2 * np.pi * self.outer_radii[j] * l))  # [cm3/day] -> [cm /day]
                try:
                    cyl.solve(dt)
                except:
                    str = "RhizoMappedSegments.solve: dumux exception with boundaries in flow {:g} cm3/day, out flow {:g} cm3/day, segment radii [{:g}-{:g}] cm"
                    str = str.format(proposed_inner_fluxes[j] / (2 * np.pi * self.radii[j] * l), proposed_outer_fluxes[j] / (2 * np.pi * self.radii[j] * l), self.radii[j], self.outer_radii[j])
                    raise Exception(str)
        else:
            print(self.mode)
            raise Exception("RhizoMappedSegments.initialize: unknown solver {}".format(self.mode))

    def get_water_volume(self):
        """ returns the water volume of the cylindrical models [cm3] """
        volumes = np.zeros((len(self.cyls),))
        if self.mode.startswith("dumux"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                cyl_water = 0.
                cyl_water_content = cyl.getWaterContent()  # getWaterContent() in cpp/pyhton_binding/richards.hh
                nodes = cyl.getPoints()
                for k, wc in enumerate(cyl_water_content):
                    r1 = nodes[k]
                    r2 = nodes[k + 1]
                    cyl_water += np.pi * (r2 * r2 - r1 * r1) * self.seg_length[j] * wc
                volumes[i] = cyl_water
        elif self.mode.startswith("python"):
            for i, cyl in enumerate(self.cyls):  # run cylindrical models
                j = self.eidx[i]  # for one process j == i
                cyl_water = 0.
                cyl_water_content = cyl.get_water_content()  # segment 0
                for k, wc in enumerate(cyl_water_content):
                    r1 = cyl.grid.nodes[k]
                    r2 = cyl.grid.nodes[k + 1]
                    cyl_water += np.pi * (r2 * r2 - r1 * r1) * self.seg_length[j] * wc
                volumes[i] = cyl_water
        else:
            raise Exception("RhizoMappedSegments.get_water_volume: unknown solver {}".format(self.mode))
        return self._map(self._flat0(comm.gather(volumes, root = 0)))  # gathers and maps correctly

    def _map(self, x):
        """Converts @param x to a numpy array and maps it to the right indices                 """
        indices = self._flat0(comm.gather(self.eidx, root = 0))  # segment indices
        if indices:  # only for rank 0 it is not empty
            assert len(indices) == len(x), "RhizoMappedSegments._map: indices and values have different length"
            p = np.zeros((len(x),), dtype = np.float64)
            for i in range(0, len(indices)):  #
                p[indices[i]] = x[i]
            return p
        else:
            return 0

    def _flat0(self, xx):
        """flattens the gathered list in rank 0, empty list for other ranks """
        if rank == 0:
            return [item for sublist in xx for item in sublist]
        else:
            return []

    def plot_cylinder(self, i):
        """ plots a specific cylinder (DUMUX only, TODO) """
        cyl = self.cyls[i]
        x_ = cyl.getDofCoordinates()
        y_ = cyl.getSolutionHead()

        SMALL_SIZE = 22
        MEDIUM_SIZE = 22
        BIGGER_SIZE = 22
        plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
        plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
        plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
        plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
        plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
        plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
        plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title
        # plt.xlim([0., 0.2 ])
        plt.plot(x_, y_)
        plt.xlabel("distance [cm]")
        plt.ylabel("matric potential [cm]")
        plt.show()

    def plot_cylinders(self):
        """ plots a specific cylinder (DUMUX only, TODO) """
        inner, outer = [], []
        zz = -self.minBound.z
        for i, cyl in enumerate(self.cyls):
            x_ = cyl.getDofCoordinates()
            y_ = cyl.getSolutionHead()
            inner.append(y_[0])
            outer.append(y_[-1])
            j = self.segments[i].y
            z = self.nodes[j].z
            col_i = int(-z / zz * 255.)
            c_ = '#%02x%02x%02x' % (col_i, col_i, 64)
            plt.plot(x_, y_, alpha = 0.1, c = c_)
        plt.xlabel("distance [cm], deeper roots are yellow")
        plt.ylabel("matric potential [cm]")
        # plt.xlim([0.05, 0.6])
        # plt.ylim([-8500, 0. ])
        plt.show()
        return  np.argmin(inner), np.argmax(inner), np.argmin(outer), np.argmax(inner)


def plot_transpiration(t, soil_uptake, root_uptake, potential_trans, title = ""):
    """ plots potential and actual transpiration over time 
    
    depending on discretisation soil_uptake and root_uptake might differ  
    
    @param t                  times [day]
    @param soil_uptake        actual transpiration [cm3/day] of soil 
    @param root_uptake        actual transpiration [cm3/day] according to root model 
    @param potential_trans    function in t stating the potential transpiration [cm3/day]
    """
    fig, ax1 = plt.subplots()
    ax1.plot(t, [potential_trans(t_) for t_ in t], 'k', label = "potential transpiration")  # potential transpiration
    ax1.plot(t, -np.array(soil_uptake), 'g', label = "soil uptake")  # actual transpiration  according to soil model
    ax1.plot(t, -np.array(root_uptake), 'r:', label = "root system uptake")  # actual transpiration according root model
    ax1.set_xlabel("Time [d]")
    ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
    ax2 = ax1.twinx()
    dt = np.diff(t)
    so = np.array(soil_uptake)
    cum_transpiration = np.cumsum(-np.multiply(so[:-1], dt))
    ax2.plot(t[1:], cum_transpiration, 'c--')  # cumulative transpiration (neumann)
    ax2.set_ylabel("Cumulative soil uptake $[cm^3]$")
    print("Cumulative soil uptake", cum_transpiration[-1], "[cm^3]")
    fig.legend()
    plt.title(title)
    if title:
        plt.savefig(title + ".png")
    plt.show()


def plot_info(x_, water_collar_cell, water_cyl, collar_sx, min_sx, min_rx, min_rsx, water_uptake, water_domain):
    """ 2x2 plot with additional information """
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    ax1.set_title("Water amount")
    ax1.plot(x_, np.array(water_collar_cell), label = "water cell")
    ax1.plot(x_, np.array(water_cyl), label = "water cylindric")
    ax1.legend()
    ax1.set_xlabel("Time (days)")
    ax1.set_ylabel("(cm3)")
    ax2.set_title("Pressure")
    ax2.plot(x_, np.array(collar_sx), label = "soil at root collar")
    ax2.plot(x_, np.array(min_sx), label = "min soil")
    ax2.plot(x_, np.array(min_rx), label = "min xylem")
    ax2.plot(x_, np.array(min_rsx), label = "min 1d at root surface")
    ax2.set_ylim([-15000, 0])
    ax2.legend()
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Matric potential (cm)")
    ax3.set_title("Water uptake")
    ax3.plot(x_, -np.array(water_uptake))
    ax3.set_xlabel("Time (days)")
    ax3.set_ylabel("Uptake (cm/day)")
    ax4.set_title("Water in domain")
    ax4.plot(x_, np.array(water_domain))
    ax4.set_xlabel("Time (days)")
    ax4.set_ylabel("cm3")
    plt.show()

# print(len(rs.nodeCTs), len(rs.segments))
# ana2 = pb.SegmentAnalyser(r.rs.nodes, r.rs.segments, r.rs.nodeCTs[1:], r.rs.radii)
# types = np.array(r.rs.subTypes, dtype=np.float64)
# ana2.addData("subType", types)
# ana2.addData("age", r.get_ages())
# pd = vp.segs_to_polydata(ana2, 1., ["radius", "subType", "creationTime", "age"])
# vp.plot_roots(pd, "creationTime")

