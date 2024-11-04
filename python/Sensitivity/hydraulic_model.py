""" 
             TODO TODO TODO
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); size = comm.Get_size()
import numpy as np
import timeit
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from rosi_richardsnc import RichardsNCSP  # C++ part (Dumux binding), macroscopic soil model
from rosi_richards import RichardsSP
# from rosi_richards import RichardsSPnum as  RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model

import plantbox as pb  # CPlantBox
import functional.van_genuchten as vg
from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.PlantHydraulicModel import HydraulicModel_Doussan
from functional.PlantHydraulicModel import HydraulicModel_Meunier
import evapotranspiration as evap
from datetime import *


def set_all_sd(rs, s):
    """ # sets all standard deviation to a percantage, i.e. value*s """
    for p in rs.getOrganRandomParameter(pb.OrganTypes.root):
        p.a_s = p.a * s
        p.lbs = p.lb * s
        p.las = p.la * s
        p.lns = p.ln * s
        p.lmaxs = p.lmax * s
        p.rs = p.r * s
        p.thetas = p.theta * s
        p.rlts = p.rlt * s  # no used
        p.ldelays = p.ldelay * s
    seed = rs.getOrganRandomParameter(pb.OrganTypes.seed)[0]  # SeedRandomParameter
    seed.firstBs = seed.firstB * s
    seed.delayBs = seed.delayB * s
    seed.maxBs = seed.maxB * s
    seed.firstSBs = seed.firstSB * s
    seed.delaySBs = seed.delaySB * s
    seed.delayRCs = seed.delayRC * s
    seed.nCs = seed.nCs * s
    seed.nzs = seed.nzs * s
    # todo seed position s


def create_mapped_rootsystem(min_b , max_b , cell_number, soil_model, fname, stochastic = False, mods = None, model = "Meunier"):
    """ loads a rmsl file, or creates a rootsystem opening an xml parameter set,  
        and maps it to the soil_model """

    global picker  # make sure it is not garbage collected away...

    params = PlantHydraulicParameters()

    if fname.endswith(".rsml"):  # opens root geometry from RSML file
        if model == "Meunier":
            r = HydraulicModel_Meunier(fname, params)
        elif model == "Doussan":
            r = HydraulicModel_Doussan(fname, params)
        else:
            raise "create_mapped_rootsystem(): unknown model"
    elif fname.endswith(".xml"):  # functional plant model with parameters given in xml file
        if rank == 0:
            if stochastic:
                seed = np.random.randint(0, 1e6)
            else:
                seed = 1  # always the same random seed
        else:
            seed = None
        seed = comm.bcast(seed, root = 0)  # random seed must be the same for each process; TODO check test for that ...
        # print("create_mapped_rootsystem(): Seed rank {:g}:".format(rank), seed)

        rs = pb.MappedRootSystem()  ####################################################################################
        rs.setSeed(seed)
        rs.readParameters(fname)

        srp = rs.getOrganRandomParameter(pb.OrganTypes.seed)
        src = srp[0].maxB  # 14

        if not stochastic:
            set_all_sd(rs, 0.)

        if mods is not None:  # apply modifications
            rrp = rs.getOrganRandomParameter(pb.OrganTypes.root)
            srp = rs.getOrganRandomParameter(pb.OrganTypes.seed)
            if "lmax" in mods:  # all types
                for i in range(0, len(rrp)):
                    rrp[i].lmax *= mods["lmax"]
                mods.pop("lmax")
            if "lmax145" in mods:
                rrp[1].lmax *= mods["lmax145"]
                rrp[4].lmax *= mods["lmax145"]
                if len(rrp) > 5:
                    rrp[5].lmax *= mods["lmax145"]
                mods.pop("lmax145")
            if "lmax2" in mods:
                rrp[2].lmax *= mods["lmax2"]
                mods.pop("lmax2")
            if "lmax3" in mods:
                rrp[3].lmax *= mods["lmax3"]
                mods.pop("lmax3")
            if "lmax4" in mods:
                rrp[4].lmax *= mods["lmax4"]
                mods.pop("lmax4")
            if "theta45" in mods:
                if len(rrp) > 5:
                    print("shootbore (theta45)")
                    rrp[5].theta = mods["theta45"]
                else:
                    print("seminal (theta45)")
                    rrp[4].theta = mods["theta45"]
                mods.pop("theta45")
            if "r145" in mods:
                rrp[1].r *= mods["r145"]
                rrp[4].r *= mods["r145"]
                if len(rrp) > 5:
                    rrp[5].r *= mods["r145"]
                mods.pop("r145")
            if "r2" in mods:
                rrp[2].r *= mods["r2"]
                mods.pop("r2")
            if "r3" in mods:
                rrp[3].r *= mods["r3"]
                mods.pop("r3")
            if "r" in mods:  # all types
                for i in range(0, len(rrp)):
                    rrp[i].r *= mods["r"]
                mods.pop("r")
            if "ln" in mods:  # all types
                for i in range(0, len(rrp)):
                    rrp[i].ln *= mods["ln"]
                mods.pop("ln")
            if "ln1" in mods:
                rrp[1].ln *= mods["ln1"]
                mods.pop("ln1")
            if "a" in mods:  # all types
                for i in range(0, len(rrp)):
                    rrp[i].a *= mods["a"]
                mods.pop("a")
            if "src" in mods:
                srp[0].maxB = mods["src"]
                print("maxB changed to", mods["src"])
                print(srp[0].firstB)
                print(srp[0].delayB)
                mods.pop("src")
            if "delaySB" in mods:
                srp[0].delaySB = mods["delaySB"]
                mods.pop("delaySB")
            if mods:  # something unused in mods
                print("\nscenario_setup.create_mapped_rootsystem() WARNING mods have unused parameters:")
                for k, v in mods.items():
                    print("key:", k)
                print()
                raise
            # print("theta45", rrp[1].theta, rrp[4].theta)  # 1.5707999
            # print("src", srp[0].maxB)  # 14

        rs.setGeometry(pb.SDF_PlantBox(1.e6, 1.e6, np.abs(min_b[2])))
        rs.initializeLB(4, 5)
        rs.simulate(1., True)

        if model == "Meunier":
            r = HydraulicModel_Meunier(rs, params)
        elif model == "Doussan":
            r = HydraulicModel_Doussan(rs, params)
        else:
            raise "create_mapped_rootsystem(): unknown model"

    r.ms.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), cut = False)

    picker = lambda x, y, z: soil_model.pick([x, y, z])
    r.ms.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions

    return r, params
