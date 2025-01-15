""" 
    Helps to set up the root architecture or a single root (and applies possible modifications), and 
    the root hydraulic model (RootHydraulicModel) 
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
    seed.seedPoss.x = seed.seedPos.x * s
    seed.seedPoss.y = seed.seedPos.y * s
    seed.seedPoss.z = seed.seedPos.z * s
    # print(seed.seedPos, seed.seedPoss)
    # print(seed.firstB)
    # print(seed.delayB)


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

        rs = pb.MappedPlant()
        rs.setSeed(seed)
        rs.readParameters(fname)

        if not stochastic:
            set_all_sd(rs, 0.)

        initial_age = 1
        if mods is not None:  # apply modifications
            if "initial_age" in mods:
                initial_age = mods["initial_age"]
                mods.pop("initial_age")
                
            apply_mods(mods, rs)

        rs.setGeometry(pb.SDF_PlantBox(1.e6, 1.e6, np.abs(min_b[2])))
        rs.initializeLB(4, 5)

        rs.simulate(initial_age, True)

        if model == "Meunier":
            r = HydraulicModel_Meunier(rs, params)
        elif model == "Doussan":
            r = HydraulicModel_Doussan(rs, params)
        else:
            raise "create_mapped_rootsystem(): unknown model"

    r.ms.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), cut = False)

    picker = lambda x, y, z: soil_model.pick([0., 0., z])
    r.ms.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions

    return r, params


def apply_mods(mods, plant):
    """
    applies changes to RootRandomParameters @param rrp and SeedRandomParameters @parm srp
    
    maximal root length:             lmax, lmax145, lmax2, lmax3, lmax4          scaling
    insertion angle from base root:  theta45, theta2, theta3                     angles in radians
    initial growth rate:             r, r145, r2, r3                             scaling
    interlateral spacing:            ln, ln145, ln2                              scaling
    root radius                      a, a145, a2, a3                             absolute
    seminal roots:                   src [number of], src_first, src_delay       scaling 
    tropism:                         tropismN, tropismN145, tropismN2, tropismN3 scaling 
                                     tropismS, tropismS145, tropismS2, tropismS3 scaling 
    root hairs:                      hairsZone, hairsLength, hairsElongation     absolute
    """
    rrp = plant.getOrganRandomParameter(pb.OrganTypes.root)
    srp = plant.getOrganRandomParameter(pb.OrganTypes.seed)

    if "lmax" in mods:
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
    if "theta2" in mods:
        rrp[2].theta = mods["theta2"]
    if "theta3" in mods:
        rrp[3].theta = mods["theta2"]

    if "r" in mods:  # all types
        for i in range(0, len(rrp)):
            rrp[i].r *= mods["r"]
        mods.pop("r")
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

    if "ln" in mods:  # all types
        for i in range(0, len(rrp)):
            rrp[i].ln *= mods["ln"]
    if "ln145" in mods:
        rrp[1].ln *= mods["ln145"]
        rrp[4].ln *= mods["ln145"]
        if len(rrp) > 5:
            rrp[5].ln *= mods["ln145"]
        mods.pop("ln")
    if "ln1" in mods:
        rrp[1].ln *= mods["ln1"]
        mods.pop("ln1")
    if "ln2" in mods:
        rrp[2].ln *= mods["ln2"]
        mods.pop("ln2")

    if "a" in mods:  # all types
        for i in range(0, len(rrp)):
            rrp[i].a = mods["a"]
        mods.pop("a")
    if "a145" in mods:
        rrp[1].a = mods["a145"]
        rrp[4].a = mods["a145"]
        if len(rrp) > 5:
            rrp[5].a = mods["a145"]
        mods.pop("a145")
    if "a2" in mods:
        rrp[2].a = mods["a2"]
        mods.pop("a2")
    if "a3" in mods:
        rrp[3].a = mods["a3"]
        mods.pop("a3")

    if "src" in mods:  # seminal root count, called basal roots in cplantbox
        srp[0].maxB = mods["src"]
        mods.pop("src")
    if "src_first" in mods:
        srp[0].firstB *= mods["src_first"]
        mods.pop("src_first")
    if "src_delay" in mods:
        srp[0].delayB *= mods["src_delay"]
        mods.pop("src_delay")

    if "tropismN" in mods:  # all types
        for i in range(0, len(rrp)):
            rrp[i].tropismN *= mods["tropismN"]
        mods.pop("tropismN")
    if "tropismN145" in mods:
        print(rrp[1].tropismN)
        print(type(rrp[1].tropismN))
        print(mods["tropismN145"])
        print(type(mods["tropismN145"]))
        rrp[1].tropismN *= mods["tropismN145"]
        rrp[4].tropismN *= mods["tropismN145"]
        if len(rrp) > 5:
            rrp[5].tropismN *= mods["tropismN145"]
        mods.pop("tropismN145")
    if "tropismN2" in mods:
        rrp[2].tropismN *= mods["tropismN2"]
        mods.pop("tropismN2")
    if "tropismN3" in mods:
        rrp[3].tropismN *= mods["tropismN3"]
        mods.pop("tropismN3")

    if "tropismS" in mods:  # all types
        for i in range(0, len(rrp)):
            rrp[i].tropismS *= mods["tropismS"]
        mods.pop("tropismS")
    if "tropismS145" in mods:
        rrp[1].tropismS *= mods["tropismS145"]
        rrp[4].tropismS *= mods["tropismS145"]
        if len(rrp) > 5:
            rrp[5].tropismS *= mods["tropismS145"]
        mods.pop("tropismS145")
    if "tropismS2" in mods:
        rrp[2].tropismS *= mods["tropismS2"]
        mods.pop("tropismS2")
    if "tropismS3" in mods:
        rrp[3].tropismS *= mods["tropismS3"]
        mods.pop("tropismS3")

    if  "hairsZone" in mods:
        for i in range(0, len(rrp)):
            rrp[i].hairsZone = mods["hairsZone"]
        mods.pop("hairsZone")
    if "hairsZone145" in mods:
        rrp[1].hairsZone = mods["hairsZone145"]
        rrp[4].hairsZone = mods["hairsZone145"]
        if len(rrp) > 5:
            rrp[5].hairsZone = mods["hairsZone145"]
        mods.pop("hairsZone145")
    if "hairsZone2" in mods:
        rrp[2].hairsZone = mods["hairsZone2"]
        mods.pop("hairsZone2")
    if "hairsZone3" in mods:
        rrp[3].hairsZone = mods["hairsZone3"]
        mods.pop("hairsZone3")

    if "hairsLength" in mods:
        for i in range(0, len(rrp)):
            rrp[i].hairsLength = mods["hairsLength"]
        mods.pop("hairsLength")
    if "hairsLength145" in mods:
        rrp[1].hairsLength = mods["hairsLength145"]
        rrp[4].hairsLength = mods["hairsLength145"]
        if len(rrp) > 5:
            rrp[5].hairsLength = mods["hairsLength145"]
        mods.pop("hairsLength145")
    if "hairsLength2" in mods:
        rrp[2].hairsLength = mods["hairsLength2"]
        mods.pop("hairsLength2")
    if "hairsLength3" in mods:
        rrp[3].hairsLength = mods["hairsLength3"]
        mods.pop("hairsLength3")

    if "hairsElongation" in mods:
        for i in range(0, len(rrp)):
            rrp[i].hairsElongation = mods["hairsElongation"]
        mods.pop("hairsElongation")

    if "dx" in mods:
        for i in range(0, len(rrp)):
            rrp[i].dx = mods["dx"]
        mods.pop("dx")

    if mods:  # something unused in mods
        print("\nscenario_setup.create_mapped_rootsystem() WARNING mods have unused parameters:")
        for k, v in mods.items():
            print("key:", k)
        print()
        raise


def create_mapped_singleroot(min_b , max_b , cell_number, soil_model, stochastic = False, mods = None, model = "Meunier"):
    """ creates a mapped sinlge root"""

    global picker  # make sure it is not garbage collected away...

    params = PlantHydraulicParameters()

    rs = create_singleroot(ns = 100, l = 50 , a = 0.05)

    if mods is not None:  # apply modifications
        apply_mods(mods, rs)

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

    picker = lambda x, y, z: soil_model.pick([0., 0., z])
    r.ms.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions

    return r, params


def create_singleroot(ns = 100, l = 50 , a = 0.05):
    """ creates a single root with @param ns segments, length l, and radius a """
    radii = np.array([a] * ns)
    nodes = [pb.Vector3d(0, 0, 0)]
    segs = []
    dx = l / ns
    z_ = np.linspace(-dx, -l , ns)
    for i in range(0, ns):
        nodes.append(pb.Vector3d(0, 0, z_[i]))
        segs.append(pb.Vector2i(i, i + 1))
    rs = pb.MappedSegments(nodes, segs, radii)
    return rs

