""" 
    Helps to set up the root architecture or a single root (and applies possible modifications), and 
    the root hydraulic model (RootHydraulicModel) 
    
    TODO check if single root is still working 
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); size = comm.Get_size()
import numpy as np
import timeit
import copy
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
        print(fname)

        if not stochastic:
            set_all_sd(rs, 0.)

        initial_age = 1.
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

    if soil_model is not None:
        picker = lambda x, y, z: soil_model.pick([0., 0., z])
        r.ms.setSoilGrid(picker)  # maps segments, maps root segments and soil grid indices to each other in both directions

    return r, params


def get_indices(key, max_st):
    """ parses names for subTypes, and '_a' for absolute value or '_s' for scaling (default is scaling)"""

    ind_indices = [list(range(0, max_st)), [1, 4], [2, 3], [1], [2], [3], [4]]
    if max_st < 5:
        ind_indices[1] = [1, 4]
    ind_names = ["", "145", "23", "1", "2", "3", "4"]

    abs_ = False  # default
    if key.endswith("_a"):
        abs_ = True
        key = key[0:-2]
    if key.endswith("_s"):
        abs_ = False
        key = key[0:-2]

    for i in range(1, len(ind_indices)):
        if key.endswith(ind_names[i]):
            return ind_indices[i], abs_, key[0:-len(ind_names[i])]

    return ind_indices[0], abs_, key


def apply_mods(mods, plant):
    """
    applies changes to RootRandomParameters @param rrp and SeedRandomParameters @parm srp
    
    maximal root length:             lmax, lmax145, lmax2, lmax3, lmax4              scaling
    insertion angle from base root:  theta45, theta2, theta3                         angles in radians
    initial growth rate:             r, r145, r2, r3                                 scaling
    interlateral spacing:            ln, ln145, ln2                                  scaling
    root radius                      a, a145, a2, a3                                 scaling
    seminal roots:                   src [number of], src_first, src_delay           scaling 
    tropism:                         tropismN, tropismN145, tropismN2, tropismN3     scaling 
                                     tropismS, tropismS145, tropismS2, tropismS3     scaling 
    root hairs:                      hairsZone, hairsLength, hairsElongation         scaling
    aixal resolution                 dx                                              always absolute  
    """
    rrp = plant.getOrganRandomParameter(pb.OrganTypes.root)
    srp = plant.getOrganRandomParameter(pb.OrganTypes.seed)

    # print("rrp")
    # print(len(rrp))
    # print(rrp)
    # for i, r in enumerate(rrp):
    #     print("Index", i)
    #     print("name", r.name)
    mods_ = copy.deepcopy(mods)
    for key in mods_.keys():

        ind_, abs_, key_ = get_indices(key, len(rrp))

        if key_ == "lmax":
            if abs_:  # absolute value
                for i in ind_:
                    rrp[i].lmax = mods[key]
            else:  # scaling
                for i in ind_:
                    rrp[i].lmax *= mods[key]
            mods.pop(key)

        if key_ == "theta":
            if abs_:  # absolute value
                for i in ind_:
                    rrp[i].theta = mods[key]
            else:  # scaling
                for i in ind_:
                    rrp[i].theta *= mods[key]
            mods.pop(key)

        if key_ == "r":
            if abs_:  # absolute value
                for i in ind_:
                    rrp[i].r = mods[key]
            else:  # scaling
                for i in ind_:
                    rrp[i].r *= mods[key]
            mods.pop(key)

        if key_ == "ln":
            if abs_:  # absolute value
                for i in ind_:
                    rrp[i].ln = mods[key]
            else:  # scaling
                for i in ind_:
                    rrp[i].ln *= mods[key]
            mods.pop(key)

        if key_ == "a":
            if abs_:  # absolute value
                for i in ind_:
                    rrp[i].a = mods[key]
            else:  # scaling
                for i in ind_:
                    rrp[i].a *= mods[key]
            mods.pop(key)

        if key_ == "hairsZone":
            if abs_:  # absolute value
                for i in ind_:
                    rrp[i].hairsZone = mods[key]
            else:  # scaling
                for i in ind_:
                    rrp[i].hairsZone *= mods[key]
            mods.pop(key)

        if key_ == "hairsLength":
            if abs_:  # absolute value
                for i in ind_:
                    rrp[i].hairsLength = mods[key]
            else:  # scaling
                for i in ind_:
                    rrp[i].hairsLength *= mods[key]
            mods.pop(key)

        if key_ == "hairsElongation":
            if abs_:  # absolute value
                for i in ind_:
                    rrp[i].hairsElongation = mods[key]
            else:  # scaling
                for i in ind_:
                    rrp[i].hairsElongation *= mods[key]
            mods.pop(key)

        if key_ == "tropismN":
            if abs_:  # absolute value
                for i in ind_:
                    rrp[i].tropismN = mods[key]
            else:  # scaling
                for i in ind_:
                    rrp[i].tropismN *= mods[key]
            mods.pop(key)

        if key_ == "tropismS":
            if abs_:  # absolute value
                for i in ind_:
                    rrp[i].tropismS = mods[key]
            else:  # scaling
                for i in ind_:
                    rrp[i].tropismS *= mods[key]
            mods.pop(key)

        if key_ == "src":  # seminal root count, called basal roots in cplantbox
            if abs_:  # absolute value
                srp[0].maxB = int(mods[key] + 0.5)  # round
            else:  # scaling
                srp[0].maxB *= mods[key]
            mods.pop(key)

        if key_ == "src_first":
            if abs_:  # absolute value
                srp[0].firstB = mods[key]
            else:  # scaling
                srp[0].firstB *= mods[key]
            mods.pop(key)

        if key_ == "src_delay":
            if abs_:  # absolute value
                srp[0].delayB = mods[key]
            else:  # scaling
                srp[0].delayB *= mods[key]
            mods.pop(key)

        if key_ == "dx":
            for i in ind_:
                rrp[i].dx = mods[key]
            mods.pop(key)


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

