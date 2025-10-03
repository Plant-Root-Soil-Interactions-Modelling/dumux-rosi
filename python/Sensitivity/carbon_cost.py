""" 
    Carbon cost
    
    Daniel Leitner, 2025  
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import numpy as np

import plantbox as pb
import scenario_setup


def carbon_(anatomy, a):
    """ carbon per length [g/cm] of model 'anatomical' """
    fresh_density = 1.03  # g / cm3 (1.02–1.05 g/cm3)
    dry_matter_content = 0.08  #

    bd_cortex = 0.3  # carbon_per_drymass
    bd_stele = 0.5
    bd_epidermis = 0.35

    epidermis_area = (a * a - (a - 0.02) * (a - 0.02)) * np.pi  # two cell layers with 0.01 cm diameter
    stele_area = anatomy[4] * anatomy[4] * np.pi
    cortex_area = a * a * np.pi - epidermis_area - stele_area
    cortex_area *= (1. - anatomy[0])  # aerenchyma percentage

    return fresh_density * dry_matter_content * (bd_stele * stele_area + bd_cortex * cortex_area + bd_epidermis * epidermis_area)


def carbon_cost(plant, params, model):
    """ calculates the carbon cost of the plant model = ["""

    fresh_density = 1.03  # g / cm3 (1.02–1.05 g/cm3)
    dry_matter_content = 0.08  #
    carbon_per_drymass = 0.45  # 45% (40-45%)

    if not params["conductivity_mode"] == "from_mecha":  # carbon cost can only calculated with an underlying anatomical model
        model = "volume"

    if model == "volume":
        print(plant.getSummed("volume") , "cm3")
        dry_mass = plant.getSummed("volume") * fresh_density * dry_matter_content

        return dry_mass * carbon_per_drymass

    elif model == "simple":  # only arenchymna percentage is substracted from total volume
        r145 = scenario_setup.get_anatomy(params["conductivity_index1"])
        r2 = scenario_setup.get_anatomy(params["conductivity_index2"])
        r3 = scenario_setup.get_anatomy(params["conductivity_index3"])
        ana = pb.SegmentAnalyser(plant)
        subType = ana.getParameter("subType")
        vol = ana.getParameter("volume")
        v = np.sum(np.multiply((1. - r145[0]) * np.equal(subType, 1), vol))
        v += np.sum(np.multiply((1. - r145[0]) * np.equal(subType, 4), vol))
        v += np.sum(np.multiply((1. - r145[0]) * np.equal(subType, 5), vol))
        v += np.sum(np.multiply((1. - r2[0]) * np.equal(subType, 2), vol))
        v += np.sum(np.multiply((1. - r3[0]) * np.equal(subType, 3), vol))
        dry_mass = v * fresh_density * dry_matter_content

        return dry_mass * carbon_per_drymass

    elif model == "anatomical":  # stele, cortex and epidermis
        scenario_setup.attach_conductivitis(params)
        r145 = scenario_setup.get_anatomy(params["conductivity_index1"])
        r2 = scenario_setup.get_anatomy(params["conductivity_index2"])
        r3 = scenario_setup.get_anatomy(params["conductivity_index3"])
        c145 = carbon_(r145, params["a145_a"])
        c2 = carbon_(r145, params["a2_a"])
        c3 = carbon_(r145, params["a3_a"])
        ana = pb.SegmentAnalyser(plant)
        subType = ana.getParameter("subType")
        len = ana.getParameter("length")
        c = np.sum(np.multiply(c145 * np.equal(subType, 1), len))
        c += np.sum(np.multiply(c145 * np.equal(subType, 4), len))
        c += np.sum(np.multiply(c145 * np.equal(subType, 5), len))
        c += np.sum(np.multiply(c2 * np.equal(subType, 2), len))
        c += np.sum(np.multiply(c3 * np.equal(subType, 3), len))

        return c

    raise "carbon_cost() unknown model " + model
