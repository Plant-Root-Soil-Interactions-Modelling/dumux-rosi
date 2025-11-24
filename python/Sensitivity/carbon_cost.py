""" 
    Carbon cost models
        
    carbon_cost(plant, params, model) calls one of the models below:
        
    carbon(...)                    proportional to volume
    carbon_cost_simple(...)        carbon cost proportional to the plant volume subtracting arenchymna percentage
    carbon_cost_anatomical(...)    carbon cost as weighted sum with varying carbon dry masses 
        for cortex, stele, and epidermis. Aarenchyma percentage is subtracted from the cortex, see also carbon_()
    
    
    Daniel Leitner, 2025  
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import numpy as np

import plantbox as pb
import scenario_setup


def carbon_cost_volume(vol):
    """ returns the carbon cost propotional to the plant volume """
    fresh_density = 1.03  # g / cm3 (1.02–1.05 g/cm3)
    dry_matter_content = 0.08  #
    carbon_per_drymass = 0.45  # 45% (40-45%)

    vol_ = np.sum(vol)
    dry_mass = vol_ * fresh_density * dry_matter_content
    return dry_mass * carbon_per_drymass


def carbon_cost_simple(vol, subType, ap145, ap2, ap3):
    """ returns the carbon cost proportional to the plant volume subtracting arenchymna percentage """
    fresh_density = 1.03  # g / cm3 (1.02–1.05 g/cm3)
    dry_matter_content = 0.08  #
    carbon_per_drymass = 0.45  # 45% (40-45%)

    v = np.sum(np.multiply((1. - ap145) * np.equal(subType, 1), vol))
    v += np.sum(np.multiply((1. - ap145) * np.equal(subType, 4), vol))
    v += np.sum(np.multiply((1. - ap145) * np.equal(subType, 5), vol))
    v += np.sum(np.multiply((1. - ap2) * np.equal(subType, 2), vol))
    v += np.sum(np.multiply((1. - ap3) * np.equal(subType, 3), vol))
    dry_mass = v * fresh_density * dry_matter_content

    return dry_mass * carbon_per_drymass


def carbon_cost_anatomical(vol, subType, segLen, ana145, ana2, ana3, a145, a2, a3):
    """ returns the carbon cost as weighted sum with varying carbon dry masses for
        cortex, stele, and epidermis. Aarenchyma percentage is subtracted from the cortex
        see also carbon_()
    """
    c145 = carbon_(ana145, a145)
    c2 = carbon_(ana2, a2)
    c3 = carbon_(ana3, a3)
    c = np.sum(np.multiply(c145 * np.equal(subType, 1), segLen))
    c += np.sum(np.multiply(c145 * np.equal(subType, 4), segLen))
    c += np.sum(np.multiply(c145 * np.equal(subType, 5), segLen))
    c += np.sum(np.multiply(c2 * np.equal(subType, 2), segLen))
    c += np.sum(np.multiply(c3 * np.equal(subType, 3), segLen))
    return c


def carbon_(anatomy, a):
    """ carbon per length [g/cm] of model 'anatomical' """
    fresh_density = 1.03  # g / cm3 (1.02–1.05 g/cm3)
    dry_matter_content = 0.08  #

    c_cortex = 0.3  # carbon_per_drymass
    c_stele = 0.5
    c_epidermis = 0.35

    epidermis_area = (a * a - (a - 0.02) * (a - 0.02)) * np.pi  # two cell layers with 0.01 cm diameter
    stele_area = anatomy[4] * anatomy[4] * np.pi
    cortex_area = a * a * np.pi - epidermis_area - stele_area
    cortex_area *= (1. - anatomy[0])  # aerenchyma percentage

    return fresh_density * dry_matter_content * (c_stele * stele_area + c_cortex * cortex_area + c_epidermis * epidermis_area)


def carbon_cost(plant, params, model):
    """ calculates the carbon cost of the plant model = [volume, simple, or anatomical] """

    carbon_per_drymass = 0.45  # 45% (40-45%)

    if not params["conductivity_mode"] == "from_mecha":  # carbon cost can only calculated with an underlying anatomical model
        if not model == "volume":
            print("carbon_cost() warning: switching carbon model to proportional to volume since conductivity_mode is not 'from_mecha'")
        model = "volume"

    ana = pb.SegmentAnalyser(plant)
    vol = ana.getParameter("volume")

    if model == "volume":
        return carbon_cost_volume(vol)

    elif model == "simple":  # only arenchymna percentage is substracted from total volume
        ind = scenario_setup.sort_indices(params)
        ana145 = scenario_setup.get_anatomy(ind[0])
        ana2 = scenario_setup.get_anatomy(ind[1])
        ana3 = scenario_setup.get_anatomy(ind[2])
        subType = ana.getParameter("subType")
        return carbon_cost_simple(vol, subType, ana145[0], ana2[0], ana3[0])

    elif model == "anatomical":  # stele, cortex and epidermis
        scenario_setup.attach_conductivities(params)
        ind = scenario_setup.sort_indices(params)
        ana145 = scenario_setup.get_anatomy(ind[0])
        ana2 = scenario_setup.get_anatomy(ind[1])
        ana3 = scenario_setup.get_anatomy(ind[2])
        subType = ana.getParameter("subType")
        len = ana.getParameter("length")
        a145 = params["a145_a"]
        a2 = params["a2_a"]
        a3 = params["a3_a"]
        return carbon_cost_anatomical(vol, subType, len, ana145, ana2, ana3, a145, a2, a3)

    raise "carbon_cost() unknown model " + model
