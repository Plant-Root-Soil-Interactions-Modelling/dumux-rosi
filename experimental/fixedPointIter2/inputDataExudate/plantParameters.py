


import numpy as np
import pandas as pd
import timeit
import pandas as pd
import matplotlib; matplotlib.use('agg')
import sys;
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import timeit
from helpfull import write_file_array, write_file_float, div0, div0f
from functional.xylem_flux import sinusoidal2
from scipy import interpolate
import functional.van_genuchten as vg

# def init_conductivities(r, kr_const = 1.8e-4, kx_const = 0.1):
    # """ Hydraulic conductivities  kr [1/day], kx [cm3/day] """
    # r.setKr([kr_const, kr_const, kr_const, kr_const, kr_const, kr_const])
    # r.setKx([1.e3, kx_const, kx_const, kx_const, kx_const, kx_const])

def init_conductivities(r, skr = 1., skx = 1.):
    """ Hydraulic conductivities for maize modified from Couvreur et al. (2012), originally from Doussan et al. (1998) """ #[age, value]

    #kr00 = np.array([[0., 0.]])  # artificial shoot
    #kx00 = np.array([[0., 1.e3]])  # artificial shoot
    
    # kr00 = np.array([[-1e4, 0.]])  # artificial shoot
    # kx00 = np.array([[0., 0.000864]])  # artificial shoot
    # kr001 = np.minimum(skr * kr00[:, 1], 1.)
    # kx001 = np.minimum(skx * kx00[:, 1], 1.)
    
    #artificial shoot does not seem to be needed? 

    kr0 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.000181], [8., 0.000181], [10, 0.0000648], [18, 0.0000648], [25, 0.0000173], [300, 0.0000173]])
    kr1 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.000181], [10., 0.000181], [16, 0.0000173], [300, 0.0000173]])

    kx0 = np.array([[0., 0.000864], [5., 0.00173], [12., 0.0295], [15., 0.0295], [20., 0.432], [300., 0.432]])
    kx1 = np.array([[0., 0.000864], [5., 0.0000864], [10., 0.0000864], [12., 0.0006048], [20., 0.0006048], [23., 0.00173], [300., 0.00173]])

    kr01 = np.minimum(skr * kr0[:, 1], 1.)
    kr11 = np.minimum(skr * kr1[:, 1], 1.)
    
    # r.setKr([kr001, kr01, kr11, kr11, kr01, kr01],
                  # [kr00[:, 0], kr0[:, 0], kr1[:, 0], kr1[:, 0], kr0[:, 0], kr0[:, 0]])
    r.setKr([kr01, kr11, kr11, kr01, kr01],
                  [kr0[:, 0], kr1[:, 0], kr1[:, 0], kr0[:, 0], kr0[:, 0]])
    kx01 = np.minimum(skx * kx0[:, 1], 1.)
    kx11 = np.minimum(skx * kx1[:, 1], 1.)

    r.setKx([kx01, kx11, kx11, kx01, kx01],
                  [kx0[:, 0], kx1[:, 0], kx1[:, 0], kx0[:, 0], kx0[:, 0]])    
    # r.setKx([kx001, kx01, kx11, kx11, kx01, kx01],
                  # [kx00[:, 0], kx0[:, 0], kx1[:, 0], kx1[:, 0], kx0[:, 0], kx0[:, 0]])


    return r
    
    
def prescribed_exudation(soil_type, ifexu):
    
    times = [0, 42, 63, 98, 154]
    
    if ifexu == "True": 
        if soil_type == 'loam':
            exu_prop = np.array([0.001,0.001,0.00055,0.00039,0.00045])#[kg C/(m2 root surface  day)] 
        elif soil_type == 'sand':
            exu_prop = np.array([0.0011,0.0011,0.0005,0.0003,0.00033]) #[kg C/(m2 root surface day)] 
        else:
            print('No exudate properties found')
    else: 
        exu_prop = np.array([0.,0.,0.,0.,0.]) #[kg/(m2 day)]

    f = interpolate.interp1d(times, exu_prop)  
    
    return f     

def exudation_rates(f, t):   
    
    kex = np.array([[0., 3.5], [f(t),f(t)/2]])
        
    return kex
