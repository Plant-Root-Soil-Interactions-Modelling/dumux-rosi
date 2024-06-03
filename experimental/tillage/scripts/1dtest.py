import sys;

sys.path.append("../modules/");
sys.path.append("../inputData/");
sys.path.append("../../../../CPlantBox/");
sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../../build-cmake/cpp/python_binding/");

import matplotlib; matplotlib.use('agg')

from rosi_richards10c_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)

from richards import RichardsWrapper  # Python part
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
# import smallTest_ads_functions as stf
import scenario_setup as stf
from scenario_setup import write_file_array, write_file_float
import scenario_setup
import functional.van_genuchten as vg

""" 
Cylindrical 1D model, diffusion only (DuMux), Michaelis Menten

everything scripted, no input file needed, also works parallel with mpiexec
"""

# directory where the results will be printed
results_dir="./results/1dtest/"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    test = os.listdir(results_dir)
    for item in test:
        try:
            os.remove(results_dir+item)
        except:
            pass

# outer time step (outside of fixed-point iteration loop)
dt = 20/60/24
dt_inner_init = dt#1/60/60/24 # dt
# min, max, objective number of iteration for the fixed-point iteration
minIter = 4 # empirical minimum number of loop to reduce error
k_iter = 30
targetIter= 5
# which functional modules to implement
doSoluteFlow = False # only water (False) or with solutes (True)
noAds = True # stop adsorption?
doPhloemFlow = False
doPhotosynthesis = False # photosynthesis-Transpiration (True) or just xylem flow (False)?
# when there is no transpiration, we use the plant wat. pot.
# at the beginning of the time step. Otherwise does not converge
beforeAtNight = True  
# static or growing organism
static_plant = False
# print debug messages in cyl3plant faile
mpiVerbose = False
# print debug messages in rhizo_mappedPlant
mpiVerboseInner = False
# how many files are printed. use 'False' in debug mode
# ATT: for short ismulations only
doMinimumPrint =  True
# use moles (mol) and not mass (g) in dumux
usemoles = True

paramIndx_ = 0
spellData = {'scenario':"none",'spellStart':0,'spellEnd':11, 'condition':"wet"}
weatherInit = scenario_setup.weather(1.,dt, spellData)
s = scenario_setup.create_soil_model(usemoles = usemoles, 
                                           results_dir = results_dir, 
                                           p_mean_ = -weatherInit['p_mean'], 
                                    paramIndx=paramIndx_,
                                    noAds = noAds, doSoluteFlow = doSoluteFlow)

                          
"""
start solve, water in  -0.2031812559287983 water out 0.0 array([0.06236577, 0.06236577, 0.06236577, 0.06236577, 0.06236577,
       0.06236577, 0.06236577, 0.06236577, 0.06236577]) array([-117., -117., -117., -117., -117., -117., -117., -117., -117.]) shape array([0.02      , 0.02947389, 0.04343551, 0.06401067, 0.09433218,
       0.13901681, 0.20486831, 0.3019133 , 0.44492796, 0.65568789]) array([0.00036812, 0.00079948, 0.0017363 , 0.00377085, 0.00818944,
       0.01778562, 0.03862636, 0.08388779, 0.18218545]) 0.02 0.6556878931143485

"""
def initialize_dumux_nc_( soilModel, gId=0, a_in=0.02, a_out=0.6556878931143485,seg_length=0.25,
                                    x=-117.,   # cm
                                    cAll = [0.,0.,0.,0.,0.,0., 0.,0.],             # mol/mol scv
                                         Cells = [],NC = 10,
                                         logbase = 0.5):                                   # cm
    verbose = False
    lId =gId
    
    if a_in < a_out:
    
        cyl = RichardsNoMPIWrapper(RichardsNCCylFoam(), usemoles)  # only works for RichardsCylFoam compiled without MPI
        cyl.initialize(verbose = False)
        cyl.setVGParameters([soilModel.soil])
        lb =  logbase
        
        points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), 
                              NC, base = lb)
        
        cyl.createGrid1d(points)# cm
            
        cyl.setParameter("SpatialParams.Temperature","293.15") # todo: redefine at each time step
        cyl.setParameter("Soil.BC.dzScaling", "1")
        cyl.setParameter( "Soil.css1Function", str(soilModel.css1Function))
        cyl.setParameter("Problem.verbose", "0")
        cyl.setParameter("Problem.reactionExclusive", "0")    
        cyl.seg_length = seg_length
        cyl.setParameter("Problem.segLength", str(seg_length))   # cm
        cyl.l = seg_length   
        cyl.setParameter( "Soil.Grid.Cells",str( NC-1)) # -1 to go from vertices to cell (dof)
        if verbose:
            print("Soil.IC.P", cyl.dumux_str(x), Cells)
        cyl.setParameter("Soil.IC.P", cyl.dumux_str(x))# cm
        
        #default: no flux
        cyl.setInnerBC("fluxCyl", 0.)  # [cm/day] #Y 0 pressure?
        cyl.setOuterBC("fluxCyl", 0.)
        
        cyl.setParameter("Soil.MolarMass", str(soilModel.solidMolarMass))
        cyl.setParameter("Soil.solidDensity", str(soilModel.solidDensity))
        cyl.setParameter("Flux.UpwindWeight", "0.5")
        cyl.setParameter("Soil.betaC", str(soilModel.betaC ))
        cyl.setParameter("Soil.betaO", str(soilModel.betaO))
        cyl.setParameter("Soil.C_S_W_thresC", str(soilModel.C_S_W_thresC )) #mol/cm3
        cyl.setParameter("Soil.C_S_W_thresO", str(soilModel.C_S_W_thresO )) #mol/cm3
        cyl.setParameter("Soil.k_decay", str(soilModel.k_decay))
        cyl.setParameter("Soil.k_decay2", str(soilModel.k_decay2 ))
        cyl.setParameter("Soil.k_DC", str(soilModel.k_DC  )) # 1/d
        cyl.setParameter("Soil.k_DO", str(soilModel.k_DO  )) # 1/d
        cyl.setParameter("Soil.k_growthC", str(soilModel.k_growthC))
        cyl.setParameter("Soil.k_growthO", str(soilModel.k_growthO))
        cyl.setParameter("Soil.K_L", str(soilModel.K_L))#[mol/cm3]
        cyl.setParameter("Soil.k_phi", str(soilModel.k_phi ))
        cyl.setParameter("Soil.k_RC", str(soilModel.k_RC))
        cyl.setParameter("Soil.k_RO", str(soilModel.k_RO ))

        cyl.setParameter("Soil.k_SC", str(soilModel.k_SC )) #cm^3/mol/d
        cyl.setParameter("Soil.k_SO", str(soilModel.k_SO )) #cm^3/mol/d
        
        cyl.setParameter("Soil.m_maxC", str(soilModel.m_maxC  ))# 1/d
        cyl.setParameter("Soil.m_maxO", str(soilModel.m_maxO  ))# 1/d
        cyl.setParameter("Soil.micro_maxC", str(soilModel.micro_maxC ))# 1/d
        cyl.setParameter("Soil.micro_maxO", str(soilModel.micro_maxO ))# 1/d
        cyl.setParameter("Soil.v_maxL", str(soilModel.v_maxL))#[d-1]

        cyl.setParameter("Soil.k_sorp", str(soilModel.k_sorp)) # mol / cm3 or mol
        cyl.setParameter("Soil.f_sorp", str(soilModel.f_sorp)) #[-]

        cyl.setParameter("Soil.kads", str(soilModel.kads)) #[cm3/mol/d]
        cyl.setParameter("Soil.kdes", str(soilModel.kdes)) #[1/d]            
        cyl.setParameter("Soil.CSSmax", str(soilModel.CSSmax)) #[mol/cm3 scv zone 1] or mol
        cyl.setParameter("Soil.alpha", str(soilModel.alpha)) #[1/d]


        cyl.setParameter("Soil.C_aOLim", str(soilModel.C_aOLim)) #[molC/cm3 scv]
        cyl.setParameter("Soil.C_aCLim", str(soilModel.C_aCLim)) #[molC/cm3 scv]
        cyl.setParameter("1.Component.LiquidDiffusionCoefficient", str(soilModel.Ds)) #m^2/s

        cyl.setParameter("2.Component.LiquidDiffusionCoefficient", str(soilModel.Dl)) #m^2/s

        
        for j in range( 1, soilModel.numComp):
            cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Type", str(3))
            cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Type", str(3))
            cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Value", str(0)) 
            cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Value", str(0 )) 
            
        for j in range( 1, soilModel.numComp):       
            cyl.setParameter("Soil.IC.C"+str(j), cyl.dumux_str(cAll[j-1]) ) 

        if len(Cells) > 0:#in case we update a cylinder, to allow dumux to do interpolation
            assert(len(cAll[j-1])==len(Cells))
            CellsStr = cyl.dumux_str(Cells/100)#cm -> m
            cyl.setParameter("Soil.IC.Z",CellsStr)# m
            if len(Cells)!= len(x):
                print("Cells, x",Cells, x, len(Cells), len(x))
                raise Exception
            for j in range( 1, soilModel.numComp):
                cyl.setParameter("Soil.IC.C"+str(j)+"Z",CellsStr) # m
                if len(Cells)!= len( cAll[j-1]):
                    print("Cells,  cAll[j-1]",Cells,  cAll[j-1], 
                            len(Cells), len(cAll[j-1]), j)
                    raise Exception
        if False:
            cyl.setParameter("Newton.MaxSteps", "200")
            cyl.setParameter("Newton.MaxTimeStepDivisions", "100")
            cyl.setParameter("Newton.EnableResidualCriterion", "true") # sometimes helps, sometimes makes things worse
            cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
            cyl.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")            
            cyl.setParameter("Newton.MaxRelativeShift","1e-15")#str(soilModel.MaxRelativeShift))
        cyl.setParameter("Newton.Verbosity", "0") 
        cyl.setParameter("Problem.verbose", "1");
        cyl.setParameter("Soil.CriticalPressure", "-15000")
        cyl.initializeProblem(maxDt = 250/(3600*24))
        cyl.setCriticalPressure(soilModel.wilting_point)  # cm pressure head
        cyl.bulkDensity_m3 = soilModel.bulkDensity_m3
        cyl.solidDensity =soilModel.solidDensity 
        cyl.solidMolarMass =soilModel.solidMolarMass
        cyl.solidMolDensity =soilModel.solidMolDensity         
        cyl.k_sorp =soilModel.k_sorp               
        cyl.CSSmax = soilModel.CSSmax               
        cyl.f_sorp =soilModel.f_sorp    
        cyl.css1Function = soilModel.css1Function
        cyl.ddt = 1.e-5
        cyl.gId = gId    
        cyl.theta_wilting_point = vg.water_content( soilModel.wilting_point, soilModel.vg_soil)
        cyl.vg_soil = soilModel.vg_soil   
        cyl.a_in = a_in    
        cyl.a_out = a_out
        ThetaCyl = cyl.getWaterContent()
        scenario_setup.setDefault(cyl)
        try:
            assert (ThetaCyl >= soilModel.vg_soil.theta_R).all()
            assert (ThetaCyl <= soilModel.vg_soil.theta_S).all()
        except:
            print('issue thetaCyl',rank,ThetaCyl, soilModel.vg_soil.theta_R, soilModel.vg_soil.theta_S )
            raise Exception
        pHeadcyl = cyl.getSolutionHead()
        
        if verbose:
            print('end initialize_',gId,'wat vol?',sum(cyl.getWaterVolumesCyl()),seg_length,
              pHeadcyl , x,pHeadcyl - x,'Cells',Cells)
        try:
            assert len(pHeadcyl) == (NC - 1)
            x_divide = np.where(x!=0,x,1)
            assert (np.logical_or( (abs((pHeadcyl - x)/x_divide)*100 < 1e-5) , 
                                   (abs(pHeadcyl - x) < 1e-9) )).all()
        except:
            print('error: issue with cylinder creations', rank)
            print('len(pHeadcyl) == (NC - 1)?', len(pHeadcyl), (NC - 1))
            print('(abs((pHeadcyl - x)/x_divide)*100 > 1e-5).all()?',abs((pHeadcyl - x)/x_divide)*100, pHeadcyl,x,'x_divide',x_divide)
            raise Exception
        
        return cyl
    else:
        print("RhizoMappedSegments.initialize_dumux_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))
        return []
            
cyl = initialize_dumux_nc_(s)
QIN = -0.2031812559287983 # cm3/day
qIn =QIN / (2 * np.pi * cyl.a_in * cyl.l) 
cyl.setInnerFluxCyl(qIn)
cyl.ddt =min( 1.e-5,cyl.ddt)
cyl.solve(dt)
inner_fluxes_real, outer_fluxes_real = cyl.getSavedBC(cyl.a_in, cyl.a_out)     
print('inner_fluxes_real', inner_fluxes_real,'QIN', QIN, 'theta', cyl.getWaterContent(), cyl.getSolutionHead() )