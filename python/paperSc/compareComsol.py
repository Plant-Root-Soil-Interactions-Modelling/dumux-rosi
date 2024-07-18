import sys; sys.path.append("../modules_fpit/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../../CPlantBox/");
sys.path.append("../../../CPlantBox/src")

import matplotlib; matplotlib.use('agg')

from rosi_richards10c_cyl import Richards10CCylFoam  # C++ part (Dumux binding)

from richards import RichardsWrapper  # Python part
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import scenario_setup as stf #
#import smallTest_ads_functions as stf2
import scenario_setup
from scenario_setup import write_file_array, write_file_float

""" 
Cylindrical 1D model, diffusion only (DuMux), Michaelis Menten

everything scripted, no input file needed, also works parallel with mpiexec
"""

results_dir="./results/compareComsol/"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    test = os.listdir(results_dir)
    for item in test:
        try:
            os.remove(results_dir+item)
        except:
            pass
        
dt = 20/60/24
k_iter = 100 #50 
targetIter= 40 #
l_ks =  "dx_2" 
organism = "plant" 
weightBefore = False
SRIBefore = False
beforeAtNight = True
adaptRSI_  = False
static_plant = False
useOuterFluxCyl_w = False
useOuterFluxCyl_sol = False
css1Function_ = 9
lightType ="" 
mpiVerbose = False
mpiVerboseInner = False
noAds = False
doSimple =False
doMinimumPrint =  True
doMean = True
usemoles = True
soil_ = scenario_setup.vg_SPP(0)


weatherInit = scenario_setup.weather(1.,dt, spellData)
s, soil = scenario_setup.create_soil_model(soil_type, year, soil_,#comp, 
            min_b, max_b, cell_number, demoType = mode, times = None, net_inf = None,
            usemoles = usemoles, dirResults = results_dir, p_mean_ = weatherInit['p_mean'], 
                                     css1Function = css1Function_,
                                    paramIndx=paramIndx_,
                                    noAds = noAds)

def initialize_dumux_nc_(soilModel, gId,a_in,a_out,seg_length, x,                                                # cm
                                cAll = [0,0,0,0,0,0,0,0],             # mol/mol scv
                                     Cells = [],NC = 10,
                                         logbase = 0.5):                                   # cm
    verbose = False
    a_in = self.radii[gId]
    a_out = self.outer_radii[gId]
    lId = 0
    #print('to wrap', rank,gId, a_in,a_out )
    if a_in < a_out:
        cyl = RichardsNoMPIWrapper(Richards10CCylFoam(), self.useMoles)  # only works for RichardsCylFoam compiled without MPI
        cyl.initialize(verbose = False)
        cyl.setVGParameters([soilModel.soil])
        lb =  logbase

        if False:
            nCells = self.NC
            cyl.createGrid([0.02], [0.6], [nCells])# cm
        else:
            points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), 
                                 NC, base = lb)
            cyl.createGrid1d(points)# cm
        l_ks = "dx_2"    
        if l_ks == "dx_2":
            cyl.setParameter("Soil.BC.dzScaling", "1")
        elif l_ks == "dx":
            cyl.setParameter("Soil.BC.dzScaling", "2")
        else:
            raise Exception
        cyl.seg_length = seg_length
        cyl.setParameter( "Soil.css1Function", str(soilModel.css1Function))
        cyl.setParameter("Problem.verbose", "0")
        cyl.setParameter("Problem.reactionExclusive", "0")    
        cyl.setParameter("Problem.segLength", str(seg_length))   # cm
        cyl.setParameter( "Soil.Grid.Cells",str( NC-1)) # -1 to go from vertices to cell (dof)
        if self.recreateComsol:
            cyl.setHomogeneousIC(-100.)  # cm pressure head
        else:
            # cyl.setHomogeneousIC(x)  # cm pressure head
            if verbose:
                print("Soil.IC.P", cyl.dumux_str(x), Cells)
            cyl.setParameter("Soil.IC.P", cyl.dumux_str(x))# cm
        # cyl.setICZ_solute(c)  # [kg/m2] 

        #default: no flux
        cyl.setInnerBC("fluxCyl", 0.)  # [cm/day] #Y 0 pressure?
        #cyl.setInnerBC_solute("fluxCyl", 0.)  # [kg/m2], that s fair
        cyl.setOuterBC("fluxCyl", 0.)
        #cyl.setOuterBC_solute("fluxCyl", 0.)


        cyl.setParameter( "Soil.MolarMass", str(soilModel.solidMolarMass))
        cyl.setParameter( "Soil.solidDensity", str(soilModel.solidDensity))
        cyl.setParameter("Flux.UpwindWeight", "1")#very important because we get high solute gradient.
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
        if (soilModel.css1Function == 8):
            cyl_vol = (a_out**2 - a_in**2) * lId * np.pi
            cssmax_ = soilModel.CSSmax * (cyl_vol / soilModel.cell_size) # fraction of cssmax in this cylinder
            #print('cssmax_',soilModel.CSSmax ,cyl_vol / soilModel.cell_size)
        else:
            cssmax_ = soilModel.CSSmax

        cyl.setParameter("Soil.CSSmax", str(cssmax_)) #[mol/cm3 scv zone 1] or mol
        cyl.setParameter("Soil.alpha", str(soilModel.alpha)) #[1/d]


        cyl.setParameter("Soil.C_aOLim", str(soilModel.C_aOLim)) #[molC/cm3 scv]
        cyl.setParameter("Soil.C_aCLim", str(soilModel.C_aCLim)) #[molC/cm3 scv]
        cyl.setParameter("1.Component.LiquidDiffusionCoefficient", str(soilModel.Ds)) #m^2/s

        cyl.setParameter("2.Component.LiquidDiffusionCoefficient", str(soilModel.Dl)) #m^2/s


        for j in range( 1, soilModel.numComp+1):
            cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Type", str(3))
            cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Type", str(3))
            cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Value", str(0)) 
            cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Value", str(0 )) 

        for j in range( 1, soilModel.numComp+1):       
            #print("cAll[j-1]",j-1, cAll[j-1],cyl.dumux_str(cAll[j-1])  )
            cyl.setParameter("Soil.IC.C"+str(j), cyl.dumux_str(cAll[j-1]) ) 

        if len(Cells) > 0:#in case we update a cylinder, to allow dumux to do interpolation
            assert(len(cAll[j-1])==len(Cells))
            CellsStr = cyl.dumux_str(Cells/100)#cm -> m
            cyl.setParameter("Soil.IC.Z",CellsStr)# m
            if len(Cells)!= len(x):
                print("Cells, x",Cells, x, len(Cells), len(x))
                raise Exception
            for j in range( 1, soilModel.numComp+1):
                cyl.setParameter("Soil.IC.C"+str(j)+"Z",CellsStr) # m
                if len(Cells)!= len( cAll[j-1]):
                    print("Cells,  cAll[j-1]",Cells,  cAll[j-1], 
                            len(Cells), len(cAll[j-1]), j)
                    raise Exception

        cyl.setParameter("Newton.MaxRelativeShift",str(soilModel.MaxRelativeShift))
        # cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "true") #<= also helps reach convergence
        cyl.initializeProblem()
        cyl.setCriticalPressure(soilModel.wilting_point)  # cm pressure head
        cyl.bulkDensity_m3 = soilModel.bulkDensity_m3
        cyl.solidDensity =soilModel.solidDensity 
        cyl.solidMolarMass =soilModel.solidMolarMass
        cyl.solidMolDensity =soilModel.solidMolDensity         
        cyl.k_sorp =soilModel.k_sorp               
        cyl.CSSmax = cssmax_ #soilModel.CSSmax               
        cyl.f_sorp =soilModel.f_sorp    
        cyl.css1Function = soilModel.css1Function
        cyl.ddt = 1.e-5
        cyl.gId = gId    
        cyl.l = seg_length    
        cyl.a_in = a_in    
        cyl.a_out = a_out
        cyl.DtCSS2 = soilModel.DtCSS2
        #cyl.base.setComputeDtCSS2(cyl.DtCSS2)  # maps segments
        ThetaCyl = cyl.getWaterContent()
        setDefault(cyl)
        try:
            assert (ThetaCyl >= soilModel.vg_soil.theta_R).all()
            assert (ThetaCyl <= soilModel.vg_soil.theta_S).all()
        except:
            print('issue thetaCyl',rank,ThetaCyl, soilModel.vg_soil.theta_R, soilModel.vg_soil.theta_S )
            raise Exception
        pHeadcyl = cyl.getSolutionHead()
        #print(cyl.getSolutionHead(),x, gId, self.seg2cell[gId])#cyl.getSolutionHead_(),
        #raise Exception
        
        try:
            assert len(pHeadcyl) == (self.NC - 1)
            x_divide = np.where(x!=0,x,1)
            assert (np.logical_or( (abs((pHeadcyl - x)/x_divide)*100 < 1e-5) , 
                                   (abs(pHeadcyl - x) < 1e-9) )).all()
        except:
            print('error: issue with cylinder creations', rank)
            print('len(pHeadcyl) == (self.NC - 1)?', len(pHeadcyl), (self.NC - 1))
            print('(abs((pHeadcyl - x)/x_divide)*100 > 1e-5).all()?',abs((pHeadcyl - x)/x_divide)*100, pHeadcyl,x,'x_divide',x_divide)
            raise Exception

        return cyl
    else:
        print("RhizoMappedSegments.initialize_dumux_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))
        return []
    
    
cyl = initialize_dumux_nc_(s)

s.dirResults = dirResults
cyl.solve(dt, maxDt = 250./(24.*3600.))
  