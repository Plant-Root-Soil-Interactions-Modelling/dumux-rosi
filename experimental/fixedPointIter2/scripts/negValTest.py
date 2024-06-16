import sys;

sys.path.append("../modules/");
sys.path.append("../inputDataTraiRhizo/");
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
import weatherFunctions
from numpy import array
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

np.set_printoptions(precision=20)
# outer time step (outside of fixed-point iteration loop)
dt = 20/60/24
dt_inner_init = dt#1/60/60/24 # dt
# min, max, objective number of iteration for the fixed-point iteration
minIter = 4 # empirical minimum number of loop to reduce error
k_iter = 30
targetIter= 5
# which functional modules to implement
doSoluteFlow = True # only water (False) or with solutes (True)
noAds = False # stop adsorption?
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

paramIndx_ = 96
spellData = {'scenario':"none",'spellStart':0,'spellEnd':11, 'condition':"wet"}
weatherInit = weatherFunctions.weather(1.,dt, spellData)
s = scenario_setup.create_soil_model(usemoles = usemoles, 
                                           results_dir = results_dir, 
                                           p_mean_ = -weatherInit['p_mean'], 
                                    paramIndx=paramIndx_,
                                    noAds = noAds, doSoluteFlow = doSoluteFlow)

def _reset_newton_solver(cyl,MaxRelativeShift=None,
                                EnableResidualCriterion=None, 
                                EnableAbsoluteResidualCriterion=None,
                                SatisfyResidualAndShiftCriterion=None, 
                                max_steps=None, 
                                max_divisions=None):
    if MaxRelativeShift is not None:
        cyl.setParameter("Newton.MaxRelativeShift", str(MaxRelativeShift))
    if EnableResidualCriterion is not None:
        cyl.setParameter("Newton.EnableResidualCriterion", str(EnableResidualCriterion)) 
    if EnableAbsoluteResidualCriterion is not None:
        cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", str(EnableAbsoluteResidualCriterion))
    if SatisfyResidualAndShiftCriterion is not None:
        cyl.setParameter("Newton.SatisfyResidualAndShiftCriterion", str(SatisfyResidualAndShiftCriterion))
    if max_steps is not None:
        cyl.setParameter("Newton.MaxSteps", str(max_steps))
    if max_divisions is not None:
        cyl.setParameter("Newton.MaxTimeStepDivisions", str(max_divisions))
    cyl.createNewtonSolver()
                
def initialize_dumux_nc_( soilModel,gId=0, a_in=0.02, 
                a_out=0.14469319849003554,seg_length=0.9999999999999999,
                x=[-114.49999999999997 for i in range(9)] ,   # cm
                cAll = [
                    [1.2196341356893417e-04,1.0239102121168215e-04,8.4765239392816990e-05,
6.8698956213938920e-05,5.3244443766544481e-05,3.8059847063775529e-05,
2.4206706500210844e-05,1.3629725904762179e-05,8.0572097658702614e-06],
[0.0007611838346166332,0.0007400622533179592,0.0007167918695770283,
0.0006946311647690286,0.0006757417785640214,0.0006615185911821278,
0.0006527162177727044,0.0006485700413368975,0.0006472163897630698],
[7.7842451102772223e-06,7.7839693574834915e-06,7.7828223470564676e-06,
7.7791765348360603e-06,7.7733792724899486e-06,7.7453090774562467e-06,
7.6132414144994961e-06,7.1317592364065047e-06,5.8228672299783495e-06],
[1.9724522459356758e-06,1.9726292326912743e-06,1.9734817601654045e-06,
1.9763013009014333e-06,1.9807695863746551e-06,2.0031793301229514e-06,
2.1109919686639225e-06,2.5200751910490590e-06,3.6967823544079324e-06],
[6.1643157718718706e-05,6.1192254523506295e-05,6.0391513637403367e-05,
5.8804717148817058e-05,5.5823224098803740e-05,4.9192552774855050e-05,
3.8279158115516276e-05,2.4051305839915572e-05,1.2011078520340050e-05],
[7.4493747190930748e-167,1.0285631241230136e-166,2.0223336113348035e-166,
1.7669149255113308e-151,4.0267093888921536e-106,2.1521049277820653e-019,
3.2617092217595747e-013,1.1614700679545641e-008,1.9258941859702000e-006],
[0.002668461405895332,0.0024617761313399378,0.002251736713488391,
0.0020203118660787316,0.0017381743627321654,0.001386349491973654,
0.0009833521748241472,0.0006054024671620996,0.0003762621516120395],
[7.218639931148065e-05,7.156780439794207e-05,7.046904813089088e-05,
6.829111612211548e-05,6.419836251557750e-05,5.509308320765652e-05,
4.010115530922683e-05,2.055822482932874e-05,6.616728599831435e-06]],             # mol/mol scv
                                         #Cells = np.array([0.0289878 , 0.03302684, 0.03888096, 0.04736588, 0.05966385,0.07748838, 0.1033231 , 0.14076767, 0.19503948]),
                                         NC = 10,
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
        Cells = np.array(cyl.base.getCellCenters()).flatten()*100
        
            
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
        cyl.MaxRelativeShift = soilModel.MaxRelativeShift
        cyl.setParameter("Newton.MaxRelativeShift",str(cyl.MaxRelativeShift))
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
            cyl.setParameter("Newton.MaxRelativeShift","1e-15")#str(soilModel.MaxRelativeShift))
        cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
        cyl.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")     
        cyl.setParameter("Newton.MaxSteps", "50")
        cyl.setParameter("Newton.MaxTimeStepDivisions", "20")
        cyl.setParameter("Newton.EnableResidualCriterion", "true") # sometimes helps, sometimes makes things worse 
            
        cyl.setParameter("Newton.MaxRelativeShift","1e-8")#str(soilModel.MaxRelativeShift))
        cyl.setParameter("Newton.Verbosity", "0") 
        cyl.setParameter("Problem.verbose", "0");
        cyl.setParameter("Soil.CriticalPressure", "-15000")
        cyl.initializeProblem(maxDt = 250/(24*3600))
        cyl.eps_regularization = soilModel.eps_regularization
        if cyl.eps_regularization is not None:
            cyl.setRegularisation(cyl.eps_regularization, cyl.eps_regularization) # needs to be low when using sand parameters. 
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
            print('(abs((pHeadcyl - x)/x_divide)*100 > 1e-5).all()?\n',abs((pHeadcyl - x)/x_divide)*100, '\n pHeadcyl',pHeadcyl,
            '\n x',x,
            '\n x_divide',x_divide)
            raise Exception
        
        return cyl
    else:
        print("RhizoMappedSegments.initialize_dumux_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))
        return []

def _check_Ccontent( cyl):
    for ncomp in range(cyl.numSoluteComp):
        try:
            assert (np.array(cyl.getSolution(ncomp + 1)).flatten() >= 0).all()
        except:
            print(rank, '(np.array(cyl.getSolution(ncomp + 1)).flatten() < 0).any()', 
                  ncomp, np.array(cyl.getSolution(ncomp + 1)).flatten(), 
                  (np.array(cyl.getSolution(ncomp + 1)).flatten() >= 0).all())
            return  False
    return True
def dothesolve(cyl):            
    
    QIN = 0.# cm3/day 
    qIn =QIN / (2 * np.pi * cyl.a_in * cyl.l) 
    cyl.setInnerFluxCyl(qIn)

    Qflowout = 0.00017406685183970837
    cyl.distributeSource(dt, Qflowout,QIN*dt, eqIdx=0)

    inner_fluxes_solMucil_temp= np.array([4.3799793455387955e-06 ,1.4413903510145875e-08,0.,0.,0.,0.,0.,0.,0.] )
    qIn_solMucil = inner_fluxes_solMucil_temp / (2 * np.pi * cyl.a_in  * cyl.l)
    typeBC = np.full(9, 3) 
    cyl.setSoluteBotBC(typeBC, qIn_solMucil)

    outer_fluxes_solMucil =np.array([ 2.81851246e-08 ,-5.38128930e-10]) *0.
    QflowOutCellLim = cyl.distributeSources(dt, source = outer_fluxes_solMucil,
                                   inner_fluxs=inner_fluxes_solMucil_temp * dt, 
                                   eqIdx =  np.array([nc+1 for nc in range(2)]))
    print('\n\nQflowOutCellLim',QflowOutCellLim)
    #raise Exception
    cyl.ddt =min( 1.e-5,cyl.ddt)
    cyl.solve(dt)
    #assert _check_Ccontent(cyl)
    print('phead',cyl.getSolutionHead())
    inner_fluxes_real, outer_fluxes_real = cyl.getSavedBC(cyl.a_in, cyl.a_out)     
    print('inner_fluxes_real', inner_fluxes_real,'QIN', QIN, 'theta', cyl.getWaterContent())
    print('solutes')
    for nc in range(cyl.numComp-1):
        print(cyl.getSolution(nc+1), sum(cyl.getSolution(nc+1)))
        #assert min(cyl.getSolution(nc+1)) >= 0
    if False:
        cyl.reset()
        print('reset')

        print('phead',cyl.getSolutionHead())
        print('solutes')
        for nc in range(cyl.numComp-1):
            print(cyl.getSolution(nc+1), sum(cyl.getSolution(nc+1)))
            assert min(cyl.getSolution(nc+1)) >= 0
        #cyl.base.printParams() # s.base.printParams()
        _reset_newton_solver( cyl,
                                        MaxRelativeShift = s.MaxRelativeShift,
                                        EnableResidualCriterion="false", 
                                        EnableAbsoluteResidualCriterion="false",
                                        SatisfyResidualAndShiftCriterion="false", 
                                        max_steps=18, max_divisions=10)
cyl = initialize_dumux_nc_(s)
for i in range(1):
    dothesolve(cyl)

print('plant vol',repr(cyl.getCellVolumes()))
print('plant point',repr(cyl.getPoints()))
print('setsources', cyl.setSourceBu)