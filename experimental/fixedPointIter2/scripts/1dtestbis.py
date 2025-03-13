import sys;
import os
sys.path.append("../modules/");
sys.path.append("../inputDataPuptake/");
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
import helpfull

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
dt = 0.013888888888888888 #20/60/24
dt_inner_init = dt#1/60/60/24 # dt
# min, max, objective number of iteration for the fixed-point iteration
minIter = 4 # empirical minimum number of loop to reduce error
k_iter = 30
targetIter= 5
# which functional modules to implement
doSoluteUptake = True # active plant solute uptake
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

path = "../inputDataPuptake/"
paramIndx_ = 44
spellData = {'scenario':"none",'spellStart':0,'spellEnd':11, 'condition':"wet"}
weatherInit = weatherFunctions.weather(1.,dt, spellData)
s = scenario_setup.create_soil_model(usemoles = usemoles, 
                                           results_dir = results_dir, 
                                           p_mean_ = -weatherInit['p_mean'], 
                                    paramIndx=paramIndx_,
                                    noAds = noAds, doSoluteFlow = doSoluteFlow,
                                         doBioChemicalReaction = False,
                                        MaxRelativeShift = 1e-8)


# all thread need a plant object, but only thread 0 will make it grow
perirhizalModel, plantModel = scenario_setup.create_mapped_plant(8., s, "P3.xml" ,
                                        path, 
                                        doPhloemFlow = False,
                                        static_plant = False,
                                        usemoles = usemoles,
                                        limErr1d3d = 5e-12, spellData = spellData)  



def initialize_dumux_nc_( soilModel, gId=0, a_in=0.02, 
                a_out=0.25584836252689325,seg_length=0.25,
                x=[-135.00735278137205, -134.9588223416131 , -134.91140713271668,
 -134.86560970364644, -134.8221731248102  ,-134.78219633367598,
 -134.74730510699206 ,-134.71990585412667 ,-134.70356315666214] ,   # cm
                cAll = [
                    [1.6241395869693224e-10, 4.4301388549752135e-08,
       9.7236021083269233e-09, 2.1996686065104217e-08,
       2.0165130948001519e-08, 2.0279337046633592e-08,
       2.0275276068003719e-08, 2.0275046865587812e-08,
       2.0274959297533126e-08],
        [0., 0., 0., 0., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0., 0., 0.]],             # mol/mol scv
                                         Cells = [],NC = 10,
                                         logbase = 0.5):                                   # cm
    verbose = False
    lId =gId
    
    if a_in < a_out:

        cyl = RichardsNoMPIWrapper(RichardsNCCylFoam(), soilModel.useMoles)  # only works for RichardsCylFoam compiled without MPI
        cyl.pindx = soilModel.pindx
        cyl.results_dir = soilModel.results_dir
        cyl.setParameter("Newton.Verbosity", "0") 
        cyl.initialize(verbose = False)
        cyl.setVGParameters([soilModel.soil])
        lb = logbase

        points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), 
                             NC, base = lb)

        cyl.createGrid1d(points)# cm
        Cells = np.array(cyl.base.getCellCenters()).flatten()*100

        cyl.setParameter("Problem.segLength", str(seg_length))  # cm
        cyl.setParameter("SpatialParams.Temperature","293.15") # todo: redefine at each time step
        cyl.setParameter("Soil.BC.dzScaling", "1")
        cyl.setParameter( "Soil.css1Function", str(soilModel.css1Function))
        cyl.setParameter("Problem.verbose", "1")
        cyl.setParameter("Problem.reactionExclusive", "0")
        cyl.setParameter("Soil.CriticalPressure", str(soilModel.wilting_point))
        if verbose:
            print("Soil.IC.P", cyl.dumux_str(x), Cells)
        cyl.setParameter("Soil.IC.P", cyl.dumux_str(x))# cm

        #default: no flux
        cyl.setInnerBC("fluxCyl", 0.)  # [cm/day] #Y 0 pressure?
        cyl.setOuterBC("fluxCyl", 0.)

        cyl.setParameter("Soil.MolarMass", str(soilModel.solidMolarMass))
        cyl.setParameter("Soil.solidDensity", str(soilModel.solidDensity))
        cyl.setParameter("Flux.UpwindWeight", "1")
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

        if doSoluteUptake:
            # for maiz:
            # http://dx.doi.org/10.1590/S0103-90162004000100012
            #Vmax = 45,6 µmol / g /h
            #Km = 33,7µmol / l
            #average root surafce area: 818,33 cm²
            #average root weight: 1,51 g


            RS_Uptake_Vmax = 6.56 *  1e-12 * 24 * 3600 # pmol / cm² / s =>  mol / cm² / d
            RS_Uptake_km =  7.1 *1e-3*1e-6# mumol/L=> mol/cm3
            #Vmax = 3.844e-10*(24*3600)/1e4# (kg m–2 s–1) * (s/d) * (m2/cm2) => (
            # kg cm–2 d–1)
            #Vmax = soilModel.Vmax#Vmax * 1000. / soilModel.molarMassC # (kg cm–2 d–1) * (g/kg) / (g/mol) => mol cm-2 d-1
            cyl.setParameter("RootSystem.Uptake.Vmax", cyl.dumux_str(RS_Uptake_Vmax))  # mol /cm^2 / s - > mol /cm^2 / day 
            #km = 1.054e-4 * 1e-6 # (kg m–3) => kg cm-3
            #km = km * 1000. / soilModel.molarMassC # (kg cm-2 * (g/kg) / (g/mol)) 
            cyl.setParameter("RootSystem.Uptake.Km", cyl.dumux_str(RS_Uptake_km))  # mol / cm3                
            cyl.setParameter( "Soil.BC.Bot.C1Type", str(8))

        # Maximum uptake rate (kg m–2 s–1)	3.844e-10	(Teo et al. (1992a)), from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7489101/
        # Michaelis constant (kg m–3)	1.054e-4	(Teo et al. (1992b)), from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7489101/

        cyl.doBioChemicalReaction = soilModel.doBioChemicalReaction
        cyl.molarMassC = soilModel.molarMassC
        cyl.MaxRelativeShift = soilModel.MaxRelativeShift_1DS
        cyl.EnableResidualCriterion = soilModel.EnableResidualCriterion
        cyl.EnableAbsoluteResidualCriterion = soilModel.EnableAbsoluteResidualCriterion
        cyl.SatisfyResidualAndShiftCriterion = soilModel.SatisfyResidualAndShiftCriterion
        cyl.MaxTimeStepDivisions = soilModel.MaxTimeStepDivisions
        cyl.MaxSteps = soilModel.MaxSteps
        cyl.initializeProblem(maxDt = soilModel.maxDt_1DS)#soilModel.maxDt)
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
        cyl.ddt = 1.e-3
        cyl.gId = gId    
        #cyl.theta_wilting_point = soilModel.theta_wilting_point    
        cyl.vg_soil = soilModel.vg_soil
        cyl.a_in = a_in    
        cyl.a_out = a_out
        ThetaCyl = cyl.getWaterContent()
        scenario_setup.setDefault(cyl)
        cyl.createNewtonSolver() # make sure solver parameters are implemented.
        try:
            assert (ThetaCyl >= soilModel.vg_soil.theta_R).all()
            assert (ThetaCyl <= soilModel.vg_soil.theta_S).all()
        except:
            print('issue thetaCyl',rank,ThetaCyl, soilModel.vg_soil.theta_R, soilModel.vg_soil.theta_S )
            raise Exception
        pHeadcyl = cyl.getSolutionHead()

        if verbose:
            print('end initialize_',gId,soilModel.seg2cell[gId],'wat vol?',sum(cyl.getWaterVolumesCyl()),soilModel.seg_length[gId],
              pHeadcyl , x,pHeadcyl - x,'Cells',Cells, a_in, a_out, cyl.segLength)


        try:
            x_divide = np.where(np.array(x)!=0,x,1)
            thetainit = cyl.getWaterContent()
            thetainitth = np.array([vg.water_content( p_mean_, soilModel.vg_soil) for
                                    p_mean_ in pHeadcyl])
            theta_divide = np.where(thetainitth!=0,thetainitth,1)

            assert abs(seg_length -
                       cyl.segLength) < 1e-16 or abs((seg_length -
                       cyl.segLength)/cyl.segLength) < 1e-5
            # we could have some issue because of the water
            # retention curve regularisation in dumux
            assert len(pHeadcyl) == (NC - 1)
            assert (np.logical_or( (abs((pHeadcyl - x)/x_divide)*100 < 1e-5) , 
                                   (abs(pHeadcyl - x) < 1e-9) )).all()
            assert (np.logical_or( (abs((thetainit - thetainitth)/theta_divide)*100 < 1e-5) ,
                                   (abs(thetainit- thetainitth) < 1e-9) )).all()
        except:
            print('error: issue with cylinder creations', rank)
            print('len(pHeadcyl) == (NC - 1)?', len(pHeadcyl), (NC - 1))
            print('(abs((pHeadcyl - x)/x_divide)*100 > 1e-5).all()?',abs((pHeadcyl - x)/x_divide)*100, pHeadcyl,x,'x_divide',x_divide)
            print('theta',thetainit,thetainitth)
            raise Exception
        return cyl

    else:
        print("RhizoMappedSegments.initialize_dumux_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(gId, a_in, a_out))
        raise Exception
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
    
    QIN = -0.03499198906790653 *0# cm3/day 
    qIn =QIN / (2 * np.pi * cyl.a_in * cyl.segLength) 
    cyl.setInnerFluxCyl(qIn)

    Qflowout = 0.008768593546613732 
    qOut = cyl.distributeSource(dt, Qflowout,QIN*dt, eqIdx=0, plantM = perirhizalModel)
    print('qOut',qOut)
    print('theta', cyl.getWaterContent())
    Watvol=cyl.getWaterVolumes()
    print('Watvol',Watvol )
    cellVol = cyl.getCellVolumes()
    print('added theta',qOut/cellVol)
    print((sum(Watvol)-0.03499198906790653 *dt)/sum(cellVol))
    #raise Exception
    inner_fluxes_solMucil_temp= np.array([0.,0.,
                                          0.,0.,0.,0.,0.,0.] )
    qIn_solMucil = inner_fluxes_solMucil_temp / (2 * np.pi * cyl.a_in  * cyl.segLength)
    typeBC = np.full(8, 3) 
    cyl.setSoluteBotBC(typeBC, qIn_solMucil)

    outer_fluxes_solMucil =np.array([0.,0.]) 
    QflowOutCellLim = cyl.distributeSources(dt, source = outer_fluxes_solMucil,
                                   inner_fluxs=inner_fluxes_solMucil_temp * dt, 
                                   eqIdx =  np.array([nc+1 for nc in range(2)]), plantM = perirhizalModel)
        
    cyl.ddt =min( 1.e-5,cyl.ddt)
    #cyl.solve(dt)
    print("go to solve")
    helpfull.run_with_timeout(5., cyl.solve, dt) # after Xmn time out
    print("finished solve")
    
    assert _check_Ccontent(cyl)
    print('phead',cyl.getSolutionHead())
    inner_fluxes_real, outer_fluxes_real = cyl.getSavedBC(cyl.a_in, cyl.a_out)     
    print('inner_fluxes_real', inner_fluxes_real,'QIN', QIN, 'theta', cyl.getWaterContent(), cyl.getSolutionHead() )
    print('solutes')
    for nc in range(cyl.numComp-1):
        print(cyl.getSolution(nc+1), sum(cyl.getSolution(nc+1)))
        assert min(cyl.getSolution(nc+1)) >= 0
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
    #cyl.base.printParams()
    dothesolve(cyl)
