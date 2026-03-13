import sys;
import os

#os.chdir('experimental/fixedPointIter2/scripts')
sys.path.append("../modules/");
sys.path.append("../inputDataExudate/");
sys.path.append("../../../../CPlantBox/");
sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../../build-cmake/cpp/python_binding/");


from rosi_richards4c_cyl import Richards4CCylFoam  # C++ part (Dumux binding)

from richards import RichardsWrapper  # Python part
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
import pandas as pd
import functional.van_genuchten as vg
import os
import scenario_setup
matplotlib.use('TkAgg', force=True)

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


def getSoilTextureAndShape():  
    """ soil shape and texture data
        to adapt according to the soil represented
    """
    solidDensity = 2650 # [kg/m^3 solid] #taken from google docs TraiRhizo
    solidMolarMass = 60.08e-3 # [kg/mol] 
    # theta_r, theta_s, alpha, n, Ks
    soilVG = [0.045, 0.43, 0.04, 1.6, 50]
    soilTextureAndShape = {
                           "solidDensity":solidDensity,
                        'solidMolarMass': solidMolarMass,
                           'soilVG':soilVG}
    
    return soilTextureAndShape
    
def setSoilParam(s):    
    """ save the soil parameters
        @param: the dumux soil object
    """
    soilTexture = getSoilTextureAndShape()
    s.solidDensity = soilTexture['solidDensity'] #[kg/m^3 solid] 
    s.solidMolarMass = soilTexture['solidMolarMass']# [kg/mol] 
    s.soil =  soilTexture['soilVG'] 
    
    s.vg_soil = vg.Parameters(s.soil) 
    # [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
    s.solidMolDensity = s.solidDensity/s.solidMolarMass
    # [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
    s.bulkDensity_m3 = s.solidMolDensity*(1.- s.vg_soil.theta_S)
    s.bulkMassDensity_gpercm3 = s.solidDensity*(1.- s.vg_soil.theta_S)*1000/1e6

    s.setParameter( "Soil.MolarMass", str(s.solidMolarMass))
    s.setParameter( "Soil.solidDensity", str(s.solidDensity))
    s.setVGParameters([s.soil])
    
    return s
                          

solidDensity = 2650 # [kg/m^3 solid] #taken from google docs TraiRhizo
solidMolarMass = 60.08e-3 # [kg/mol] 
theta_S = 0.494
# [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
solidMolDensity = solidDensity/solidMolarMass
# [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
bulkDensity_m3 = solidMolDensity*(1.- theta_S)
bulkMassDensity_gpercm3 = solidDensity*(1.- theta_S)*1000/1e6


soil_type = 0
kads = 7.07e+02 # m3/kgC/yr, see 10.1016/j.soilbio.2020.107912, A.3
kdes =  1.63e+03 # [1/yr] see 10.1016/j.soilbio.2020.107912, A.3
k_clay_silt = {}
k_clay_silt[0] = 0.67
k_clay_silt[1] = 0.082
molarMassC = 12.011
mg_per_molC = molarMassC * 1000.
yr_per_d = 1/365 # [yr/d]
m3_per_cm3 = 1e-6; # m3/cm3
cm3_per_m3 = 1e6; # cm3/m3

# [kg/g] * [g/mol] = kg/mol
kgC_per_mol = (1/1000) * molarMassC
# [m3/kgC/yr] * [yr/d] * [cm3/m3] * [kgC/mol] = [cm3/mol/d]
kads = kads * yr_per_d * cm3_per_m3 * kgC_per_mol # [cm3/mol/d]
kdes = kdes * yr_per_d # [1/d]
Qmmax = k_clay_silt[soil_type] * 0.079 # max ratio gOC-gmineral soil, see 10.1016/j.soilbio.2020.107912
# [g OC / g mineral soil] * [g mineral soil/ cm3 bulk soil] *[ mol C/g C]
CSSmax_ = Qmmax * bulkMassDensity_gpercm3*(1/molarMassC)
CSSmax = CSSmax_ # mol C/cm3 bulk soil

def initialize_dumux_nc_(  gId=0, a_in=0.02,
                a_out=0.03 ,seg_length=1.0 ,
                x=[-10000] ,   # cm
                NC = 10,logbase = 0.5, doCells = False,points = [0.02,0.03]):                                   # cm
    verbose = False
    lId =gId
    
    if a_in < a_out:
    
        cyl = RichardsNoMPIWrapper(Richards4CCylFoam())  # only works for RichardsCylFoam compiled without MPI
        if False:
            cyl.setParameter("Newton.EnableResidualCriterion", "true") # sometimes helps, sometimes makes things worse
            cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
            cyl.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
            cyl.setParameter("Newton.MaxRelativeShift", "1e-9")# reset value
        cyl.initialize(verbose = False) # No parameter file found. Continuing without parameter file.
        #soilVG = [0.08, 0.43, 0.04, 1.6, 50]
        #cyl.setVGParameters([soilVG])
        lb =  logbase
        
        #points = getPoints()
        
        cyl.createGrid1d(points)# cm
        if doCells:
            Cells = cyl.getCellCenters_().reshape(-1)
        else:
            Cells = []
        cyl.setParameter("Soil.kads", str(kads)) #[cm3/mol/d]
        cyl.setParameter("Soil.kdes", str(kdes)) #[1/d]            
        cyl.setParameter("Soil.CSSmax", str(CSSmax)) #[mol/cm3 scv zone 1] or mol
        cyl.setParameter("Soil.vmax_decay", str(0)) #mol C / m^3 scv / s
        cyl.setParameter("Problem.doDecay",str(0))
        cyl.wilting_point = -15000

        cyl.setParameter("Flux.UpwindWeight", "1")
        cyl.setParameter("SpatialParams.Temperature","293.15") # todo: redefine at each time step
        cyl.setParameter("Soil.BC.dzScaling", "1")
        cyl.setParameter( "Soil.css1Function", str(9))
        cyl.setParameter("Problem.verbose", "0")
        cyl.setParameter("Problem.reactionExclusive", "0")    
        cyl.setParameter("Soil.CriticalPressure", str(-15000))
        cyl.seg_length = seg_length
        cyl.setParameter("Problem.segLength", str(seg_length))   # cm
        cyl.l = seg_length   
        cyl.setParameter( "Soil.Grid.Cells",str( NC-1)) # -1 to go from vertices to cell (dof)
        print("Soil.IC.P", cyl.dumux_str(x))
        cyl.setParameter("Soil.IC.P", cyl.dumux_str(x))# cm
        
        #default: no flux
        cyl.setInnerBC("fluxCyl", -1)#-0.1)  # [cm/day] #Y 0 pressure?
        cyl.setOuterBC("fluxCyl", 0.1)
        cyl.setParameter("Newton.MaxTimeStepDivisions",
                     str( 100) )
        cyl.setParameter("Newton.MaxSteps",
                     str( 100) )
            
        
        for j in range( 1, cyl.numComp):
            cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Type", str(3))
            cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Type", str(3))
            cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Value", str(0)) 
            cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Value", str(0 ))

        cyl.setParameter( "Soil.BC.Bot.C1Value", str(1.)) 
        #cyl.setParameter( "Soil.BC.Top.C1Value", str(-0.1))
        
        
        
        if len(Cells) > 0:#in case we update a cylinder, to allow dumux to do interpolation
            #assert(len(c2)==len(Cells))
            c2 = np.zeros(len(Cells))
            CellsStr = cyl.dumux_str(Cells/100)#cm -> m
            cyl.setParameter("Soil.IC.Z",CellsStr)# m
            if len(Cells)!= len(x):
                print("Cells, x",Cells, x, len(Cells), len(x))
                raise Exception
                
            j = 2
            cyl.setParameter("Soil.IC.C"+str(j)+"Z",CellsStr) # m
            if len(Cells)!= len( c2):
                print("Cells,  cAll[j-1]",Cells,  c2, 
                        len(Cells), len(c2), j)
                raise Exception

        # cyl.setParameter("Soil.BC.Top.SType", "3")  # michaelisMenten=8 (SType = Solute Type)
        # cyl.setParameter("Soil.BC.Top.CValue", "-0.05")  # michaelisMenten=8 (SType = Solute Type)
        # cyl.setParameter("Soil.BC.Bot.SType", "3")  # michaelisMenten=8 (SType = Solute Type)
        # cyl.setParameter("Soil.BC.Bot.CValue", "0.1")
        cyl.setParameter("Soil.IC.C1", cyl.dumux_str(1) ) 

        cyl.maxDt = 250/(3600*24) # soilModel.maxDt_1DS
        cyl.setParameter("Soil.mucilCAffectsW", "true")
        cyl.setParameter("Problem.verbose", "0")
        cyl.setParameter("Newton.Verbosity", "0")
        setSoilParam(cyl)
        cyl.initializeProblem(maxDt=cyl.maxDt ) # ewton solver configured with the following options and parameters:

        
        cyl.eps_regularization = 1e-14
        cyl.setRegularisation(cyl.eps_regularization, cyl.eps_regularization)

        cyl.setCriticalPressure(-15000)  # cm pressure head
        
        #cyl.setRegularisation(1e-16, 1e-16)
        
        return cyl
    else:
        print("RhizoMappedSegments.initialize_dumux_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))
        return []
        
#cyl = initialize_dumux_nc_(c2 = mlFr, NC = NC,x=[-10000],a_in=a_in)
lb = 0.5
a_in=0.02
a_out =0.03
NC = 10        
#points = np.array([a_in + (a_out - a_in)*i/(NC-1) for i in range(NC)])
points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), 
                              NC, base = lb)
CC = np.array([(points[i] + points[i+1])/2 for i in range(len(points)-1)])#/100
print('CC1',CC[[0,-1]])
#raise Exception
cyl = initialize_dumux_nc_(x=[-100. for i in range(NC-1)],a_in=a_in,a_out=a_out,
                            NC = NC, doCells = True,points=points)
length = cyl.seg_length
cellVolumes = cyl.getCellVolumes()
dt = 20./24./60.
rank = 0
outer_fluxes_solMucil = [-0.1 ]
inner_fluxes_solMucil = [1.]

from rhizo_modelsPlant import RhizoMappedSegments 
import plantbox as pb
perirhizalModel = RhizoMappedSegments(soilModel = cyl, 
                         ms = pb.MappedRootSystem(),
                         limErr1d3dAbs = 1e-8, 
                         RichardsNCCylFoam = Richards4CCylFoam) 

QflowOutCellLim = cyl.distributeSources(dt, source = outer_fluxes_solMucil,
                       inner_fluxs=inner_fluxes_solMucil, 
                       eqIdx =  np.array([1]),plantM=perirhizalModel)
print('QflowOutCellLim',(QflowOutCellLim[0]).sum())
for i in range(1):
    print('cyl.numSoluteComp',cyl.numSoluteComp)
    cyl.ddt = 1e-3
    

    Wvolbefore = cyl.getWaterVolumes() # cellVolumes * cyl.getWaterContent() # cm3
    Smassbefore = cyl.getConcentration(1) * Wvolbefore # mol  cyl.getContent(1) #c

    print('cyl.getConcentration(1)',cyl.getContent(1),cyl.getContent(1).sum(), (cyl.getConcentration(1)  * cyl.getWaterVolumes()).sum())
    print('cyl.getConcentration(2)',cyl.getContent(2))
    cyl.solve(dt)
    print('cyl.getConcentration(1)',cyl.getContent(1),cyl.getContent(1).sum(), (cyl.getConcentration(1)  * cyl.getWaterVolumes()).sum())
    print('cyl.getConcentration(2)',cyl.getContent(2))
    print('cyl.getReactions(1)',cyl.getReactions(1)*cellVolumes*dt,(cyl.getReactions(1)*cellVolumes*dt).sum())
    Wvolafter = cyl.getWaterVolumes() # cellVolumes*s.getWaterContent() # cm3
    Smassafter =  cyl.getConcentration(1) * Wvolafter # mol cyl.getContent(1) #


    rootSoilFluxes = cyl.getInnerFlow(0, length) * dt # cm3
    rootSoilFluxesS = cyl.getInnerFlow(1, length) * dt # mol
    soilSoilFluxes = cyl.getOuterFlow(0, length) * dt # cm3
    soilSoilFluxesS = cyl.getOuterFlow(1, length) * dt # mol
    soilReactions = cyl.getReactions(1)*cellVolumes*dt
    soilSource = (QflowOutCellLim[0]).sum()*dt

    # TODO: currently, setSource not properly implemented for richards and richards 2c
    # so left out of the mass balance.
    # scvSources = cyl.getSource(0) * cellVolumes * dt # cm3
    # scvSourcesS = cyl.getSource(1) * cellVolumes * dt # kg

    if rank == 0:
        print('\tChange in water volume [cm3] per voxel:',sum(Wvolafter-Wvolbefore))
        print('\tChange in solute mass [mol] per voxel:',sum(Smassafter-Smassbefore))
        print('\n\n\tRMSE for water volume balance [cm3]:',np.mean(np.abs(rootSoilFluxes+soilSoilFluxes+sum(Wvolafter-Wvolbefore))))
        print('\tRMSE for solute mass balance [mol]:',rootSoilFluxesS+soilSoilFluxesS+sum(Smassafter-Smassbefore-soilReactions)-soilSource,'\n\n')
        print('QflowOutCellLim',(QflowOutCellLim[0]).sum(), 'dt',(QflowOutCellLim[0]).sum()*dt)


    for ii in range(cyl.numSoluteComp):
        print(cyl.getSolution(ii+1))
    print('getSolutionHead',cyl.getSolutionHead())
    print('\n\n\n')
    