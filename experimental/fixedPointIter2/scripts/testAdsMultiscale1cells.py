import sys;
import os

#os.chdir('experimental/fixedPointIter2/scripts')
sys.path.append("../modules/");
sys.path.append("../inputDataExudate/");
sys.path.append("../../../../CPlantBox/");
sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../../build-cmake/cpp/python_binding/");


from rosi_richards4c_cyl import Richards4CCylFoam as RichardsNCCylFoam  # C++ part (Dumux binding)
#from rosi_richards10c_cyl import Richards10CCylFoam as RichardsNCCylFoam
from rosi_richards4c import Richards4CSPILU as RichardsNCSP 

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
from scipy.interpolate import PchipInterpolator,  CubicSpline
from scipy import sparse
import scipy.sparse.linalg as LA
from plantbox import Perirhizal

perirhizalModel = Perirhizal()

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



def getBiochemParam(s,paramIdx):    
    """ define TraiRhizo biochemical parameters 
        @param: the dumux soil object
        @ param: index of the TraiRhizo parameter set to use
        
    """
    # file containing the TraiRhizo parameter sets
    paramSet = pd.read_csv('./output_random_rows.csv').iloc[paramIdx].to_dict() # select one specific parameter set from index
    s.molarMassC = 12.011
    s.mg_per_molC = s.molarMassC * 1000.
    Ds = paramSet['DS_W'] #cm^2/d
    Dl = 0.003456 #cm^2/d
    s.Ds = Ds /(24*3600) /10000 # m^2/s
    s.Dl = Dl /(24*3600) /10000# m^2/s
    s.BufferPower = BufferPower # random value
    if False:
        kads = 7.07e+02 # m3/kgC/yr, see 10.1016/j.soilbio.2020.107912, A.3
        yr_per_d = 1/365 # [yr/d]
        m3_per_cm3 = 1e-6; # m3/cm3
        cm3_per_m3 = 1e6; # cm3/m3
        
        # [kg/g] * [g/mol] = kg/mol
        kgC_per_mol = (1/1000) * s.molarMassC
        # [m3/kgC/yr] * [yr/d] * [cm3/m3] * [kgC/mol] = [cm3/mol/d]
        s.kads = kads * yr_per_d * cm3_per_m3 * kgC_per_mol
        
        kdes =  1.63e+03 # [1/yr] see 10.1016/j.soilbio.2020.107912, A.3
        s.kdes = kdes * yr_per_d
    s.kads = 10**5 # cm3/mol
    s.kdes = 1.# -
    s.Qmmax = 0.45 * 0.079 # max ratio gOC-gmineral soil, see 10.1016/j.soilbio.2020.107912
    # [g OC / g mineral soil] * [g mineral soil/ cm3 bulk soil] *[ mol C/g C]
    CSSmax_ = s.Qmmax * s.bulkMassDensity_gpercm3*(1/s.molarMassC)
    s.CSSmax = CSSmax_ # mol C/cm3 bulk soil
    #s.CSSmax = s.Qmmax * s.bulkDensity_m3 / 1e6 # mol OC/mol soil * [mol soil/m3] * [m3/cm3] =  mol/cm3

            
    s.css1Function = 5 # current adsorption function implemented.

    return s
    
    
def getSoilTextureAndShape(soil_= "loam"):  
    """ soil shape and texture data
        to adapt according to the soil represented
    """
    min_b = np.array([-2., -1., -1.]) # np.array( [5, 5, 0.] )
    max_b =np.array( [0., 0., 0.]) #  np.array([-5, -5, -5.])
    cell_number = np.array( [1,1,1]) #np.array( [1,1,1]) # 1cm3
    area = 3*3
    solidDensity = 2650 # [kg/m^3 solid] #taken from google docs TraiRhizo
    solidMolarMass = 60.08e-3 # [kg/mol] 
    # theta_r, theta_s, alpha, n, Ks
    
    if soil_ == "loam":
        i = 0
    else:
        i = 1
     
    soilVG = vg_SPP(i)
    
    Kc_value = np.array([1,1,1,1.2,1.2,1.2])
    Kc_days = np.array([1,42,63,98,154,288])
    
    Kc = np.zeros((Kc_days[-1]))
    dummy = 0
    for i in range(0,len(Kc)):
        if i+1 in Kc_days:
            Kc[i] = Kc_value[np.where(Kc_days == (i + 1))[0][0]]
            dummy = dummy+1
        else:
            slope = (Kc_value[dummy]-Kc_value[dummy-1])/(Kc_days[dummy]-Kc_days[dummy-1])
            Kc[i] = Kc_value[dummy-1]+slope*((i+1)-Kc_days[dummy-1])
    
    
    soilTextureAndShape = {'min_b' : min_b,'max_b' : max_b,
                            'area':area,
                           'cell_number':cell_number,
                           "solidDensity":solidDensity,
                        'solidMolarMass': solidMolarMass,
                           'soilVG':soilVG,
                           'Kc':Kc}
    
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
    
    s.setParameter( "Soil.MolarMass", str(s.solidMolarMass))#0.06008
    s.setParameter( "Soil.solidDensity", str(s.solidDensity))
    
    s.setVGParameters([s.soil])
    
    return s

def setBiochemParam(s):
    """ send the TraiRhizo biochemical parameters to dumux
        @param: the dumux soil object
    """

    s.setParameter( "Soil.css1Function", str(s.css1Function))
    
    s.setParameter("1.Component.LiquidDiffusionCoefficient", str(s.Ds)) #m^2/s
    s.setParameter("2.Component.LiquidDiffusionCoefficient", str(s.Dl)) #m^2/s
    s.setParameter("Component.BufferPower", str(s.BufferPower))
    if False:
        
        s.setParameter("Soil.f_sorp", str(s.f_sorp)) #[-]
        
        s.setParameter("Soil.alpha", str(s.alpha)) #[1/d]
    s.setParameter("Soil.kads", str(s.kads)) #[cm3/mol] 
    s.setParameter("Soil.kdes", str(s.kdes)) #[-]
    s.setParameter("Soil.CSSmax", str(s.CSSmax)) #[mol/cm3 scv] 
         
        
    if s.dimWorld == 3:
        # 1 == True
        # if we define a source or sink for the cell 
        # (from the results of the 1d models),
        # do not compute on top of biochemical reactions in dumux
        s.setParameter("Problem.reactionExclusive", "1")
    
    return s                          


def setDefault(s):
    """ Defined some usefull default parameters
    """
    molarMassWat = s.molarMassWat # [g/mol]
    densityWat = 1. #[g/cm3]
    # [mol/cm3] = [g/cm3] /  [g/mol] 
    molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 
    s.molarDensityWat = molarDensityWat

    # low MaxRelativeShift == higher precision in dumux
    s.setParameter("Problem.dobioChemicalReaction",str(s.doBioChemicalReaction))
    s.setParameter("Problem.verbose", "0")
    s.setParameter("Newton.Verbosity", "0") 
    
    # force solute mole fraction > 0 and water in possible pressure ranges
    s.setParameter("Newton.EnableChop", "true")
    
    # UpwindWeight = 1, better when we have high solute gradient.
    # UpwindWeight = 0.5, better when have high water flow and low solute gradient
    s.setParameter("Flux.UpwindWeight", "1")#very important because we get high solute gradient.
    

    s.EnableResidualCriterion = False
    s.setParameter("Newton.EnableResidualCriterion", 
                     str( s.EnableResidualCriterion ))
    s.EnableAbsoluteResidualCriterion = False
    s.setParameter("Newton.EnableAbsoluteResidualCriterion", 
                     str( s.EnableAbsoluteResidualCriterion ))
    s.SatisfyResidualAndShiftCriterion = False
    s.setParameter("Newton.SatisfyResidualAndShiftCriterion",
                     str( s.SatisfyResidualAndShiftCriterion) )  
    s.MaxTimeStepDivisions = 10
    s.setParameter("Newton.MaxTimeStepDivisions",
                     str( s.MaxTimeStepDivisions) )  
    s.MaxSteps = 18
    s.setParameter("Newton.MaxSteps",
                     str( s.MaxSteps) )  
    s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift))
    
    return s
def setIC3D(s, soil_type, ICcc = None):
    return setIC(s, soil_type, ICcc)
def setIC(s, soil_type, ICcc = None):
    """ Defined the initial concentraition of the solutes
        [mol C / cm3 water] for disolved solutes and [mol C / cm3 scv]
        for solutes in the soil phase
        @param: s the dumux soil object
        @param: soil_Type, 0 = loam, 1 = sand
        @param: ICcc (optional) predefined initial conditions
    """
    if ICcc is None:
        #paramSet = pd.read_csv('./output_random_rows.csv').loc[paramIdx]
        #C_S = paramSet['CS_init'] /s.mg_per_molC## small C solutes in mol/cm3 water
        #C_L = paramSet['CL_init'] /s.mg_per_molC## large C solutes in mol/cm3 water
        C_S = 0
        C_L = 0

        # concentraiton of adsobed C_S
        #s.CSS_init  = getCSS(s, C_S) #mol C/ cm3 scv
        
            
        #unitConversion = 1.0e6 # mol/cm3  => mol/m3 
        #addedVar = 1. * float(s.doSoluteFlow) # empirical factor
        #s.CSW_init = C_S * unitConversion
        s.ICcc = np.array([0.,#C_S *unitConversion*addedVar,
                           0,# s.CSS_init*unitConversion*addedVar,
                           ])# in mol/m3 water or mol/m3 scv
    else:
        s.ICcc = ICcc
    
    for i in range(s.numSoluteComp):
        #mol/m3 to mol/mol
        molarC = s.ICcc[i] #/ s.phaseDensity(isDissolved = (i < s.numDissolvedSoluteComp)) 
        s.setParameter( "Soil.IC.C"+str(i+1), str(molarC ))# send data to dumux
    return s


def setupOther(s, soil_type, simMax):
    """ define remaining soil parameters """ 
    
    # climate data 
    soilTextureAndShape = getSoilTextureAndShape(soil_type) 
    cell_number = soilTextureAndShape['cell_number']
    times, net_inf = None, None #evap.net_infiltration(soil_type, simMax, soilTextureAndShape['Kc'])
    
    s.setParameter("SpatialParams.Temperature","293.15") # todo: redefine at each time step?
    
    # BC
    p_mean_ = initP
    if s.dimWorld == 1: # 1d model
        s.setParameter("Soil.IC.P", s.dumux_str(p_mean_))
        s.setInnerBC("fluxCyl",  s.win)  # [cm/day]
        s.setOuterBC("fluxCyl", 0.)
    if s.dimWorld == 3:# 3d model
        s.setParameter("Soil.IC.P", s.dumux_str(p_mean_))
        if times is not None:
            s.setTopBC("atmospheric", 0.5, [times, net_inf])  # 0.5 is dummy value
        else:
            s.setTopBC("noFlux")
        s.setBotBC("freeDrainage")#"noFlux")#
    
        for i in range(1, s.numComp):# no flux
            s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(2))
            s.setParameter( "Soil.BC.Top.C"+str(i)+"Type", str(2))
            s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(0)) 
            s.setParameter( "Soil.BC.Top.C"+str(i)+"Value", str(0 ))       
        
    
    # IC  
    s.maxDt =  250./(3600.*24.)
    s.maxDt_1DS = s.maxDt # [s], lower maxDt for 1D models
    s.initializeProblem(s.maxDt)

    s.eps_regularization = 1e-14
    s.setRegularisation(s.eps_regularization, s.eps_regularization)
    #df = pd.read_csv("../inputDataExudate/data/init_pot_2019.csv")  # initial potential
    #h  = np.flip(df[soil_type].loc[:].values) #cm
    #h = np.repeat(h[:,np.newaxis],cell_number[0],axis=1) #x-axis
    #h = np.repeat(h[:,:,np.newaxis],cell_number[1],axis=2) #y-axis
    #h = h.flatten()
    #h = np.ones((20*44*75))*-100 #TODO
    #s.setInitialConditionHead(h)  # cm
    
    # for boundary conditions constantFlow, constantFlowCyl, and atmospheric
    s.wilting_point = -15000
    s.setCriticalPressure(s.wilting_point)  
    s.ddt = 1.e-3  # [day] initial Dumux time step
    s.bulkMassErrorWater_rel = 0.
    s.bulkMassErrorWater_relLim = 0.    
    s.totC3dInit = sum(s.getTotCContent()) # mol    
    # initial soil water and solute content
    cell_volumes = s.getCellVolumes()  # cm3
    s.buWSoilInit = sum(np.multiply(np.array(s.getWaterContent()), cell_volumes)) # cm3 water
    return s

def vg_SPP(i = int(1)):
    """ Van Genuchten parameter, called by maize()  """
        
    soil = {}
    # theta_r, theta_s, alpha, n, Ks
    #soil[0] = [0.08, 0.43, 0.04, 1.6, 50] #Mona
    soil[0] = [0.041, 0.494, 0.0256, 1.49, 245]
    soil[1] = [0.03, 0.414, 0.038, 2, 1864]
    return soil[i]
    
def create_soil_model(simMax, usemoles, results_dir , soil_='loam',
                     noAds = False, ICcc = None, doSoluteFlow = True,
                     doBioChemicalReaction=True, 
                     MaxRelativeShift = 1e-8):
    """
        Creates a soil domain from @param min_b to @param max_b with resolution @param cell_number
        homogeneous domain 
        initial potentials are linear and mean potential is @p_mean_
        @ param: noAds: turn off adsorption?
        @param: paramIndx index of the TraiRhizo parameter set to use
        @param: ICcc (optional) initial concentraiton values for the solute components
        @param: usemoles [bool] dumux uses moles (True) or grammes (False)
        returns soil_model (RichardsWrapper(RichardsSP())) and soil parameter (vg.Parameters)
    """
    if soil_ == "loam":
        soil_type = 0
    else:
        soil_type = 1
        
    s = RichardsWrapper(RichardsNCSP(), usemoles)  # water and N solute          
    s.results_dir = results_dir   
    s.pindx = soil_type
    
    s.MaxRelativeShift = MaxRelativeShift # 1e-10
    s.MaxRelativeShift_1DS = MaxRelativeShift
    
    soilTextureAndShape = getSoilTextureAndShape(soil_type) 
    min_b = soilTextureAndShape['min_b']
    max_b = soilTextureAndShape['max_b']
    cell_number = soilTextureAndShape['cell_number']
    s.cell_size = np.prod((max_b - min_b) / cell_number) # cm3 
    s.setParameter( "Soil.Grid.Cells", s.dumux_str(cell_number))  # send data to dumux
    s.noAds = noAds # no adsorption?
    s.doSoluteFlow = doSoluteFlow
    s.doBioChemicalReaction = doBioChemicalReaction
    
    s.setParameter("Newton.Verbosity", "0") 
    s.initialize() 
    setDefault(s)
    
    setSoilParam(s)
    getBiochemParam(s,soil_type)
    setBiochemParam(s)
    setIC3D(s, soil_type, ICcc)
    s.isPeriodic = True
    s.createGrid(min_b, max_b, cell_number, s.isPeriodic)  # [cm] 
    s = setupOther(s, soil_, simMax)
    
    return s
    
    
def initialize_dumux_nc_(  gId=0, a_in=0.02,
                a_out=0.03 ,seg_length=1.0 ,
                x=[-10000] ,   # cm
                c2 = 0.06,             # mol/mol scv
                NC = 10,logbase = 0.5, Cells = []):                                   # cm
    verbose = False
    lId =gId
    
    if a_in < a_out:
    
        cyl = RichardsNoMPIWrapper(RichardsNCCylFoam(), True)  # only works for RichardsCylFoam compiled without MPI
        if False:
            cyl.setParameter("Newton.EnableResidualCriterion", "true") # sometimes helps, sometimes makes things worse
            cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
            cyl.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
            cyl.setParameter("Newton.MaxRelativeShift", "1e-9")# reset value
        cyl.initialize(verbose = False) # No parameter file found. Continuing without parameter file.

        lb =  logbase
    
        points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), 
                             NC, base = lb)
        
        cyl.createGrid1d(points)# cm
        setSoilParam(cyl)
        cyl.doSoluteFlow = True
        cyl.noAds = False
        cyl = getBiochemParam(cyl,61)
        cyl = setBiochemParam(cyl)

        cyl.setParameter("Flux.UpwindWeight", "1")
        cyl.setParameter("SpatialParams.Temperature","293.15") # todo: redefine at each time step
        cyl.setParameter("Soil.BC.dzScaling", "1")
        cyl.setParameter( "Soil.css1Function", str(5))
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
        cyl.setInnerBC("fluxCyl", exW)  # [cm/day] #Y 0 pressure?
        cyl.setOuterBC("fluxCyl", 0.)
        #cyl.setParameter("2.Component.LiquidDiffusionCoefficient", str( 0.003456/(24*3600) /10000)) #m^2/s
        cyl.setParameter("Newton.MaxTimeStepDivisions",
                     str( 100) )
        cyl.setParameter("Newton.MaxSteps",
                     str( 100) )
            
        
        for j in range( 1, cyl.numComp):
            #cyl.setParameter(str(j)+".Component.LiquidDiffusionCoefficient", str( 0.003456/(24*3600) /10000)) #m^2/s
            cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Type", str(3))
            cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Type", str(3))
            cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Value", str(0)) 
            cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Value", str(0 ))
        #cyl.setParameter( "Soil.BC.Bot.C2Value", str(ex2)) # mol/cm2/d
        cyl.setParameter( "Soil.BC.Bot.C1Value", str(ex2)) # mol/cm2/d
        # mol/cm2/day
        cyl.setParameter("Soil.IC.C1",cyl.dumux_str(c2) ) 
        #cyl.setParameter("Soil.IC.C2",cyl.dumux_str(c2) ) 
    
        if len(Cells) > 0:#in case we update a cylinder, to allow dumux to do interpolation
            assert(len(c2)==len(Cells))
            CellsStr = cyl.dumux_str(Cells/100)#cm -> m
            cyl.setParameter("Soil.IC.Z",CellsStr)# m
            if len(Cells)!= len(x):
                print("Cells, x",Cells, x, len(Cells), len(x))
                raise Exception
                
            j = 1
            cyl.setParameter("Soil.IC.C"+str(j)+"Z",CellsStr) # m
            if len(Cells)!= len( c2):
                print("Cells,  cAll[j-1]",Cells,  c2, 
                        len(Cells), len(c2), j)
                raise Exception
                        
        #print('str(c2)',str(c2))


        cyl.maxDt = 250/(3600*24) # soilModel.maxDt_1DS
        cyl.setParameter("Soil.mucilCAffectsW", "false")
        cyl.setParameter("Problem.verbose", "0")
        cyl.setParameter("Newton.Verbosity", "0")
        cyl.initializeProblem(maxDt=cyl.maxDt ) # ewton solver configured with the following options and parameters:

        
        cyl.eps_regularization = 1e-14
        cyl.setRegularisation(cyl.eps_regularization, cyl.eps_regularization)

        cyl.setCriticalPressure(-15000)  # cm pressure head
        
        #cyl.setRegularisation(1e-16, 1e-16)
        
        return cyl
    else:
        print("RhizoMappedSegments.initialize_dumux_: Warning, segment {:g} might not be in domain, radii [{:g}, {:g}] cm".format(i, a_in, a_out))
        return []
    
def storeWData1d():#cm3    
    storeMassData1_ = np.array([0. for nn in range(s3d.numberOfCellsTot)])
    storeMassData1_[0] = np.array([ cyl.getWaterVolumes().sum() for cyl in cyls ]).sum(axis=0)
    return storeMassData1_
    
def storeMassData1d():#mol    
    #print('np.array([ cyl.getTotCContent_each().sum(axis=1) for cyl in cyls ])',
    #np.array([ cyl.getTotCContent_each().sum(axis=1) for cyl in cyls ]), 
    #'cyl.getCss1()',cyl1.getCss1())
    storeMassData1_ = np.array([np.array([0.,0.]) for nn in range(s3d.numberOfCellsTot)])
    storeMassData1_[0] =  np.array([ cyl.getTotCContent_each().sum(axis=1) for cyl in cyls ]).sum(axis=0)
    return storeMassData1_
                
    
def compute1dChangesSolute():#
    sources_Wsol_from1d = np.full( (s3d.numComp, s3d.numberOfCellsTot),0. )
    sources_Wsol_from1d[0,cellIds] = np.array(
            soil_W1d_perVoxelAfter - soil_W1d_perVoxelBefore - outer_R_bc_sol[0]
        )[cellIds]/dt
    res = {i: sources_Wsol_from1d[0][i] for i in range(len(sources_Wsol_from1d[0]))} 
    #print('compute1dChangesSolutesoil_solute1d_perVoxelAfter',soil_solute1d_perVoxelAfter, soil_solute1d_perVoxelBefore,outer_R_bc_sol)  
    #print('setSource',0, res)
    s3d.setSource(res.copy(), eq_idx = 0)  # [mol/day], in 
    for nc in range(1,s3d.numComp):
        sources_Wsol_from1d[nc,cellIds] = np.array(
            soil_solute1d_perVoxelAfter[cellIds][nc-1] - soil_solute1d_perVoxelBefore[cellIds][nc-1] - outer_R_bc_sol[nc][cellIds]
        )/dt
        #print('compute1dChangesSolute',nc,sources_Wsol_from1d[nc,cellIds],
        #        soil_solute1d_perVoxelAfter[cellIds][nc-1], soil_solute1d_perVoxelBefore[cellIds][nc-1] , outer_R_bc_sol[nc][cellIds])
        #sendSource2dumux(sources_sol_from1d[nc], nc + 1)
        res = {i: sources_Wsol_from1d[nc][i] for i in range(len(sources_Wsol_from1d[nc]))}                  
       
        # send to dumux
        #print('setSource',nc, res)
        s3d.setSource(res.copy(), eq_idx = nc)  # [mol/day], in 
    return sources_Wsol_from1d
    
def getTBalance3d(dt):
    # check that source + flow == dS for each 1d and 3d cell      
    scvFluxes = s3d.getFluxesPerCell(0) * dt  # cm3
    scvFluxesS = s3d.getFluxesPerCell(1) * dt  # kg
    scvFluxesS2 = s3d.getFluxesPerCell(2) * dt  # kg
    scvSources = s3d.getSource(0) * s3d.getCellVolumes() * dt # cm3
    scvSourcesS = s3d.getSource(1) * s3d.getCellVolumes() * dt # kg
    scvSourcesS2 = s3d.getSource(2) * s3d.getCellVolumes() * dt # kg
    print('\tChange in water volume [cm3] per voxel:',(Wvolafter-Wvolbefore).sum())#,Wvolafter,Wvolbefore)
    print('\tChange in solute mass [mol] per voxel:',(Smassafter-Smassbefore).sum())#, Smassafter,Smassbefore)
    print('\tRMSE for water volume balance [cm3]:',np.mean(np.abs(scvFluxes+(Wvolafter-Wvolbefore)-scvSources)))
    print('\tRMSE for water volume balance cell0 [cm3]:',scvFluxes[0]+(Wvolafter[0]-Wvolbefore[0])-scvSources[0])
    print('\tRMSE for solute mass balance [mol]:',np.mean(np.abs(scvFluxesS+scvFluxesS2+(Smassafter.sum(0)-Smassbefore.sum(0))-scvSourcesS-scvSourcesS2)))
    #print('getSolution',s3d.getSolution(1), s3d.getSolution(2), s3d.getCss1() )
    #print('scvFluxes,Wvolafter,Wvolbefore,scvSources,',scvFluxes,Wvolafter,Wvolbefore,scvSources,'\n\n')
    #print('Smassafter',Smassafter, Smassbefore,Smassafter.sum(),Smassafter.sum(0),Smassafter.sum(1))
    #print('scvSourceS',  scvSourcesS, 'scvSourcesS2',scvSourcesS2)
    print('scvFluxes,',scvFluxesS,scvFluxesS2,Smassafter, Smassbefore,scvSourcesS,scvSourcesS2,'\n\n')
        
def getTBalance1d(dt,cyl,Wvolafter,Wvolbefore,Smassafter,Smassbefore):
    # check that source + flow == dS for each 1d domaine  
    rootSoilFluxes = cyl.getInnerFlow(0, cyl.seg_length) * dt # cm3
    rootSoilFluxesS = cyl.getInnerFlow(1, cyl.seg_length) * dt # kg
    soilSoilFluxes = cyl.getOuterFlow(0, cyl.seg_length) * dt # cm3
    soilSoilFluxesS = cyl.getOuterFlow(1, cyl.seg_length) * dt # kg
    
    print('\tChange in water volume [cm3] 1D:',(Wvolafter-Wvolbefore))
    print('\tChange in solute mass [mol] 1D:',(Smassafter-Smassbefore))
    print('\tRMSE for water volume balance [cm3] 1D:',np.mean(np.abs(rootSoilFluxes+soilSoilFluxes+(Wvolafter-Wvolbefore))))
    print('\tRMSE for solute mass balance [mol] 1D:',rootSoilFluxesS[0]+soilSoilFluxesS[0]+(Smassafter-Smassbefore))
    print('content',cyl.getContent(1), cyl.getCss1() * cyl.getCellVolumes_() )
    print('concentration', cyl.getConcentration(1),cyl.getCss1(), getCSS(cyl, cyl.getConcentration(1)))
    print('css1 data', cyl.kads, cyl.kdes,cyl.CSSmax)
    print(rootSoilFluxesS,'rootSoilFluxesS',rootSoilFluxesS/dt/cyl.getFaceSurfaces_(cyl.seg_length)[0],
                soilSoilFluxesS,Smassafter,Smassbefore)
    print()
    print()
    
        
def get1d3dBalance(dt):
    scvFluxes = s3d.getFluxesPerCell(0) * dt  # cm3
    print('\tRMSE Delta 1d3d cell0:',scvFluxes[0]+(Wvolafter[0]-Wvolbefore[0])- ( sum(cylsWVolAfter) - sum(cylsWVolBefore)))
    print('\tRMSE 1d3d cell0:',scvFluxes[0]+(Wvolafter[0])- ( sum(cylsWVolAfter)))
    
    # for each 3d cell sum(1d) + Q_3D == 3d

def getVolumes(vertices, length):
    """ volume of each cell of a hollow cylinder mesh """
    return   np.pi * (vertices[1:] ** 2 - vertices[:-1]**2) * length

def interAndExtraPolation_(pt,pointsOld_, chip, spl):
    if (pt >= min(pointsOld_)) and (pt <= max(pointsOld_)): #interpolation
        return chip(pt)
    else:#extrapolation
        return spl(pt)
    
def interAndExtraPolation(pointsNew,pointsOld, valOld):
    chip = PchipInterpolator(pointsOld, valOld)#for interpolation
    spl =  CubicSpline(pointsOld, valOld, bc_type='not-a-knot') #extrapolation  
    return np.array([interAndExtraPolation_(pt, pointsOld, chip, spl) for pt in pointsNew])

def update_concentration( totContent, changeRatio, gradient, phaseVolOrMolFrOldold, 
                volumes, isWater, verbose=False):
    """
    Update the concentration to get specific total content and gradient.
    @param totContent: new total content of the element in the cylinder (cm3 water or mol solute)
    @param changeRatio: ratio between toContent and old total content of the element (for trouble shooting)
    @param gradient: old concentration gradient
    @param phaseVolOrMolFrOldold: old water content (cm3/cm3) or volume of bulk soil
    @param volumes: volume of each cell of the cylinder (cm3)
    @param isWater : water element (True) or solute (False)
    @param verbose
    """

    def print_verbose(*args):
        if verbose:
            print(rank, *args)

    def solve_concentration_matrix(matrix, aB):
        SPmatrix = sparse.csc_matrix(sparse.coo_matrix(matrix))
        return LA.spsolve(SPmatrix, aB, use_umfpack=True)

    def validate_concentration(val_new, totContent, gradient, volumes):
        """ check that tot content and gradient is correct """
        try:
            assert abs(sum(val_new * volumes) - totContent) < 1e-13
            assert (abs((val_new[1:] - val_new[:-1]) - gradient) < 1e-13).all()
        except:
            print_verbose('update_concentration error', 'val_new', val_new, 'totContent', totContent, 'gradient', gradient, 
            'phaseVolOrMolFrOldold', phaseVolOrMolFrOldold, 'volumes', volumes, 'changeRatio', changeRatio)
            raise Exception

    # Verbose initial information
    print_verbose('isWater', isWater, 'update_concentration error', 'totContent', totContent, 'gradient', gradient, 'phaseVolOrMolFrOldold', phaseVolOrMolFrOldold, 'volumes', volumes, 'changeRatio', changeRatio)

    # Create and solve concentration matrix
    matrix_size = NC - 1
    matrix = np.diag(np.full(matrix_size - 1, -1.), k=1) + np.diag(np.full(matrix_size, 1.), k=0)
    matrix[-1, :] = volumes
    aB = np.append(-gradient, totContent)
    val_new = solve_concentration_matrix(matrix, aB)
    
    print_verbose('val_new', val_new)
    validate_concentration(val_new, totContent, gradient, volumes)

    # Set bounds based on water or solute
    if isWater:
        maxVal = s3d.vg_soil.theta_S
        minVal = s3d.theta_wilting_point
    else:
        maxVal = np.inf
        minVal = 0.
    # check that the mean concentration respects the boundary 
    assert ( (minVal - totContent/sum(volumes)) < 1e-14) 
    assert ( (totContent/sum(volumes) - maxVal) < 1e-14) 
        

    # Verbose before changes
    print_verbose('BEFORE possible changes.val_new',val_new, 'new content', val_new * volumes, 'totContent', totContent, 'sum new content', sum(val_new * volumes), 'sum(val_new * volumes) - totContent', sum(val_new * volumes) - totContent, 'maxVal',maxVal, 'minVal',minVal)

    # Adapt values if necessary
    try:
        val_new = np.array(perirhizalModel.adapt_values(val_new_ = val_new, minVal_ = minVal, maxVal_=maxVal, volumes_=volumes, 
                                             divideEqually_=True, verbose_ = verbose))
    except:
        print("issue adapt_values val_new, minVal, maxVal, volumes",
              val_new, minVal, maxVal, volumes)
        raise Exception
    assert (val_new >= minVal).all()
    assert (val_new <= maxVal).all()
    assert abs(sum(val_new * volumes) - totContent) < 1e-14

    return val_new
    
def getCSS(s, CSW):
    """ @return concentration of adsobed carbon in the soil
        according to @param CSW [mol C / cm3 water], the concentration
        of small solutes in the soil water
        @param: s the dumux soil object
    """
    return  (s.kads * CSW * s.CSSmax)/(s.kads * CSW + s.kdes) #kd*CSW
    
def getCSWfromC_total(s, C_total):
    """
    Compute the dissolved concentration C_SW (mol/cm³ water)
    from total concentration C_total (mol/cm³ soil),
    using equilibrium adsorption.

    Parameters:
    - C_total: float, total concentration [mol/cm³ soil]
    - kads: float, adsorption rate constant [cm³ water / mol]
    - kdes: float, desorption rate constant [dimensionless]
    - CSSmax: float, max sorption site capacity [mol/cm³ soil]

    Returns:
    - C_SW: float, dissolved concentration [mol/cm³ water]
    """
    a = s.kads
    d = s.kdes
    Cmax = s.CSSmax
    Ct = C_total

    # Coefficients of the quadratic equation: a*x^2 + B*x - Ct*d = 0
    B = -Ct * a + d + a * Cmax

    discriminant = B**2 + 4 * a * Ct * d
    if discriminant < 0:
        raise ValueError("getCSWfromC_total: No real solution exists for the given parameters.")

    # Only the positive root makes physical sense
    C_SW = (-B + math.sqrt(discriminant)) / (2 * a)
    return C_SW
    
def cylinderGrowth():
    ls = [1.,3.]
    ratio1 = ls[0]/sum(ls)
    rOuts = [(Vc1*ratio1/(np.pi*ls[0])+rIns[0]*rIns[0])**(1/2),(Vc1*(1-ratio1)/(np.pi*ls[1])+rIns[1]*rIns[1])**(1/2)]
    
    thetaLeftOver = 0.
    print('Vol Old', sum(cyls[0].getCellVolumes()),sum(cyls[1].getCellVolumes()),
            sum([sum(cyls[0].getCellVolumes()),sum(cyls[1].getCellVolumes())]))
    print('Wvol Old', sum(cyls[0].getWaterVolumes()),sum(cyls[1].getWaterVolumes()),
            sum([sum(cyls[0].getWaterVolumes()),sum(cyls[1].getWaterVolumes())]))
            
    for gId, cyl in enumerate(cyls):
        oldPoints = np.array(cyl.getPoints()).flatten()
        points = np.logspace(np.log(rIns[gId]) / np.log(lb), np.log(rOuts[gId]) / np.log(lb), 10, base = lb) # cm

        volOld = cyl.getCellVolumes()
        ## new shape: 
        volNew = getVolumes(points, ls[gId] ) # cm^3
        deltaVol = sum(volNew) - sum(volOld) 
        
        ## shape:
        centersNew = (points[1:] + points[:-1]) / 2  # cm

        changeRatio = min(sum(volNew)/sum(volOld), 1.)

        ##  water:
        theta_old = cyl.getWaterContent() # cm3/cm3
        gradientOld = (theta_old[1:] - theta_old[:-1])
        gradientNew = interAndExtraPolation(points[1:-1],oldPoints[1:-1], gradientOld)

        wOld = sum(theta_old*volOld)
        maxVal = cyl.vg_soil.theta_S  # keep that or set lower higher bound?
        minVal = s3d.theta_wilting_point

        if deltaVol < 0.:
            changeRatioW = max(min(changeRatio, sum(maxVal*volNew)/wOld),sum(minVal*volNew)/wOld)
        else:
            changeRatioW = changeRatio
    
        theta_new = update_concentration(totContent = wOld*changeRatioW + max(thetaLeftOver*deltaVol,0.),
                                                             changeRatio=changeRatioW, 
                                                         gradient =gradientNew, phaseVolOrMolFrOldold = theta_old, volumes = volNew,isWater = True)
        newHead = np.array([vg.pressure_head(nt, cyl.vg_soil) for nt in theta_new])
        if  deltaVol < 0.:                                                   
            thetaLeftOver = sum(theta_old * volOld - theta_new * volNew)/(-deltaVol)
            

        ## new contents:  
        molKonzOld =np.array( [np.array(cyl.getConcentration(nC+1)) for nC in range(s3d.numSoluteComp)])   #mol/mol 
        volWatNew = theta_new *volNew
        molKonzNew = []
        molFrNew = []
        for nComp in range(1, s3d.numComp):
            if (molKonzOld[nComp -1] != 0.).any():
                isDissolved = (nComp <= s3d.numDissolvedSoluteComp)

                if isDissolved: # cm3 phase
                    molKonzOld[nComp -1] *= theta_old # to go from C_S^l in mol/cm3 wat to C_S^l + C_S^s in mol/cm3 scv
                    if nComp == 1:
                        molKonzOld[nComp -1] += cyl.getCss1() # add linear adsorption

                gradientOld = (molKonzOld[nComp -1][1:] - molKonzOld[nComp -1][:-1])   
                gradientNew = interAndExtraPolation(points[1:-1],oldPoints[1:-1], gradientOld)
                cOld = sum(molKonzOld[nComp -1] * volOld  ) 
                
                if nComp == 1:
                    assert abs(cOld - sum(cyl.getCSS()*volOld) - sum(cyl.getContent(nComp) ))< 1e-13
                else:
                    assert abs(cOld - sum(cyl.getContent(nComp)))< 1e-13
                
                #get mole fraction for each cell
                molKonzNew.append(
                    update_concentration(totContent = cOld*changeRatio+ max(konzLeftOver[nComp -1]*deltaVol,0.),
                                                         changeRatio=changeRatio,gradient =gradientNew, 
                        phaseVolOrMolFrOldold = molFrOld[nComp -1], volumes = volNew,isWater = False, verbose = False))
                if isDissolved:
                    molKonzNew[0] = getCSWfromC_total(molKonzNew[0])# from C_S^l + C_S^s in mol/cm3 scv to C_S^l in mol/cm3 wat
            else:
                molKonzNew.append(molKonzOld[nComp -1])
                
            molFrNew.append( molKonzNew[nComp -1] * s3d.phaseDensity(isDissolved)/1e6 )# go from concentration to mol fraction
            
    
        cyls[gId] = initialize_dumux_nc_( gId, seg_length=ls[gId],
                                                    x = newHead,# cm
                                                    a_in=rIns[gId],a_out=rOuts[gId],
                                                    c2 = molFrNew, # mol/mol water or mol/mol scv
                                                    Cells = centersNew)
                                                    
    print('Vol New', sum(cyls[0].getCellVolumes()),sum(cyls[1].getCellVolumes()),
            sum([sum(cyls[0].getCellVolumes()),sum(cyls[1].getCellVolumes())]))
    print('Wvol New', sum(cyls[0].getWaterVolumes()),sum(cyls[1].getWaterVolumes()),
            sum([sum(cyls[0].getWaterVolumes()),sum(cyls[1].getWaterVolumes())]))
    print()
    raise Exception
      
molMarssW = 18 # g/mol
# molMassMulC = 30 #g/mol
rhoW = 1 #g/cm3
rhoWM = rhoW/molMarssW #[g/cm3]*[g/mol] = mol/cm3
cll_i = 1e-5 *0#*30# *0# mol/cm3
mlFr = cll_i/rhoWM #cw / molMassMulC * molMarssW
initP = -10000.
BufferPower = 100.
exW = -0.1
s3d = create_soil_model(10., True, '.' , soil_='loam',
                     noAds = False, ICcc = [mlFr, 0.], doSoluteFlow = True,
                     doBioChemicalReaction=True, 
                     MaxRelativeShift = 1e-8)

s3d.theta_wilting_point = vg.water_content( s3d.wilting_point, s3d.vg_soil)
cellIds = 0
Vc1 = s3d.getCellVolumes()[cellIds]
rIns = [0.01, 0.02]
ls = [1.,1.]
ratio1 = ls[0]/sum(ls)
rOuts = [(Vc1*ratio1/(np.pi*ls[0])+rIns[0]*rIns[0])**(1/2),(Vc1*(1-ratio1)/(np.pi*ls[1])+rIns[1]*rIns[1])**(1/2)]
ex2 = 1.#e-4#/2
NC = 2        
lb = 0.5
#points = np.array([a_in + (a_out - a_in)*i/(NC-1) for i in range(NC)])

# CC = np.array([(points[i] + points[i+1])/2 for i in range(len(points)-1)])#/100
cyls = []
for nn in range(2):
    points = np.logspace(np.log(rIns[nn]) / np.log(lb), np.log(rOuts[nn]) / np.log(lb), 
                              NC, base = lb)
    cells = np.array([cc/2 for cc in np.diff(points)])+ points[:-1]
    
    cyls.append(initialize_dumux_nc_(x=[initP for i in range(NC-1)],
                            a_in=rIns[nn],a_out=rOuts[nn],seg_length=ls[nn],
                            c2 = [mlFr for i in range(NC-1)],
                            NC = NC, Cells =cells))

timeT = 0.


# Q = passage 1d to 3d ok?
# if buffer only in some cell, ok?
dt = 3600./24./3600.#10./24./60.

cylsWVolAfter = [cyl.getWaterVolumes().sum() for cyl in cyls]
cylsSAfter = [cyl.getTotCContent_each().sum() for cyl in cyls]
cylsVol = [cyl.getCellVolumes().sum() for cyl in cyls]

soil_W1d_perVoxelAfter = storeWData1d()
soil_solute1d_perVoxelAfter = storeMassData1d()
#print('soil_solute1d_perVoxelAfter',
#            soil_solute1d_perVoxelAfter)

outer_R_bc_sol = np.zeros((s3d.numComp, len(s3d.getSolution(1))))
Wvolafter = s3d.getWaterVolumes()#.sum()
Smassafter = s3d.getTotCContent_each()#.sum()

print(s3d.getSolutionHead_(),cyls[0].getSolutionHead_())
print('\tRMSE 1d3d cell0:',(Wvolafter[0]), ( sum(cylsWVolAfter)))
print('\tRMSE 1d3d cell0:',(s3d.getCellVolumes()[0]),( sum(cylsVol)))
print('\tRMSE 1d3d cell0S:',(Smassafter[0]), ( sum(cylsSAfter)))

for i in range(2):
    cylsWVolBefore = cylsWVolAfter
    cylsSBefore = cylsSAfter
    
    Wvolbefore = Wvolafter
    Smassbefore = Smassafter
    soil_W1d_perVoxelBefore = soil_W1d_perVoxelAfter
    soil_solute1d_perVoxelBefore = soil_solute1d_perVoxelAfter
    #print('soil_solute1d_perVoxelAfter, soil_solute1d_perVoxelBefore',
    #        soil_solute1d_perVoxelAfter, soil_solute1d_perVoxelBefore)
    for cyl in cyls:
        cyl.ddt = 1e-3
        cyl.solve(dt)
    cylsWVolAfter = [cyl.getWaterVolumes().sum() for cyl in cyls]
    cylsSAfter = [cyl.getTotCContent_each().sum() for cyl in cyls]    

    print('\n\n')
    for nn, cyl in enumerate(cyls):
        getTBalance1d(dt,cyl,cylsWVolAfter[nn],cylsWVolBefore[nn],cylsSAfter[nn],cylsSBefore[nn])
        raise Exception
    print('\n\n')
    
    
    soil_W1d_perVoxelAfter = storeWData1d()
    soil_solute1d_perVoxelAfter = storeMassData1d()
    
    # get cm3 and mol changes in the 1DS
    sources_Wsol_from1d = compute1dChangesSolute()
    #print('soil_solute1d_perVoxelAfter, soil_solute1d_perVoxelBefore',
    #        soil_solute1d_perVoxelAfter, soil_solute1d_perVoxelBefore)
    #print('sources_Wsol_from1d',sources_Wsol_from1d)
    
    s3d.ddt = 1e-3
    s3d.solve(dt)
    Wvolafter = s3d.getWaterVolumes()#.sum()
    Smassafter = s3d.getTotCContent_each()#.sum()
    #print('Smassbefore, Smassafter',Smassbefore, Smassafter)
    
    
    timeT += dt
    outer_R_bc_sol = np.array([s3d.getFluxesPerCell( nc) *dt for nc in range(s3d.numComp)])
    
    getTBalance3d(dt)
    print()
    raise Exception
    get1d3dBalance(dt)
    #fillOut(i+1)
    print('timeT',timeT)
    
    cylinderGrowth()
    raise Exception
    
#doPrints(cyl,a_in,timeT,0)
#cc = (cyl1.getCellCenters_().reshape(-1))#/100.
getPressureHead = [[],[],[]]
CsolutionsWatc2 = [[[] for i in range(3)] for j in range(cyl2.numSoluteComp)]
CsolutionsWatc1 = [[[] for i in range(3)] for j in range(cyl1.numSoluteComp)]
CsolutionsWatc3 = [[[] for i in range(3)] for j in range(s3d.numSoluteComp)]
theta1 = [[] for i in range(3)]
theta2 = [[] for i in range(3)]
theta3 = [[] for i in range(3)]

def fillOut(i):
    
    for j in range(cyl1.numDissolvedSoluteComp):
        CsolutionsWatc1[j][i] = cyl1.getSolution(j+1)*rhoWM 
    for j in range(cyl1.numDissolvedSoluteComp,cyl1.numSoluteComp):
        CsolutionsWatc1[j][i] = cyl1.getSolution(j+1)*cyl1.solidMolDensity /1e6 
    theta1[i] = cyl1.getWaterContent()

    for j in range(cyl2.numDissolvedSoluteComp):
        CsolutionsWatc2[j][i] = cyl2.getSolution(j+1)*rhoWM 
    for j in range(cyl2.numDissolvedSoluteComp,cyl2.numSoluteComp):
        CsolutionsWatc2[j][i] = cyl2.getSolution(j+1)*cyl2.solidMolDensity /1e6 
    theta2[i] = cyl2.getWaterContent()

    for j in range(s3d.numDissolvedSoluteComp):
        CsolutionsWatc3[j][i] = s3d.getSolution(j+1)*rhoWM 
    for j in range(s3d.numDissolvedSoluteComp,s3d.numSoluteComp):
        CsolutionsWatc3[j][i] = s3d.getSolution(j+1)*s3d.solidMolDensity /1e6 
    theta3[i] = s3d.getWaterContent()
    
outputs = {'cc':cc,'Cmucil':CsolutionsWat[0],'theta':theta}
import pickle

with open("results/dataDumuxAds.pkl", "wb") as file:
    pickle.dump(outputs, file)
    
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 8), sharex=True)

# Second subplot (top-right): Cmucil
for i in range(3):
    axes[0].plot(cc, CsolutionsWat[0][i], label=str(dt * i))
axes[0].grid(True)
axes[0].legend()
axes[0].ticklabel_format(style='plain', useOffset=False)
axes[0].set_title("Cmucil")

for i in range(3):
    axes[1].plot(cc, theta[i], label=str(dt * i))
axes[1].grid(True)
axes[1].legend()
axes[1].ticklabel_format(style='plain', useOffset=False)
axes[1].set_title("theta")


# Adjust layout
plt.grid(True)
plt.legend()
plt.tight_layout()
#plt.show()
plt.savefig("adsorption.png", dpi=300, bbox_inches='tight')#, transparent=True)
plt.close()
