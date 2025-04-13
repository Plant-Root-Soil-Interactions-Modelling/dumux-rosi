import sys;
import os

#os.chdir('experimental/fixedPointIter2/scripts')
sys.path.append("../modules/");
sys.path.append("../inputDataTraiRhizo/");
sys.path.append("../../../../CPlantBox/");
sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../../build-cmake/cpp/python_binding/");


from rosi_richards22c_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)

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



def getBiochemParam(s,paramIdx):    
    """ define TraiRhizo biochemical parameters 
        @param: the dumux soil object
        @ param: index of the TraiRhizo parameter set to use
        
    """
    # file containing the TraiRhizo parameter sets
    paramSet = pd.read_csv('./output_random_rows.csv').iloc[paramIdx].to_dict() # select one specific parameter set from index
    s.molarMassC = 12.011
    s.mg_per_molC = s.molarMassC * 1000.
    s.betaC = paramSet['beta_C'] # -
    s.betaO = paramSet['beta_O'] # -

    ## ATT! C_S_W_thres has to be > 0
    s.C_S_W_thresC =  paramSet['C_thres,C']/s.mg_per_molC# mgC/cm3 water => mol/cm3 water
    s.C_S_W_thresO =  paramSet['C_thres,O']/s.mg_per_molC# mgC/cm3 water => mol/cm3 water

    Ds = paramSet['DS_W'] #cm^2/d
    Dl = 0.003456 #cm^2/d
    s.Ds = Ds /(24*3600) /10000 # m^2/s
    s.Dl = Dl /(24*3600) /10000# m^2/s

    s.k_SC =  paramSet['k_C,S'] * s.mg_per_molC/paramSet['theta'] #[cm^3 water /mgC/d] => [cm^3 bulk soil/mol/d ]
    s.k_DC = paramSet['k_d,C'] #1/d
    s.k_DO = paramSet['k_d,O']#1/d

    s.K_L =  paramSet['K_L'] /s.mg_per_molC#[mol cm-3 soil]
    s.k_SO =  paramSet['k_O,S'] * s.mg_per_molC/paramSet['theta']
    #s.thetaInit = paramSet['theta']

    s.k_RC = paramSet['k_r,C']#1/d
    s.k_RO = paramSet['k_r,O']#1/d

    s.m_maxC =  paramSet['m_max,C']#1/d,  Maximum maintenance rate coefficient for the corresponding microbial groups
    s.m_maxO =  paramSet['m_max,O']#1/d,  Maximum maintenance rate coefficient for the corresponding microbial groups
    s.k_decay2 =  paramSet['p_L'] # - , Proportion of large polymers formed from dead microbial biomass due to maintenance

    s.micro_maxC = paramSet['u_max,C']#1/d
    s.micro_maxO =  paramSet['u_max,O']#1/d

    s.v_maxL = paramSet['v_max,L'] #1/d
    s.k_decay = paramSet['Y'] # -, Maintenance yield
    s.k_growthC = paramSet['Y_C'] # -
    s.k_growthO = paramSet['Y_O'] #- 
    # can cause issue
    s.k_sorp = paramSet['k_sorp'] /s.mg_per_molC# mg C/cm 3 soil solution =>  mol C/cm3 soil solution
    
    s.alpha = 0.1 # -
    s.f_sorp = 0# 0.5
    s.k_phi = 0.1
    s.C_aOLim=1.e-10 * float(s.doSoluteFlow) # so that microbe community can always regrow
    s.C_aCLim=1.e-10 * float(s.doSoluteFlow) # so that microbe community can always regrow
    s.f_Im = 4.1*1e-5*(24*3600)*1e-6 #[cm3/mol/d] or [1/d]
    s.f_Min = 2.3*1e-5*(24*3600)*1e-6 #[1/d]
    s.v_maxNH4 = 4.4*1e4*(24*3600) #[cm3/mol/d] or [1/d]
    s.v_maxNO3 = 4.4*1e4*(24*3600) #[1/d]
    
    if s.noAds:
        s.CSSmax = 0.
        s.alpha = 0.
        s.kadsN = 0.
        s.kdesN = 0.
    else:
        s.Qmmax = 0.45 * 0.079 # max ratio gOC-gmineral soil, see 10.1016/j.soilbio.2020.107912
        # [g OC / g mineral soil] * [g mineral soil/ cm3 bulk soil] *[ mol C/g C]
        CSSmax_ = s.Qmmax * s.bulkMassDensity_gpercm3*(1/s.molarMassC)
        s.CSSmax = CSSmax_ # mol C/cm3 bulk soil
        #s.CSSmax = s.Qmmax * s.bulkDensity_m3 / 1e6 # mol OC/mol soil * [mol soil/m3] * [m3/cm3] =  mol/cm3

        
    s.css1Function = 9 # current adsorption function implemented.
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
    return s
    
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

def setBiochemParam(s):
    """ send the TraiRhizo biochemical parameters to dumux
        @param: the dumux soil object
    """

    s.setParameter( "Soil.css1Function", str(s.css1Function))
    s.setParameter("Soil.C_aOLim", str(s.C_aOLim)) #[molC/cm3 scv]
    s.setParameter("Soil.C_aCLim", str(s.C_aCLim)) #[molC/cm3 scv]

    s.setParameter("1.Component.LiquidDiffusionCoefficient", str(s.Ds)) #m^2/s
    s.setParameter("2.Component.LiquidDiffusionCoefficient", str(s.Dl)) #m^2/s


    s.setParameter("Soil.betaC", str(s.betaC)) # -
    s.setParameter("Soil.betaO", str(s.betaO )) # -
    s.setParameter("Soil.C_S_W_thresC", str(s.C_S_W_thresC )) #mol/cm3
    s.setParameter("Soil.C_S_W_thresO", str(s.C_S_W_thresO )) #mol/cm3
    s.setParameter("Soil.k_decay", str(s.k_decay ))
    s.setParameter("Soil.k_decay2", str(s.k_decay2))
    s.setParameter("Soil.k_DC", str(s.k_DC )) # 1/d
    s.setParameter("Soil.k_DO", str(s.k_DO )) # 1/d
    s.setParameter("Soil.k_growthC", str(s.k_growthC ))
    s.setParameter("Soil.k_growthO", str(s.k_growthO ))
    s.setParameter("Soil.K_L", str(s.K_L))#[mol/cm3]
    s.setParameter("Soil.k_phi", str(s.k_phi ))
    s.setParameter("Soil.k_RC", str(s.k_RC))
    s.setParameter("Soil.k_RO", str(s.k_RO))

    s.setParameter("Soil.k_SC", str(s.k_SC )) #cm^3/mol/d
    s.setParameter("Soil.k_SO", str(s.k_SO )) #cm^3/mol/d

    s.setParameter("Soil.m_maxC", str(s.m_maxC  ))# 1/d
    s.setParameter("Soil.m_maxO", str(s.m_maxO  ))# 1/d
    s.setParameter("Soil.micro_maxC", str(s.micro_maxC ))# 1/d
    s.setParameter("Soil.micro_maxO", str(s.micro_maxO ))# 1/d
    s.setParameter("Soil.v_maxL", str(s.v_maxL))#[d-1]

    s.setParameter("Soil.k_sorp", str(s.k_sorp)) # mol / cm3 soil water or mol
    s.setParameter("Soil.f_sorp", str(s.f_sorp)) #[-]
    s.setParameter("Soil.CSSmax", str(s.CSSmax)) #[mol/cm3 scv zone 1] or mol
    s.setParameter("Soil.alpha", str(s.alpha)) #[1/d]
    s.setParameter("Soil.kads", str(s.kads)) #[cm3/mol/d] or [1/d]
    s.setParameter("Soil.kdes", str(s.kdes)) #[1/d]
    if s.noAds:
        s.setParameter("Soil.kadsN", str(s.kadsN)) #[cm3/mol/d] or [1/d]
        s.setParameter("Soil.kdesN", str(s.kdesN)) #[1/d]
    s.setParameter("Soil.f_Im", str(s.f_Im)) #[cm3/mol/d] or [1/d]
    s.setParameter("Soil.f_Min", str(s.f_Min)) #[1/d]
    s.setParameter("Soil.v_maxNH4", str(s.v_maxNH4))
    s.setParameter("Soil.v_maxNO3", str(s.v_maxNO3))
    
    
    if s.dimWorld == 3:
        # 1 == True
        # if we define a source or sink for the cell 
        # (from the results of the 1d models),
        # do not compute on top of biochemical reactions in dumux
        s.setParameter("Problem.reactionExclusive", "1")
    
    return s                          

def initialize_dumux_nc_(  gId=0, a_in=0.02,
                a_out=0.03 ,seg_length=1.0 ,
                x=[-10000] ,   # cm
                c2 = 0.06,             # mol/mol scv
                NC = 10,logbase = 0.5, doCells = False,points = [0.02,0.03]):                                   # cm
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
        #soilVG = [0.08, 0.43, 0.04, 1.6, 50]
        #cyl.setVGParameters([soilVG])
        lb =  logbase
        
        #points = getPoints()
        
        cyl.createGrid1d(points)# cm
        setSoilParam(cyl)
        cyl.doSoluteFlow = True
        cyl.noAds = False
        cyl = getBiochemParam(cyl,61)
        cyl = setBiochemParam(cyl)
        if doCells:
            Cells = cyl.getCellCenters_().reshape(-1)
        else:
            Cells = []
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
        cyl.setInnerBC("fluxCyl", -0.1)  # [cm/day] #Y 0 pressure?
        cyl.setOuterBC("fluxCyl", 0.)
        #cyl.setParameter("2.Component.LiquidDiffusionCoefficient", str( 0.003456/(24*3600) /10000)) #m^2/s
        cyl.setParameter("Newton.MaxTimeStepDivisions",
                     str( 100) )
        cyl.setParameter("Newton.MaxSteps",
                     str( 100) )
            
        
        for j in range( 1, cyl.numComp):
            cyl.setParameter(str(j)+".Component.LiquidDiffusionCoefficient", str( 0.003456/(24*3600) /10000)) #m^2/s
            cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Type", str(3))
            cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Type", str(3))
            cyl.setParameter( "Soil.BC.Bot.C"+str(j)+"Value", str(0)) 
            cyl.setParameter( "Soil.BC.Top.C"+str(j)+"Value", str(0 ))
        cyl.setParameter( "Soil.BC.Bot.C2Value", str(ex2)) # mol/cm2/d
        cyl.setParameter( "Soil.BC.Bot.C1Value", str(ex2)) # mol/cm2/d
        # mol/cm2/day
        for j in range( 1, cyl.numComp):
            cyl.setParameter("Soil.IC.C"+str(j),cyl.dumux_str(c2) ) 
        #cyl.setParameter("Soil.IC.C2",cyl.dumux_str(c2) ) 
    
        if len(Cells) > 0:#in case we update a cylinder, to allow dumux to do interpolation
            assert(len(c2)==len(Cells))
            CellsStr = cyl.dumux_str(Cells/100)#cm -> m
            cyl.setParameter("Soil.IC.Z",CellsStr)# m
            if len(Cells)!= len(x):
                print("Cells, x",Cells, x, len(Cells), len(x))
                raise Exception
                
            #j = 2
            for j in range( 1, cyl.numComp):
                cyl.setParameter("Soil.IC.C"+str(j)+"Z",CellsStr) # m
            if len(Cells)!= len( c2):
                print("Cells,  cAll[j-1]",Cells,  c2, 
                        len(Cells), len(c2), j)
                raise Exception
                        
        #print('str(c2)',str(c2))


        cyl.maxDt = 250/(3600*24) # soilModel.maxDt_1DS
        cyl.setParameter("Soil.mucilCAffectsW", "true")
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
        
def doPrints(cyl,a_in,timeT, initMolm):
    print('CC2',(cyl.getCellCenters_().reshape(-1)[[0,-1]]-a_in),'cm')#/100
    print('mucil',cyl.getSolution(2)[[0,-1]]*rhoWM)# mol/cm3
    #print('getSolution0',cyl.getSolution(0)[[0,-1]])
    print('getSolutionHead',cyl.getSolutionHead()[[0,-1]])
    #print('getPressure',cyl.base.getPressure()[[0,-1]])
    pHead_mucil = np.array(cyl.base.getPressureHead())[[0,-1]]
    print('pHead_mucil',pHead_mucil)
    print('effect mucil',cyl.getSolutionHead()[[0,-1]] - pHead_mucil)
    print('getWaterContent',cyl.getWaterContent()[[0,-1]])
    print('getViscosity',np.array(cyl.base.getViscosity())[[0,-1]])
    print('getConductivity',np.array(cyl.base.getConductivity())[[0,-1]])  # [m/s]?
    # mol/cm2/d * s * d/s * cm * cm/m = mol/m 
    molinput = ex2 * timeT /(24.*3600.) * a_in*2.*np.pi #*100
    # mol/molW * molW/cm3W * cm3W/cm3 * cm2 * cm/m = mol/m
    currentCls = sum(cyl.getSolution(2)*rhoWM*cyl.getWaterContent()*cyl.getCellSurfacesCyl().reshape(-1))#*100 # 
    print(currentCls,currentCls-molinput-initMolm,'ex2',ex2,'a_in',a_in,'timeT',timeT)
ex2 = 1e-6/2#*0
a_in=0.02
a_out =0.03
NC = 10        
cw = 0.03 #77859 #0.4762 # g m /gW
molMarssW = 18 # g/mol
molMassMulC = 30 #g/mol
rhoW = 1 #g/cm3
rhoWM = rhoW/molMarssW #[g/cm3]*[g/mol] = mol/cm3
cll_i = 1e-3 #*30# *0# mol/cm3
mlFr = cll_i/rhoWM #cw / molMassMulC * molMarssW
print('mlFr',mlFr)
#cyl = initialize_dumux_nc_(c2 = mlFr, NC = NC,x=[-10000],a_in=a_in)
lb = 0.5
#points = np.array([a_in + (a_out - a_in)*i/(NC-1) for i in range(NC)])
points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), 
                              NC, base = lb)
CC = np.array([(points[i] + points[i+1])/2 for i in range(len(points)-1)])#/100
print('CC1',CC[[0,-1]])
print(len(points), len(CC),len([mlFr * CC[i] * 3e3 for i in range(NC-1)]))
#raise Exception
cyl = initialize_dumux_nc_(x=[-10000. for i in range(NC-1)],a_in=a_in,a_out=a_out,
                            c2 = [mlFr * CC[i] for i in range(NC-1)], #
                            NC = NC, doCells = True,points=points)
timeT = 0.
#doPrints(cyl,a_in,timeT,0)
cc = (cyl.getCellCenters_().reshape(-1))#/100.
getPressureHead = [[],[],[]]

CsolutionsWat = [[[] for i in range(3)] for j in range(cyl.numSoluteComp)]

#Cmucil = [[],[],[]]
theta = [[],[],[]]
Viscosity = [[],[],[]]
initpressure = [[],[],[]]
Ks = [[],[],[]]
getPressureHead[0] = np.array(cyl.base.getPressureHead())#*100.
for j in range(cyl.numDissolvedSoluteComp):
    CsolutionsWat[j][0] = cyl.getSolution(j+1)*rhoWM 
for j in range(cyl.numDissolvedSoluteComp,cyl.numSoluteComp):
    CsolutionsWat[j][0] = cyl.getSolution(j+1)*cyl.solidMolDensity /1e6 
#Cmucil[0] = cyl.getSolution(2)*rhoWM
theta[0] = cyl.getWaterContent()
Viscosity[0] = np.array(cyl.base.getViscosity())
initpressure[0] = cyl.getSolutionHead()#/100.
Ks[0] = np.array(cyl.base.getConductivity())*100.*24.*3600.
initCls = sum(cyl.getSolution(2)*rhoWM*cyl.getWaterContent()*cyl.getCellSurfacesCyl().reshape(-1))#*100 # 

for j in range(cyl.numSoluteComp):
    print(CsolutionsWat[j][0])
dt = 3600./2./24./3600.#10./24./60.
#raise Exception
for i in range(2):
    print(cyl.numSoluteComp)
    cyl.ddt = 1e-3
    cyl.solve(dt)#36000/24/3600/2, saveInnerDumuxValues_ =True)
    timeT += dt
    #doPrints(cyl,a_in,timeT, initCls)
    print('timeT',timeT)
    if False:
        print('molFr_mucil',cyl.getSolution(2))
        print('mucil',cyl.getSolution(2)*rhoWM)# mol/cm3
        print('getSolution0',cyl.getSolution(0))
        print('getSolutionHead',cyl.getSolutionHead())
        print('getPressure',cyl.base.getPressure())
        print('getPressureHead',cyl.base.getPressureHead())
        print('effect mucil',cyl.getSolutionHead() - cyl.base.getPressureHead())
        print('getWaterContent',cyl.getWaterContent())
        print('Viscosity',cyl.base.getViscosity())# [Pa*s]
        print('getConductivity',cyl.base.getConductivity())  
        print('getConductivity_cmd',cyl.base.getConductivity()[0]*100.*24.*3600.) # [m/s]?
        print('mobility',cyl.base.getMobility(0))
    print('end')
    #for j in range(cyl.numSoluteComp):
    #    CsolutionsWat[j][i+1] = cyl.getSolution(j+1)*rhoWM 
    for j in range(cyl.numDissolvedSoluteComp):
        CsolutionsWat[j][i+1] = cyl.getSolution(j+1)*rhoWM 
    for j in range(cyl.numDissolvedSoluteComp,cyl.numSoluteComp):
        CsolutionsWat[j][i+1] = cyl.getSolution(j+1)*cyl.solidMolDensity /1e6 
        
    #Cmucil[i+1] = cyl.getSolution(2)*rhoWM 
    theta[i+1] = cyl.getWaterContent()
    Viscosity[i+1] = np.array(cyl.base.getViscosity())
    getPressureHead[i+1] = np.array(cyl.base.getPressureHead())
    initpressure[i+1] = cyl.getSolutionHead()#/100.
    Ks[i+1] = np.array(cyl.base.getConductivity())*100.*24.*3600.
    

outputs = {'cc':cc,'getPressureHead':getPressureHead,'Cmucil':CsolutionsWat[1],'theta':theta,
            'Viscosity':Viscosity,'initpressure':initpressure,'Ks':Ks}
import pickle

with open("results/dataDumux10N.pkl", "wb") as file:
    pickle.dump(outputs, file)
    
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(10, 8), sharex=True)

# First subplot (top-left): getPressureHead
for i in range(3):
    axes[0, 0].plot(cc, getPressureHead[i], label=str(dt * i))
axes[0, 0].grid(True)
axes[0, 0].ticklabel_format(style='plain', useOffset=False)
axes[0, 0].legend()
axes[0, 0].set_title("Pressure Head")

# Second subplot (top-right): Cmucil
for i in range(3):
    axes[0, 1].plot(cc, CsolutionsWat[1][i], label=str(dt * i))
axes[0, 1].grid(True)
axes[0, 1].legend()
axes[0, 1].ticklabel_format(style='plain', useOffset=False)
axes[0, 1].set_title("Cmucil")

for i in range(3):
    axes[1, 0].plot(cc, theta[i], label=str(dt * i))
axes[1, 0].grid(True)
axes[1, 0].legend()
axes[1, 0].ticklabel_format(style='plain', useOffset=False)
axes[1, 0].set_title("theta")
for i in range(3):
    axes[1, 1].plot(cc, Viscosity[i], label=str(dt * i))
axes[1, 1].grid(True)
axes[1, 1].legend()
axes[1, 1].ticklabel_format(style='plain', useOffset=False)
axes[1, 1].set_title("Viscosity")
for i in range(3):
    axes[2, 0].plot(cc, initpressure[i], label=str(dt * i))
axes[2, 0].grid(True)
axes[2, 0].legend()
axes[2, 0].ticklabel_format(style='plain', useOffset=False)
axes[2, 0].set_title("initpressure")
for i in range(3):
    axes[2, 1].plot(cc, Ks[i], label=str(dt * i))
axes[2, 1].grid(True)
axes[2, 1].legend()
axes[2, 1].ticklabel_format(style='plain', useOffset=False)
axes[2, 1].set_title("Ks")


# Adjust layout
plt.grid(True)
plt.legend()
plt.tight_layout()
#plt.show()
plt.savefig("mucilN1.png", dpi=300, bbox_inches='tight')#, transparent=True)
plt.close()

fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(12, 8), sharex=True)

# First subplot (top-left): getPressureHead
for i in range(3):
    axes[0, 0].plot(cc, CsolutionsWat[0][i], label=str(dt * i))
axes[0, 0].grid(True)
axes[0, 0].legend()
axes[0, 0].ticklabel_format(style='plain', useOffset=False)
axes[0, 0].set_title("CL^l")
for i in range(3):
    axes[1, 0].plot(cc, CsolutionsWat[2][i], label=str(dt * i))
axes[1, 0].grid(True)
axes[1, 0].legend()
axes[1, 0].ticklabel_format(style='plain', useOffset=False)
axes[1, 0].set_title("CH^l")

for i in range(3):
    axes[0, 1].plot(cc, CsolutionsWat[1][i], label=str(dt * i))
axes[0, 1].grid(True)
axes[0, 1].legend()
axes[0, 1].ticklabel_format(style='plain', useOffset=False)
axes[0, 1].set_title("mucilage")

for i in range(3):
    axes[1, 1].plot(cc, CsolutionsWat[3][i], label=str(dt * i))
axes[1, 1].grid(True)
axes[1, 1].legend()
axes[1, 1].ticklabel_format(style='plain', useOffset=False)
axes[1, 1].set_title("NL^l")

for i in range(3):
    axes[2, 0].plot(cc, CsolutionsWat[4][i], label=str(dt * i))
axes[2, 0].grid(True)
axes[2, 0].legend()
axes[2, 0].ticklabel_format(style='plain', useOffset=False)
axes[2, 0].set_title("NH^l")
for i in range(3):
    axes[2, 1].plot(cc, CsolutionsWat[5][i], label=str(dt * i))
axes[2, 1].grid(True)
axes[2, 1].legend()
axes[2, 1].ticklabel_format(style='plain', useOffset=False)
axes[2, 1].set_title("NH4^l")
for i in range(3):
    axes[3, 0].plot(cc, CsolutionsWat[6][i], label=str(dt * i))
axes[3, 0].grid(True)
axes[3, 0].legend()
axes[3, 0].ticklabel_format(style='plain', useOffset=False)
axes[3, 0].set_title("NO3")
# Adjust layout
plt.grid(True)
plt.legend()
plt.tight_layout()
#plt.show()
plt.savefig("mucilN2.png", dpi=300, bbox_inches='tight')#, transparent=True)
plt.close()

fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(12, 8), sharex=True)
for i in range(3):
    axes[0, 0].plot(cc, CsolutionsWat[9][i], label=str(dt * i))
axes[0, 0].grid(True)
axes[0, 0].legend()
axes[0, 0].ticklabel_format(style='plain', useOffset=False)
axes[0, 0].set_title("Cca")
for i in range(3):
    axes[0, 1].plot(cc, CsolutionsWat[10][i], label=str(dt * i))
axes[0, 1].grid(True)
axes[0, 1].legend()
axes[0, 1].ticklabel_format(style='plain', useOffset=False)
axes[0, 1].set_title("Ccd")
for i in range(3):
    axes[1, 0].plot(cc, CsolutionsWat[15][i], label=str(dt * i))
axes[1, 0].grid(True)
axes[1, 0].legend()
axes[1, 0].ticklabel_format(style='plain', useOffset=False)
axes[1, 0].set_title("Nca")
for i in range(3):
    axes[1, 1].plot(cc, CsolutionsWat[16][i], label=str(dt * i))
axes[1, 1].grid(True)
axes[1, 1].legend()
axes[1, 1].ticklabel_format(style='plain', useOffset=False)
axes[1, 1].set_title("Ncd")

for i in range(3):
    axes[2, 0].plot(cc, CsolutionsWat[7][i], label=str(dt * i))
axes[2, 0].grid(True)
axes[2, 0].legend()
axes[2, 0].ticklabel_format(style='plain', useOffset=False)
axes[2, 0].set_title("Coa")
for i in range(3):
    axes[2, 1].plot(cc, CsolutionsWat[8][i], label=str(dt * i))
axes[2, 1].grid(True)
axes[2, 1].legend()
axes[2, 1].ticklabel_format(style='plain', useOffset=False)
axes[2, 1].set_title("Cod")
for i in range(3):
    axes[3, 0].plot(cc, CsolutionsWat[13][i], label=str(dt * i))
axes[3, 0].grid(True)
axes[3, 0].legend()
axes[3, 0].ticklabel_format(style='plain', useOffset=False)
axes[3, 0].set_title("Noa")
for i in range(3):
    axes[3, 1].plot(cc, CsolutionsWat[14][i], label=str(dt * i))
axes[3, 1].grid(True)
axes[3, 1].legend()
axes[3, 1].ticklabel_format(style='plain', useOffset=False)
axes[3, 1].set_title("Nod")
plt.grid(True)
plt.legend()
plt.tight_layout()
#plt.show()
plt.savefig("mucilN3.png", dpi=300, bbox_inches='tight')#, transparent=True)
plt.close()


fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(10, 8), sharex=True)
for i in range(3):
    axes[0, 0].plot(cc, CsolutionsWat[11][i], label=str(dt * i))
axes[0, 0].grid(True)
axes[0, 0].legend()
axes[0, 0].ticklabel_format(style='plain', useOffset=False)
axes[0, 0].set_title("CL^s")
for i in range(3):
    axes[0, 1].plot(cc, CsolutionsWat[18][i], label=str(dt * i))
axes[0, 1].grid(True)
axes[0, 1].legend()
axes[0, 1].ticklabel_format(style='plain', useOffset=False)
axes[0, 1].set_title("NL^s")
for i in range(3):
    axes[1, 0].plot(cc, CsolutionsWat[12][i], label=str(dt * i))
axes[1, 0].grid(True)
axes[1, 0].legend()
axes[1, 0].ticklabel_format(style='plain', useOffset=False)
axes[1, 0].set_title("CO2")
for i in range(3):
    axes[1, 1].plot(cc, CsolutionsWat[19][i], label=str(dt * i))
axes[1, 1].grid(True)
axes[1, 1].legend()
axes[1, 1].ticklabel_format(style='plain', useOffset=False)
axes[1, 1].set_title("NN2")
for i in range(3):
    axes[2, 0].plot(cc, CsolutionsWat[17][i], label=str(dt * i))
axes[2, 0].grid(True)
axes[2, 0].legend()
axes[2, 0].ticklabel_format(style='plain', useOffset=False)
axes[2, 0].set_title("NH4^s")
plt.grid(True)
plt.legend()
plt.tight_layout()
#plt.show()
plt.savefig("mucilN4.png", dpi=300, bbox_inches='tight')#, transparent=True)
plt.close()
'''
maxCggS = 10*1e-3 g/g soil
rho_r = 1.26# g cm-3
rho_w = 1# g cm-3
theta = # cm3 cm-3
maxCggW = 10 * rho_r / (rho_w + theta)

##
cw = 0.077859 #0.4762 # g m /gW

'''