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
    Ds = paramSet['DS_W'] #cm^2/d
    Dl = 0.003456 #cm^2/d
    s.Ds = Ds /(24*3600) /10000 # m^2/s
    s.Dl = Dl /(24*3600) /10000# m^2/s
    s.BufferPower = 0.#140 # random value
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
    
    s.setParameter("1.Component.LiquidDiffusionCoefficient", str(s.Ds)) #m^2/s
    s.setParameter("2.Component.LiquidDiffusionCoefficient", str(s.Dl)) #m^2/s
    s.setParameter("Component.BufferPower", str(s.BufferPower))
    if False:
        
        s.setParameter("Soil.f_sorp", str(s.f_sorp)) #[-]
        s.setParameter("Soil.CSSmax", str(s.CSSmax)) #[mol/cm3 scv zone 1] or mol
        s.setParameter("Soil.alpha", str(s.alpha)) #[1/d]
        s.setParameter("Soil.kads", str(s.kads)) #[cm3/mol/d] or [1/d]
        s.setParameter("Soil.kdes", str(s.kdes)) #[1/d]
        if s.noAds:
            s.setParameter("Soil.kadsN", str(s.kadsN)) #[cm3/mol/d] or [1/d]
            s.setParameter("Soil.kdesN", str(s.kdesN)) #[1/d]
         
        
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
for j in range(cyl.numDissolvedSoluteComp):
    CsolutionsWat[j][0] = cyl.getSolution(j+1)*rhoWM 
for j in range(cyl.numDissolvedSoluteComp,cyl.numSoluteComp):
    CsolutionsWat[j][0] = cyl.getSolution(j+1)*cyl.solidMolDensity /1e6 
#Cmucil[0] = cyl.getSolution(2)*rhoWM

theta[0] = cyl.getWaterContent()

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
    print('end')
    #for j in range(cyl.numSoluteComp):
    #    CsolutionsWat[j][i+1] = cyl.getSolution(j+1)*rhoWM 
    for j in range(cyl.numDissolvedSoluteComp):
        CsolutionsWat[j][i+1] = cyl.getSolution(j+1)*rhoWM 
    for j in range(cyl.numDissolvedSoluteComp,cyl.numSoluteComp):
        CsolutionsWat[j][i+1] = cyl.getSolution(j+1)*cyl.solidMolDensity /1e6 
        
    #Cmucil[i+1] = cyl.getSolution(2)*rhoWM 
    theta[i+1] = cyl.getWaterContent()
    initpressure[i+1] = cyl.getSolutionHead()#/100.
    

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
