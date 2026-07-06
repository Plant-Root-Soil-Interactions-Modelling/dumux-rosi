import sys; 
sys.path.append("../../experimental/fixedPointIter2/modules"); 
sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richards_cyl import RichardsCylFoam #as RichardsCylFoamInit # C++ part (Dumux binding)
from rosi_richards_cyl import RichardsCylFoamAna
from rosi_richards5c_cyl import Richards5CCylFoam #as Richards5CCylFoamAna
from rosi_richards5c_cyl import Richards5CCylFoamAna
#from rosi_richards_cylExtrusion import  RichardsCylFoamExtrusion # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
from richards_no_mpi import RichardsNoMPIWrapper  
import functional.van_genuchten as vg

import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import numpy as np

""" 
Get the inner and outer boundary flows for the 1d axisymmetric soil
"""


def SorptionParams(soil_type, sorp): 
    # kads = [[0.16,0.73,15.1],[0.03,0.12, 2.51]] #cm^3/mol/d
    # kdes = [[205.7,14.4,96],[205.7, 14.4, 96]]  #1/d
    kads = [[5100,581,185],[867,98.77, 31.45]] #cm^3/mol/d
    kdes = [[74.5, 6.79, 0.12],[74.5, 6.79, 0.12]]  #1/d
    kads_ = kads[soil_type][sorp]
    kdes_ = kdes[soil_type][sorp]
    kd_values = {'kads': kads_, 'kdes': kdes_}
    
    return kd_values
    
def DecayParams(): 
    Vmax = 2.28e-4 #mol C / m^3 water/ s,0.73,15.1]
    Km = 20.6 #mol C / m^3 water
    decay_params =  {'Vmax': Vmax, 'Km': Km}
    
    return decay_params   
    
def DiffusionParams(Diff): 
    Dl = [1e-11, 1e-10, 1e-9, 1e-8]
    return Dl[Diff]
    
def getBiochemParam(s,soil_type, sorp, diff):    
    """ define TraiRhizo biochemical parameters 
        @param: the dumux soil object
        @ param: index of the TraiRhizo parameter set to use
        
    """
    s.doSimpleReaction = 1 #only diffusion, decay, sorption and not Mona's complete model 
    s.molarMassC = 12.011
    s.mg_per_molC = s.molarMassC * 1000.
    s.Ds = DiffusionParams(diff)  #m^2/s
    
    decay_params = DecayParams()
    s.vmax_decay = decay_params['Vmax'] *0 #mol C / m^3 scv / s #max decay rate from Nideggen et al. 
    s.km_decay = decay_params['Km'] #mol C / m^3 scv #michaelis constant from Nideggen et al. 
    
    if s.doAds:
        # kads = 1e1 #2.86e-04 #7.07e+02 # m3/kgC/yr, see 10.1016/j.soilbio.2020.107912, A.3
        # kdes =  1e3 #1.67e-01 #1.63e+03 # [1/yr] see 10.1016/j.soilbio.2020.107912, A.3
        kd_values = SorptionParams(soil_type, sorp)
        kads = kd_values['kads']*s.molarMassC/s.bulkMassDensity_gpercm3 #[cm^3/mol/d] optimized from doi: 10.1111/j.1365-2389.2010.01244.x
        kdes = kd_values['kdes']  #[1/d]             
        k_clay_silt = {}
        k_clay_silt[0] = 0.67
        k_clay_silt[1] = 0.082
        
        yr_per_d = 1/365 # [yr/d]
        m3_per_cm3 = 1e-6; # m3/cm3
        cm3_per_m3 = 1e6; # cm3/m3
        
        # [kg/g] * [g/mol] = kg/mol
        kgC_per_mol = (1/1000) * s.molarMassC
        # [m3/kgC/yr] * [yr/d] * [cm3/m3] * [kgC/mol] = [cm3/mol/d]
        s.kads = kads * yr_per_d * cm3_per_m3 * kgC_per_mol # [cm3/mol/d]
        s.kdes = kdes * yr_per_d # [1/d]
        s.Qmmax = k_clay_silt[soil_type] * 0.079 # max ratio gOC-gmineral soil, see 10.1016/j.soilbio.2020.107912
        # [g OC / g mineral soil] * [g mineral soil/ cm3 bulk soil] *[ mol C/g C]
        CSSmax_ = s.Qmmax * s.bulkMassDensity_gpercm3*(1/s.molarMassC)
        s.CSSmax = CSSmax_ # mol C/cm3 bulk soil
        #s.CSSmax = s.Qmmax * s.bulkDensity_m3 / 1e6 # mol OC/mol soil * [mol soil/m3] * [m3/cm3] =  mol/cm3
        
    else: 
        s.CSSmax = 0.
        s.alpha = 0.
        #s.css1Function = 0
        s.kads = 1
        s.kdes = 1
    
    
    s.css1Function = 9 # 5: instantaneous, 9: PDE adsorption
    assert (( s.css1Function == 5) or ( s.css1Function == 9))# or ( s.css1Function == 0))
    return s
    
def setBiochemParam(s):
    """ send the TraiRhizo biochemical parameters to dumux
        @param: the dumux soil object
    """
    s.setParameter( "Soil.doSimpleReaction", str( s.doSimpleReaction))
    s.setParameter( "Soil.css1Function", str(s.css1Function))

    #diffusion
    s.setParameter("1.Component.LiquidDiffusionCoefficient", str(s.Ds)) #m^2/s
    
    #decay
    s.setParameter("Soil.vmax_decay", str(s.vmax_decay)) #mol C / m^3 scv / s 
    s.setParameter("Soil.km_decay", str(s.km_decay)) #mol C / m^3 scv
   

    #sorption
    s.setParameter("Soil.CSSmax", str(s.CSSmax)) #[mol/cm3] 
    s.setParameter("Soil.kads", str(s.kads)) #[cm3/mol/d]
    s.setParameter("Soil.kdes", str(s.kdes))  #[1/d]
    
    for i in range(s.numSoluteComp):
        #mol/m3 to mol/mol
        s.setParameter( "Soil.IC.C"+str(i+1), str(0))
    
    if s.dimWorld == 3:
        # 1 == True
        # if we define a source or sink for the cell 
        # (from the results of the 1d models),
        # do not compute on top of biochemical reactions in dumux
        s.setParameter("Problem.reactionExclusive", "1")
    
    return s



def getSoilTextureAndShape(res = 1, soil_= 'loam'):  
    """ soil shape and texture data
        to adapt according to the soil represented
    """
    # min_b = np.array([-3./2, -3./2, -5.]) # np.array( [5, 5, 0.] )
    # max_b =np.array( [3./2, 3./2, 0.]) #  np.array([-5, -5, -5.])
    # cell_number = np.array([3,3,5]) #  [2,1,1])#np.array( [1,1,1]) # 1cm3
    # area = 3*3

    if res == 1: 
        min_b = np.array([-20/2, -45/2, -74.]) # np.array( [5, 5, 0.] )
        max_b =np.array( [20/2, 45/2, 0.]) #  np.array([-5, -5, -5.])
        cell_number = np.array([20,45,74]) # np.array( [3,12,40]) #np.array( [1,1,1]) # 1cm3
        area = 20 * 45  # cm2 45
    elif res == -1: 
        min_b = np.array([-1e3, -1e3, -2.]) # np.array( [5, 5, 0.] )
        max_b =np.array( [1e3, 1e3, 0.]) #  np.array([-5, -5, -5.])
        cell_number = np.array([1,1,1]) # np.array( [3,12,40]) #np.array( [1,1,1]) # 1cm3
        area = (max_b[0] - min_b[0]) * (max_b[1] - min_b[1])   # cm2 
    elif res == 2: 
        min_b = np.array([-20/2, -44/2, -74.])
        max_b =np.array( [20/2, 44/2, 0.]) 
        cell_number = np.array([10,22,37]) 
        area = 20 * 44  # cm2 45 
    elif res == 4: 
        min_b = np.array([-20/2, -44/2, -72.])
        max_b =np.array( [20/2, 44/2, 0.]) 
        cell_number = np.array([5,11,18]) 
        area = 20 * 44  # cm2 45        
    elif res == 5: 
        min_b = np.array([-20/2, -45/2, -75.])
        max_b =np.array( [20/2, 45/2, 0.]) 
        cell_number = np.array([4,9,15])
        area = 20 * 45  # cm2 45
    else: 
        print('Wrong resolution chosen') 
    
    solidDensity = 2650 # [kg/m^3 solid] #taken from google docs TraiRhizo
    solidMolarMass = 60.08e-3 # [kg/mol] 
    # theta_r, theta_s, alpha, n, Ks
   
    
    if soil_ == "loam":
        i = 0
    else:
        i = 1
        
    soilVG = [0.062, 0.337, 0.0182, 2.733, 1174] # vg_SPP(i)
    
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
def setSoilParam(s, soilTexture):    
    """ save the soil parameters
        @param: the dumux soil object
    """
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
    
def setDefault(s):
    """ Defined some usefull default parameters
    """
    molarMassWat = 18. # [g/mol]
    densityWat = 1. #[g/cm3]
    # [mol/cm3] = [g/cm3] /  [g/mol] 
    molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 
    s.molarDensityWat = molarDensityWat

    # low MaxRelativeShift == higher precision in dumux
    # s.setParameter("Problem.dobioChemicalReaction",str(s.doBioChemicalReaction))
    # s.setParameter("Problem.doDecay",str(s.doDecay))
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
    s.MaxTimeStepDivisions = 100
    s.setParameter("Newton.MaxTimeStepDivisions",
                     str( s.MaxTimeStepDivisions) )  
    s.MaxSteps = 50
    s.setParameter("Newton.MaxSteps",
                     str( s.MaxSteps) )  
    s.MaxRelativeShift = 1e-8
    s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift))
    
    return s


r_in = 0.02  # cm
r_out = 3.1916009083031205
length = 0.9999999999999999
    
def solve(simtimes, N, analytical = False):

    #loam = [0.08, 0.43, 0.04, 1.6, 5]  # K = 5 !
    if analytical:
        s = RichardsNoMPIWrapper(Richards5CCylFoamAna()) 
    else:
        s = RichardsNoMPIWrapper(Richards5CCylFoam()) 
        
        
    #s.setParameter( "Assembly.NumericDifference.BaseEpsilon", str(1e-10));
    #s.setParameter("Assembly.NumericDifference.PriVarMagnitude", s.dumux_str([1e-10, 1e-2, 1e-10, 1e-10]));
    s.initialize()
    s = setDefault(s)
    s.doAds = False
    s.doDecay = False
    s.doBioChemicalReaction = False
    s.setParameter("Problem.segLength", str(length))#
    s.setParameter("Soil.MolarMass", str(length))#
    s.createGrid1d(np.linspace(r_in, r_out, N))  # [cm]
    #s.setVGParameters([loam])
    soilTextureAndShape = getSoilTextureAndShape()
    setSoilParam(s, soilTextureAndShape)
    getBiochemParam(s,soil_type = 0,sorp = 0, diff = 0)
    setBiochemParam(s) 

    s.setHomogeneousIC(-1000.)  # cm pressure head
    s.setTopBC("constantFlux", 0.1)  #  [cm/day] "noFlux")#
    s.setBotBC("constantFlux", -0.2*0) # "noFlux")#
    cin = 9.8703454e-07 # mol /day
    for i in range(1, 2):#s.numComp):# no flux
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(2))
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(cin / (2 * np.pi * r_in * length))) # mol / cm -2 /day
    s.setParameter( "Soil.IC.C1", str(0 ))
        
    if analytical:
        s.setParameter("Problem.useExtrusion", "true")
    else:
        s.setParameter("Problem.useExtrusion", "false")
        
    s.initializeProblem(maxDt = 0.01)
    #s.eps_regularization = 1e-10 
    #s.setRegularisation(s.eps_regularization, s.eps_regularization) 
    cellVolumes = s.getCellSurfacesCyl() * length # cm3
    if rank == 0:
        print(s)

    s.ddt = 1.e-3  # days

    simtimes.insert(0, 0)
    dt_ = np.diff(simtimes)
    
    watpots = []
    #x = s.get
    
    for r, dt in enumerate(dt_):
        

        time = simtimes[r] 
        print('time',time, s.useMoles())
            
        #if time >= 5:
        #    s.setSoluteTopBC([1], [0.])

        if rank == 0:
            print("*****", "#", r, "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

        Wvolbefore = cellVolumes * s.getWaterContent() # cm3
        Smassbefore = s.getContent(1) # mol  * Wvolbefore
        #print("Wvolbefore",Wvolbefore,sum(Wvolbefore))
        print("Smassbefore",Smassbefore,sum(Smassbefore))
        s.solve(dt, saveInnerDumuxValues_ = True)
        
        Wvolafter = cellVolumes*s.getWaterContent() # cm3
        #print("Wvolafter",Wvolafter,sum(Wvolafter))
        Smassafter = s.getContent(1) #* Wvolafter # mol
        print("Smassafter",Smassafter,sum(Smassafter))
        
            
        rootSoilFluxes = s.getInnerFlow(0, length) * dt # cm3
        rootSoilFluxesS = s.getInnerFlow(1, length) * dt # mol
        soilSoilFluxes = s.getOuterFlow(0, length) * dt # cm3
        soilSoilFluxesS = s.getOuterFlow(1, length) * dt # mol
        print("rootSoilFluxesS",rootSoilFluxesS, s.getInnerFlow(1, length), dt)
        print("soilSoilFluxesS",soilSoilFluxesS, s.getOuterFlow(1, length), dt)
        print('getBoundaryFluxesPerFace_',s.getBoundaryFluxesPerFace_(1, length)) # mol/day
        print('getFaceSurfaces_',s.getFaceSurfaces_(length)) # cm2
        print("base.getScvfBoundaryFluxes()",s.base.getScvfBoundaryFluxes()[1][0] *24 * 3600 / 1e4) # mol/cm2/day
        print('neumann', s.base.getNeumann(0,1))
        
        # TODO: currently, setSource not properly implemented for richards and richards 2c
        # so left out of the mass balance.
        # scvSources = s.getSource(0) * cellVolumes * dt # cm3
        # scvSourcesS = s.getSource(1) * cellVolumes * dt # kg
        
        if rank == 0:
            print('\tChange in water volume [cm3] per voxel:',sum(Wvolafter-Wvolbefore))
            print('\tChange in solute mass [g] per voxel:',sum(Smassafter-Smassbefore))
            print('\tRMSE for water volume balance [cm3]:',np.mean(np.abs(rootSoilFluxes+soilSoilFluxes+sum(Wvolafter-Wvolbefore))))
            print('\tRMSE for solute mass balance [g]:',np.mean(np.abs(rootSoilFluxesS+soilSoilFluxesS+sum(Smassafter-Smassbefore))),'\n\n')
        print('getContent',s.getContent(1), s.getContent(2), s.getContent(3))   

       


if __name__ == "__main__":


    simTimes = [0.5]#,0.78,1]  # days
    simdataAna = solve(simTimes, 3, True)
    #simdata = solve(simTimes, 3, False)

"""
when we have no BC flow, seems to work ok, just inner flow.

todo:
- Y results so different for solutes in ana and numeric?
- error when i have more than one cell
- solutes works on its own but not when water is moving
    - jacobian not well filled?
    - cach not well rolled back or updated?
    - ussie that carbon advection according to numercial diff vs water flow from analytical diff?
    - ok so now inner BC is wrong for the numerical
"""