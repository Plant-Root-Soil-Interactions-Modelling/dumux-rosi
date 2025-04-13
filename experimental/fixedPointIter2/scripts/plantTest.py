import sys;
import os

#os.chdir('experimental/fixedPointIter2/scripts')
sys.path.append("../modules/");
sys.path.append("../inputDataTraiRhizo/");
sys.path.append("../../../../CPlantBox/");
sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../../build-cmake/cpp/python_binding/");


import plantbox as pb 
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
results_dir="./results/plantN/"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    test = os.listdir(results_dir)
    for item in test:
        try:
            os.remove(results_dir+item)
        except:
            pass
path = "../../../../CPlantBox/modelparameter/structural/plant/"
fname = path+"small_2020.xml"
from rhizo_modelsPlant import RhizoMappedSegments  # Helper class for cylindrical rhizosphere models
from functional.phloem_flux import PhloemFluxPython  # root system Python hybrid solver


def weather(simDuration, dt):
        Qnigh = 0; Qday = 960e-6 
        hp = 1.
        Tnigh = 15.8; Tday = 22
        RHday = 0.6; RHnigh = 0.88
        Pair = 1010.00 #hPa
        pmean = -100.
        cs = 350e-6
            
        coefhours = 1.
        RH_ = RHnigh + (RHday - RHnigh) * coefhours
        TairC_ = Tnigh + (Tday - Tnigh) * coefhours
        Q_ = Qnigh + (Qday - Qnigh) * coefhours
        
        es =  6.112 * np.exp((17.67 * TairC_)/(TairC_ + 243.5))
        ea = es*RH_
        assert ea < es
        
        assert ((RH_ > 0) and(RH_ < 1))
        bl_thickness = 1/1000 #1mm * m_per_mm
        diffusivity= 2.5e-5#m2/sfor 25*C
        rbl =bl_thickness/diffusivity #s/m 13
        
        Kcanopymean = 1e-1 # m2/s
        meanCanopyL = (2/3) * hp /2
        rcanopy = meanCanopyL/Kcanopymean
        windSpeed = 2 #m/s
        zmzh = 2 #m
        karman = 0.41 #[-]
        
        rair = 1
        if hp > 0:
            rair = np.log((zmzh - (2/3)*hp)/(0.123*hp)) * np.log((zmzh - (2/3)*hp)/(0.1*hp)) / (karman*karman*windSpeed)

        weatherVar = {'TairC' : TairC_,'TairK' : TairC_ + 273.15,'Pair':Pair,"es":es,
                        'Qlight': Q_,'rbl':rbl,'rcanopy':rcanopy,'rair':rair,"ea":ea,
                        'cs':cs, 'RH':RH_, 'p_mean':pmean
                     }
        return weatherVar
def main():
    initSim = 10.
    weatherInit = weather(initSim,0)
    seed = 1
    ms = pb.MappedPlant(seed)
    ms.setSeed(seed)
    ms.readParameters(path + fname)
    ms.setGeometry(pb.SDF_PlantBox(np.inf, np.inf, 500))
    ms.initialize(verbose = False)
    ms.simulate(10,verbose= False)
    plantModel = PhloemFluxPython(ms,psiXylInit = -659.8 -10,ciInit = 350e-6*0.5) 
    ms.constantLoc = True
    picker = lambda x, y, z: int(-10*z)
    plantModel.rs.setSoilGrid(picker)

    plantModel.wilting_point = -15000.
    # set kr and kx for root system or plant

    import plantParameters
    plantParameters.init_conductivities(r = plantModel)
    plantModel = plantParameters.phloemParam(plantModel, weatherInit)
    plantModel.setKrm2([[0.]], False)
    plantModel.setKrm1([[0.]], False) 

    ms.simulate(initSim)
    print(ms.organTypes)
    sxx = np.array([-1000 for i in range(len(ms.organTypes))])
    sxx[ms.organTypes==2] = -100
    plantModel.solve_photosynthesis(sim_time_ = 10, 
                        sxx_=sxx, 
                        cells_ = False,
                        ea_ = weatherInit["ea"],#not used
                        es_= weatherInit["es"],#not used
                        verbose_ = False, doLog_ = False,
                        TairC_= weatherInit["TairC"],#not used
                                    soil_k_ = [], # [day-1]
                        outputDir_= "./results/rhizoplantExud")
    print(plantModel.psiXyl,np.array(plantModel.outputFlux))
    plantModel.withInitVal = True
    plantModel.initValST = 2.
    Qout = [[[] for i in range(3)] for j in range(10,15)]
    for i in range(3):
        dt = 1.
        plantModel.startPM(10., 10.+dt, 1, ( weatherInit["TairC"]  +273.15) , False, ".")
        Nt = len( ms.nodes)
        #print(plantModel.Q_out)
        print(plantModel.Q_Grmax)
        print(plantModel.Fpsi)
        for idc, coef in enumerate(range(10,15)):
            Qout[idc][i] = plantModel.Q_out[(Nt*coef):(Nt*(coef+1))]
        ms.simulate(dt,verbose= False)
    cc = np.array([i for i in range(Nt)])
    if True:
        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(10, 8), sharex=True)

        # First subplot (top-left): getPressureHead
        for i in range(3):
            axes[0, 0].plot(cc, Qout[0][i], label=str(dt * i))
        axes[0, 0].grid(True)
        axes[0, 0].ticklabel_format(style='plain', useOffset=False)
        axes[0, 0].legend()
        axes[0, 0].set_title("QN_ST")

        # Second subplot (top-right): Cmucil
        for i in range(3):
            axes[0, 1].plot(cc, Qout[1][i], label=str(dt * i))
        axes[0, 1].grid(True)
        axes[0, 1].legend()
        axes[0, 1].ticklabel_format(style='plain', useOffset=False)
        axes[0, 1].set_title("QN_Xyl")

        for i in range(3):
            axes[1, 0].plot(cc, Qout[2][i], label=str(dt * i))
        axes[1, 0].grid(True)
        axes[1, 0].legend()
        axes[1, 0].ticklabel_format(style='plain', useOffset=False)
        axes[1, 0].set_title("QN_Cell")
        for i in range(3):
            axes[1, 1].plot(cc, Qout[3][i], label=str(dt * i))
        axes[1, 1].grid(True)
        axes[1, 1].legend()
        axes[1, 1].ticklabel_format(style='plain', useOffset=False)
        axes[1, 1].set_title("QN_Struct")
        for i in range(3):
            axes[2, 0].plot(cc, Qout[4][i], label=str(dt * i))
        axes[2, 0].grid(True)
        axes[2, 0].legend()
        axes[2, 0].ticklabel_format(style='plain', useOffset=False)
        axes[2, 0].set_title("QN_Store")
        
        ot = list(ms.organTypes)
        ot = [2] + ot
        
        
        for i in range(3):
            axes[2, 1].plot(cc,np.array(ot), label=str(dt * i))
        axes[2, 1].grid(True)
        axes[2, 1].legend()
        axes[2, 1].ticklabel_format(style='plain', useOffset=False)
        axes[2, 1].set_title("organTypes")
                
        # Adjust layout
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        #plt.show()
        plt.savefig("plantN.png", dpi=300, bbox_inches='tight')#, transparent=True)
        plt.close()
        
main()