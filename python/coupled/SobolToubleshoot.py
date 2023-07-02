import sys; 
CPBdir = "../../../../cpb3101/CPlantBox"
sys.path.append(CPBdir+"/src");
sys.path.append(CPBdir);
sys.path.append("../../..");sys.path.append(".."); 
sys.path.append(CPBdir+"/src/python_modules");
sys.path.append("../build-cmake/cpp/python_binding/") # dumux python binding
sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../modules/") # python wrappers 
import numpy as np
import json
import plantbox as pb
import importlib
importlib.reload(pb)
import plantbox
importlib.reload(plantbox)
from phloem_flux import PhloemFluxPython 
import pickle
from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.test_functions import Ishigami
from joblib import Parallel, delayed
import os 
from helpUqr_seb import *


simDuration=10; condition="wet";idtorun=151
weatherX = weather(simDuration, condition, 1)
assert weatherX['Qlight']>0




namesVars = ["maxTil" ,"firstTil","delayTil" ,"maxB","firstB","delayB",
             "delayLat","delayNGStart","delayNGEnd",
    "lar0", "lbr0", "lnr0", "lmaxr0", "rr0", "ar0", "tropismNr0", "tropismSr0", "thetar0",
    "lar1", "lbr1", "lnr1", "lmaxr1", "rr1", "ar1", "tropismNr1", "tropismSr1", "thetar1",
    "lar2", "lbr2", "lnr2", "lmaxr2", "rr2", "ar2", "tropismNr2", "tropismSr2", "thetar2",
    "last", "lbst", "lnst", "lmaxst", "rst", "ast", "tropismNst", "tropismSst", "thetast",
     "lmaxle", "rle", "ale", "tropismNle", "tropismSle", "thetale","Width_petiole","Width_blade",
            "areaMax"]   


limsVars = [[0.5,1.5] for ii in namesVars]
problem = {
    'num_vars': len(namesVars),
    'names': namesVars,
    'bounds': limsVars
}


##### TO CHANGE ######
Nreps = 1#8
reps =2**Nreps
maxcore =  os.cpu_count()
######################

param_values = saltelli.sample(problem, reps)



Yall = [np.array([]) for i in range(7)]#Yexuds,Ygrs,Yrms,Ytrans,Yassi,Ywue]
namesY = ["Yexuds","Ygrs","Yrms","Ytrans","Yassi","Ywue","YIwue"]

myid=idtorun;varNames=namesVars; varLims=list(param_values[myid]); testType="SEB"


assert testType=="Xylem" or testType == "Phloem" or testType == "SEB"
test_values = varLims
DictVal = {}
for key in varNames:
    for value in test_values:
        DictVal[key] = value
        test_values.remove(value)
        break
DictVal={'maxTil': 1.09375, 'firstTil': 1.46875, 'delayTil': 1.46875, 'maxB': 0.65625, 'firstB': 1.28125, 
         'delayB': 0.96875, 'delayLat': 0.53125, 'delayNGStart': 0.84375, 'delayNGEnd': 1.46875, 'lar0': 1.15625, 
         'lbr0': 1.09375, 'lnr0': 1.40625, 'lmaxr0': 0.65625, 'rr0': 0.65625, 'ar0': 1.34375, 'tropismNr0': 1.03125, 
         'tropismSr0': 0.78125, 'thetar0': 0.59375, 'lar1': 0.78125, 'lbr1': 1.03125, 'lnr1': 0.71875, 'lmaxr1': 0.90625, 
         'rr1': 0.59375, 'ar1': 0.59375, 'tropismNr1': 0.53125, 'tropismSr1': 1.15625, 'thetar1': 1.15625, 'lar2': 0.84375, 
         'lbr2': 1.21875, 'lnr2': 1.34375, 'lmaxr2': 1.40625, 'rr2': 1.34375, 'ar2': 1.09375, 'tropismNr2': 0.90625, 
         'tropismSr2': 1.46875, 'thetar2': 0.96875, 'last': 1.28125, 'lbst': 1.03125, 'lnst': 1.03125, 'lmaxst': 0.90625, 
         'rst': 0.78125, 'ast': 1.34375, 'tropismNst': 0.96875, 'tropismSst': 0.53125, 'thetast': 1.03125, 'lmaxle': 0.71875, 'rle': 1.28125, 
         'ale': 0.78125, 'tropismNle': 1.28125, 'tropismSle': 0.53125, 'thetale': 0.96875, 
         'Width_petiole': 0.53125, 'Width_blade': 1.03125, 'areaMax': 0.71875}
pl = pb.MappedPlant(seednum = 2) 
path = CPBdir+"/modelparameter/plant/"
name = "Triticum_aestivum_adapted_2023"

pl.readParameters(path + name + ".xml")
depth = 60
sdf = pb.SDF_PlantBox(np.Inf, np.Inf, depth )
pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil



for p in pl.getOrganRandomParameter(pb.leaf):
    p.lmax *= DictVal['lmaxle']
    p.areaMax *= DictVal['areaMax']
    p.tropismN *= DictVal['tropismNle']
    p.tropismS *= DictVal['tropismSle']
    p.r *= DictVal["rle"]
    p.a *= DictVal["ale"]
    p.theta *= DictVal["thetale"]
    p.Width_petiole *= DictVal["Width_petiole"]
    p.Width_blade *= DictVal["Width_blade"]

for p in pl.getOrganRandomParameter(pb.stem):
    p.lmax *= DictVal['lmaxst']
    p.tropismN *= DictVal['tropismNst']
    p.tropismS *= DictVal['tropismSst']
    p.r *= DictVal["rst"]
    p.a *= DictVal["ast"]
    p.la *= DictVal["last"]
    p.lb *= DictVal["lbst"]
    p.ln *= DictVal["lnst"]
    p.delayLat *= DictVal["delayLat"]
    p.delayNGStart *= DictVal["delayNGStart"]
    p.delayNGEnd *= DictVal["delayNGEnd"]

for p in pl.getOrganRandomParameter(pb.root):
    if (p.subType ==0):
        pass
    elif (p.subType ==1)or(p.subType >3):
        p.lmax *= DictVal['lmaxr0']
        p.tropismN *= DictVal['tropismNr0']
        p.tropismS *= DictVal['tropismSr0']
        p.r *= DictVal["rr0"]
        p.a *= DictVal["ar0"]
        p.la *= DictVal["lar0"]
        p.lb *= DictVal["lbr0"]
        p.ln *= DictVal["lnr0"]
        p.theta *= DictVal["thetar0"]
    elif (p.subType ==2):
        p.lmax *= DictVal['lmaxr1']
        p.tropismN *= DictVal['tropismNr1']
        p.tropismS *= DictVal['tropismSr1']
        p.r *= DictVal["rr1"]
        p.a *= DictVal["ar1"]
        p.la *= DictVal["lar1"]
        p.lb *= DictVal["lbr1"]
        p.ln *= DictVal["lnr1"]
        p.theta *= DictVal["thetar1"]
    elif (p.subType ==3):
        p.lmax *= DictVal['lmaxr2']
        p.tropismN *= DictVal['tropismNr2']
        p.tropismS *= DictVal['tropismSr2']
        p.r *= DictVal["rr2"]
        p.a *= DictVal["ar2"]
        p.theta *= DictVal["thetar2"]
    else:
        print("root subtype not recognized",p.subType)
        raise Exception("root subtype not recognized")

for p in pl.getOrganRandomParameter(pb.seed):
    p.maxTil =int(round(p.maxTil* DictVal['maxTil']))
    p.firstTil *= DictVal['firstTil']
    p.delayTil *= DictVal['delayTil']
    p.firstB *= DictVal['firstB']
    p.delayB *= DictVal['delayB']
    p.maxB = int(round(p.maxB*DictVal["maxB"]))

pl.initialize(verbose = True)#, stochastic = False)
pl.simulate(simDuration, False)#, "outputpm15.txt")

picker = lambda x,y,z : max(int(np.floor(-z)),-1)   
pl.setSoilGrid(picker)  # maps segment


r = PhloemFluxPython(pl,psiXylInit =-1000,ciInit = 0.5 *weatherX['cs'])


hp = max([tempnode[2] for tempnode in r.get_nodes()]) /100
weatherX = weather(simDuration, condition, hp)

plant =setDefaultVals(r,weatherX,DictVal)

#xXx = np.full((len(r.plant.nodes)),weatherX["p_mean"])
p_bot = weatherX["p_mean"] + depth/2
p_top = weatherX["p_mean"] - depth/2
xXx = np.linspace(p_top, p_bot, depth)



plant.Patm = weatherX["Pair"]
##resistances
plant.g_bl = resistance2conductance(weatherX["rbl"],weatherX,r) / r.a2_bl
plant.g_canopy = resistance2conductance(weatherX["rcanopy"],weatherX,r) / r.a2_canopy
plant.g_air = resistance2conductance(weatherX["rair"],weatherX,r) / r.a2_air

verbose_phloem = True


plant.solve_photosynthesis(sim_time_ = simDuration, sxx_=xXx, 
                           cells_ = True,ea_ = weatherX["ea"],es_ = weatherX["es"],
        verbose_ = False, doLog_ = False,TairC_= weatherX["TairC"] ,outputDir_= "./results/")

filename ="./inPMerror.txt"
plant.startPM(simDuration, simDuration + 3/(24), 1, ( weatherX["TairC"]  +273.15) , True, filename)