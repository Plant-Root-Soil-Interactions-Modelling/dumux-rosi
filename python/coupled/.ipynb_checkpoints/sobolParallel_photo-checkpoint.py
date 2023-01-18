directoryN = "/ipynbmorris/"
import sys; 
CPBdir = "/home/m.giraud/DUMUX/CPlantBox"
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



from helpUqr import *
Tmin = 20.7; Tmax = 30.27
specificHumidity = 0.0111
Pair = 1070.00 #hPa
thetaInit = 20/100
weatherX = weather(0, Tmin,Tmax,specificHumidity,Pair,thetaInit)




#			default value of env variable, to re-set at runtime
Patm = 1013.15;#default [hPa]
cs = 350e-6; #example from Dewar2002 [mol mol-1]
TairC = 20; #[°C]
Qlight = 900e-6;#mean absorbed photon irradiance per leaf segment [mol photons m-2 s-1]  
Chl = 55,; 
oi = 210e-3;#leaf internal [O2] [mol mol-1]

SPAD= 100.0
chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10

#			parameter to re-parametrise , put in phloem files
#water stress factor, parametrised from data of Corso2020
fwr = 9.308e-2; #residual opening when water stress parametrised with data from corso2020 [-]
sh = 3.765e-4;#sensibility to water stress
p_lcrit = -0.869;#min psiXil for stomatal opening [Mpa]
#influence of N contant, to reparametrise!, 
VcmaxrefChl1 = 1.28/2;#otherwise value too high, original: 1.31
VcmaxrefChl2 = 8.33/2; #otherwise value too high, original: 8,52
#for Vc, Vj, Ag, to reparametrise!
a1=4.; #g0+ fw[i] * a1 *( An[i] + Rd)/(ci[i] - deltagco2[i]);#tuzet2003
a3 = 1.7;#Jrefmax = Vcrefmax * a3 ;#Eq 25
alpha = 0.2; #or 0.44 , coefb = -(alpha * Qlight + Jmax);, alpha * Qlight * Jmax;
theta = 0.9;#or 0.67 coefa = theta;
gamma0 = 28e-6;  gamma1 = 0.0509;  gamma2 = 0.001;
g0 = 0.3e-3;#residual stomatal opening to CO2, Tuzet 2003 [mol CO2 m-2 s-1]


#Rostamza2020
alphaMin = 0.02
alphaMax = 0.95
thetaMin = 0.04
thetaMax = 1.0

namesVars = ['Chl','oi','fwr','sh','p_lcrit',
              'VcmaxrefChl1', 'VcmaxrefChl2',
             'a1','a3', 'alpha','theta',
             'gamma0', 'gamma1', 'gamma2','g0',
             'kr','kx']
namesLim = [[chl_*0.1, chl_],[oi*0.5,oi*2],[fwr*0.5,fwr*2],[sh*0.5,sh*2],
            [-1.25,0],
             [0, 10], [0, 10],
             [a1*0.5,a1*2],[a3*0.5,a3*2], [alphaMin,alphaMax],[thetaMin,thetaMax],
             [gamma0*0.5,gamma0*2], [gamma1*0.5,gamma1*2], [gamma2*0.5,gamma2*2],
            [g0*0.5,g0*2],
             [0.5,2],[0.5,2]]
problem = {
    'num_vars': len(namesVars),
    'names': namesVars,
    'bounds': namesLim
}

##### TO CHANGE ######
Nreps = 14
reps =2**Nreps
maxcore =  os.cpu_count()
######################

param_values = saltelli.sample(problem, reps)




Yvcs = np.array([])
Yvjs = np.array([])
Yevs = np.array([])


def SobolTest(myid,weatherX, X):
    pl = pb.MappedPlant(seednum = 2) 
    path = CPBdir+"/modelparameter/plant/"
    name = "Triticum_aestivum_adapted_2021"#"wheat_uqr15" #"manyleaves"## "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010

    pl.readParameters(path + name + ".xml")
    depth = 60
    sdf = pb.SDF_PlantBox(np.Inf, np.Inf, depth )
    pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil
    pl.initialize(verbose = True)#, stochastic = False)
    pl.simulate(21, False)#, "outputpm15.txt")

    picker = lambda x,y,z : max(int(np.floor(-z)),-1)   
    pl.setSoilGrid(picker)  # maps segment
    
    r = PhloemFluxPython(pl,psiXylInit =-1000,ciInit = weatherX["cs"]*0.5)
    plant =setDefaultVals(r,weatherX,1,1)
    krr = np.array([np.array(xi) for xi in plant.kr], dtype=object)
    kxx = np.array([np.array(xi) for xi in plant.kx], dtype=object)

    xXx = np.full((len(r.plant.nodes)),-1000.0)
    plant.Chl = np.array([X[0]])
    plant.oi= X[1]
    plant.fwr = X[2]
    plant.sh = X[3]
    plant.p_lcrit = X[4]
    plant.VcmaxrefChl1 = X[5]
    plant.VcmaxrefChl2 = X[6]
    plant.a1= X[7]
    plant.a3= X[8]
    plant.alpha= X[9]
    plant.theta= X[10]
    plant.gamma0 = X[11]
    plant.gamma1 = X[12]
    plant.gamma2 = X[13]
    plant.g0 = X[14]
    plant.kr = krr*X[15];plant.kx = kxx*X[16]
    
    plant.solve_photosynthesis(sim_time_ = 21, sxx_=xXx, 
                               cells_ = True,RH_ = weatherX["RH"],
            verbose_ = False, doLog_ = False,TairC_= weatherX["TairC"] )
    #print("in SobolTest",myid, X)
    npLB = np.array(plant.ci )
    Yvc = sum(np.array(plant.Vc)[np.where(npLB >0)])
    Yvj = sum(np.array(plant.Vj)[np.where(npLB >0)])
    Ev = np.mean(np.array(plant.Ev)[np.where(npLB >0)])
    #print("finish SobolTest",myid)
    return Yvc,Yvj,Ev

didDo = 0
parallelizer = Parallel(n_jobs= max(1,maxcore - 1))
while didDo < len(param_values):
    n_jobs_ = min(maxcore - 1,len(param_values) -  didDo)
    print("didDo",didDo, ", to do",len(param_values), ", will do",n_jobs_,"n_job",max(1,maxcore - 1))
    tasks_iterator = (delayed(SobolTest)
                            (i,weatherX, param_values[i])
                        for i in range(n_jobs_))
    results = parallelizer(tasks_iterator)
    Yvcs = np.concatenate((Yvcs,[smalldic[0] for smalldic in results ]))
    Yvjs = np.concatenate((Yvjs,[smalldic[1] for smalldic in results ]))
    Yevs = np.concatenate((Yevs,[smalldic[2] for smalldic in results ]))
    didDo += n_jobs_
    
#print("all outputs model")    
print("testLengths",reps,param_values.shape,len(Yvcs),len(Yvjs),len(Yevs))

with open('YAlls'+'.pkl','wb') as f:
     pickle.dump({'Yvcs':Yvcs,
                 'Yvjs':Yvjs,
                 'Yevs':Yevs},f, protocol=pickle.HIGHEST_PROTOCOL)
        
Sic = sobol.analyze(problem, Yvcs)
Sij = sobol.analyze(problem, Yvjs)
Siev = sobol.analyze(problem, Yevs)

print(Sic['S1'],Sij['S1'],Siev['S1'])
print(Sic['ST'],Sij['ST'],Siev['ST'])
print();print("ST_conf/ST")
print(Sic['ST_conf']/Sic['ST'],Sij['ST_conf']/Sij['ST'],Siev['ST_conf']/Siev['ST'])
print();print()

with open('Sic_R'+repr(Nreps)+'.pkl','wb') as f:
     pickle.dump(Sic,f, protocol=pickle.HIGHEST_PROTOCOL)
with open("Sic_R"+repr(Nreps)+".txt", 'w') as f: 
    for key, value in Sic.items(): 
        f.write('%s:%s\n' % (key, value))
        
with open('Sij_R'+repr(Nreps)+'.pkl','wb') as f:
     pickle.dump(Sij,f, protocol=pickle.HIGHEST_PROTOCOL)
with open("Sij_R"+repr(Nreps)+".txt", 'w') as f: 
    for key, value in Sij.items(): 
        f.write('%s:%s\n' % (key, value))
        
with open('Siev_R'+repr(Nreps)+'.pkl','wb') as f:
     pickle.dump(Siev,f, protocol=pickle.HIGHEST_PROTOCOL)
with open("Siev_R"+repr(Nreps)+".txt", 'w') as f: 
    for key, value in Siev.items(): 
        f.write('%s:%s\n' % (key, value))