directoryN = "/ipynbmorris/"
import sys; 
CPBdir = "../../../../cpb2705/StarchDynamicNew/CPlantBox"
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
from functional.phloem_flux import PhloemFluxPython 
import pickle
from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.test_functions import Ishigami
from joblib import Parallel, delayed
import os 
from helpUqr import *

def runallSobol(simDuration, condition,directoryN):
    
    weatherX = weather(simDuration, condition, 1)
    assert weatherX['Qlight']>0

    #			default value of env variable, to re-set at runtime
    cs = 390e-6; 

    Chl = 55,; 
    oi = 210e-3;#leaf internal [O2] [mol mol-1]

    SPAD= 100.0
    chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10


    Q10 = 2.; #double TrefQ10 = 20;//to compute effect of T on growth (see CN-wheat, residual respiration @Barillot 2016, appendix)
    #psiMax = 0;  
    psiMin = -2000*(1/0.9806806);#limit wat. pot. in xylem for water-limited growth, [cm]
    KMfu = 0.2; #@see C_fluxes,Michaelis menten coef for active sucrose usage
    Vmaxloading = 0.019872;#mmol cm-1 d-1 for leaf blade 1cm wide
    CSTimin = 0.4;#minimum CST value below which there is no sink of sucrose
    #ALSO target value for starch
    beta_loading = 1;#@see C_fluxes, feedback effect of C_ST on Q_FL
    Mloading = 0.2;#@see C_fluxes,Michaelis menten coef for Fl
    Gr_Y = 0.75;#growth efficiency
    psi_osmo = -4000*1.0197
    delta_osmo_min = -psi_osmo - psiMin
    
    #kg1 shouls stay between 0.9/(1-0.9) and 0.1/(1-0.1)



    namesVars_phlo = ['Q10','delta_osmo_min','psiMin','KMfu','Vmaxloading','CSTimin',
                  'beta_loading', 'Mloading','Gr_Y','rhoSucrose',
                 'kx_st','kr_st','Across','krm2','krm1','rmax',"k_S_", "leafGrowthZone"]
    namesLim_phlo = [Q10,delta_osmo_min,psiMin,KMfu,Vmaxloading,
                     CSTimin,beta_loading,Mloading,Gr_Y,1,
                 1,1,1,1,1,1,1,1 ]
    
    namesVars_xyl = ['Chl','oi','fw1r','k_fw1','psi_t,crit,1',
                  'kchl1', 'kchl2','kg1','kjmax', 'alpha',
                     'omega',#'gamma0', 'gamma1', 'gamma2',
                     'gm',
                 'kx_x','kr_x', "l_kr"]
    namesLim_xyl = [[0.5,2],[0.5,2],[0.5,2],[0.5,2],[0.5,2],
                     [0.5,2],[0.5,2],[0.5,2],[0.5,2],[0.02/0.4,0.95/0.4],
                 [0.04/0.6,0.98/0.6],#[0.5,2],[0.5,2],[0.5,2],
                    [0.5,2],[0.5,2] ,[0.5,2],[0.5,2]]
    problem = {
        'num_vars': len(namesVars_xyl),
        'names': namesVars_xyl,
        'bounds': namesLim_xyl
    }
    namesVars = namesVars_xyl + namesVars_phlo

    ##### TO CHANGE ######
    Nreps = 12
    if (simDuration > 10) and (condition == "dry"):
        Nreps = 12
    reps =2**Nreps
    maxcore = os.cpu_count()
    print('nreps, cpu',Nreps, maxcore)
    ######################

    param_values = saltelli.sample(problem, reps)
    #print(len(param_values))
    
    #print(param_values[13])
    #assert False




    Yev  = np.array([])
    Yag = np.array([])


    didDo = 0
    parallelizer = Parallel(n_jobs= max(1,maxcore - 1))
    while didDo < len(param_values):
        n_jobs_ = min(maxcore - 1,len(param_values) -  didDo)
        print("didDo",didDo, ", to do",len(param_values), ", will do",n_jobs_,"n_job",max(1,maxcore - 1),len(param_values))
        tasks_iterator = (delayed(SobolTest)
                                (i,weatherX, 
                                 namesVars, #all the var names
                                 list(param_values[didDo + i])+list(namesLim_phlo), #all the var values
                                 simDuration, condition,"Xylem",directoryN)
                            for i in range(n_jobs_))
        results = parallelizer(tasks_iterator)
        #Yexuds = np.concatenate((Yexuds,[smalldic[0] for smalldic in results ]))
        #Ygrs = np.concatenate((Ygrs,[smalldic[1] for smalldic in results ]))
        #Yrms = np.concatenate((Yrms,[smalldic[2] for smalldic in results ]))
        Yev = np.concatenate((Yev,[smalldic[0] for smalldic in results ]))
        Yag = np.concatenate((Yag,[smalldic[1] for smalldic in results ]))
        didDo += n_jobs_

    #print("all outputs model")    
    
    addToName = repr(Nreps)+"_"+repr(simDuration)+"_"+condition

    with open(directoryN+'YAlls_R'+addToName+'.pkl','wb') as f:
         pickle.dump({'Yev':Yev,
                     'Yag':Yag},f, protocol=pickle.HIGHEST_PROTOCOL)

    Siev = sobol.analyze(problem, Yev)
    Siag = sobol.analyze(problem, Yag)
    #Sirm = sobol.analyze(problem, Yrms)

    print();print("S1 and ST")
    print(Siev['S1'],Siag['S1'])
    print(Siev['ST'],Siag['ST'])
    print();print("ST_conf/ST")
    print(Siev['ST_conf']/Siev['ST'],Siag['ST_conf']/Siag['ST'])
    print();print()

    with open(directoryN+'Siev_R'+addToName+'.pkl','wb') as f:
         pickle.dump(Siev,f, protocol=pickle.HIGHEST_PROTOCOL)
    with open(directoryN+"Siev_R"+addToName+".txt", 'w') as f: 
        for key, value in Siev.items(): 
            f.write('%s:%s\n' % (key, value))

    with open(directoryN+'Siag_R'+addToName+'.pkl','wb') as f:
         pickle.dump(Siag,f, protocol=pickle.HIGHEST_PROTOCOL)
    with open(directoryN+"Siag_R"+addToName+".txt", 'w') as f: 
        for key, value in Siag.items(): 
            f.write('%s:%s\n' % (key, value))




if __name__ == '__main__':
    
    simDuration = int(sys.argv[1])
    condition = sys.argv[2]
    jobID = sys.argv[3]
    assert condition=="dry" or condition == "wet"
    print(simDuration, condition, jobID)
    directoryN = "results/"+os.path.basename(__file__)[:-3]+"/"+str(simDuration)+condition+jobID+"/"
    
    main_dir=os.environ['PWD']#dir of the file
    results_dir = main_dir +"/"+directoryN
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    else:
        import shutil
        shutil.rmtree(results_dir)
        os.makedirs(results_dir)
    runallSobol(simDuration, condition,directoryN)
    
    #simDuration = int(sys.argv[1])
    #condition = sys.argv[2]
    #assert condition=="dry" or condition == "wet"
    #assert simDuration== 7 or simDuration == 21
    #print(simDuration, condition)
    #runallSobol(simDuration, condition)