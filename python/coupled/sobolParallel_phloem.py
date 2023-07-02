
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
    beta_loading = 1;#@see C_fluxes, feedback effect of C_ST on Q_FL
    Mloading = 0.2;#@see C_fluxes,Michaelis menten coef for Fl
    Gr_Y = 0.75;#growth efficiency
    psi_osmo = -10000*1.0197
    delta_osmo_min = -psi_osmo - psiMin



    namesVars_phlo = ["maxTil","maxB", "a",
                     "lmaxr", "rr", "lnr", 
                     "lmaxst", "rst", "lnst", #"ast",
                     "lmaxl", "rl", "Width_blade", 
                        'Q10',
                      'delta_osmo_min','psiMin',
                      'KMfu','Vmaxloading',
                  'beta_loading', 'Mloading','Gr_Y',
                      'rhoSucrose',
                     'k_st',
                      'Across',
                      'krm', 
                      "k_starch"]
    namesLim_phlo = [[0.5,2],[0.5,2],[0.5,2],
                     [0.5,2],[0.5,2],[0.5,2],#[0.5,2],
                     [0.5,2],[0.5,2],[0.5,2],
                     [0.5,2],[0.5,2],[0.5,2],
                    [Q10*0.5, Q10*2],
                     [delta_osmo_min/2,delta_osmo_min*2],[psiMin*2,psiMin/2],
                     [KMfu*0.5,KMfu*2], [Vmaxloading*0.5,Vmaxloading*2],
                     [beta_loading*0.5, beta_loading*2], [Mloading/2,Mloading*2],[Gr_Y*0.5,Gr_Y*2],
                     [0.5,2],
                     [0.5,2],
                     [0.5,2],
                     [0.5,2],
                     [0.5,2] ]
    
    
    
    namesVars_xyl = [ 'k_fw1','psi_t,crit,1',              #'oi', 'Chl','fw1r',
                     'kchl', 'kg1','kjmax', 'alpha',         #'kchl1', 'kchl2',
                     'omega', 'gm', 'k_x']
    namesLim_xyl = [1 for ii in namesVars_xyl]
    problem = {
        'num_vars': len(namesVars_phlo),
        'names': namesVars_phlo,
        'bounds': namesLim_phlo
    }
    namesVars = namesVars_xyl + namesVars_phlo

    ##### TO CHANGE ######
    Nreps = 11
    #if (simDuration > 12):# and (condition == "dry"):
    #    Nreps = 8
    reps =2**Nreps
    maxcore =  os.cpu_count()
    print('nreps, cpu',Nreps, maxcore)
    ######################

    param_values_phlo = saltelli.sample(problem, reps)




    Yexuds = np.array([])
    Ygrs  = np.array([])
    Yrms = np.array([])


    didDo = 0
    parallelizer = Parallel(n_jobs= max(1,maxcore - 1))
    while didDo < len(param_values_phlo):
        n_jobs_ = min(maxcore - 1,len(param_values_phlo) -  didDo)
        print("didDo",didDo, ", to do",len(param_values_phlo), ", will do",n_jobs_,"n_job",max(1,maxcore - 1),len(param_values_phlo))
        tasks_iterator = (delayed(SobolTest)
                                (i,weatherX, 
                                 namesVars, #all the var names
                                 list(namesLim_xyl)+ list(param_values_phlo[didDo + i]), #all the var values
                                 simDuration,condition,"Phloem",directoryN)
                            for i in range(n_jobs_))
        results = parallelizer(tasks_iterator)
        Yexuds = np.concatenate((Yexuds,[smalldic[0] for smalldic in results ]))
        Ygrs = np.concatenate((Ygrs,[smalldic[1] for smalldic in results ]))
        Yrms = np.concatenate((Yrms,[smalldic[2] for smalldic in results ]))
        didDo += n_jobs_

    #print("all outputs model")    
    
    addToName = repr(Nreps)+"_"+repr(simDuration)+"_"+condition

    #with open('sobolPhlo4/YAlls_R'+addToName+'.txt','w') as f:
     #    pickle.dump({'Yexuds':Yexuds,
      #               'Ygrs':Ygrs,
       #              'Yrms':Yrms},f, protocol=pickle.HIGHEST_PROTOCOL)
    with open(directoryN+'YAlls_R'+addToName+'.pkl','wb') as f:
         pickle.dump({'Yexuds':Yexuds,
                     'Ygrs':Ygrs,
                     'Yrms':Yrms},f, protocol=pickle.HIGHEST_PROTOCOL)

    Siex = sobol.analyze(problem, Yexuds)
    Sigr = sobol.analyze(problem, Ygrs)
    Sirm = sobol.analyze(problem, Yrms)

    print(Siex['S1'],Sigr['S1'],Sirm['S1'])
    print(Siex['ST'],Sigr['ST'],Sirm['ST'])
    print();print("ST_conf/ST")
    print(Siex['ST_conf']/Siex['ST'],Sigr['ST_conf']/Sigr['ST'],Sirm['ST_conf']/Sirm['ST'])
    print();print()

    with open(directoryN+'Siex_R'+addToName+'.pkl','wb') as f:
         pickle.dump(Siex,f, protocol=pickle.HIGHEST_PROTOCOL)
    with open(directoryN+"Siex_R"+addToName+".txt", 'w') as f: 
        for key, value in Siex.items(): 
            f.write('%s:%s\n' % (key, value))

    with open(directoryN+'Sigr_R'+addToName+'.pkl','wb') as f:
         pickle.dump(Sigr,f, protocol=pickle.HIGHEST_PROTOCOL)
    with open(directoryN+"Sigr_R"+addToName+".txt", 'w') as f: 
        for key, value in Sigr.items(): 
            f.write('%s:%s\n' % (key, value))

    with open(directoryN+'Sirm_R'+addToName+'.pkl','wb') as f:
         pickle.dump(Sirm,f, protocol=pickle.HIGHEST_PROTOCOL)
    with open(directoryN+"Sirm_R"+addToName+".txt", 'w') as f: 
        for key, value in Sirm.items(): 
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