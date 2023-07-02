directoryN = "/ipynbmorris/"
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

def runallSobol(simDuration, condition, power2):
    
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
    #power2 = 6
    reps =2**power2
    maxcore =  os.cpu_count()
    print('power2, cpu',power2, maxcore)
    ######################

    param_values = saltelli.sample(problem, reps)




    #Yexuds = np.array([])
    #Ygrs  = np.array([])
    #Yrms = np.array([])
    #Ytrans = np.array([])
    #Yassi  = np.array([])
    #Ywue = np.array([])
    
    Yall = [np.array([]) for i in range(7)]#Yexuds,Ygrs,Yrms,Ytrans,Yassi,Ywue]
    namesY = ["Yexuds","Ygrs","Yrms","Ytrans","Yassi","Ywue","YIwue"]
    allLoops = np.array([])


    didDo = 0
    parallelizer = Parallel(n_jobs= max(1,maxcore - 1))
        
        
    while didDo < len(param_values):
        n_jobs_ = min(maxcore - 1,len(param_values) -  didDo)
        print("didDo",didDo, ", to do",len(param_values), ", will do",n_jobs_,"n_job",max(1,maxcore - 1),len(param_values))
        if True:#idtorun == -1:
            tasks_iterator = (delayed(SobolTest)
                                    (i,weatherX, 
                                     namesVars,  list(param_values[didDo + i]), #all the var values
                                     simDuration,condition,"SEB",power2)
                                for i in range(n_jobs_))
            results = parallelizer(tasks_iterator)
        else:
            SobolTest(idtorun,weatherX, 
                     namesVars,  list(param_values[idtorun]), #all the var values
                     simDuration,condition,"SEB")
        #Yexuds = np.concatenate((Yexuds,[smalldic[0] for smalldic in results ]))
        #Ygrs = np.concatenate((Ygrs,[smalldic[1] for smalldic in results ]))
        #Yrms = np.concatenate((Yrms,[smalldic[2] for smalldic in results ]))
        #Ytrans = np.concatenate((Ytrans,[smalldic[3] for smalldic in results ]))
        #Yassi = np.concatenate((Yrms,[smalldic[4] for smalldic in results ]))
        #Ywue = np.concatenate((Ywue,[smalldic[5] for smalldic in results ]))
        
        
        for j in range(len(Yall)):
            Yall[j] = np.concatenate((Yall[j] ,[smalldic[j] for smalldic in results ]))
        allLoops   = np.concatenate((allLoops,[smalldic[7] for smalldic in results ])) 
        didDo += n_jobs_

    
    addToName = repr(power2)+"_"+repr(simDuration)+"_"+condition
    
    test_values = Yall.copy()
    DictVal = {}
    for key in namesY:
        for value in test_values:
            DictVal[key] = value
            test_values.remove(value)
            break

    print("finished sobol")
    with open('results/sobolSEB/'+str(power2)+condition+str(simDuration)+"/"+'YAlls_R'+addToName+'.pkl','wb') as f:
         pickle.dump(DictVal,f, protocol=pickle.HIGHEST_PROTOCOL)
    print("finished sobolWRITING")
    
            
    for j in range(7):#, Yout in enumerate(Yall):
        
        Yout = Yall[j]
        print(namesY[j])
        Si = sobol.analyze(problem, Yout)
        print(Si['S1'])
        print(Si['ST'])
        print("ST_conf/ST")
        print(Si['ST_conf']/Si['ST'])
        print();print()

        with open('results/sobolSEB/'+str(power2)+condition+str(simDuration)+"/"+namesY[j]+'_R'+addToName+'.pkl','wb') as f:
             pickle.dump(Si,f, protocol=pickle.HIGHEST_PROTOCOL)
        with open('results/sobolSEB/'+str(power2)+condition+str(simDuration)+"/"+namesY[j]+'_R'+addToName+".txt", 'w') as f: 
            for key, value in Si.items(): 
                f.write('%s:%s\n' % (key, value))
    with open('results/sobolSEB/'+str(power2)+condition+str(simDuration)+'/loops.txt', 'w') as log:
        log.write(','.join([num for num in map(str, allLoops)])  +'\n')


simDuration = int(sys.argv[1])
condition = sys.argv[2]
if len(sys.argv) > 3:
    power2= int( sys.argv[3])
else:
    power2= 6
    
assert condition=="dry" or condition == "wet"
print(simDuration, condition, power2)

main_dir=os.environ['PWD']#dir of the file
results_dir = main_dir +"/results/sobolSEB/"+str(power2)+condition+str(simDuration)+"/"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    import shutil
    shutil.rmtree(results_dir)
    os.makedirs(results_dir)
runallSobol(simDuration, condition, power2)