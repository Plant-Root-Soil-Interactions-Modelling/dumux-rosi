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

def runallSobol(simDuration, condition):
    if condition =="dry":
        Tmin = 20.7; Tmax = 30.27
        specificHumidity = 0.0111
        Pair = 1070.00 #hPa
        thetaInit = 20/100
    else:
        Tmin = 15.8; Tmax = 22
        specificHumidity = 0.0097
        Pair = 1010.00 #hPa
        thetaInit = 30/100
    weatherX = weather(0, Tmin,Tmax,specificHumidity,Pair,thetaInit)
    assert weatherX['Qlight']>0

    #			default value of env variable, to re-set at runtime
    cs = 350e-6; #example from Dewar2002 [mol mol-1] 

    Chl = 55,; 
    oi = 210e-3;#leaf internal [O2] [mol mol-1]

    SPAD= 100.0
    chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10


    Q10 = 2.; #double TrefQ10 = 20;//to compute effect of T on growth (see CN-wheat, residual respiration @Barillot 2016, appendix)
    psiMax = 0;  psiMin = -2000*(1/0.9806806);#limit wat. pot. in xylem for water-limited growth, [cm]
    KMfu = 0.2; #@see C_fluxes,Michaelis menten coef for active sucrose usage
    Vmaxloading = 0.019872;#mmol cm-1 d-1 for leaf blade 1cm wide
    CSTimin = 0.4;#minimum CST value below which there is no sink of sucrose
    beta_loading = 1;#@see C_fluxes, feedback effect of C_ST on Q_FL
    Mloading = 0.2;#@see C_fluxes,Michaelis menten coef for Fl
    Gr_Y = 0.75;#growth efficiency



    namesVars = ['Q10','psiMax','psiMin','KMfu','Vmaxloading','CSTimin',
                  'beta_loading', 'Mloading','Gr_Y','rhoSucrose',
                 'kx','kr','Across','krm2','krm1','rmax']
    namesLim = [[Q10*0.5, Q10*2],[psiMax + psiMin/2+100,psiMax],[psiMin*2,psiMin/2],[KMfu*0.5,KMfu*2],[Vmaxloading*0.5,Vmaxloading*2],
                [CSTimin*0.5,CSTimin*2],
                 [beta_loading*0.5, beta_loading*2], [Mloading/2,Mloading*2],
                 [Gr_Y*0.5,Gr_Y*2],[0.5,2],
                 [0.5,2],[0.5,2],[0.5,2],[0.5,2],[0.5,2],[0.5,2]]
    problem = {
        'num_vars': len(namesVars),
        'names': namesVars,
        'bounds': namesLim
    }

    ##### TO CHANGE ######
    Nreps = 9
    reps =2**Nreps
    maxcore =  os.cpu_count()
    print('nreps, cpu',Nreps, maxcore)
    ######################

    param_values = saltelli.sample(problem, reps)




    Yexuds = np.array([])
    Ygrs  = np.array([])
    Yrms = np.array([])


    def SobolTest(myid,weatherX, X, simDuration):
        pl = pb.MappedPlant(seednum = 2) 
        path = CPBdir+"/modelparameter/plant/"
        name = "Triticum_aestivum_adapted_2021"#"wheat_uqr15" #"manyleaves"## "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010

        pl.readParameters(path + name + ".xml")
        depth = 60
        sdf = pb.SDF_PlantBox(np.Inf, np.Inf, depth )
        pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil
        pl.initialize(verbose = True)#, stochastic = False)
        pl.simulate(simDuration, False)#, "outputpm15.txt")

        picker = lambda x,y,z : max(int(np.floor(-z)),-1)   
        pl.setSoilGrid(picker)  # maps segment

        r = PhloemFluxPython(pl,psiXylInit =-1000,ciInit = weatherX["cs"]*0.5)
        plant =setDefaultVals(r,weatherX,1,1,X[10],X[11],X[12],X[9], X[13],X[14],X[15])

        xXx = np.full((len(r.plant.nodes)),weatherX["p_mean"])
        plant.Q10 = X[0]
        plant.psiMax= X[1]
        plant.psiMin = X[2]
        plant.Vmaxloading = X[3]
        plant.CSTimin = X[4]
        plant.beta_loading = X[5]
        plant.Mloading = X[6]
        plant.Gr_Y= X[7]
        
        verbose_phloem = True
        filename = "sobolPhoto4/inPM"+repr(myid)+".txt"
        #plant.doTroubleshooting=True
        try:

            plant.solve_photosynthesis(sim_time_ = 21, sxx_=xXx, 
                                       cells_ = True,RH_ = weatherX["RH"],
                    verbose_ = False, doLog_ = False,TairC_= weatherX["TairC"] )

            plant.startPM(simDuration, simDuration + 1/(24), 1, ( weatherX["TairC"]  +273.15) , verbose_phloem, filename)
            Nt = len(plant.plant.nodes) 
            Q_Rm    = np.array(plant.Q_out[(Nt*2):(Nt*3)])
            Q_Exud  = np.array(plant.Q_out[(Nt*3):(Nt*4)])
            Q_Gr    = np.array(plant.Q_out[(Nt*4):(Nt*5)])
        except:
            print("error at my id ",myid)
            raise Exception
        del plant
        del pl
        Yexud = sum(Q_Exud)
        Ygr = sum(Q_Gr)
        Ev = sum(Q_Rm)
        
        return Yexud,Ygr,Ev

    didDo = 0
    parallelizer = Parallel(n_jobs= max(1,maxcore - 1))
    while didDo < len(param_values):
        n_jobs_ = min(maxcore - 1,len(param_values) -  didDo)
        print("didDo",didDo, ", to do",len(param_values), ", will do",n_jobs_,"n_job",max(1,maxcore - 1),len(param_values))
        tasks_iterator = (delayed(SobolTest)
                                (i,weatherX, param_values[didDo + i],simDuration)
                            for i in range(n_jobs_))
        results = parallelizer(tasks_iterator)
        Yexuds = np.concatenate((Yexuds,[smalldic[0] for smalldic in results ]))
        Ygrs = np.concatenate((Ygrs,[smalldic[1] for smalldic in results ]))
        Yrms = np.concatenate((Yrms,[smalldic[2] for smalldic in results ]))
        didDo += n_jobs_

    #print("all outputs model")    
    
    addToName = repr(Nreps)+"_"+repr(simDuration)+"_"+condition

    with open('sobolPhlo4/YAlls_R'+addToName+'.pkl','wb') as f:
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

    with open('sobolPhlo4/Siex_R'+addToName+'.pkl','wb') as f:
         pickle.dump(Siex,f, protocol=pickle.HIGHEST_PROTOCOL)
    with open("sobolPhlo4/Siex_R"+addToName+".txt", 'w') as f: 
        for key, value in Siex.items(): 
            f.write('%s:%s\n' % (key, value))

    with open('sobolPhlo4/Sigr_R'+addToName+'.pkl','wb') as f:
         pickle.dump(Sigr,f, protocol=pickle.HIGHEST_PROTOCOL)
    with open("sobolPhlo4/Sigr_R"+addToName+".txt", 'w') as f: 
        for key, value in Sigr.items(): 
            f.write('%s:%s\n' % (key, value))

    with open('sobolPhlo4/Sirm_R'+addToName+'.pkl','wb') as f:
         pickle.dump(Sirm,f, protocol=pickle.HIGHEST_PROTOCOL)
    with open("sobolPhlo4/Sirm_R"+addToName+".txt", 'w') as f: 
        for key, value in Sirm.items(): 
            f.write('%s:%s\n' % (key, value))

simDuration = int(sys.argv[1])
condition = sys.argv[2]
assert condition=="dry" or condition == "wet"
print(simDuration, condition)
runallSobol(simDuration, condition)