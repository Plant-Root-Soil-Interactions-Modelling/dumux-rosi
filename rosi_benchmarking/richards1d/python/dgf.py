import numpy as np

def createDGF_1D(filename, N, depth, top, bot, domainId, vertex3d = False):

    z_ = np.linspace(0,-depth,N)
    initial = np.linspace(top,bot,N) # per node
    initialC = np.linspace(top,bot,N-1) # per cell    
    id = range(0,N)               
        
    file = open(filename,"w") 
    
    file.write("DGF\nVertex\n") 
    file.write('parameters 2\n') # initial data, domain index 
    for i in range(0,N):
        if vertex3d: 
            file.write('{:g} {:g} {:g} {:g} {:g}\n'.format(np.zeros(N,), np.zeros(N,), z_[i], initialC[i], domainId[i]))
        else:        
            file.write('{:g} {:g} {:g}\n'.format(z_[i], initial[i], domainId[i]))
             
    file.write('#\nSimplex\n') 
    file.write('parameters 2\n') # initial data, domain index 
    for i in range(0,N-1):
        file.write('{:g} {:g} {:g} {:g}\n'.format(id[i], id[i+1], initialC[i], domainId[i]));
                   
    file.write('#\nBOUNDARYSEGMENTS\n') # not used... 
    file.write('2 0\n')            
    file.write('3 {:g}\n'.format(N-1)) # vertex id, index starts with 0
    file.write('#\nBOUNDARYDOMAIN\ndefault 1\n#\n')
    file.close() 





if __name__ == "__main__":
    
    # Jan 1 (Figure 2abc)
    domain_jan1 = np.hstack((np.ones(50,), 2*np.ones(151,)))
    createDGF_1D("jan1.dgf",201,2.,-200,-200.,domain_jan1)
    
    # Jan 2 (Figure 3)
    domain_jan2 = np.ones(55,)
    createDGF_1D("jan2.dgf",55,.54,-0,-54,domain_jan2)
    domain_jan2 = np.ones(109,)  # same with resolution of 0.5 cm
    createDGF_1D("jan2b.dgf",109,.54,-0,-54,domain_jan2)
    domain_jan2 = np.ones(217,)  # same with resolution of 0.25 cm
    createDGF_1D("jan2c.dgf",217,.54,-0,-54,domain_jan2)
    
    # Jan 3 (Figure 4)
    domain_jan1 = np.ones(201,)
    createDGF_1D("jan3.dgf",201,2.,-400,-400.,domain_jan1)

    print("its done.")