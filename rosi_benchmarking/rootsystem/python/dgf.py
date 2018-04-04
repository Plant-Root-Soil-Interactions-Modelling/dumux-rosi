import numpy as np

def createDGF_1D(filename, N, depth, top, bot, domainId, vertex3d = False):

    # prepare arrays
    z_ = np.linspace(0,-depth,N)
    initial = np.linspace(top,bot,N) # per node
    initialC = np.linspace(top,bot,N-1) # per cell    
    id = range(0,N)               
    
    # write file     
    file = open(filename,"w") 
 
    file.write("DGF\n") 
    file.write('Vertex\n')
    file.write('parameters 2\n'); 
    for i in range(0,N):
        if vertex3d: 
            file.write('{:g} {:g} {:g} {:g} {:g}\n'.format(0., 0., z_[i], initialC[i], domainId[i]))
        else:        
            file.write('{:g} {:g} {:g}\n'.format(z_[i], initial[i], domainId[i]))
             
    file.write('#\n');
    file.write('Simplex\n'); 
    file.write('parameters 2\n'); 
    for i in range(0,N-1):
        file.write('{:g} {:g} {:g} {:g}\n'.format(id[i], id[i+1], initialC[i], domainId[i]));
        
    # not used...        
    file.write('#\n')
    file.write('BOUNDARYSEGMENTS\n') # how do i get the boundary segments into DUMUX ?
    file.write('2 0\n')            
    file.write('3 {:g}\n'.format(N-1)) # vertex id, but index starts with 0
    file.write('#\n');
    file.write('BOUNDARYDOMAIN\n')
    file.write('default 1\n');
    file.write('#\n')

    file.close() 


def createDGF_1Droots(filename, nodes, seg):
        
    file = open(filename,"w") # write file 
 
    file.write("DGF\n") 
    file.write('Vertex\n')
    # file.write('parameters 2\n'); 
    for i in range(0,len(nodes)):
        file.write('{:g} {:g} {:g} \n'.format(nodes[i,0], nodes[i,1], nodes[i,2]))
             
    file.write('#\n');
    file.write('Simplex\n'); 
    # file.write('parameters 2\n'); 
    for i in range(0,len(seg)):
        file.write('{:g} {:g} \n'.format(seg[i,0], seg[i,1]));
        
    # not used...        
    file.write('#\n')
    file.write('BOUNDARYSEGMENTS\n') # how do i get the boundary segments into DUMUX ?
    file.write('2 0\n')            
    file.write('3 {:g}\n'.format(len(seg))) # vertex id, but index starts with 0
    file.write('#\n');
    file.write('BOUNDARYDOMAIN\n')
    file.write('default 1\n');
    file.write('#\n')

    file.close() 



if __name__ == "__main__":
    
#     # Jan 1 (Figure 2abc)
#     domain_jan1 = np.hstack((np.ones(50,), 2*np.ones(151,)))
#     createDGF_1D("jan1.dgf",201,2.,-200,-200.,domain_jan1)
#     #
#     domain_jan1hd = np.hstack((np.ones(500,), 2*np.ones(1501,)))
#     createDGF_1D("jan1hd.dgf",2001,2.,-200,-200.,domain_jan1hd)    
#     
#     # Jan 2 (Figure 3)
#     domain_jan2 = np.ones(55,)
#     createDGF_1D("jan2.dgf",55,.54,-0,-54,domain_jan2)
#     domain_jan2 = np.ones(109,)  # same with resolution of 0.5 cm
#     createDGF_1D("jan2b.dgf",109,.54,-0,-54,domain_jan2)
#     domain_jan2 = np.ones(217,)  # same with resolution of 0.25 cm
#     createDGF_1D("jan2c.dgf",217,.54,-0,-54,domain_jan2)
#     
#     # Jan 3 (Figure 4)
#     domain_jan1 = np.ones(201,)
#     createDGF_1D("jan3.dgf",201,2.,-400,-400.,domain_jan1)

  
    nnz = 100
    L = 0.5 # length of single straight root (m)

    # create grid
    nodes = np.zeros((nnz,3))
    seg = np.zeros(((nnz-1),2), dtype=int) 
    c = 0
    for i in range(1, nnz):
        seg[c,0] = i-1
        seg[c,1] = i
        c += 1    
        nodes[i,:] = [0.,0.,-i*L/(nnz-1)]  
        
    sn = len(seg)
    nn = len(nodes) 

    createDGF_1Droots("roots.dgf", nodes, seg)

    print("its done.")
    
    
    
    
    