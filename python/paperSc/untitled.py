
import sys;
sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../modules/");
sys.path.append("/data");
sys.path.append("../../../CPlantBox/src/python_modules");
sys.path.append("../../../CPlantBox/src/functional/");
sys.path.append("../../../CPlantBox/src/rsml/");
sys.path.append("../../../CPlantBox/src/visualisation/")
sys.path.append("../../../CPlantBox/src/structural/")
sys.path.append("../../../CPlantBox/src/external/")
sys.path.append("../../../CPlantBox/");
sys.path.append("../../../CPlantBox/src");
import os
import numpy as np
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

import plantbox as pb  # CPlantBox
import van_genuchten as vg
#import evapotranspiration as evap
from functional.xylem_flux import *
from datetime import *
from functional.plant_conductivities import init_conductivities
from functional.phloem_flux import PhloemFluxPython  # root system Python hybrid solver

def write_file_float(name, data, directory_, allranks = False):
    if (rank == 0) or allranks:
        name2 = directory_+ name+ '.txt'
        with open(name2, 'a') as log:
            log.write(repr( data)  +'\n')
        
def write_file_array(name, data, space =",", directory_ ="./results/", fileType = '.txt', allranks = False ):
    np.array(data).reshape(-1)
    try:
        if (rank == 0) or allranks:
            name2 = directory_+ name+ fileType
            # print('write_file_array',name2)
            with open(name2, 'a') as log:
                log.write(space.join([num for num in map(str, data)])  +'\n')
    except:
        print(name, data,data.shape)
        raise Exception

results_dir="./results/testFiles/"

if not os.path.exists(results_dir):
    os.makedirs(results_dir)
    os.makedirs(results_dir+'cyl_val/')
else:
    test = os.listdir(results_dir)
    for item in test:
        try:
            os.remove(results_dir+item)
        except:
            pass


write_file_array("other_content", 
             np.array([1,1,1])*2, 
             directory_ =results_dir, 
                 allranks = True)
write_file_array("Cyl_content", 
             np.array([1,2,5])*2, 
             directory_ =results_dir+'cyl_val/', 
                 allranks = True)

rs = pb.MappedPlant()

r = PhloemFluxPython(rs,psiXylInit = -659.8,ciInit = 350e-6*0.5) 

#number of vascular bundles
VascBundle_leaf = 32
VascBundle_stem = 52
VascBundle_root = 1 #valid for all root type

#numPerBundle
numL = 18
numS = 21
numr0 = 33
numr1 = 25
numr2 = 25
numr3 = 1

#radius of phloem type^4 * number per bundle
rad_s_l   = numL* (0.00025 **4)# * 2; rad_x_l_2   = (0.0005 **4) * 2   
rad_s_s   = numS *(0.00019 **4) #* 3; rad_x_s_2   = (0.0008 **4) * 1     
rad_s_r0  = numr0 *(0.00039 **4) #* 4    
rad_s_r12 = numr1*(0.00035**4) #* 4; rad_x_r12_2 = (0.00087**4) * 1
rad_s_r3  = numr3 *(0.00068**4) #* 1      

# axial conductivity [cm^3/day] , mu is added later as it evolves with CST  
beta = 0.9 #Thompson 2003a
kz_l   = VascBundle_leaf * rad_s_l   * np.pi /8 * beta  
kz_s   = VascBundle_stem * rad_s_s   * np.pi /8 * beta
kz_r0  = VascBundle_root * rad_s_r0  * np.pi /8 * beta
kz_r12 = VascBundle_root * rad_s_r12 * np.pi /8 * beta
kz_r3  = VascBundle_root * rad_s_r3  * np.pi /8 * beta

#print([[kz_r0,kz_r12,kz_r12,kz_r3],[kz_s,kz_s ],[kz_l]])
#raise Exception
#radial conductivity [1/day],
kr_l  = 0.#3.83e-4 * hPa2cm# init: 3.83e-4 cm/d/hPa
kr_s  = 0.#1.e-20  * hPa2cm # set to almost 0
# kr_r0 = 1e-1
# kr_r1 = 1e-1
# kr_r2 = 1e-1
# kr_r3 = 1e-1
kr_r0 = 5e-2
kr_r1 = 5e-2
kr_r2 = 5e-2
kr_r3 = 5e-2
l_kr =  0.8 #cm

r.setKr_st([[kr_r0,kr_r1 ,kr_r2 ,kr_r0],[kr_s,kr_s ],[kr_l]] , kr_length_= l_kr, verbose = False)
r.setKx_st([[kz_r0,kz_r12,kz_r12,kz_r0],[kz_s,kz_s ],[kz_l]],  False)
