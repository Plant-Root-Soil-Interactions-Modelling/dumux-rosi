"""
    print simulation timings 
"""
import sys; sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");

from scenario_setup import *

def make_list():
    jobs = []

    method = ['sra', "agg", "par"]  
    plant = ['springbarley', 'maize']  
    dim = ["1D", "3D"] 
    soil = ['hydrus_loam', 'hydrus_clay', 'hydrus_sandyloam']  
    outer_radius = ['voronoi'] # 'length', 'surface', 'volume', 

    print("Creating", len(method) * len(plant) * len(dim) * len(soil) * len(outer_radius), "simulations")
    print()

    for m in method:
        for p in plant:
            for d in dim:
                for s in soil:
                    for o in outer_radius:
                        jobs.append([m, p, d, s, o])
    return jobs

if __name__ == "__main__":
    
    job_list = make_list()
    names, timings = print_timings(job_list)
    print()
    c = 0
    for l in range(0,3):
        for i, p in enumerate(['springbarley', 'maize']):
            for j, d in enumerate(["1D", "3D"]): # dimesnions
                for k in range(0,3): # soils                   
                        base_index = 6*i+3*j+k
                        print(names[c+base_index], ": ", 100*(1/(timings[c+base_index]/timings[base_index]) -1) )
            print() 
        c+=12
        
    

