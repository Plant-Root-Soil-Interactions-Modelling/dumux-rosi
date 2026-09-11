import sys
import numpy as np
from os import walk
import re
import glob
import os
import csv
from pathlib import Path
import pickle


path2file = "../results_field_mixedmodel_reduced/"


def get_num(filename): 
    regex = re.compile(r'\d+')
    return [int(x) for x in regex.findall(filename)]
    

def extract_number(filename):
    name = os.path.basename(filename)
    numbers = re.findall(r'\d+', name)
    return int(numbers[0]) if numbers else -1  
    
def pad_rows_top(arr, target_rows):
    if arr.shape[0] < target_rows:
        # Create a NaN array of missing rows
        padding = np.full((target_rows - arr.shape[0], arr.shape[1]), np.nan)
        return np.vstack([padding, arr])  # NaNs on top
    return arr
                           


soiltype = ['loam','sand']
diffusion = ['low', 'medium', 'mediumhigh', 'high']
sorption = ['low', 'medium', 'high']
target_day = 30

gradients = {}
for i in range(0, len(soiltype)): 
    for m in range(0, len(sorption)): 
        for n in range(0, len(diffusion)):  

            scenario = soiltype[i]+'_diffusion'+diffusion[n]+'_sorption'+sorption[m]+'/'
            time = np.loadtxt(path2file + scenario + "time.txt", delimiter=",")[:-1,0]
            time_ = target_day

            #micro
            path_cyl = path2file+"/"+scenario+"cyl_val/"
            folder = Path(path_cyl)
            files = glob.glob(os.path.join(folder, "*time*"))
            files_sorted = sorted(files, key=extract_number,reverse=True)
            filenames = [os.path.basename(f) for f in files_sorted]


            concentration = []
            age = []
            distance = []
            for o, file in enumerate(filenames): 
                num = extract_number(file)
                
                with open(path_cyl + 'Cyl_cellVol_'+str(num)+".txt") as f:
                    cellvol_ = [list(map(float, row)) for row in csv.reader(f)]
                max_len = max(len(row) for row in cellvol_)
                cellvol = np.array([row + [np.nan]*(max_len-len(row)) for row in cellvol_])    
                
                with open(path_cyl + 'Cyl_watercontent_'+str(num)+".txt") as f:
                    wc_ = [list(map(float, row)) for row in csv.reader(f)]
                wc = np.array([row + [np.nan]*(max_len-len(row)) for row in wc_])  
                watvol_micro_ = np.multiply(wc[:len(cellvol)],cellvol[:len(wc)])
                
                with open(path_cyl + 'Cyl_content1_'+str(num)+".txt") as f:
                    totC_ = [list(map(float, row)) for row in csv.reader(f)]
                totC = np.array([row + [np.nan]*(max_len-len(row)) for row in totC_])  
                conc = np.divide(totC[:len(watvol_micro_)],watvol_micro_[:len(totC)]) #per cm^3 water
                
                with open(path_cyl + 'Cyl_coord_'+str(num)+".txt") as f:
                    coord_ = [list(map(float, row)) for row in csv.reader(f)]
                dist = np.array([row + [np.nan]*(max_len-len(row)) for row in coord_])  
                
                time_cyl = np.array(np.loadtxt(path_cyl + 'Cyl_time_'+str(num)+".txt"))
                if isinstance(time_cyl, np.ndarray):
                    if not time_cyl.ndim > 0 and time_cyl.size > 0:
                        time_cyl = np.array([time_cyl])
                        
                        
                age_single = time_-time_cyl[0]
                if age_single<=0: 
                    continue
                else: 
                    target_idx = np.abs(time_cyl - time_).argmin()
                    distance_ = (dist[target_idx,:]-dist[target_idx,0])*10 #from cm to mm 
                    concentration_ = conc[target_idx,:] #mol
                    age_ = (np.ones((len(distance_)))*age_single).T
                
                
                    distance.append(np.array(distance_))
                    concentration.append(np.array(concentration_)) 
                    age.append(np.array(age_)) 

            distance = np.array(distance) 
            concentration = np.array(concentration) 
            age = np.array(age) 
            gradients[(i, m, n)] = {
                "distance": distance,
                "concentration": concentration,
                "age": age
            }
with open("gradients.pkl", "wb") as f:
    pickle.dump(gradients, f)   
    
    