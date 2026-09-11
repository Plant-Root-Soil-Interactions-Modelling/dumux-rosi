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

conc_volume = {}
for i in range(0, len(soiltype)): 
    for m in range(0, len(sorption)): 
        for n in range(0, len(diffusion)):  

            scenario = soiltype[i]+'_diffusion'+diffusion[n]+'_sorption'+sorption[m]+'/'
            
            #macro 
            res = 4
            
            with open(path2file + scenario + "WaterC_macro.csv") as f:
                wc = [list(map(float, row)) for row in csv.reader(f)]
            max_len = max(len(row) for row in wc)
            wc_ = np.array([row + [np.nan]*(max_len-len(row)) for row in wc])  
            watvol_macro = wc_ * (res**3)        
            
            with open(path2file + scenario + "TotC_macro.csv") as f:
                totC = [list(map(float, row)) for row in csv.reader(f)]
            totC_ = np.array([row + [np.nan]*(max_len-len(row)) for row in totC])    
            conc_macro = np.divide(totC_ , watvol_macro)

            time = np.loadtxt(path2file + scenario + "time.txt", delimiter=",")[:-1,0]


            #micro
            path_cyl = path2file+"/"+scenario+"cyl_val/"
            folder = Path(path_cyl)
            files = list(folder.glob("*time*.txt"))
            counts = np.sort([int(re.search(r'time_(\d+)', str(f)).group(1)) for f in files])
            
            conc_list = []
            watvol_list = []
            for o in counts: 
                with open(path_cyl + 'Cyl_cellVol_'+str(o)+".txt") as f:
                    cellvol_ = [list(map(float, row)) for row in csv.reader(f)]
                # cellvol_ = np.genfromtxt(path_cyl + 'Cyl_cellVol_'+str(o)+".txt", delimiter=',')
                max_len = max(len(row) for row in cellvol_)
                cellvol = np.array([row + [np.nan]*(max_len-len(row)) for row in cellvol_])    
                
                with open(path_cyl + 'Cyl_watercontent_'+str(o)+".txt") as f:
                    wc_ = [list(map(float, row)) for row in csv.reader(f)]
                # wc_ = np.genfromtxt(path_cyl + 'Cyl_watercontent_'+str(o)+".txt", delimiter=',')
                wc = np.array([row + [np.nan]*(max_len-len(row)) for row in wc_])  
                wc = wc[:len(cellvol)]
                cellvol = cellvol[:len(wc)]
                watvol_micro_ = np.multiply(wc,cellvol)
                
                with open(path_cyl + 'Cyl_content1_'+str(o)+".txt") as f:
                    totC_ = [list(map(float, row)) for row in csv.reader(f)]
                # totC_ = np.genfromtxt(path_cyl + 'Cyl_content1_'+str(o)+".txt", delimiter=',')
                totC = np.array([row + [np.nan]*(max_len-len(row)) for row in totC_])  
                conc_micro_ = np.divide(totC[:len(watvol_micro_)],watvol_micro_[:len(totC)])
                
                if o == 0: 
                    max_rows = np.shape(cellvol)[0]
                else: 
                    conc_micro_ = pad_rows_top(conc_micro_, max_rows)
                    watvol_micro_ = pad_rows_top(watvol_micro_, max_rows)
                conc_list.append(conc_micro_)
                watvol_list.append(watvol_micro_)
                

            conc_micro = np.hstack(conc_list)
            watvol_micro = np.hstack(watvol_list)
            print(scenario)     
            print("conc_macro", conc_macro.shape)
            print("conc_micro", conc_micro.shape)
            print("watvol_micro", watvol_micro.shape)
            print(conc_micro.nbytes / 1024**3, "GB")
            if np.shape(conc_macro)[0]>np.shape(conc_micro)[0]: 
                conc_macro = conc_macro[:np.shape(conc_micro)[0]]
                watvol_macro = watvol_macro[:np.shape(watvol_micro)[0]]
            conc_ = np.hstack([conc_macro, conc_micro])
            watvol_ = np.hstack([watvol_macro, watvol_micro])
            idx = np.argsort(-conc_, axis = 1) 
            conc = np.take_along_axis(conc_, idx, axis = 1)
            watvol = np.take_along_axis(watvol_, idx, axis = 1)
            watvol_cum = np.cumsum(watvol, axis=1)

            time_ = []
            conc_ = []
            watvol_cum_ = []
            for p in range(0, np.shape(conc)[0]): 
                if np.around(int(time[p] *1000)/1000-int(time[p]),2) == 0.5 : #only plot
                    time_.append(time[p])
                    conc_.append(np.array(conc[p,:])) 
                    watvol_cum_.append(np.array(watvol_cum[p,:])) 

            time_ = np.array(time_) 
            concentration = np.array(conc_) 
            cum_watvolume = np.array(watvol_cum_) 
            conc_volume[(i, m, n)] = {
                "time": time_,
                "concentration": concentration,
                "cum_watvolume": cum_watvolume
            }
with open("concentration_volume.pkl", "wb") as f:
    pickle.dump(conc_volume, f)   
    
    
