import random
import numpy as np

for j in range(0,200):
    
    res = np.zeros(6,)
    pickable = list(range(1,10))
    for i in range(0,6):
        rn = random.randint(0,len(pickable)-1)
        res[i] = pickable[rn]
        pickable.pop(rn)
        
    print(res)
