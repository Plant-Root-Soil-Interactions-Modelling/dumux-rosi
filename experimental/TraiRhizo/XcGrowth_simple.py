import sys; sys.path.append("./modules_fpit"); 
sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")


from rosi_richards10c import RichardsNCSPnum as RichardsNCSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import numpy as np
import scenario_setup as scenario
from cyl3plant_simple import simulate_const

if __name__ == '__main__':

    min_b = np.array([-3./2, -12./2, -40.])
    cell_number =np.array( [3,12,40]) # 1cm3 
    max_b =np.array( [3./2, 12./2, 0.])
       
    results_dir="./results/"
    
    """
     Cylindric rhizosphere models, C exudate from finite volumes
    """

    fname = 'L_WT_phenolics'

    """scenario"""
    year = 2013
    soil_type = "loam"
    genotype = "WT"
    comp = "phenolics"
    usemoles = True
    """ parameters   """
    soil_ = scenario.vg_SPP(0)
    
    """ initialize """
    s = RichardsWrapper(RichardsNCSP(), usemoles)  # water and N solute          
    s.dirResults = "./results/"
    
    #@see dumux-rosi\cpp\python_binding\solverbase.hh
    s.noAds = True
    paramIndx=5
    ICcc=None
    p_mean_=-100
    css1Function=9
    scenario.setDefault(s)
    scenario.setSoilParam(s,paramIndx)
    s.theta_init = 0.4# vg.water_content(p_mean_, s.vg_soil)
    
    s.css1Function = css1Function
    s.setParameter( "Soil.css1Function", str(s.css1Function))
    scenario.getBiochemParam(s,paramIndx,s.noAds )
    
    scenario.setBiochemParam(s)
    
    s.initialize(doMPI_=False)
    def getCSS2Init(cs):
        return 0, 0
    s.getCSS2Init = getCSS2Init
    scenario.setIC3D(s, paramIndx, ICcc)
    #s.createGrid(min_b, max_b, cell_number, False)  # [cm] 
    #s.createGrid([-1.5, -6., -40.], [1.5, 6., 0.], [3,12,40])  # [cm]
    s.createGrid([-3./2, -12./2, -50.], [3./2, 12./2, 0.], [3,12,20])  # [cm]
    cell_number_ = cell_number
    cell_number= s.dumux_str(cell_number)#.replace("[", "");cell_number=cell_number.replace("]", "");cell_number=cell_number.replace(",", "");
    s.setParameter( "Soil.Grid.Cells", cell_number)   
    
    s.setParameter("Problem.reactionExclusive", "1")

    s, s.vg_soil = scenario.setupOther(s, p_mean_)
    

    
    s.base.printParams()

    while True:
        print('enter')
        s.solve(1)
        print('leave')
        #cyl3.simulate_const(s,1)
        
          
    print("fin", rank)
    #return results_dir

        

#if __name__ == '__main__':
#    results_dir = XcGrowth()