"""
functions for the cylindrical coupling approach
"""
import numpy as np
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
import sys
from xylem_flux import *
import vtk_plot as vtk
import plantbox as pb
#import evapotranspiration as evap
import timeit
import visualisation.vtk_plot as vp
import functional.van_genuchten as vg
from scenario_setup import weather, resistance2conductance
from decimal import *
from functional.xylem_flux import sinusoidal2, sinusoidal

from scenario_setup import write_file_array, write_file_float


        
def simulate_const(s, dt):
    
    k_soil_solve = 0
    redoSolve = True
    maxRelShift = s.MaxRelativeShift
    while redoSolve:
        s.ddt =min( 1.e-5,s.ddt)#or just reset to 1e-5?
        try:
            decreaseMaxRelShift = False
            print("3d solve")
            s.solve(dt, doMPIsolve_=False)  # in modules/solverbase.py
            print("done")

            solComp = [s.getSolution(ncom+1) for ncom in range(s.numComp)]
            whereError = None
            if rank == 0:
                whereError = [np.where(SC <0.) for SC in solComp]
                solComp = [min(SC) for SC in solComp]
            solComp = comm.bcast(solComp, root = 0)
            whereError = comm.bcast(whereError, root = 0)
            if min(solComp) <0.:
                print("min(solComp) <0.", rank, solComp, whereError)
                decreaseMaxRelShift = True
                raise Exception
            redoSolve = False
            # newton parameters are re-read at each 'solve()' calls
            s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift))# reset value
            s.setParameter("Newton.EnableResidualCriterion", "false") # sometimes helps, sometimes makes things worse
            s.setParameter("Newton.EnableAbsoluteResidualCriterion", "false")
            s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "false")
            
            s.setParameter("Newton.MaxSteps", "18")
            s.setParameter("Newton.MaxTimeStepDivisions", "10")
        except Exception as err:
            s.setParameter("Newton.EnableResidualCriterion", "false") # sometimes helps, sometimes makes things worse
            s.setParameter("Newton.EnableAbsoluteResidualCriterion", "false")
            s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "false")
            
            print(rank, f"Unexpected {err=}, {type(err)=}", 'k_soil_solve',k_soil_solve)
            if k_soil_solve > 6:
                raise Exception
            if k_soil_solve == 0:
                s.setParameter("Newton.MaxSteps", "200")
                s.setParameter("Newton.MaxTimeStepDivisions", "100")
            elif k_soil_solve == 1: # worth making the computation more precise?
                print(rank,
                      'soil.solve() failed. making the computation more precise')
                s.setParameter("Newton.EnableResidualCriterion", "true") # sometimes helps, sometimes makes things worse
                s.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
                s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
                s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift/10.))# reset value
            else:
                if decreaseMaxRelShift:
                    change = 0.1
                else:
                    change = 10
                print(rank,
                      'soil.solve() failed. NewtonMaxRelativeShift from',
                      maxRelShift,'to',maxRelShift*change)
                maxRelShift *= change
                # newton parameters are re-read at each 'solve()' calls
                s.setParameter("Newton.MaxRelativeShift", str(maxRelShift))
            s.reset()
            for ncomp in range(s.numComp):
                try:
                    assert (np.array(s.getSolution(ncomp + 1)).flatten() >= 0).all()
                except:
                    raise Exception
            k_soil_solve += 1
    
    