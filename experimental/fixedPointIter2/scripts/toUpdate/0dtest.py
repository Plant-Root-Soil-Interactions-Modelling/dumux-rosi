import os
import sys
import numpy as np

"""
def run_with_timeout(timeout, func, *args, **kwargs):
    # Start the timer
    signal.alarm(int(timeout))
    try:
        result = func(*args, **kwargs)
    finally:
        # Disable the alarm
        signal.alarm(0)
    return result


def example_func(arg):
    import time
    time.sleep(arg)
    return f"Function completed in {arg} seconds"


try:
    result = run_with_timeout(2, example_func, 3)  # This should timeout
except TimeoutException as e:
    print(e)
except Exception as e:
    print(f"An error occurred: {e}")
else:
    print(result)
"""

import threading
import time

class TimeoutException(Exception):
    pass

def run_with_timeout(timeout, func, *args, **kwargs):
    class FuncThread(threading.Thread):
        def __init__(self):
            super().__init__()
            self.result = None
            self.exception = None

        def run(self):
            try:
                self.result = func(*args, **kwargs)
            except Exception as e:
                self.exception = e

    thread = FuncThread()
    thread.start()
    thread.join(timeout)

    if thread.is_alive():
        raise TimeoutException("Function call timed out")
    elif thread.exception:
        raise thread.exception
    else:
        return thread.result

# Example function that takes a specific amount of time to run
def example_func(arg):
    time.sleep(arg)
    return f"Function completed in {arg} seconds"


while redoSolve:
    # s.ddt =min( 1.e-5,s.ddt)#or just reset to 1e-5?
    s.setMaxTimeStepSize(s.maxDt)  # reset each time
    try:
        decreaseMaxRelShift = False
        if rank == 0:
            print("solve 3d soil")
        # s.solve(dt)  # in modules/solverbase.py
        helpfull.run_with_timeout(60., s.solve, dt)  # time out error after Xmn
        if rank == 0:
            print("solve 3d soil finished")

        # if we got solute content < 0, throw exception
        solComp = [s.getSolution(ncom + 1) for ncom in range(s.numSoluteComp)]
        whereError = None
        if rank == 0:
            whereError = [np.where(SC < 0.) for SC in solComp]
            solComp = [min(SC) for SC in solComp]
        solComp = comm.bcast(solComp, root=0)
        whereError = comm.bcast(whereError, root=0)
        if min(solComp) < 0.:
            print("min(solComp) <0.", rank, solComp, whereError)
            decreaseMaxRelShift = True
            raise Exception

        redoSolve = False
        ## solving succeded, reset solving parameters (in case  they were changed)
        # newton parameters are re-read at each 'solve()' calls
        s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift))  # reset value
        s.setParameter("Newton.EnableResidualCriterion",
                       "false")  # sometimes helps, sometimes makes things worse
        s.setParameter("Newton.EnableAbsoluteResidualCriterion", "false")
        s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "false")

        s.setParameter("Newton.MaxSteps", "18")
        s.setParameter("Newton.MaxTimeStepDivisions", "10")
        s.createNewtonSolver()  # re-create Newton solver to implement the new newton parameters

    except Exception as err:
        s.setParameter("Newton.EnableResidualCriterion",
                       "false")  # sometimes helps, sometimes makes things worse
        s.setParameter("Newton.EnableAbsoluteResidualCriterion", "false")
        s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "false")

        print(rank, f"Unexpected {err=}, {type(err)=}", 'k_soil_solve', k_soil_solve)

        if k_soil_solve > 6:  # to avoid going to infinit loop
            raise Exception
        if k_soil_solve == 0:
            # for the 1st fail, simply increase number of steps allowed
            s.setParameter("Newton.MaxSteps", "200")
            s.setParameter("Newton.MaxTimeStepDivisions", "100")
        elif k_soil_solve == 1:
            # 2nd fail: try making the computation more precise
            print(rank,
                  'soil.solve() failed. making the computation more precise')
            s.setParameter("Newton.EnableResidualCriterion",
                           "true")  # sometimes helps, sometimes makes things worse
            s.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
            s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
            s.setParameter("Newton.MaxRelativeShift",
                           str(s.MaxRelativeShift / 10.))  # reset value
        else:
            # try decreasing MaxRelShift (if got solute content < 0)
            # try increasing maxrelshift (if got a non-convergence error)
            if decreaseMaxRelShift:
                change = 0.1
            else:
                change = 10
            print(rank,
                  'soil.solve() failed. NewtonMaxRelativeShift from',
                  maxRelShift, 'to', maxRelShift * change)
            maxRelShift *= change
            # newton parameters are re-read at each 'solve()' calls
            s.setParameter("Newton.MaxRelativeShift", str(maxRelShift))
        s.reset()  # reset solution vector
        s.createNewtonSolver()  # re-create Newton solver to implement the new newton parameters
        k_soil_solve += 1
