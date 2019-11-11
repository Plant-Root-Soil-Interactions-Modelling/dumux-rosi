import os
import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/soil_richards/")
import richards_yasp_solver as solver

print(os.getcwd())

s = solver.RichardsYaspSolver()

s.initialize([""])
s.initialize(["-input","../input/b1a_3d.input"])
s.createGrid() # wrong param group (todo)
# s.createGrid([-1, -1, -1], [1,1,1], [2, 2, 2], "false, false, false")
#s.createGrid("../grids/b1.dgf") # does only take intervals

s.initializeProblem()

print()
print("Number of points", len(s.getPoints()))
print("Number of cells", len(s.getCells()))

print(); 
print(s)

day = 3600*24
s.ddt = day
for i in range(0, 10):
    print("***************time step****************", i, " internal initial time step", s.ddt)
    s.simulate(day, -1)
    #
    # do wild stuff 
    #
 
print(len(s.initialValues[0]))
print(s.initialValues[0])
print(len(s.solution[0]))
print(s.solution)
 
print("done")
