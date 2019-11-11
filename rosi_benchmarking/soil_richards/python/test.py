import sys
sys.path.append("../../../build-cmake/rosi_benchmarking/soil_richards/")
import richards_yasp_solver as solver

print("alive")

s = solver.RichardsYaspSolver()

print(s.initialValues)
print(s.solution)

# s.createGrid()
# solver.initialize(["input/b1a_3d.input"])
