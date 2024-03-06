import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

""" CPLantBox tutorial example 6c (see CPlantBox/tutorial/latex/PlantBox_RootSystem/tutorial.tex) """
""" coupling with DuMux as solver for the soil, run in dumux-rosi """
from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
from functional.root_conductivities import *  # hard coded conductivities
import plantbox as pb
#import visualisation.vtk_plot as vp
from richards import RichardsWrapper  # Python part

import rsml.rsml_reader as rsml
import numpy as np
#import matplotlib.pyplot as plt
import timeit
import pandas as pd

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()



def testSolver(solverType, name, SLURM_NTASKS,NODELIST,MEM ):      
    def sinusoidal(t):
        return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.

    """ Parameters """
    min_b = [-4., -4., -25.]
    max_b = [4., 4., 0.]
    cell_number = [8, 8, 25]  # [16, 16, 30]  # [32, 32, 60]
    periodic = True

    path = "../../../CPlantBox/modelparameter/structural/rootsystem/"
    fname = "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
    loam = [0.08, 0.43, 0.04, 1.6, 50]
    initial = -659.8 + 12.5  # -659.8

    trans = 6.4  # cm3 /day (sinusoidal)
    wilting_point = -15000  # cm

    sim_time = 7  # [day] for task b
    rs_age = 10  # root system initial age
    age_dependent = False  # conductivities
    dt = 120. / (24 * 3600)  # [days] Time step must be very small
    
    extraprint = False
    if solverType == 0:
        from rosi_richards import RichardsSP 
    elif solverType == 1:
        from rosi_richards import RichardsSPSeq as RichardsSP
    elif solverType == 2:
        from rosi_richards import RichardsSPILU as RichardsSP
    elif solverType == 3:
        from rosi_richards import RichardsSPILURes as RichardsSP
    elif solverType == 4:
        from rosi_richards import RichardsSPSSORC as RichardsSP
    else: 
        raise Exception   
    # from rosi_richards import RichardsSPnum as RichardsSP

    """ Initialize macroscopic soil model """
    s = RichardsWrapper(RichardsSP())
    s.initialize()
    s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
    s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
    s.setTopBC("noFlux")
    s.setBotBC("noFlux")
    s.setVGParameters([loam])
    if False:
        s.setParameter("Newton.MaxSteps", "200")
        s.setParameter("Newton.MaxTimeStepDivisions", "100")
        s.setParameter("Newton.EnableResidualCriterion", "true") 
        s.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
        s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
    #s.setParameter("Newton.MaxRelativeShift", "1e-9")# reset value
    if extraprint:
        print("created s", flush=True)
    s.initializeProblem(maxDt = 2500/(24*3600))
    if extraprint:
        print("created s", flush=True)
    s.setCriticalPressure(wilting_point)

    """ Initialize xylem model """
    if False:
        fname = "../../grids/RootSystem8.rsml"
        r = XylemFluxPython(fname)
        r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), 
                                pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                                pb.Vector3d(cell_number[0], cell_number[1], 
                                            cell_number[2]), True)  # cutting
        init_conductivities(r, False)
    else:

        rs = pb.MappedRootSystem()
        rs.setSeed(1)
        rs.readParameters(path + fname + ".xml")
        if not periodic:
            sdf = pb.SDF_PlantBox(0.99 * (max_b[0] - min_b[0]), 0.99 * (max_b[1] - min_b[1]), max_b[2] - min_b[2])
        else:
            sdf = pb.SDF_PlantBox(np.Inf, np.Inf, max_b[2] - min_b[2])
        rs.setGeometry(sdf)
        rs.initialize()
        rs.simulate(rs_age, False)
        r = XylemFluxPython(rs)
        init_conductivities(r, age_dependent)

    """ Coupling (map indices) """
    picker = lambda x, y, z: s.pick([x, y, z])
    if extraprint:
        print("created s", flush=True)
    r.rs.setSoilGrid(picker)  # maps segments
    if extraprint:
        print("created s", flush=True)
    r.rs.setRectangularGrid(pb.Vector3d(min_b), pb.Vector3d(max_b), pb.Vector3d(cell_number), True)
    if extraprint:
        print("created s", flush=True)
    r.test()  # sanity checks
    nodes = r.get_nodes()
    if extraprint:
        print("created s", flush=True)
    cci = picker(nodes[0, 0], nodes[0, 1], nodes[0, 2])  # collar cell index
    if extraprint:
        print("created s", flush=True)

    """ Numerical solution """
    start_time = timeit.default_timer()
    x_, y_, x0, rx0 = [], [], [], []
    if extraprint:
        print("created s", flush=True)
    sx = comm.bcast(s.getSolutionHead() ,root=0) # inital condition, solverbase.py
    if extraprint:
        print("created s", flush=True)
    N = round(sim_time / dt)
    t = 0.

    for i in range(0, N):

        if extraprint:
            print("created s", flush=True)

        rx = r.solve(rs_age + t, -trans * sinusoidal(t), sx[cci], sx, True, wilting_point)  # xylem_flux.py
        if extraprint:
            print("created s", flush=True)
        x_.append(t)
        y_.append(float(r.collar_flux(rs_age + t, rx, sx)))

        fluxes = r.soilFluxes(rs_age + t, rx, sx, False)
        s.setSource(fluxes)  # richards.py
        s.solve(dt)
        sx =  comm.bcast(s.getSolutionHead()  ,root=0) # richards.py
        x0.append(sx[cci])
        rx0.append(rx[0])

        min_sx, min_rx, max_sx, max_rx = np.min(sx), np.min(rx), np.max(sx), np.max(rx)
        n = round(float(i) / float(N) * 100.)
        if extraprint:
            print(solverType, "[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g} days {:g}"
                .format(min_sx, max_sx, min_rx, max_rx, s.simTime, rx[0]))
        t += dt
    solvetime =  timeit.default_timer() - start_time
    print ("Coupled benchmark solved in ",solvetime, " s")
    with open("./results/"+name+"example6c_time.txt", 'a') as log:
        log.write(", ".join([num for num in map(str, np.array([solverType,solvetime, name, SLURM_NTASKS,NODELIST ,MEM]))])  +'\n')

    with open("./results/"+name+"example6c_y_"+str(solverType)+".txt", 'w') as log:
        log.write(", ".join([num for num in map(str, y_)])  +'\n')
    with open("./results/"+name+"example6c_cci_"+str(solverType)+".txt", 'w') as log:
        log.write(", ".join([num for num in map(str, x0)])  +'\n')
    with open("./results/"+name+"example6c_rx0_"+str(solverType)+".txt", 'w') as log:
        log.write(", ".join([num for num in map(str, rx0)])  +'\n')


    return x0, y_, rx0


if __name__ == '__main__':
    name = "" 
    if len(sys.argv)>1:
        name = sys.argv[1]
    SLURM_NTASKS=1
    if len(sys.argv)>2:
        SLURM_NTASKS = sys.argv[2]
    NODELIST="node01"
    if len(sys.argv)>3:
        NODELIST = sys.argv[3]
    MEM="10"
    if len(sys.argv)>4:
        MEM = sys.argv[4]
    print('extra name',name)   
    x0s, ys, rx0 = [], [], []
    for solverType in range(5):
        x0, y_, rx0_ = testSolver(solverType, name, SLURM_NTASKS,NODELIST,MEM )
        x0s.append(x0)
        ys.append(y_)
        rx0.append(rx0_)


    newnames = ['AMG_analytic','AMG_seq','ILU', 'ILURes', 'SORC']# 'AMG_num',
    if True:
        df = pd.DataFrame(x0s, index = newnames ,dtype=float)
        df = df.T
        df.to_csv('./results/'+name+'6dallSolvers_x.csv', index=False)
        df = pd.DataFrame(ys, index = newnames ,dtype=float)
        df = df.T
        df.to_csv('./results/'+name+'6dallSolvers_y.csv', index=False)
        df = pd.DataFrame(rx0, index = newnames ,dtype=float)
        df = df.T
        df.to_csv('./results/'+name+'6dallSolvers_rx.csv', index=False)

    if False:
        import numpy as np
        import matplotlib.pyplot as plt

        sim_time = 7 
        dt = 120. / (24 * 3600)  # [days] Time step must be very small
        N = round(sim_time / dt)
        x = np.array([i *dt for i in range(0, N)])


        dfrx = pd.read_csv('./results/6dallSolvers_rx.csv')
        dfsx = pd.read_csv('./results/6dallSolvers_x.csv')
        dfy  = pd.read_csv('./results/6dallSolvers_y.csv')

        for i, column in enumerate(dfrx.columns):
            plt.plot(x, dfrx[column], label=newnames[i])

        plt.xlabel('day')
        plt.ylabel('psi_x at root colar (cm)')
        #plt.title('Line Plot of DataFrame Columns')
        plt.legend()
        plt.grid(True)
        plt.show()

        for i, column in enumerate(dfsx.columns):
            plt.plot(x, dfsx[column], label=newnames[i])

        plt.xlabel('day')
        plt.ylabel('psi_soil around root colar (cm)')
        #plt.title('Line Plot of DataFrame Columns')
        plt.legend()
        plt.grid(True)
        plt.show()

        for i, column in enumerate(dfy.columns):
            plt.plot(x, -dfy[column], label=newnames[i])

        plt.xlabel('day')
        plt.ylabel("Transpiration $[cm^3 d^{-1}]$")
        #plt.title('Line Plot of DataFrame Columns')
        plt.legend()
        plt.grid(True)
        plt.show()

        for i, column in enumerate(dfy.columns):
            plt.plot(x, np.cumsum(-dfy[column]*dt), label=newnames[i])

        plt.xlabel('day')
        plt.ylabel("Cumulative Transpiration $[cm^3 d^{-1}]$")
        #plt.title('Line Plot of DataFrame Columns')
        plt.legend()
        plt.grid(True)
        plt.show()