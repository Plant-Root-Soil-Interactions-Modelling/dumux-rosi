#include "solverbase.hh"

// initialize
#include <dune/common/parallel/mpihelper.hh> // in dune parallelization is realized with MPI
int DUMUX_VERSION = 100; // the future (needed in dumuxmessage.hh)
#include <dumux/common/dumuxmessage.hh> // for fun (a static class)
#include <dumux/common/parameters.hh> // global parameter tree with defaults and parsed from args and .input file

// createGrid
#include <dumux/io/grid/gridmanager.hh>

// simulate
#include <dumux/common/timeloop.hh> // timeloop is need for step size control
// #include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/porousmediumflow/richards/newtonsolver.hh>
#include <dumux/common/timeloop.hh>


///**
// * Writes the Dumux welcome message, and creates the global Dumux parameter tree
// *
// * Normally you state an input file, that contains all parameters that are needed for the simulation.
// * SolverBase will try to set some of them dynamically.
// */
//template<class Problem, class Assembler, class LinearSolver>
//void SolverBase<Problem, Assembler, LinearSolver>::initialize(std::vector<std::string> args)
//{
//    std::vector<char*> cargs;
//    cargs.reserve(args.size());
//    for(size_t i = 0; i < args.size(); ++i) {
//        cargs.push_back(const_cast<char*>(args[i].c_str()));
//    } // its a beautiful language
//
//	auto& mpiHelper = Dune::MPIHelper::instance(cargs.size(), &cargs[0]); // of type MPIHelper, or FakeMPIHelper (in mpihelper.hh)
//
//    if (mpiHelper.rank() == 0) { // rank is the process number
//    	Dumux::DumuxMessage::print(/*firstCall=*/true); // print dumux start message
//    }
//
//    Dumux::Parameters::init(cargs.size(), &cargs[0]); // parse command line arguments and input file
//}

/**
 * Creates a grid from the parameter tree
 */
template<class Problem, class Assembler, class LinearSolver>
void SolverBase<Problem, Assembler, LinearSolver>::createGrid()
{
//    GridManagerFix<Grid> gridManager;
//    gridManager.init();
//    grid = gridManager.gridPtr();
//    try {
//        gridData = gridManager.getGridData(); // test for null by getter
//    } catch(...) { }
}

///**
// * Creates a rectangular grid with a given resolution
// *
// * depending on the Grid you choose at compile time this will work, or not
// * e.g. the method does not make a lot of sense for unstructured grids
// */
//template<class Problem, class Assembler, class LinearSolver>
//void SolverBase<Problem, Assembler, LinearSolver>::createGrid(VectorType boundsMin, VectorType boundsMax, VectorType numberOfCells, std::string periodic)
//{
//    auto p = Dumux::Parameters::paramTree(); // had to modify parameters.hh, its private an no way I can pull it out
//    std::ostringstream bmin;
//    std::ostringstream bmax;
//    std::ostringstream cells;
//    bmin << boundsMin[0] << ", " << boundsMin[1]<< ", " << boundsMin[2];
//    bmax << boundsMax[0] << ", " << boundsMax[1]<< ", " << boundsMax[2];
//    cells << numberOfCells[0] << ", " << numberOfCells[1]<< ", " << numberOfCells[2];
//    p["Soil.Grid.UpperRight"] =  bmin.str();
//    p["Soil.Grid.LowerLeft"] = bmax.str();
//    p["Soil.Grid.Cells"] = cells.str();
//    p["Soil.Grid.Periodid"] = periodic;
//    createGrid();
//}
//
///**
// * Creates a grid from a file
// *
// * depending on the Grid you choose at compile time it will accept the file type, or not.
// */
//template<class Problem, class Assembler, class LinearSolver>
//void SolverBase<Problem, Assembler, LinearSolver>::createGrid(std::string file)
//{
//    auto p = Dumux::Parameters::paramTree();
//    p["Grid.File"] = file;
//    createGrid();
//}
//
///**
// * After the grid is created, the problem can be initialized.
// */
//template<class Problem, class Assembler, class LinearSolver>
//void SolverBase<Problem, Assembler, LinearSolver>::initializeProblem()
//{
//    gridGeometry = std::make_shared<FVGridGeometry>(grid->leafGridView());
//    gridGeometry->update();
//
//    problem = std::make_shared<Problem>(gridGeometry);
//
//    int dof = gridGeometry->numDofs();
//    x = SolutionVector(dof);
//
//    problem->applyInitialSolution(x); // Dumux way of saying x = problem->applyInitialSolution()
//
//    for (int i=0; i<numberOfEquations; i++) { // copy to Python binding
//        initialValues[i].resize(dof);
//        for (int j=0; j<dof; j++) {
//            initialValues[i][j] = x[j];
//        }
//    }
//
//    gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
//    gridVariables->init(x); // initialize all variables , updates volume variables to the current solution, and updates the flux variable cache
//
//}
//
//template<class Problem, class Assembler, class LinearSolver>
//void SolverBase<Problem, Assembler, LinearSolver>::simulate(double dt)
//{
//    using namespace Dumux;
//
//    if (ddt<0) {
//        ddt = dt/10; // initialize
//    }
//
//    double maxDt = getParam<double>("TimeLoop.MaxTimeStepSize", dt); // maximal time step size
//    double initialDt = getParam<double>("TimeLoop.DtInitial", ddt); // initial time step
//    std::shared_ptr<CheckPointTimeLoop<double>> timeLoop =
//        std::make_shared<CheckPointTimeLoop<double>>(/*start time*/0., initialDt, /*final time*/ dt); // the main time loop is moved to Python
//    timeLoop->setMaxTimeStepSize(maxDt);
//
//    auto assembler= std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop); // dynamic
//    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
//    using NonLinearSolver = RichardsNewtonSolver<Assembler, LinearSolver>;
//    auto nonLinearSolver = std::make_shared<NonLinearSolver>(assembler, linearSolver);
//
//    timeLoop->start();
//    auto xOld = x;
//    do {
//        assembler->setPreviousSolution(xOld); // set previous solution for storage evaluations
//
//        nonLinearSolver->solve(x, *timeLoop); // solve the non-linear system with time step control
//
//        xOld = x; // make the new solution the old solution
//
//        gridVariables->advanceTimeStep();
//        timeLoop->advanceTimeStep(); // advance to the time loop to the next step
//        timeLoop->reportTimeStep(); // report statistics of this time step
//
//        ddt = nonLinearSolver->suggestTimeStepSize(timeLoop->timeStepSize());
//        timeLoop->setTimeStepSize(ddt); // set new dt as suggested by the newton solver
//
//        problem->setTime(simTime + timeLoop->time()); // pass current time to the problem ddt?
//
//    } while (!timeLoop->finished());
//
//}
//
//template<class Problem, class Assembler, class LinearSolver>
//std::vector<VectorType> SolverBase<Problem, Assembler, LinearSolver>::getPoints()
//{
//// todo
//	return { VectorType({0,0,0}) };
//}
//
//template<class Problem, class Assembler, class LinearSolver>
//std::vector<VectorType> SolverBase<Problem, Assembler, LinearSolver>::getCells()
//{
//// todo
//	return { VectorType({0,0,0}) };
//}
//
//template<class Problem, class Assembler, class LinearSolver>
//int SolverBase<Problem, Assembler, LinearSolver>::pickIndex(VectorType pos)
//{
//	// todo
//	return 0;
//}
//
//
//template<class Problem, class Assembler, class LinearSolver>
//double SolverBase<Problem, Assembler, LinearSolver>::solutionAt(VectorType pos)
//{
//	return solution[pickIndex(pos)];
//}
