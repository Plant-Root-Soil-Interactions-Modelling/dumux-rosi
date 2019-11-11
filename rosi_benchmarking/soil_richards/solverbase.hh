#ifndef SOLVERBASE_H_
#define SOLVERBASE_H_

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/pybindxi/numpy.h>
namespace py = pybind11;

#include <sstream>
#include <stdexcept>

#include <dumux/io/grid/gridmanager.hh> // grid types
#include <dune/grid/common/mcmgmapper.hh> // the element and vertex mappers
#include <dumux/common/defaultmappertraits.hh> // nothingness

// initialize
#include <dune/common/parallel/mpihelper.hh> // in dune parallelization is realized with MPI
#include <dumux/common/dumuxmessage.hh> // for fun (a static class)
#include <dumux/common/parameters.hh> // global parameter tree with defaults and parsed from args and .input file

// createGrid
#include <dumux/io/grid/gridmanager.hh>

// initialize
#include <dune/common/parallel/mpihelper.hh> // in dune parallelization is realized with MPI
#include <dumux/common/dumuxmessage.hh> // for fun (a static class)
#include <dumux/common/parameters.hh> // global parameter tree with defaults and parsed from args and .input file

// simulate
#include <dumux/common/timeloop.hh> // timeloop is need for step size control
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/porousmediumflow/richards/newtonsolver.hh>

using VectorType = std::array<double,3>;

/**
 * Derived class will pass ownership
 */
template<class Grid>
class GridManagerFix :public Dumux::GridManager<Grid> {
public:
	std::shared_ptr<Grid> gridPtr() {
		return this->gridPtr_;
	}
};

/**
 * Dumux as a solver with a simple Python interface
 *
 * first we want to run a simulation with .input file parameters.
 * in future these parameters will be replaced step by step.
 */
template<class Problem, class Assembler, class LinearSolver>
class SolverBase {
public:

	std::string gridType = "YaspGrid"; // <- for better description and warnings
	int dim = Problem::dimWorld;
	bool isBox = Problem::isBox; // numerical method
	const std::vector<std::string> primNames = { "Matric potential [Pa]" };

	SolverBase()
	{ };

	/**
	 * Writes the Dumux welcome message, and creates the global Dumux parameter tree
	 *
	 * Normally you state an input file, that contains all parameters that are needed for the simulation.
	 * SolverBase will try to set some of them dynamically.
	 */
	void initialize(std::vector<std::string> args)
	{
		std::vector<char*> cargs;
		cargs.reserve(args.size());
		for(size_t i = 0; i < args.size(); ++i) {
			cargs.push_back(const_cast<char*>(args[i].c_str()));
		} // its a beautiful language

		int argc = cargs.size();
		char** argv  =  &cargs[0];

		auto& mpiHelper = Dune::MPIHelper::instance(argc, argv); // of type MPIHelper, or FakeMPIHelper (in mpihelper.hh)

		if (mpiHelper.rank() == 0) { // rank is the process number
			Dumux::DumuxMessage::print(/*firstCall=*/true); // print dumux start message
		}

		Dumux::Parameters::init(argc, argv); // parse command line arguments and input file
	}

	/**
	 * Creates a grid from the (global) Dumux parameter tree
	 */
	void createGrid()
	{
		GridManagerFix<Grid> gridManager;
		gridManager.init();
		grid = gridManager.gridPtr();
		try {
			gridData = gridManager.getGridData(); // test for null by getter
		} catch(...) { }
	}

	/**
	 * Creates a rectangular grid with a given resolution
	 *
	 * depending on the Grid you choose at compile time this will work, or not
	 * e.g. the method does not make a lot of sense for unstructured grids
	 */
	virtual void createGrid(VectorType boundsMin, VectorType boundsMax,
			VectorType numberOfCells, std::string periodic = "false false false")
	{
		auto& p = Dumux::Parameters::paramTree(); // had to modify parameters.hh, its private an no way I can pull it out
		std::ostringstream bmin;
		std::ostringstream bmax;
		std::ostringstream cells;
		bmin << boundsMin[0] << " " << boundsMin[1]<< " " << boundsMin[2];
		bmax << boundsMax[0] << " " << boundsMax[1]<< " " << boundsMax[2];
		cells << numberOfCells[0] << " " << numberOfCells[1]<< " " << numberOfCells[2];
		p["Grid.UpperRight"] =  bmin.str();
		p["Grid.LowerLeft"] = bmax.str();
		p["Grid.Cells"] = cells.str();
		p["Grid.Periodid"] = periodic;
		createGrid();
	}


	/**
	 * Creates a grid from a file
	 *
	 * depending on the Grid you choose at compile time it will accept the file type, or not.
	 */
	virtual void createGrid(std::string file)
	{
		auto& p = Dumux::Parameters::paramTree();
		p["Grid.File"] = file;
		createGrid();
	}

	/**
	 * Writes a parameter into the global Dumux parameter map
	 */
	virtual void setParameter(std::string key, std::string value)
	{
		auto& p = Dumux::Parameters::paramTree();
		p[key] = value;
	}

	////    virtual void initialConditions(); // TODO
	////    virtual void boundaryConditions(); // TODO
	////    virtual void sourceTerm(); // TODO

	/**
	 * After the grid is created, the problem can be initialized
	 *
	 * creates (a) the GridGeometry, (b) the Problem, (c) applies initial conditions (todo),
	 * (d) GridVariables (e) resets simtime, ddt
	 */
	virtual void initializeProblem()
	{
		gridGeometry = std::make_shared<FVGridGeometry>(grid->leafGridView());
		gridGeometry->update();
		problem = std::make_shared<Problem>(gridGeometry);
		int dof = gridGeometry->numDofs();
		x = SolutionVector(dof);
		problem->applyInitialSolution(x); // Dumux way of saying x = problem->applyInitialSolution()

		for (int i=0; i<numberOfEquations; i++) { // copy to Python binding
			initialValues[i].resize(dof);
			solution[i].resize(dof);
			for (int j=0; j<dof; j++) {
				initialValues[i][j] = x[j];
				solution[i][j] = x[j];
			}
		}

		gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
		gridVariables->init(x); // initialize all variables , updates volume variables to the current solution, and updates the flux variable cache
		simTime = 0; // reset
		ddt = -1;
	}

	/**
	 * Returns the Dune vertices (vtk points) of the grid
	 */
	virtual std::vector<VectorType> getPoints() // vtk naming
	{
		std::vector<VectorType> points;
		points.reserve(gridGeometry->numDofs()); // depending on method this is correct, or a guess
		for (const auto& v : vertices(grid->leafGridView())) {
			auto p = v.geometry().center();
			switch(dim) {
			case 1: // 1D
				points.push_back({0., 0., p[0]});
				break;
			case 2: // 2D
				points.push_back({p[0], p[1], 0.}); // ignore z
				break;
			case 3: // 3D
				points.push_back({p[0], p[1], p[2]});
				break;
			default:
				throw std::invalid_argument("SolverBase::getPoints: Dimension "+ std::to_string(dim) + "D is not supported");
			}

		}
		return points;
	}

	/**
	 * return the Dune elements (vtk cells) of the grid
	 */
	virtual std::vector<VectorType> getCells()
	{
		std::vector<VectorType> cells;
		cells.reserve(gridGeometry->numDofs()); // depending on method this is correct, or a guess
		for (const auto& e : elements(grid->leafGridView())) {
			auto c = e.geometry().center();
			cells.push_back({c[0], c[1], c[2]});
		}
		return cells;
	}

	/**
	 * Returns the points, where the DOF sit, in the same order like the solution values.
	 *
	 * DOF sit either at the vertices (points) for box method or
	 * element centers (cell centers) for CCTpfa.
	 */
	virtual std::vector<VectorType> getDof() {
		if (isBox) {
			return getPoints();
		} else {
			return getCells();
		}
	}

	/**
	 * simulates the problem for time spand dt, with initial time step ddt.
	 *
	 * Assembler needs a TimeLoop, so i have to create it in each simulate call.
	 * (could be improved, but overhead is likely to be small)
	 */
	virtual void simulate(double dt, double maxDt = -1)
	{
		using namespace Dumux;

		if (ddt<1.e-6) { // happens at the first call
			ddt = getParam<double>("TimeLoop.DtInitial", dt/10); // from params, or guess something
		}

		std::shared_ptr<CheckPointTimeLoop<double>> timeLoop =
				std::make_shared<CheckPointTimeLoop<double>>(/*start time*/0., ddt, /*final time*/ dt); // the main time loop is moved to Python
		if (maxDt<0) { // per default value take from parameter tree
			maxDt = getParam<double>("TimeLoop.MaxTimeStepSize", dt); // maximal time step size
		}
		timeLoop->setMaxTimeStepSize(maxDt);

		auto assembler= std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop); // dynamic
		auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
		using NonLinearSolver = RichardsNewtonSolver<Assembler, LinearSolver>;
		auto nonLinearSolver = std::make_shared<NonLinearSolver>(assembler, linearSolver);

		timeLoop->start();
		auto xOld = x;
		do {
			assembler->setPreviousSolution(xOld); // set previous solution for storage evaluations

			nonLinearSolver->solve(x, *timeLoop); // solve the non-linear system with time step control

			xOld = x; // make the new solution the old solution

			gridVariables->advanceTimeStep();
			timeLoop->advanceTimeStep(); // advance to the time loop to the next step
			timeLoop->reportTimeStep(); // report statistics of this time step

			ddt = nonLinearSolver->suggestTimeStepSize(timeLoop->timeStepSize());
			timeLoop->setTimeStepSize(ddt); // set new dt as suggested by the newton solver

			problem->setTime(simTime + timeLoop->time()); // pass current time to the problem ddt?

		} while (!timeLoop->finished());

		simTime += dt;

	}

	//    int pick(VectorType pos); // vertex or element index
	//    double interpolate(VectorType pos, int i = 0);

	/**
	 * Quick overview
	 */
	std::string toString()
	{
		std::ostringstream msg;
		msg << "DuMux Solver using " << gridType << " in " << dim << "D ";
		if (isBox) {
			msg << "(box method)";
		} else {
			msg << "(CCTpfa method)";
		}
		if (gridGeometry) {
			msg << " with " << gridGeometry->numDofs() << " DOF\n";
		} else {
			msg << " uninitialized ";
			return msg.str();
		}
		msg << "Simulation time is "<< simTime/3600/24 << " days, current internal time step is "<< ddt/3600/24 << " days";
		return msg.str();
	}

	double simTime = 0;
	double ddt = -1; // internal time step, minus indicates that it needs to be initialized

	static const bool numberOfEquations = 1;
	std::array<std::vector<double>,numberOfEquations> initialValues;
	std::array<std::vector<double>, numberOfEquations> solution;

protected:

	using Grid = typename Problem::Grid;
	using GridData = Dumux::GridData<Grid>;
	using FVGridGeometry = typename Problem::FVGridGeometry;
	using SolutionVector = typename Problem::SolutionVector;
	using GridVariables = typename Problem::GridVariables;

	std::shared_ptr<Grid> grid;
	std::shared_ptr<GridData> gridData;
	std::shared_ptr<FVGridGeometry> gridGeometry;
	std::shared_ptr<Problem> problem;
	std::shared_ptr<GridVariables> gridVariables;
	SolutionVector x;

};

#endif

///**
// * lacking polymorphism I have not found a way to make the gird dynamic.
// * you need to choose the grid at compile time,
// * and I don't know how to pass it via CMakeLists.txt building the Python binding
// */
//using Grid = Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>;
//using GridView = typename Dune::YaspGridFamily<3,Dune::EquidistantOffsetCoordinates<double,3>>::Traits::LeafGridView;
//using Scalar = double;
//using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
//using VertexMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
//using MapperTraits = Dumux::DefaultMapperTraits<GridView, ElementMapper, VertexMapper>;
//using FVGridGeometry = Dumux::BoxFVGridGeometry<double, GridView, /*enableCache*/ true, Dumux::BoxDefaultGridGeometryTraits<GridView, MapperTraits>>;
