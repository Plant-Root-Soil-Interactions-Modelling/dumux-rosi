#ifndef PYTHON_SOLVER_BASE_H_
#define PYTHON_SOLVER_BASE_H_

#include <type_traits>
#include <dumux/linear/istlsolvers.hh>


// initialize
#include <stdlib.h>
#include <dumux/common/initialize.hh>
#include <dune/common/parallel/mpihelper.hh> // in dune parallelization is realized with MPI
#include <dumux/common/dumuxmessage.hh> // for fun (a static class)
#include <dumux/common/parameters.hh> // global parameter tree with defaults and parsed from args and .input file

// createGrid
#include <dumux/io/grid/gridmanager.hh>
#include <dune/grid/common/gridfactory.hh>

// simulate
#include <dumux/common/timeloop.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include "../../dumux/porousmediumflow/richards/newtonsolver.hh"

// getDofIndices, getPointIndices, getCellIndices
#include <dune/grid/utility/globalindexset.hh>

// pick
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/discretization/localview.hh>

#include <ostream>
#include <iostream>
#include <limits>
#include <array>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
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
 * Dumux as a solver with a simple Python interface.
 *
 * No input file is needed. The .input file parameters are passed to
 * the solver before the problem is created.
 *
 * There is no need to write a VTK, geometry and results can be easily
 * accessed.
 *
 * This class works in combination with methods from the Python side
 * in solverbase.py. For simplicity most MPI communication is done in
 * Python.
 *
 * Examples are given in the python directory
 */
template<class Problem, class Assembler, class LinearSolver, int dim = 3 /*Problem::dimWorld */>
class SolverBase {
public:

    using VectorType = std::array<double, dim>;
    using NumEqVector = typename Problem::NumEqVector;	
    using GlobalPosition = typename Problem::GlobalPosition; 
	
    int numComp(){return nEV.size();}
    int numFluidComp(){return Problem::FluidSystem::numComponents;}
    
    NumEqVector nEV;
    bool isBox = Problem::isBox; // numerical method
    int dimWorld = dim;
	
    double simTime = 0;
    double ddt = -1; // internal time step, minus indicates that its uninitialized
    int maxRank = -1; // max mpi rank
    int rank = -1; // mpi rank
	bool problemInitialized = false;
    bool periodic = false; // periodic domain
    std::array<int, dim> numberOfCells;

    //
    std::shared_ptr<Dumux::CheckPointTimeLoop<double>> timeLoop;
    using NonLinearSolver =      Dumux::RichardsNewtonSolver<Assembler, LinearSolver>;
    using NonLinearSolverNoMPI = Dumux::RichardsNewtonSolver<Assembler, LinearSolver,
								 Dumux::PartialReassembler<Assembler>,
								Dune::Communication<Dune::FakeMPIHelper::MPICommunicator> >;
    std::shared_ptr<Assembler> assembler;
    std::shared_ptr<LinearSolver> linearSolver;
    std::shared_ptr<NonLinearSolver> nonLinearSolver;
    std::shared_ptr<NonLinearSolverNoMPI> nonLinearSolverNoMPI;


    SolverBase() {
        for (int i=0; i<dim; i++) { // initialize numberOfCells
            numberOfCells[i] = 0;
        }
        setParameter("Grid.Overlap","1");
    }

    virtual ~SolverBase() { }

    //! prints all used and unused parameters
	void printParams()
	{
		Dumux::Parameters::print();
	}
    

    /**
     * Writes the Dumux welcome message, and creates the global Dumux parameter tree from defaults and the .input file
     *
     * Normally you state an input file, that contains all parameters that are needed for the simulation.
     * SolverBase will optionally set most of them dynamically.
     */
    virtual void initialize(std::vector<std::string> args_ = std::vector<std::string>(0), bool verbose = true, bool doMPI = true) {
    
		setenv("DUMUX_NUM_THREADS","1",0);// will not overwrite the DUMUX_NUM_THREADS value if it already exists.
		
		if (const char* dumuxNumThreads = std::getenv("DUMUX_NUM_THREADS"))
		{
			double num_threads= std::max(1, std::stoi(std::string{ dumuxNumThreads }));
			if(num_threads==1)
			{
				setParameter("Assembly.Multithreading","false");
			}
		}
        
		std::vector<char*> cargs;
        cargs.reserve(args_.size());
        for(size_t i = 0; i < args_.size(); i++) {
            cargs.push_back(const_cast<char*>(args_[i].c_str()));
        } // its a beautiful language
        int argc = cargs.size();
        char** argv  =  &cargs[0];

		Dumux::initialize(argc, argv);
        auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);
        maxRank = mpiHelper.size();
        rank = mpiHelper.rank();

        if (isBox) { // add DuMux peculiarities
            setParameter("Grid.Overlap","0");
        } else {
            setParameter("Grid.Overlap","1");
        }
        setParameter("Flux.UpwindWeight", "0.5"); // Timo's advice for flows that are not strongly convection dominated, Dumux default = 1

        if ((rank == 0) && verbose) { // rank is the process number
            std::cout << "\n" << toString() << "\n" << std::flush; // add my story
            Dumux::DumuxMessage::print(/*firstCall=*/true); // print dumux start message
        }
		if(doMPI)
		{
			mpiHelper.getCommunication().barrier(); // no one is allowed to mess up the message
		}
        setParameter("Problem.Name","noname");
        Dumux::Parameters::init(argc, argv); // parse command line arguments and input file
    }

    /**
     * Creates the Grid and gridGeometry from the (global) DuMux parameter tree.
     *
     * Parameters are (there might be more):
     * Grid.UpperRight
     * Grid.LowerLeft
     * Grid.Cells
     * Grid.Periodic
     * Grid.File
     * Grid.Overlap (should = 0 for box, = 1 for CCTpfa), automatically set in SolverBase::initialize
     */
    virtual void createGrid(std::string modelParamGroup = "") {
        std::string pstr =  Dumux::getParam<std::string>("Grid.Periodic", " ");
        periodic = ((pstr.at(0)=='t') || (pstr.at(0)=='T')); // always x,y, not z
        GridManagerFix<Grid> gridManager;
        gridManager.init(modelParamGroup);
        grid = gridManager.gridPtr();
        try {
            gridData = gridManager.getGridData(); // test for null by getter (todo)
            // missing concept to add grid data dynamically
        } catch(...) { }
        gridGeometry = std::make_shared<FVGridGeometry>(grid->leafGridView());
        gridGeometry->update(gridManager.grid().leafGridView());
    }


    /**
     * Creates a rectangular grid with a given resolution
     *
     * depending on the Grid you choose at compile time this will work, or not
     * e.g. the method does not make a lot of sense for unstructured grids
     */
    virtual void createGrid(VectorType boundsMin, VectorType boundsMax,
        std::array<int, dim> numberOfCells, bool periodic = false) {
        this->numberOfCells = numberOfCells;
        auto& p = Dumux::Parameters::paramTree_(); // had to modify parameters.hh, its private an no way I can pull it out
        std::ostringstream bmin;
        std::ostringstream bmax;
        std::ostringstream cells;
        for (int i=0; i<dim; i++) {
            if (boundsMin[i] >= boundsMax[i]) {
                throw std::invalid_argument("SolverBase::createGrid: bounds min >= bounds max");
            }
            bmin << boundsMin[i] << " ";
            bmax << boundsMax[i] << " ";
            cells << numberOfCells[i] << " ";
        }
        p["Grid.LowerLeft"] = bmin.str();
        p["Grid.UpperRight"] = bmax.str();
        p["Grid.Cells"] = cells.str();
        if (dim == 1) { // todo dim = 2
            if (periodic) {
                p["Grid.Periodic"] = "true";
            } else {
                p["Grid.Periodic"] = "false";
            }
        }
        if (dim == 3) {
            if (periodic) {
                std::string xx = "true ";
                std::string yy = "true ";
                if (numberOfCells[0]==1) {
                    xx = "false ";
                }
                if (numberOfCells[1]==1) {
                    yy = "false ";
                }
                p["Grid.Periodic"] = xx+yy+"false";
            } else {
                p["Grid.Periodic"] = "false false false";
            }
        }
        createGrid();
    }

    /**
     * Creates a 1D grid from points for 1d or 1d/3d
     */
    virtual void createGrid1d(const std::vector<VectorType>& points) {
        Dune::GridFactory<Grid> factory;
        constexpr auto line = Dune::GeometryTypes::line;
        Dune::FieldVector<double, dim> p;
        for (int k = 0; k<points.size(); k++) {
            for (int i=0; i<dim; i++) {
                p[i] = points[k][i];
            }
            factory.insertVertex(p);
        }
        for (int k = 0; k<points.size()-1; k++) {
            factory.insertElement(line,  {static_cast<unsigned int>(k), static_cast<unsigned int>(k+1)});
        }
        grid = std::shared_ptr<Grid>(factory.createGrid());
        gridGeometry = std::make_shared<FVGridGeometry>(grid->leafGridView());
        gridGeometry->update(grid->leafGridView());
    }

    //    /**
    //     * Creates a 1D grid with number of cells = [1,1,points.size()-1] where points are the cell centers for 3d/3d
    //     */
    //    virtual void createGrid3d(const std::vector<VectorType>& points, const std::vector<VectorType>& p0) {
    //        assert(dim==3 && "SolverBase::createGrid3d: works only for dim==3");
    //        Dune::GridFactory<Grid> factory;
    //        constexpr auto cube = Dune::GeometryTypes::hexahedron;
    //        std::array<Dune::FieldVector<double, dim>,4> p0_;
    //        for (int j=0; j<4; j++) {
    //            for (int i=0; i<dim; i++) {
    //                p0_[j][i] = p0.at(j).at(i);
    //            }
    //        }
    //        Dune::FieldVector<double, dim> p;
    //        for (int i=0; i<dim; i++) {
    //            p[i] = points.at(0).at(i);
    //        }
    //        for (int i=0; i<4; i++) {
    //            factory.insertVertex(p+p0_[i]);
    //        }
    //        for (int k = 1; k<points.size(); k++) {
    //            for (int i=0; i<3; i++) {
    //                p[i] = points.at(k).at(i);
    //            }
    //            for (int i=0; i<4; i++) {
    //                factory.insertVertex(p+p0_[i]);
    //            }
    //        }
    //        for (int k = 0; k<points.size()-1; k++) {
    //            factory.insertElement(cube,  {
    //                static_cast<unsigned int>(4*k), static_cast<unsigned int>(4*k+1),
    //                static_cast<unsigned int>(4*k+2), static_cast<unsigned int>(4*k+3),
    //                static_cast<unsigned int>(4*k+4), static_cast<unsigned int>(4*k+5),
    //                static_cast<unsigned int>(4*k+6), static_cast<unsigned int>(4*k+7)});
    //        }
    //        grid = std::shared_ptr<Grid>(factory.createGrid());
    //        gridGeometry = std::make_shared<FVGridGeometry>(grid->leafGridView());
    //        gridGeometry->update();
    //    }

    /**
     * Creates a grid from a file
     *
     * depending on the Grid you choose at compile time it will accept the file type, or not.
     */
    virtual void readGrid(std::string file) {
        auto& p = Dumux::Parameters::paramTree_();
        p["Grid.File"] = file;
        createGrid();
    }

    /**
     * Returns a rectangular bounding box around the grid geometry
     *
     * [minx, miny, minz, maxx, maxy, maxz]
     */
    virtual std::vector<double> getGridBounds() {
        auto bMax = gridGeometry->bBoxMax();
        auto bMin = gridGeometry->bBoxMin();
        std::vector<double> bnds = std::vector<double>(2*dim);
        for (int i=0; i<dim; i++) {
            bnds[i] = bMin[i];
            bnds[dim+i] = bMax[i];
        }
        return bnds;
    }

    /**
     * Writes a parameter into the global Dumux parameter map
     */
    virtual void setParameter(std::string key, std::string value) {
        auto& p = Dumux::Parameters::paramTree_();
        p[key] = value;
    }

    /**
     * Reads a parameter from the global Dumux parameter map,
     * returns an empty string if value is not set.
     */
    virtual std::string getParameter(std::string key) {
        return Dumux::getParam<std::string>(key, "");
    }
	
	
    /**
     * After the grid is created, the problem can be initialized
     *
     * initializeProblem() creates
     * (a) the Problem (initial values, bc, and source terms determined from the input file parameters)
     * (b) the solution vector
     * (c) GridVariables
     * (d) resets simtime, and internal time step ddt
     * (e) creates global index maps
     * (d) initializes the timeLoop and solvers
     *
     * The initialize values are set to the current solution,
     * i.e. can be analyzed using getSolution().
     */
    virtual void initializeProblem(double maxDt = -1) {
		maxDt_ = maxDt;
        problem = std::make_shared<Problem>(gridGeometry);
        int dof = gridGeometry->numDofs();
        x = SolutionVector(dof);

        problem->applyInitialSolution(x); // Dumux way of saying x = problem->applyInitialSolution()

        gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        gridVariables->init(x); // initialize all variables , updates volume variables to the current solution, and updates the flux variable cache

        simTime = 0; // reset
        ddt = -1;

        pointIdx = std::make_shared<Dune::GlobalIndexSet<GridView>>(grid->leafGridView(), dim); // global index mappers
        cellIdx = std::make_shared<Dune::GlobalIndexSet<GridView>>(grid->leafGridView(), 0);

        localCellIdx.clear();
        for (const auto& e : Dune::elements(gridGeometry->gridView())) {
            int eIdx = gridGeometry->elementMapper().index(e);
            int gIdx = cellIdx->index(e);
            localCellIdx[gIdx] = eIdx;
        }

		
        globalPointIdx.resize(gridGeometry->gridView().size(dim)); // number of vertices
        for (const auto& v : Dune::vertices(gridGeometry->gridView())) {
            int vIdx = gridGeometry->vertexMapper().index(v);
            int gIdx = pointIdx->index(v);
            globalPointIdx[vIdx] = gIdx;
        }
		if (maxDt_<0) { // per default value take from parameter tree
            maxDt_ = Dumux::getParam<double>("TimeLoop.MaxTimeStepSize", 3600.); // if none, default is 1h
        }
        if (ddt<1.e-6) { // happens at the first call
            ddt = Dumux::getParam<double>("TimeLoop.DtInitial",maxDt_/10); // from params, or guess something
        }

        timeLoop = std::make_shared<Dumux::CheckPointTimeLoop<double>>(/*start time*/0., ddt, /*final time*/ 3600., false); // the main time loop is moved to Python
        
        timeLoop->setMaxTimeStepSize(maxDt_);

        assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, x); // dynamic
        
		linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
		
        nonLinearSolver = std::make_shared<NonLinearSolver>(assembler, linearSolver);
        nonLinearSolver->setVerbosity(false);
        nonLinearSolverNoMPI = std::make_shared<NonLinearSolverNoMPI>(assembler, linearSolver,
								Dune::FakeMPIHelper::getCommunication());
        nonLinearSolverNoMPI->setVerbosity(false);

		problemInitialized = true;
    }

    /**
     * Sets the initial conditions, for a MPI process
     *
     *  @param init         globally shared initial data, sorted by global index [Pa]
     */
    virtual void setInitialCondition(std::vector<double> init, int eqIdx = 0) {
        if (isBox) {
            throw std::invalid_argument("SolverBase::setInitialCondition: Not implemented yet (sorry)");
        } else {
            for (const auto& e : Dune::elements(gridGeometry->gridView())) {
                int eIdx = gridGeometry->elementMapper().index(e);
                int gIdx = cellIdx->index(e);
                x[eIdx][eqIdx] = init[gIdx];
            }
        }
    }

    /**
     * Simulates the problem for time span dt, with maximal time step maxDt.
     *
     * Assembler needs a TimeLoop, so i have to create it in each solve call.
     * (could be improved, but overhead is likely to be small)
     */
    virtual void solve(double dt, bool doMPIsolve = true, bool saveInnerDumuxValues = false) {
        checkGridInitialized();
        using namespace Dumux;
        clearSaveBC();
        clearInnerFluxesAndSources();

        auto xOld = x;
		xBackUp = x; saveInnerVals();  simTimeBackUp = simTime ;
        timeLoop->reset(/*start time*/0., ddt, /*final time*/ dt, /*verbose*/ false);
        timeLoop->setMaxTimeStepSize(maxDt_);

        timeLoop->start();
		double minddt = std::min(1.,dt);//in case we have a time step < 1s																	
        do {
			// because suggestTimeStepSize() is used after timeStepSize(), the results could be > maxDt_
			// manually limit again according to maxDt_?
			if (doMPIsolve) {
				ddt = nonLinearSolver->suggestTimeStepSize(timeLoop->timeStepSize());
			} else {
				ddt = nonLinearSolverNoMPI->suggestTimeStepSize(timeLoop->timeStepSize());
			}
            ddt = std::max(ddt, minddt); // limit minimal suggestion
            
            timeLoop->setTimeStepSize(ddt); // set new dt as suggested by the newton solver, limit according to simTime
			ddt = timeLoop->timeStepSize(); // time step to use
            problem->setTime(simTime + timeLoop->time(), ddt); // pass current time to the problem ddt?

            assembler->setPreviousSolution(xOld); // set previous solution for storage evaluations

			if (doMPIsolve) {
				nonLinearSolver->solve(x, *timeLoop); // solve the non-linear system with time step control
			} else {
				nonLinearSolverNoMPI->solve(x, *timeLoop); // solve the non-linear system with time step control
			}

            xOld = x; // make the new solution the old solution
			
            if(saveInnerDumuxValues)
            {
                doSaveBC(timeLoop->timeStepSize() );
				doSaveInnerFluxesAndSources(timeLoop->timeStepSize() );
            }
			
            gridVariables->advanceTimeStep();

            timeLoop->advanceTimeStep(); // advance to the time loop to the next step
            timeLoop->reportTimeStep(); // report statistics of this time step

        } while (!timeLoop->finished());

        simTime += dt;
    }

    /**
	 * Simulates the problem for time span dt, with maximal time step maxDt.
     *
     * Assembler needs a TimeLoop, so i have to create it in each solve call.
     * (could be improved, but overhead is likely to be small)
     */
    void solveNoMPI(double dt, bool saveInnerDumuxValues = false) {
		solve(dt, false, saveInnerDumuxValues);
    }
	
	
    /**
	 * Change the maximum step size of the time loop
     */
	void setMaxTimeStepSize(double maxDt)
	{
        checkGridInitialized();
		maxDt_ = maxDt;
		if(maxDt_ > 0.)
		{
			timeLoop->setMaxTimeStepSize(maxDt_);
		}
	}
	
    /**
	 * manually (re)create linear solver. 
	 * useful to implement new solver parameters
     */
	void createLinearSolver()
	{
		linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
	}
	
    /**
	 * manually (re)create nonlinear solver. 
	 * useful to implement new solver parameters
     */
	void createNewtonSolver()
	{
		// also reset assembler?
        nonLinearSolver = std::make_shared<NonLinearSolver>(assembler, linearSolver);
        nonLinearSolver->setVerbosity(false);
        nonLinearSolverNoMPI = std::make_shared<NonLinearSolverNoMPI>(assembler, linearSolver,
								Dune::FakeMPIHelper::getCommunication());
        nonLinearSolverNoMPI->setVerbosity(false);
	}

    /**
     * Finds the steady state of the problem.
     *
     * Optionally, solve for a time span first, to get a good initial guess.
     */
    virtual void solveSteadyState() {
        checkGridInitialized();
        using namespace Dumux;

        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables); // steady state
        auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
        using NonLinearSolver = RichardsNewtonSolver<Assembler, LinearSolver>;
        auto nonLinearSolver = std::make_shared<NonLinearSolver>(assembler, linearSolver);
        nonLinearSolver->setVerbosity(false);

        assembler->setPreviousSolution(x);
        nonLinearSolver->solve(x); // solve the non-linear system

        simTime = std::numeric_limits<double>::infinity();
    }

    /**
     * Returns the Dune vertices (vtk points) of the grid for a single mpi process.
     * Gathering and mapping is done in Python.
     */
    virtual std::vector<VectorType> getPoints() {
        checkGridInitialized();
        std::vector<VectorType> points;
        points.reserve(gridGeometry->gridView().size(dim));
        for (const auto& v : vertices(gridGeometry->gridView())) {
            auto p = v.geometry().center();
            VectorType vp;
            for (int i=0; i<dim; i++) { // found no better way
                vp[i] = p[i];
            }
            points.push_back(vp);
        }
        return points;
    }

    /**
     * return the Dune element (vtk cell) centers of the grid for a single mpi process.
     * Gathering and mapping is done in Python.
     */
    virtual std::vector<VectorType> getCellCenters() {
        checkGridInitialized();
        std::vector<VectorType> cells;
        cells.reserve(gridGeometry->gridView().size(0));
        for (const auto& e : elements(gridGeometry->gridView())) {
            auto p = e.geometry().center();
            VectorType vp;
            for (int i=0; i<dim; i++) { // found no better way
                vp[i] = p[i];
            }
            cells.push_back(vp);
        }
        return cells;
    }

    /**
     * The Dune elements (vtk cell) of the grid as global vertex indices
     *
     * This is done for a single process, gathering and mapping is done in Python.
     */
    virtual std::vector<std::vector<int>> getCells() {
        checkGridInitialized();
        std::vector<std::vector<int>> cells;
        cells.reserve(gridGeometry->gridView().size(0));
        for (const auto& e : elements(gridGeometry->gridView())) {
            std::vector<int> cell;
            for (int i =0; i<e.geometry().corners(); i++) {
                int j = gridGeometry->vertexMapper().subIndex(e, i, 3);
                cell.push_back(globalPointIdx[j]);
            }
            cells.push_back(cell);
        }
        return cells;
    }

    /**
     * The volume [m3] of each element (vtk cell)
     *
     * This is done for a single process, gathering and mapping is done in Python.
     */
    virtual std::vector<double> getCellVolumes() {
        std::vector<double> vols;
        for (const auto& e : elements(gridGeometry->gridView())) {
            vols.push_back(e.geometry().volume());
        }
        return vols;
    }

    /**
     * The volume [m3] of each element (vtk cell)
     *
     * This is done for a single process, gathering and mapping is done in Python.
     */
    virtual std::vector<double> getCellVolumesCyl() {
		DUNE_THROW(Dune::NotImplemented, "Do NOT use getCellVolumesCyl, use getCellSurfacesCyl instead and multiply by the length");

		assert(false&&"Do NOT use getCellVolumesCyl, use getCellSurfacesCyl instead");//do I need to use dune-assert here?
        std::vector<double> vols;
        for (const auto& e : elements(gridGeometry->gridView())) {
            vols.push_back(e.geometry().volume()*2*e.geometry().center()[0]*3.1415);
        }
        return vols;
    }

/**
     * The Surface [m2] of each element (vtk cell)
     *
     * This is done for a single process, gathering and mapping is done in Python.
     */
    virtual std::vector<double> getCellSurfacesCyl() {
        std::vector<double> vols;
		auto points = getPoints();//get the vertices == faces of 1D domain
        for (int i = 0; i < (points.size()-1); i++) {
			double rIn = points.at(i).at(0);
			double rOut = points.at(i + 1).at(0);
            vols.push_back((rOut*rOut - rIn*rIn)* M_PI);
        }
        return vols;
    }

    /**
     * Returns the coordinate, where the DOF sit, in the same order like the solution values.
     *
     * DOF sit either at the vertices (points) for box method or
     * element centers (cell centers) for CCTpfa.
     *
     * For a single mpi process. Gathering and mapping is done in Python
     */
    virtual std::vector<VectorType> getDofCoordinates() {
        if (isBox) {
            return getPoints();
        } else {
            return getCellCenters();
        }
    }

    /**
     * Return the indices of the grid vertices.
     * Used to map the coordinates when gathered from the processes,
     * makes little sense to call directly.
     *
     * For a single mpi process. Gathering is done in Python
     */
    virtual std::vector<int> getPointIndices() {
        std::vector<int> indices;
        indices.reserve(gridGeometry->gridView().size(dim));
        for (const auto& v : vertices(gridGeometry->gridView())) {
            indices.push_back(pointIdx->index(v));
        }
        return indices;
    }

    /**
     * Return the indices of the grid elements for a single mpi process.
     * Used to map the coordinates when gathered from the processes,
     * makes little sense to call directly.
     *
     * Gathering is done in Python
     */
    virtual std::vector<int> getCellIndices() {
        std::vector<int> indices;
        indices.reserve(gridGeometry->gridView().size(0));
        for (const auto& e : elements(gridGeometry->gridView())) {
            indices.push_back(cellIdx->index(e));
        }
        return indices;
    }

    /**
     * Return the indices of the grid elements or vertices where the DOF sit.
     * Used to map the coordinates when gathered from the processes,
     * makes little sense to call directly.
     *
     * For a single mpi process. Gathering is done in Python
     */
    virtual std::vector<int> getDofIndices() {
        if (isBox) {
            return getPointIndices();
        } else {
            return getCellIndices();
        }
    }

	/**
     * Set the current solution manualy
     */
    virtual void setSolution( std::vector<double> sol, int eqIdx = 0) {
        int n = checkGridInitialized();
        std::vector<int> dofIndices = getDofIndices() ;
        for (int c = 0; c<n; c++) {
            x[c][eqIdx] = sol[dofIndices[c]] ;
        }
    }

    /**
     * Get the last saved solution. 
	 * useful when doing fixed-point iteration
	 * 
     */
    virtual std::vector<double> getOldSolution(int eqIdx = 0) {
        int n = checkGridInitialized();
        std::vector<double> sol;
        sol.resize(n);
        for (int c = 0; c<n; c++) {
            sol[c] = xBackUp[c][eqIdx];
        }
        return sol;
    }
	
    /**
     * Returns the current solution for a single mpi process.
     * Gathering and mapping is done in Python
     */
    virtual std::vector<double> getSolution(int eqIdx = 0) {
        int n = checkGridInitialized();
        std::vector<double> sol;
        sol.resize(n);
        for (int c = 0; c<n; c++) {
            sol[c] = x[c][eqIdx];
        }
        return sol;
    }

    /**
     * Returns the current solution at a cell index
     * for all mpi processes
     *
     * TODO currently works only for CCTpfa
     */
    double getSolutionAt(int gIdx, int eqIdx = 0) {
        if (isBox) {
            throw std::invalid_argument("SolverBase::getSolutionAt: Not implemented yet (sorry)");
        } else {
            double y = 0.;
            if (localCellIdx.count(gIdx)>0) {
                int eIdx = localCellIdx[gIdx];
                y = x[eIdx][eqIdx];
            }
            return gridGeometry->gridView().comm().min(y);
            // I would prefer sum(y), BUT more than one process can have this cell (borders overlap)
        }
    }

    /**
     * Returns the maximal flux (over the boundary scvfs) of an element, given by its global element index,
     * for all mpi processes
	 * ATT: this only gives the flux AT THE END of the solve() function
     */
    virtual double getNeumann(int gIdx, int eqIdx = 0) {
        double f = 0.;
        if (localCellIdx.count(gIdx)>0) {
            int eIdx = localCellIdx[gIdx];
            auto e = gridGeometry->element(eIdx);
            auto fvGeometry = Dumux::localView(*gridGeometry); // soil solution -> volume variable
            fvGeometry.bindElement(e);

            auto elemVolVars = Dumux::localView(gridVariables->curGridVolVars());
            elemVolVars.bindElement(e, fvGeometry, x);
			auto elemFluxVars = Dumux::localView(gridVariables->gridFluxVarsCache());
			//bindElement throws error. neumann() does not use elemFluxVars anyway
			//elemFluxVars.bindElement(e, fvGeometry, elemVolVars);

            for (const auto& scvf : scvfs(fvGeometry)) {
                if (scvf.boundary()) {
                    double n = problem->neumann(e, fvGeometry, elemVolVars, elemFluxVars, scvf)[eqIdx];  // [ kg / (m2 s)]
                    f = (std::abs(n) > std::abs(f)) ? n : f;
                }
            }
        }
        return gridGeometry->gridView().comm().sum(f); // so clever
    }

    /**
     * Return the neumann fluxes of the current solution for each boundary element (as mean over all scvfs, TODO why, net-flux would make more sense).
     *
     * For a single mpi process. Gathering is done in Python
     */
    virtual std::map<int, double> getAllNeumann(int eqIdx = 0) {
        std::map<int, double> fluxes;
        for (const auto& e : Dune::elements(gridGeometry->gridView())) { // soil elements
            double f = 0.;
            auto fvGeometry = Dumux::localView(*gridGeometry); // soil solution -> volume variable
            fvGeometry.bindElement(e);
            int c = 0;
            for (const auto& scvf : scvfs(fvGeometry)) {
                if (scvf.boundary()) {
                    c++;
                    auto elemVolVars  = Dumux::localView(gridVariables->curGridVolVars());
                    elemVolVars.bindElement(e, fvGeometry, x);
					auto elemFluxVars = Dumux::localView(gridVariables->gridFluxVarsCache());
					//elemFluxVars.bindElement(e, fvGeometry, elemVolVars);

                    f += problem->neumann(e, fvGeometry, elemVolVars, elemFluxVars, scvf)[eqIdx]; // [kg / (m2 s)]
                }
            }
            if (c>0) {
                fluxes[cellIdx->index(e)] = f/c; // mean value
            }
        }
        return fluxes;
    }

    /**
     * Returns the net flux [kg/s]. TODO crashes, no idea why
     *
     * partly from velocityoutput.hh
     *
     * For a single mpi process. Gathering is done in Python
     */
    virtual std::vector<double> getNetFlux(int eqIdx = 0) {
        int n = checkGridInitialized();
        std::vector<double> fluxes;
        fluxes.resize(n);

        auto elemVolVars = Dumux::localView(gridVariables->curGridVolVars());
        auto fvGeometry = Dumux::localView(*gridGeometry); // soil solution -> volume variable
        auto elemFluxVarsCache = Dumux::localView(gridVariables->gridFluxVarsCache());

        // the upwind term to be used for the volume flux evaluation, currently not needed
        // auto upwindTerm = [eqIdx](const auto& volVars) { return volVars.mobility(eqIdx); };

        for (const auto& e : Dune::elements(gridGeometry->gridView())) { // soil elements

            fvGeometry.bindElement(e);
            elemVolVars.bindElement(e, fvGeometry, x);

            double f = 0.;
            for (const auto& scvf : scvfs(fvGeometry)) {

                if (!scvf.boundary()) {
                    // bind the element flux variables cache
                    elemFluxVarsCache.bindElement(e, fvGeometry, elemVolVars);
                    //                    FluxVariables fluxVars;
                    //                    fluxVars.init(*problem, e, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
                    //                    f += fluxVars.advectiveFlux(eqIdx, upwindTerm);
                }
                else { }

            }
            fluxes[cellIdx->index(e)] = f;
        }
        return fluxes;
    }

    /**
     * Picks a cell and returns its global element cell index
     *
     * The lucky rank who found it, maps the local index to a global one,
     * and broadcasts to the others
     */
    virtual int pickCell(VectorType pos) {//TODO: do a nonMPI version?
        checkGridInitialized();
        if (periodic) {
            auto b = getGridBounds();
            for (int i = 0; i < 2; i++) { // for x and y, not z
                double minx = b[i];
                double xx = b[i+3]-minx;
                if (!std::isinf(xx)) { // periodic in x
                    pos[i] -= minx; // start at 0
                    if (pos[i]>=0) {
                        pos[i] = pos[i] - int(pos[i]/xx)*xx;
                    } else {
                        pos[i] = pos[i] + int((xx-pos[i])/xx)*xx;
                    }
                    pos[i] += minx;
                }
            }
        }
        auto& bBoxTree = gridGeometry->boundingBoxTree();
        Dune::FieldVector<double, dim> p;
        for (int i=0; i<dim; i++) {
            p[i] = pos[i];
        }
        // std::cout << "point: " << pos[0]<<", "<< pos[1] <<","<< pos[2] << " in box "; // <<  getGridBounds();
        auto entities = Dumux::intersectingEntities(p, bBoxTree);
        int gIdx = -1;
        if (!entities.empty()) {
            auto element = bBoxTree.entitySet().entity(entities[0]);
            gIdx = cellIdx->index(element);
        }
        gIdx = gridGeometry->gridView().comm().max(gIdx); // so clever
        return gIdx;
    }

    /**
     * Picks a cell and returns its global element cell index @see pickCell
     */
    virtual int pick(VectorType x) {
        return pickCell(x);
    }

    /**
     * Quick overview
     */
    virtual std::string toString() {
        std::ostringstream msg;
        msg << "DuMux Solver " << dim << "D ";
        if (isBox) {
            msg << "(box method)";
        } else {
            msg << "(CCTpfa method)";
        }
        if (gridGeometry) {
            msg << " with " << gridGeometry->numDofs() << " DOF";
        } else {
            msg << " no problem initialized";
        }
        if (maxRank>=0) {
            msg << " on process rank "<< rank+1 << "/" << maxRank;
        }
        if (simTime>0) {
            msg << "\nSimulation time is "<< simTime/3600/24 << " days, current internal time step is "<< ddt/3600/24 << " days";
        }
        return msg.str();
    }

    /**
     * Checks if the grid was initialized, and returns number of local dof
     * i.e. createGrid() or createGrid1d() was called
     */
    virtual int checkGridInitialized() {
        if (!gridGeometry) {
            throw std::invalid_argument("SolverBase::checkGridInitialized: Grid not initialized, call createGrid() or createGrid1d() first");
        }
        if (this->isBox) {
            return this->gridGeometry->gridView().size(dim);
        } else {
            return this->gridGeometry->gridView().size(0);
        }
    }

    /**
     * Checks if the problem was initialized
     */
    virtual bool checkProblemInitialized() {
        return problemInitialized;
    }

    virtual std::map<int, int> getGlobal2localCellIdx() {
        return localCellIdx;
    }
	
    virtual std::vector<int> getLocal2globalPointIdx() {
        return globalPointIdx;
    }
    std::vector<std::vector<std::vector<double>>> inSources;//[time][cellidx][eqIdx]
    std::vector<std::vector<std::vector<double>>> inFluxes;//[time][cellidx][eqIdx]
    std::vector<std::vector<int>> face2CellIds;
	std::vector<double> inFluxes_time;
    std::vector<double> inFluxes_ddt;
    
    // std::vector<std::vector<double>> getFlux_10c_() {
        // // Compute the total flux [mol] per face during the last solve() call
        
        // // Create a 2D vector to store the total flux for each cell and component
        // std::vector<std::vector<double>> totalFlux(inFluxes[0].size(), std::vector<double>(inFluxes[0][0].size(), 0.0));//[cellidx][eqIdx]

        // // Iterate over each sub-timestep
        // for (size_t i = 0; i < inFluxes.size(); ++i) {
            // // Iterate over each cell
            // for (size_t j = 0; j < inFluxes[i].size(); ++j) {
                // // Iterate over each component within the cell
                // for (size_t k = 0; k < inFluxes[i][j].size(); ++k) {
                    // // Calculate total flux per component for each cell
                    // // Return the negated total flux
                    // totalFlux[j][k] -= inFluxes[i][j][k] * inFluxes_ddt[i];
                // }
            // }
        // }
        // return totalFlux;
    // }
    
    
    // std::vector<std::vector<double>> getFlux_10c() {
        // assert(dimWorld == 3); // only works for 3d grids

        // auto ff10c_ = getFlux_10c_();                  // Local thread flux for each cell
        // auto f2cidx_ = face2CellIds.at(0);             // Corresponding cell for each face

        // std::unordered_map<int, double> cellFluxMap;   // Map to store the sum of fluxes per cell
        // std::unordered_map<int, int> cellCountMap;     // Map to count the number of faces per cell

        // // Sum values per cell and count the faces
        // for (size_t i = 0; i < f2cidx_.size(); ++i) {
            // int cellId = f2cidx_[i];
            // cellFluxMap[cellId] += ff10c_[i];
            // cellCountMap[cellId]++;
        // }

        // std::vector<double> ff10c;
        // std::vector<int> f2cidx;

        // // Only keep cells where all 6 faces belong to this thread
        // for (const auto& [cellId, count] : cellCountMap) {
            // if (count == 6) {
                // ff10c.push_back(cellFluxMap[cellId]);
                // f2cidx.push_back(cellId);
            // }
        // }

        // // Get global indices
        // auto dofind = getDofIndics();
        // std::vector<int> f2cidx_g;
        // f2cidx_g.reserve(f2cidx.size());
        // int maxlid = -1;
        // for (int idx : f2cidx) {
            // f2cidx_g.push_back(dofind[idx]);
            // maxlid = maxlid < dofind[idx] ? dofind[idx]:maxlid;
        // }

        // // Gather global indices and fluxes
        // std::vector<int> f2cidx_gAll;
        // f2cidx_gAll.reserve(maxgCid);
        // maxgCid = gridGeometry->gridView().comm().max(maxlid);
        // for (size_t gCid = 0; gCid < maxlid ; gCid ++) {
        
            // if (localCellIdx.count(gIdx)>0) {
                // int eIdx = localCellIdx[gIdx];
                // y = x[eIdx][eqIdx];
            // }
            
            // f2cidx_gAll.at(gCid) = 
        // }
        // auto f2cidx_gAll = gather(f2cidx_g);
        // auto ff10c_All = gather(ff10c);
        

        // std::vector<std::vector<double>> flux10cCell;

        // if (rank == 0) {
            // // Flatten arrays and remove duplicates
            // std::unordered_set<int> f2cidx_gAll_unique(f2cidx_gAll.begin(), f2cidx_gAll.end());
            // std::vector<double> ff10c_All_unique;

            // for (int idx_ : f2cidx_gAll_unique) {
                // auto it = std::find(f2cidx_gAll.begin(), f2cidx_gAll.end(), idx_);
                // if (it != f2cidx_gAll.end()) {
                    // size_t index = std::distance(f2cidx_gAll.begin(), it);
                    // ff10c_All_unique.push_back(ff10c_All[index]);
                // }
            // }

            // // Transpose for [comp][cell]
            // flux10cCell = transpose(ff10c_All_unique, numComp, numberOfCellsTot);

            // //assert(flux10cCell.size() == static_cast<size_t>(numComp));
            // //assert(flux10cCell[0].size() == static_cast<size_t>(numberOfCellsTot));

            // // Convert water flux from mol to cm3
            // double molarDensityWat = densityWat / molarMassWat; // [mol/cm3]
            // for (double& flux : flux10cCell[0]) {
                // flux /= molarDensityWat;
            // }
        // }

        // return flux10cCell;
    // }
    
	// reset to value before the last call to solve()
    virtual void reset() {
        checkGridInitialized();
		x = xBackUp;
		simTime = simTimeBackUp;
		resetInnerVals();
	}
	
	// reset to values at the last call to saveManual()
    virtual void resetManual() {
        checkGridInitialized();
		x = xBackUpManual;
		resetInnerValsManual();
	}
	
    virtual void save() {
        checkGridInitialized();
		xBackUp = x;
	}
    
	// save/store value to reset them with resetManual()   
    virtual void saveManual() {
        checkGridInitialized();
		xBackUpManual = x;
		simTimeBackUpManual = simTime ;
		saveInnerValsManual();
	}
	
    virtual void resetInnerVals() {}
	virtual void saveInnerVals() {}
    virtual void resetInnerValsManual() {}
    virtual void saveInnerValsManual() {}
	
	
protected:

    using Grid = typename Problem::Grid;
    using FVGridGeometry = typename Problem::GridGeometry;
    using SolutionVector = typename Problem::SolutionVector;
    using GridVariables = typename Problem::GridVariables;
    using FluxVariables = typename Problem::FluxVariables;
	using ElementFluxVariablesCache = typename Problem::ElementFluxVariablesCache;

    using GridData = Dumux::GridData<Grid>;
    using GridView = typename Grid::Traits::LeafGridView;

    std::shared_ptr<Grid> grid;
    std::shared_ptr<GridData> gridData;
    std::shared_ptr<FVGridGeometry> gridGeometry;
    std::shared_ptr<Problem> problem;
    std::shared_ptr<GridVariables> gridVariables;

    std::shared_ptr<Dune::GlobalIndexSet<GridView>> pointIdx; // global index mappers
    std::shared_ptr<Dune::GlobalIndexSet<GridView>> cellIdx; // global index mappers
    std::map<int, int> localCellIdx; // global to local index mapper
    std::vector<int> globalPointIdx; // local to global index mapper

    SolutionVector x;
	SolutionVector xBackUp;
    SolutionVector xBackUpManual;
    double simTimeBackUp;
    double simTimeBackUpManual;
	double maxDt_;
	
	
    /**
     * @see richards_cyl::clearSaveBC
     */	 
    virtual void clearSaveBC() {}
	
    /**
     * @see richards_cyl::doSaveBC
     */	 
    virtual void doSaveBC(double ddt_current) {}
    
	
    /**
     * for 1d3d coupling (and to compute the error),
	 * useful to get from dumux 
	 * the inter-cell flow and sources
	 * only implemented for richards10c(_cyl)
     */
    void doSaveInnerFluxesAndSources(double ddt_current) {
		//need to save it per face and not per cell (for when we have mpi)
        std::vector<std::vector<double>> inFluxes_i = getFlux();//per face
        std::vector<int> face2CellId_i = getFluxScvf_idx();//per face
        std::vector<std::vector<double>> inSources_i = getSource();//per cell
		
		inSources.push_back(inSources_i);
        face2CellIds.push_back(face2CellId_i);
        inFluxes.push_back(inFluxes_i);
        inFluxes_ddt.push_back(ddt_current);
		resetSetFaceFlux();
    }
	
	
    std::vector<std::vector<double>> getFlux() {		
    	std::vector<NumEqVector> flux_ = getProblemFlux();
		int numComp_ = numComp();
		std::vector<double> flx_row(numComp_);
		std::vector<std::vector<double>> flux(flux_.size(), flx_row);
		
		for(int faceIdx_ = 0; faceIdx_ < flux.size(); faceIdx_ ++)
		{
			for(int eqIdx = 0; eqIdx < numComp_; eqIdx ++)
			{
				flux.at(faceIdx_).at(eqIdx) = flux_.at(faceIdx_)[eqIdx];
			}
		}
		return flux;
    }
	
    std::vector<std::vector<double>> getSource() {		
    	std::vector<NumEqVector> source_ = getProblemSource();
		int numComp_ = numComp();
		std::vector<double> src_row(numComp_);
		std::vector<std::vector<double>> source(source_.size(), src_row);
		
		for(int cellIdx = 0; cellIdx < source.size(); cellIdx ++)
		{
			for(int eqIdx = 0; eqIdx < numComp_; eqIdx ++)
			{
				source.at(cellIdx).at(eqIdx) = source_.at(cellIdx)[eqIdx];
			}
		}
		return source;
    }
	
	virtual std::vector<NumEqVector> getProblemFlux()
	{
		std::vector<NumEqVector> flux_;
		return flux_;
	}
	virtual std::vector<NumEqVector> getProblemSource()
	{
		std::vector<NumEqVector> source_;
		return source_;
	}
	
	virtual std::vector<int> getFluxScvf_idx()
	{
		std::vector<int> face2CellId;
		return face2CellId;
	}
	
    void clearInnerFluxesAndSources() {    
        inSources.clear();
        inFluxes.clear();
		face2CellIds.clear();
        inFluxes_time.clear();
        inFluxes_ddt.clear();
    }
	
	void virtual resetSetFaceFlux(){}
	
	
	//! true if on the point lies on the upper boundary
	bool onUpperBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[dimWorld - 1] > this->gridGeometry->bBoxMax()[dimWorld - 1] - eps_;
	}

	//! true if on the point lies on the upper boundary
	bool onLowerBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[dimWorld - 1] < this->gridGeometry->bBoxMin()[dimWorld - 1] + eps_;
	}
	static constexpr double eps_ = 1.e-7;
	
};

/**
 * pybind11
 */
template<class Problem, class Assembler, class LinearSolver, int dim = 3>
void init_solverbase(py::module &m, std::string name) {
    using Solver = SolverBase<Problem, Assembler, LinearSolver, dim>; // choose your destiny
    py::class_<Solver>(m, name.c_str())
            .def(py::init<>()) // initialization
			.def("numComp",  &Solver::numComp)
			.def("numFluidComp",  &Solver::numFluidComp)
            .def("initialize", &Solver::initialize, py::arg("args_") = std::vector<std::string>(0),
													py::arg("verbose") = true,py::arg("doMPI") = true)
            .def("createGrid", (void (Solver::*)(std::string)) &Solver::createGrid) // overloads, defaults , py::arg("modelParamGroup") = ""
            .def("createGrid", (void (Solver::*)(std::array<double, dim>, std::array<double, dim>, std::array<int, dim>, bool)) &Solver::createGrid) // overloads, defaults , py::arg("boundsMin"), py::arg("boundsMax"), py::arg("numberOfCells"), py::arg("periodic") = false
            .def("createGrid1d", &Solver::createGrid1d)
            // .def("createGrid3d", &Solver::createGrid3d)
            .def("readGrid", &Solver::readGrid)
            .def("getGridBounds", &Solver::getGridBounds)
            .def("setParameter", &Solver::setParameter)
            .def("getParameter", &Solver::getParameter)
            .def("printParams", &Solver::printParams)
			.def("setMaxTimeStepSize",&Solver::setMaxTimeStepSize)
			.def("createLinearSolver",&Solver::createLinearSolver)
			.def("createNewtonSolver",&Solver::createNewtonSolver)
            .def("initializeProblem", &Solver::initializeProblem, py::arg("maxDt") = -1)
            .def("setInitialCondition", &Solver::setInitialCondition, py::arg("init"), py::arg("eqIdx") = 0)
			.def("reset", &Solver::reset)
			.def("resetManual", &Solver::resetManual)
			.def("saveManual", &Solver::saveManual)
			.def("save", &Solver::save)
            // simulation
            .def("solve", &Solver::solve, py::arg("dt"), py::arg("doMPIsolve")=true, 
					py::arg("saveInnerDumuxValues")=false)
			.def("solveNoMPI", &Solver::solveNoMPI, py::arg("dt"), 
					py::arg("saveInnerDumuxValues")=false)
            .def("solveSteadyState", &Solver::solveSteadyState)
            // post processing (vtk naming)
            .def("getPoints", &Solver::getPoints) //
            .def("getCellCenters", &Solver::getCellCenters)
            .def("getCells", &Solver::getCells)
            .def("getCellVolumes", &Solver::getCellVolumes)
            .def("getCellVolumesCyl", &Solver::getCellVolumesCyl)
            .def("getCellSurfacesCyl", &Solver::getCellSurfacesCyl)
            .def("getDofCoordinates", &Solver::getDofCoordinates)
            .def("getPointIndices", &Solver::getPointIndices)
            .def("getCellIndices", &Solver::getCellIndices)
			.def("getGlobal2localCellIdx", &Solver::getGlobal2localCellIdx) // localCellIdx
            .def("getLocal2globalPointIdx", &Solver::getLocal2globalPointIdx) // globalPointIdx
            .def("getDofIndices", &Solver::getDofIndices)
            .def("getOldSolution", &Solver::getOldSolution, py::arg("eqIdx") = 0)
            .def("getSolution", &Solver::getSolution, py::arg("eqIdx") = 0)
            .def("getSolutionAt", &Solver::getSolutionAt, py::arg("gIdx"), py::arg("eqIdx") = 0)            
			.def("setSolution", &Solver::setSolution)
            .def("getNeumann", &Solver::getNeumann, py::arg("gIdx"), py::arg("eqIdx") = 0)
            .def("getAllNeumann", &Solver::getAllNeumann, py::arg("eqIdx") = 0)
            .def("getNetFlux", &Solver::getNetFlux, py::arg("eqIdx") = 0)
            .def("pickCell", &Solver::pickCell)
            .def("pick", &Solver::pick)
            // members
			.def_readonly("face2CellIds", &Solver::face2CellIds)
			.def_readonly("inSources", &Solver::inSources) // read only
            .def_readonly("inFluxes", &Solver::inFluxes) // read only
            .def_readonly("inFluxes_time", &Solver::inFluxes_time) // read only
            .def_readonly("inFluxes_ddt", &Solver::inFluxes_ddt) // read only
            .def_readonly("dimWorld", &Solver::dimWorld) // read only
            .def_readonly("simTime", &Solver::simTime) // read only
            .def_readwrite("ddt", &Solver::ddt) // initial internal time step
            .def_readonly("rank", &Solver::rank) // read only
            .def_readonly("maxRank", &Solver::maxRank) // read only
            .def_readonly("numberOfCells", &Solver::numberOfCells) // read only
            .def_readonly("periodic", &Solver::periodic) // read only
            // useful
            .def("__str__",&Solver::toString)
            .def("checkGridInitialized", &Solver::checkGridInitialized)
            .def("checkProblemInitialized", &Solver::checkProblemInitialized);
}

#endif

/**
 * lacking polymorphism I have not found a way to make the grid dynamic.
 * you need to choose the grid at compile time,
 * and I don't know how to pass it via CMakeLists.txt building the Python binding
 */

//using Grid = Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>;
//using GridView = typename Dune::YaspGridFamily<3,Dune::EquidistantOffsetCoordinates<double,3>>::Traits::LeafGridView;
//using Scalar = double;
//using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
//using VertexMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
//using MapperTraits = Dumux::DefaultMapperTraits<GridView, ElementMapper, VertexMapper>;
//using FVGridGeometry = Dumux::BoxFVGridGeometry<double, GridView, /*enableCache*/ true, Dumux::BoxDefaultGridGeometryTraits<GridView, MapperTraits>>;
//using SolutionVector =  Dune::BlockVector<GetPropType<TypeTag, Properties::PrimaryVariables>>; // in fvproperties
