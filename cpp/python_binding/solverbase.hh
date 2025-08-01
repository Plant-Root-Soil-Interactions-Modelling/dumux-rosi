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
// #include <dumux/porousmediumflow/richards/newtonsolver.hh>
#include "../../dumux/porousmediumflow/richards/newtonsolver.hh"
#include <dumux/porousmediumflow/velocity.hh>


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
    using VelocityVector = std::vector<Dune::FieldVector<double, dim>>;
	using GlobalPosition = typename Problem::GlobalPosition; 
	
    int numComp(){return nEV.size();}    
	int numFluidComp(){return Problem::FluidSystem::numComponents;}
    NumEqVector nEV;
    bool isBox = Problem::isBox; // numerical method
    int dimWorld = dim;

    double simTime = 0;
    double ddt = -1; // internal time step, minus indicates that its uninitialized
    double maxDt_ = -1; // max time step
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
        periodic = ((pstr.at(0)=='t') || (pstr.at(0)=='T')); // always x,y, not z: pstr = "true ..."
        if ((!periodic) and (pstr.size()>4)) { // y = true
            periodic =  (pstr.at(6)=='t') || (pstr.at(6)=='T'); // pstr = "false true false"
        }
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
    //     * Creates a 1D grid with number of cells = [1,1,points.size()-1] where points are the cem
    // 		* TODO needs testing
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
        problem = std::make_shared<Problem>(gridGeometry);
        int dof = gridGeometry->numDofs();
        x = SolutionVector(dof);

        problem->applyInitialSolution(x); // Dumux way of saying x = problem->applyInitialSolution()

        gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        gridVariables->init(x); // initialize all variables , updates volume variables to the current solution, and updates the flux variable cache
        flowVelocities = std::make_shared<PorousMediumFlowVelocity>(*gridVariables);

        simTime = 0; // reset
        ddt = -1;

        pointIdx = std::make_shared<Dune::GlobalIndexSet<GridView>>(grid->leafGridView(), dim); // global index mappers
        faceIdx = std::make_shared<Dune::GlobalIndexSet<GridView>>(grid->leafGridView(), 1);
        cellIdx = std::make_shared<Dune::GlobalIndexSet<GridView>>(grid->leafGridView(), 0);

        localCellIdx.clear();
        for (const auto& e : Dune::elements(gridGeometry->gridView())) {
            int eIdx = gridGeometry->elementMapper().index(e);
            int gIdx = cellIdx->index(e);
            localCellIdx[gIdx] = eIdx;
        }
		int nFaces = gridGeometry->numScvf();
        auto fvGeometry = Dumux::localView(*gridGeometry);
		face2CellIdx = std::vector<int>(nFaces,-1);	// to sum face flow per cell
		 for (const auto& e : Dune::elements(gridGeometry->gridView(),Dune::Partitions::interior))
        {
            fvGeometry.bindElement(e);
			//int i = 0;
            for (const auto& scvf : scvfs(fvGeometry))
            {
				// local face to local cell indx
				face2CellIdx.at(scvf.index()) = fvGeometry.scv(scvf.insideScvIdx()).dofIndex();
            }
        }
		scvfInnerFluxes.assign(nEV.size(), std::vector<double>(nFaces, 0.));
		scvfBoundaryFluxes.assign(nEV.size(), std::vector<double>(nFaces, 0.));
		scvSources.assign(nEV.size(),std::vector<double>(gridGeometry->numScv(),0.));


        globalPointIdx.resize(gridGeometry->gridView().size(dim)); // number of vertices
        for (const auto& v : Dune::vertices(gridGeometry->gridView())) {
            int vIdx = gridGeometry->vertexMapper().index(v);
            int gIdx = pointIdx->index(v);
            globalPointIdx[vIdx] = gIdx;
        }
		if (maxDt<0) { // per default value take from parameter tree
            maxDt = Dumux::getParam<double>("TimeLoop.MaxTimeStepSize", 3600.); // if none, default is 1h
        }
        if (ddt<1.e-6) { // happens at the first call
            ddt = Dumux::getParam<double>("TimeLoop.DtInitial",maxDt/10); // from params, or guess something
        }

        timeLoop = std::make_shared<Dumux::CheckPointTimeLoop<double>>(/*start time*/0., ddt, /*final time*/ 3600., false); // the main time loop is moved to Python

        timeLoop->setMaxTimeStepSize(maxDt);
    	maxDt_ = maxDt; // save to re-set at each solve() function

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

		scvSources.assign(nEV.size(),std::vector<double>(gridGeometry->numScv(),0.));
		scvfBoundaryFluxes.assign(nEV.size(), std::vector<double>(gridGeometry->numScvf(), 0.));
		scvfInnerFluxes.assign(nEV.size(), std::vector<double>(gridGeometry->numScvf(), 0.));

        auto xOld = x;
		xBackUp = x; saveInnerVals();  simTimeBackUp = simTime ;

        timeLoop->reset(/*start time*/0., ddt, /*final time*/ dt, /*verbose*/ false);
        timeLoop->setMaxTimeStepSize(maxDt_); // set maxDt after reset()

        timeLoop->start();
		double minddt = std::min(1.,dt);//in case we have a time step < 1s
        do {
			// because suggestTimeStepSize() is used after timeStepSize(), the results could be > maxDt
			// manually limit again according to maxDt?
			if (doMPIsolve) {
				ddt = nonLinearSolver->suggestTimeStepSize(timeLoop->timeStepSize());
			} else {
				ddt = nonLinearSolverNoMPI->suggestTimeStepSize(timeLoop->timeStepSize());
			}
            ddt = std::max(ddt, minddt); // limit minimal suggestion
            timeLoop->setTimeStepSize(ddt); // set new dt as suggested by the newton solver

            ddt = timeLoop->timeStepSize(); // selected ddt
            problem->setTime(simTime + timeLoop->time(), ddt); // pass current time to the problem ddt?

            assembler->setPreviousSolution(xOld); // set previous solution for storage evaluations

			if (doMPIsolve) {
				nonLinearSolver->solve(x, *timeLoop); // solve the non-linear system with time step control
			} else {
				nonLinearSolverNoMPI->solve(x, *timeLoop); // solve the non-linear system with time step control
			}

			// auto dummy = this->getVelocities(); /////////////////////////////////

            xOld = x; // make the new solution the old solution

            if(saveInnerDumuxValues)
            {
				getScvfFluxesAtT(timeLoop->timeStepSize(),dt) ;
				getScvSourcesAtT(timeLoop->timeStepSize(),dt) ;
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
        std::vector<double> vols;
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
     * TODO replace min
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
	 *
	 * TODO revise
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
					const auto& bcTypes = problem->boundaryTypes(e, scvf);
					double n;
					if (bcTypes.hasNeumann()){
						n = problem->neumann(e, fvGeometry, elemVolVars, elemFluxVars, scvf)[eqIdx];  // [ kg / (m2 s)]
					}
					
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
     *
     * Todo revise
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
     * For a single mpi process. Gathering is done in Python
	 * [ kg /m^3/s] or [ mol / m^3/s]
     */
    void getScvSourcesAtT(double dt, double outer_dt) {
		checkGridInitialized();

        auto fvGeometry = Dumux::localView(*gridGeometry);

        auto elemVolVars = Dumux::localView(gridVariables->curGridVolVars());


        for (const auto& e : elements(gridGeometry->gridView())) { //, Dune::Partitions::interior

            fvGeometry.bind(e);

            elemVolVars.bind(e, fvGeometry, x);

			for (const auto& scv : scvs(fvGeometry))
            {
				double pos0 = 1;
				if(dimWorld == 1){
					pos0 = scv.center()[0];
				}
				NumEqVector scvfSource_(0.0);
				scvfSource_ = problem->source(e, fvGeometry, elemVolVars, scv); // [ kg / (m^3 \cdot s)] or [ mol / (m^3 \cdot s)]

				for(int eqIdx = 0; eqIdx < nEV.size(); eqIdx ++)
				{
					scvSources[eqIdx][scv.dofIndex()] += scvfSource_[eqIdx]/pos0*dt/outer_dt;
				}
			}
		}

	}
    /**
     * For a single mpi process. Gathering is done in Python
	 * [ kg/s] or [ mol/s]
     */
    void getScvfFluxesAtT(double dt, double outer_dt) {

        checkGridInitialized();

        auto fvGeometry = Dumux::localView(*gridGeometry);

        auto elemVolVars = Dumux::localView(gridVariables->curGridVolVars());

		auto elemFluxVarsCache = Dumux::localView(gridVariables->gridFluxVarsCache());

		//localResidual = std::make_shared<assembler->localResidual()>;


        for (const auto& e : elements(gridGeometry->gridView(), Dune::Partitions::interior)) {

            fvGeometry.bind(e);

            elemVolVars.bind(e, fvGeometry, x);
            // std::cout << "f\n" << std::flush;


            elemFluxVarsCache.bind(e, fvGeometry, elemVolVars); // check workspace
            //std::cout << "g" << std::flush;

			//LocalAssembler localAssembler(assembler, e, x);
			//std::make_shared<Assembler>

            for (const auto& scvf : scvfs(fvGeometry))
            {
				double pos0 = 1;
				//double scvf_area = scvf.area();
				if(dimWorld == 1){
					pos0 = scvf.center()[0];
					//scvf_area = 2 * M_PI * pos0 * segLength;//m2 TODO handle 1D area?
				}
				NumEqVector scvfFlux_(0.0);
				if(scvf.boundary())
				{
					const auto& bcTypes = problem->boundaryTypes(e, scvf);
					if (bcTypes.hasNeumann()){
						scvfFlux_ = problem->neumann(e, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
					}
					for(int eqIdx = 0; eqIdx < nEV.size(); eqIdx ++)
					{
						scvfBoundaryFluxes[eqIdx][scvf.index()] += scvfFlux_[eqIdx]/pos0*dt/outer_dt;
					}
				}
				else{
					scvfFlux_ = assembler->localResidual().computeFlux(*problem, e, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
					for(int eqIdx = 0; eqIdx < nEV.size(); eqIdx ++)
					{
						scvfInnerFluxes[eqIdx][scvf.index()] += scvfFlux_[eqIdx]/pos0*dt/outer_dt;
					}
				}
            }
        }
    }

	std::vector<double> getFaceSurfaces() {

        checkGridInitialized();

        auto fvGeometry = Dumux::localView(*gridGeometry);
		std::vector<double> scvfSurface(gridGeometry->numScvf());
        for (const auto& e : elements(gridGeometry->gridView(), Dune::Partitions::interior)) {

            fvGeometry.bind(e);


            for (const auto& scvf : scvfs(fvGeometry))
            {
				scvfSurface[scvf.index()] = scvf.area();
            }
        }
		return scvfSurface;
    }

	std::vector<double> getCylFaceCoordinates() {

        checkGridInitialized();

        auto fvGeometry = Dumux::localView(*gridGeometry);
		std::vector<double> scvfSurface(gridGeometry->numScvf());
        for (const auto& e : elements(gridGeometry->gridView(), Dune::Partitions::interior)) {

            fvGeometry.bind(e);


            for (const auto& scvf : scvfs(fvGeometry))
            {
				scvfSurface[scvf.index()] = scvf.center()[0];
            }
        }
		return scvfSurface;
    }

        /**
     * For a single mpi process. Gathering is done in Python
     */
    virtual std::vector<std::array<double,3>> getVelocities(int eqIdx = 0) {

        // std::cout << "a\n" << std::flush;
        int n = checkGridInitialized();
        VelocityVector velocity = VelocityVector();
        velocity.resize(n);
        // std::cout << "b\n" << std::flush;

        auto fvGeometry = Dumux::localView(*gridGeometry);
        // std::cout << "c\n" << std::flush;

        auto elemVolVars = Dumux::localView(gridVariables->curGridVolVars());
		auto elemFluxVarsCache = Dumux::localView(gridVariables->gridFluxVarsCache());
        //std::cout << "d\n" << std::flush;

        for (const auto& e : elements(gridGeometry->gridView(), Dune::Partitions::interior)) {

            fvGeometry.bind(e);
            // std::cout << "e\n" << std::flush;

            elemVolVars.bind(e, fvGeometry, x);
            // std::cout << "f\n" << std::flush;


            elemFluxVarsCache.bind(e, fvGeometry, elemVolVars); // check workspace
            //std::cout << "g" << std::flush;

            flowVelocities->calculateVelocity(velocity, e, fvGeometry, elemVolVars, elemFluxVarsCache, 0 ); // veloctiy.hh
            //std::cout << "h\n" << std::flush;

//            // where does velocity go?
//            int idx = cellIdx->index(e);
//            fluxes[idx] = velocity;

        }

        std::vector<std::array<double,3>> v;
        for (auto v_ : velocity) {
            std::array<double,3> a = {v_[0], v_[1], v_[2]};
            v.push_back(a);
        }

        return v;
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
        // std::cout << "point: " << pos[0]<<", "<< pos[1] <<","<< pos[2] << " in box "<<  ", " << periodic << "\n";
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


	bool useMoles()
	{
		return problem->useMoles;
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

    virtual std::vector<int> getFace2CellIdx() {
        return face2CellIdx;
    }
	/*
	 * [ kg/s] or [ mol/s]
	 */
    virtual std::vector<std::vector<double>> getScvfInnerFluxes() {
        return scvfInnerFluxes;
    }
	/*
	 * [ kg /m^2/s] or [ mol / m^2/s]
	 */
    virtual std::vector<std::vector<double>> getScvfBoundaryFluxes() {
        return scvfBoundaryFluxes;
    }
	/*
	 * [ kg /m^3/s] or [ mol / m^3/s]
	 */
    virtual std::vector<std::vector<double>> getScvSources() {
        return scvSources;
    }

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

	//! true if on the point lies on the upper boundary
	// TODO: use to get specific fluxes
	bool onUpperBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[dimWorld - 1] > this->gridGeometry->bBoxMax()[dimWorld - 1] - eps_;
	}

	//! true if on the point lies on the upper boundary
	bool onLowerBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[dimWorld - 1] < this->gridGeometry->bBoxMin()[dimWorld - 1] + eps_;
	}
	static constexpr double eps_ = 1.e-7;

protected:

    using Grid = typename Problem::Grid;
    using FVGridGeometry = typename Problem::GridGeometry;
    using SolutionVector = typename Problem::SolutionVector;
    using GridVariables = typename Problem::GridVariables;
    using FluxVariables = typename Problem::FluxVariables;
	using ElementFluxVariablesCache = typename Problem::ElementFluxVariablesCache;
	using PorousMediumFlowVelocity = Dumux::PorousMediumFlowVelocity<GridVariables, FluxVariables>;

    using GridData = Dumux::GridData<Grid>;
    using GridView = typename Grid::Traits::LeafGridView;

    std::shared_ptr<Grid> grid;
    std::shared_ptr<GridData> gridData;
    std::shared_ptr<FVGridGeometry> gridGeometry;
    std::shared_ptr<Problem> problem;
    std::shared_ptr<GridVariables> gridVariables;
    std::shared_ptr<PorousMediumFlowVelocity> flowVelocities;

    std::shared_ptr<Dune::GlobalIndexSet<GridView>> pointIdx; // global index mappers
    std::shared_ptr<Dune::GlobalIndexSet<GridView>> faceIdx; // global index mappers
    std::shared_ptr<Dune::GlobalIndexSet<GridView>> cellIdx; // global index mappers
    std::map<int, int> localCellIdx; // global to local index mapper
    std::vector<int> globalPointIdx; // local to global index mapper
    std::vector<int> face2CellIdx;
    std::vector<std::vector<double>> scvfBoundaryFluxes;
    std::vector<std::vector<double>> scvfInnerFluxes;
    std::vector<std::vector<double>> scvSources;

    SolutionVector x;
	SolutionVector xBackUp;
    SolutionVector xBackUpManual;
    double simTimeBackUp;
    double simTimeBackUpManual;



};

/**
 * pybind11
 */
template<class Problem, class Assembler, class LinearSolver, int dim = 3>
void init_solverbase(py::module &m, std::string name) {
    using Solver = SolverBase<Problem, Assembler, LinearSolver, dim>; // choose your destiny
    py::class_<Solver>(m, name.c_str())
            .def(py::init<>()) // initialization
			.def("numFluidComp",  &Solver::numFluidComp)
			.def("numComp",  &Solver::numComp)
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
			.def("getFaceSurfaces", &Solver::getFaceSurfaces)
			.def("getCylFaceCoordinates", &Solver::getCylFaceCoordinates)
            .def("getCellSurfacesCyl", &Solver::getCellSurfacesCyl)
            .def("getDofCoordinates", &Solver::getDofCoordinates)
            .def("getPointIndices", &Solver::getPointIndices)
            .def("getCellIndices", &Solver::getCellIndices)
			.def("getGlobal2localCellIdx", &Solver::getGlobal2localCellIdx) // localCellIdx
            .def("getLocal2globalPointIdx", &Solver::getLocal2globalPointIdx) // globalPointIdx
			.def("getFace2CellIdx", &Solver::getFace2CellIdx)
            .def("getDofIndices", &Solver::getDofIndices)
            .def("getOldSolution", &Solver::getOldSolution, py::arg("eqIdx") = 0)
            .def("getSolution", &Solver::getSolution, py::arg("eqIdx") = 0)
            .def("getSolutionAt", &Solver::getSolutionAt, py::arg("gIdx"), py::arg("eqIdx") = 0)
			.def("setSolution", &Solver::setSolution)
            .def("getNeumann", &Solver::getNeumann, py::arg("gIdx"), py::arg("eqIdx") = 0)
            .def("getAllNeumann", &Solver::getAllNeumann, py::arg("eqIdx") = 0)
            .def("getVelocities", &Solver::getVelocities, py::arg("eqIdx") = 0)
            .def("getScvfBoundaryFluxes", &Solver::getScvfBoundaryFluxes)
            .def("getScvfInnerFluxes", &Solver::getScvfInnerFluxes)
            .def("getScvSources", &Solver::getScvSources)
            .def("pickCell", &Solver::pickCell)
            .def("pick", &Solver::pick)
			.def("useMoles",&Solver::useMoles)
            // members
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
