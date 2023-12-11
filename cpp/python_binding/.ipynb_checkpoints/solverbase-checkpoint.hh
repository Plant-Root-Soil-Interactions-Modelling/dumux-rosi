#ifndef PYTHON_SOLVER_BASE_H_
#define PYTHON_SOLVER_BASE_H_

// initialize
#include <dune/common/parallel/mpihelper.hh> // in dune parallelization is realized with MPI
#include <dumux/common/dumuxmessage.hh> // for fun (a static class)
#include <dumux/common/parameters.hh> // global parameter tree with defaults and parsed from args and .input file

// createGrid
#include <dumux/io/grid/gridmanager.hh>
#include <dune/grid/common/gridfactory.hh>

// simulate
#include <dumux/common/timeloop.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/porousmediumflow/richards/newtonsolver.hh>

// getDofIndices, getPointIndices, getCellIndices
#include <dune/grid/utility/globalindexset.hh>

// pick
#include <dumux/common/geometry/intersectingentities.hh>
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
    bool isBox = Problem::isBox; // numerical method
    int dimWorld = dim;
    double simTime = 0;
    double ddt = -1; // internal time step, minus indicates that its uninitialized
    int maxRank = -1; // max mpi rank
    int rank = -1; // mpi rank

    bool periodic = false; // periodic domain
    std::array<int, dim> numberOfCells;
    int verbose = 0;

    SolverBase() {
        for (int i=0; i<dim; i++) { // initialize numberOfCells
            numberOfCells[i] = 0;
        }
        setParameter("Grid.Overlap","1");
    }

    virtual ~SolverBase() { }

    /**
     * Writes the Dumux welcome message, and creates the global Dumux parameter tree from defaults and the .input file
     *
     * Normally you state an input file, that contains all parameters that are needed for the simulation.
     * SolverBase will optionally set most of them dynamically.
     */
    virtual void initialize(std::vector<std::string> args_ = std::vector<std::string>(0), bool verbose = true,
							bool doMPI = true) {
		std::vector<char*> cargs;
        cargs.reserve(args_.size());
        for(size_t i = 0; i < args_.size(); i++) {
            cargs.push_back(const_cast<char*>(args_[i].c_str()));
        } // its a beautiful language
        int argc = cargs.size();
        char** argv  =  &cargs[0];

        if (isBox) { // add DuMux peculiarities
            setParameter("Grid.Overlap","0");
        } else {
            setParameter("Grid.Overlap","1");
        }
        
        setParameter("Flux.UpwindWeight", "0.5"); // Timo's advice for flows that are not strongly convection dominated, Dumux default = 1
		
		// if(doMPI)
		// {
			// DUNE_EXPORT typedef Dune::MPIHelper MPIHelper;
		// }else{
			// DUNE_EXPORT typedef Dune::MPIHelper FakeMPIHelper; 
		// }

			//auto& mpiHelper = Dune::FakeMPIHelper::instance(argc, argv);// Dune::MPIHelper::instance(argc, argv);
			auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);
			maxRank = mpiHelper.size();
			rank = mpiHelper.rank();

			if ((rank == 0) && verbose) { // rank is the process number
				std::cout << "\n" << toString() << "\n" << std::flush; // add my story
				Dumux::DumuxMessage::print(/*firstCall=*/true); // print dumux start message
			}
				//std::cout<<"mpiHelper::isFake "<<Dune::MPIHelper::isFake<<" " //Dune::FakeMPIHelper::isFake<<" "<<Dune::MPIHelper::isFake<<" size "
				//<<Dune::MPIHelper::getCollectiveCommunication().size()<<" doMPI "<<doMPI<<std::endl;
			if(doMPI)
			{
				mpiHelper.getCollectiveCommunication().barrier(); // no one is allowed to mess up the message
			}
		// }else{
			// if (verbose) { // even if we have MPI, run model sequentially.
				// std::cout << "\n" << toString() << "\n" << std::flush; // add my story
				// Dumux::DumuxMessage::print(/*firstCall=*/true); // print dumux start message
			// }
		// }
        setParameter("Problem.Name","noname");
        Dumux::Parameters::init(argc, argv); // parse command line arguments and input file
		
		//@see <dumux/common/parameters.hh> 
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
        gridGeometry->update();
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
        auto& p = Dumux::Parameters::paramTree(); // had to modify parameters.hh, its private an no way I can pull it out
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
                p["Grid.Periodic"] = "true true false";
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
        gridGeometry->update();
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
        auto& p = Dumux::Parameters::paramTree();
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
        auto& p = Dumux::Parameters::paramTree();
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
     *
     * The initialize values are set to the current solution,
     * i.e. can be analyzed using getSolution().
     */
    virtual void initializeProblem() {
        verbose =  Dumux::getParam<int>("Problem.verbose",0);
        int verboseIndexSet =  Dumux::getParam<int>("Problem.verboseIndexSet",0);
        
        problem = std::make_shared<Problem>(gridGeometry);
        int dof = gridGeometry->numDofs();
        x = SolutionVector(dof);

        problem->applyInitialSolution(x); // Dumux way of saying x = problem->applyInitialSolution()

        gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        gridVariables->init(x); // initialize all variables , updates volume variables to the current solution, and updates the flux variable cache

        simTime = 0; // reset
        ddt = -1;
        if(verboseIndexSet){std::cout<<"getpointIdx"<<std::endl;}
        pointIdx = std::make_shared<Dune::GlobalIndexSet<GridView>>(grid->leafGridView(), dim, verboseIndexSet); // global index mappers
        
        if(verboseIndexSet){std::cout<<"getcellIdx"<<std::endl;}
        cellIdx = std::make_shared<Dune::GlobalIndexSet<GridView>>(grid->leafGridView(), 0, verboseIndexSet);
        if(verboseIndexSet){std::cout<<"GOTcellIdx"<<std::endl;}
        
        //std::cout<<"getlocalCellIdx"<<std::endl;
        localCellIdx.clear();
        for (const auto& e : Dune::elements(gridGeometry->gridView())) {
            int eIdx = gridGeometry->elementMapper().index(e);
            int gIdx = cellIdx->index(e);
            localCellIdx[gIdx] = eIdx;
        }
        //std::cout<<"getglobalPointIdx"<<std::endl;
        globalPointIdx.resize(gridGeometry->gridView().size(dim)); // number of vertices
        for (const auto& v : Dune::vertices(gridGeometry->gridView())) {
            int vIdx = gridGeometry->vertexMapper().index(v);
            int gIdx = pointIdx->index(v);
            globalPointIdx[vIdx] = gIdx;
        }
        //std::cout<<"finished filling the maps"<<std::endl;
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
    virtual void solve(double dt, double maxDt = -1, bool solverVerbose = false) {
        checkInitialized();
        using namespace Dumux;

        //problem->verbose =  getParam<int>("Problem.verbose", 0);
        if (ddt<1.e-6) { // happens at the first call
            ddt = getParam<double>("TimeLoop.DtInitial", dt/10); // from params, or guess something
        }

        std::shared_ptr<CheckPointTimeLoop<double>> timeLoop =
            std::make_shared<CheckPointTimeLoop<double>>(/*start time*/0., ddt, /*final time*/ dt, false); // the main time loop is moved to Python
        if (maxDt<0) { // per default value take from parameter tree
            maxDt = getParam<double>("TimeLoop.MaxTimeStepSize", dt); // if none, default is outer time step
        }
        timeLoop->setMaxTimeStepSize(maxDt);

        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop); // dynamic
        auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
        using NonLinearSolver = RichardsNewtonSolver<Assembler, LinearSolver>;
        auto nonLinearSolver = std::make_shared<NonLinearSolver>(assembler, linearSolver);//
        nonLinearSolver->setVerbose(solverVerbose);
        if(verbose){std::cout<<rank<<" timeLoop->start();"<<std::endl;}
        timeLoop->start();
		
		xBackUp = x;
        auto xOld = x;
        double minddt = std::min(1.,dt);//in case we have very small simulation time
        do {
			if(verbose){std::cout<<rank<< " before, nonLinearSolver->suggestTimeStepSize, current time: "<<timeLoop->time()
                            <<" simTime "<<simTime<<std::endl;}
            ddt = nonLinearSolver->suggestTimeStepSize(timeLoop->timeStepSize());
            ddt = std::max(ddt, minddt); // limit minimal and maximal suggestion
            timeLoop->setTimeStepSize(ddt); // set new dt as suggested by the newton solver
            ddt = timeLoop->timeStepSize();//limited ddt to stay below dt
            problem->setTime(simTime + timeLoop->time(), ddt); // pass current time to the problem ddt?
            if(verbose){std::cout<<rank<<" before assembler->setPreviousSolution, ddt: "<<ddt<<std::endl;}
            assembler->setPreviousSolution(xOld); // set previous solution for storage evaluations
            if(verbose){std::cout<<rank<<" nonLinearSolver->solve"<<std::endl;}
            nonLinearSolver->solve(x, *timeLoop); // solve the non-linear system with time step control

            xOld = x; // make the new solution the old solution
            if(verbose){std::cout<<rank<<" gridVariables->advanceTimeStep"<<std::endl;}
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
    void solveNoMPI(double dt, double maxDt = -1, bool solverVerbose = false, bool saveBC = false) {
		
        checkInitialized();
        using namespace Dumux;
        clearSaveBC();
        //problem->verbose =  getParam<int>("Problem.verbose", 0);
        
        if (ddt<1.e-6) { // happens at the first call
            ddt = getParam<double>("TimeLoop.DtInitial", dt/10); // from params, or guess something
        }

        std::shared_ptr<CheckPointTimeLoop<double>> timeLoop =
            std::make_shared<CheckPointTimeLoop<double>>(/*start time*/0., ddt, /*final time*/ dt, false); // the main time loop is moved to Python
        if (maxDt<0) { // per default value take from parameter tree
            maxDt = getParam<double>("TimeLoop.MaxTimeStepSize", dt); // if none, default is outer time step
        }
        timeLoop->setMaxTimeStepSize(maxDt);

        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop); // dynamic
        auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
        using NonLinearSolver = RichardsNewtonSolver<Assembler, LinearSolver,
								 PartialReassembler<Assembler>,
								Dune::CollectiveCommunication<Dune::FakeMPIHelper::MPICommunicator> >;
        auto nonLinearSolver = std::make_shared<NonLinearSolver>(assembler, linearSolver, 
								Dune::FakeMPIHelper::getCollectiveCommunication());//
        nonLinearSolver->setVerbose(solverVerbose);
		//std::cout<<"1cyl, timeLoop->start();"<<std::endl;
        timeLoop->start();
		
		xBackUp = x;
        auto xOld = x;
        double minddt = std::min(1.,dt);//in case we have very small simulation time
        do {
			//std::cout<<"1cyl, nonLinearSolver->suggestTimeStepSize"<<std::endl;
            ddt = nonLinearSolver->suggestTimeStepSize(timeLoop->timeStepSize());
            ddt = std::max(ddt, minddt); // limit minimal and maximal suggestion
            timeLoop->setTimeStepSize(ddt); // set new dt as suggested by the newton solver
            ddt = timeLoop->timeStepSize();//limited ddt to stay below dt
            problem->setTime(simTime + timeLoop->time(), ddt); // pass current time to the problem ddt?
            //std::cout<<"1cyl, assembler->setPreviousSolution"<<std::endl;
            assembler->setPreviousSolution(xOld); // set previous solution for storage evaluations
            //std::cout<<"1cyl, nonLinearSolver->solve"<<std::endl;
            nonLinearSolver->solve(x, *timeLoop); // solve the non-linear system with time step control

            xOld = x; // make the new solution the old solution
            //std::cout<<"1cyl, gridVariables->advanceTimeStep"<<std::endl;
            if(saveBC)
            {
                doSaveBC(timeLoop->time() );
            }
            gridVariables->advanceTimeStep();

            timeLoop->advanceTimeStep(); // advance to the time loop to the next step
            timeLoop->reportTimeStep(); // report statistics of this time step

        } while (!timeLoop->finished());

        simTime += dt;
    }
	
    virtual void clearSaveBC() {}
    virtual void doSaveBC(double currentTime) {}
    
    virtual void reset() {
        checkInitialized();
		x = xBackUp;
	}
    /**
     * Finds the steady state of the problem.
     *
     * Optionally, solve for a time span first, to get a good initial guess.
     */
    virtual void solveSteadyState() {
        checkInitialized();
        using namespace Dumux;

        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables); // steady state
        auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
        using NonLinearSolver = RichardsNewtonSolver<Assembler, LinearSolver>;
        auto nonLinearSolver = std::make_shared<NonLinearSolver>(assembler, linearSolver);
        nonLinearSolver->setVerbose(false);

        assembler->setPreviousSolution(x);
        nonLinearSolver->solve(x); // solve the non-linear system

        simTime = std::numeric_limits<double>::infinity();
    }

    /**
     * Returns the Dune vertices (vtk points) of the grid for a single mpi process.
     * Gathering and mapping is done in Python.
     */
    virtual std::vector<VectorType> getPoints() {
        checkInitialized();
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
        checkInitialized();
        std::vector<VectorType> cells;
        cells.reserve(gridGeometry->gridView().size(0));
		
		//std::cout<<"getCellCenters rank:"<<rank<<"\n";
        for (const auto& e : elements(gridGeometry->gridView())) {
            auto p = e.geometry().center();
			
			//int gIdx = cellIdx->index(e);
			//std::cout<<"	gIdx "<<gIdx;
			
            VectorType vp;
            for (int i=0; i<dim; i++) { // found no better way
                vp[i] = p[i];
			//	std::cout<<" "<<p[i];
            }
			//std::cout<<std::endl;
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
        checkInitialized();
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
		assert(false&&"Do NOT use getCellVolumesCyl, use getCellSurfacesCyl instead");//do I need to use dune-assert here?
        std::vector<double> vols;
        for (const auto& e : elements(gridGeometry->gridView())) {
            vols.push_back(e.geometry().volume()*2*e.geometry().center()[0]*3.1415);
        }
        return vols;
    }

    /**
     * The volume [m3] of each element (vtk cell)
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
     * Returns the current solution for a single mpi process.
     * Gathering and mapping is done in Python
     * solution == mole (or mass) fraction for all elements EXCEPT the liquide phase
     * solution == pressure (absolute?) for the liquide phase
     */
    virtual std::vector<double> getSolution(int eqIdx = 0) {
        int n = checkInitialized();
        std::vector<double> sol;
        sol.resize(n);
        //std::cout<<"std::vector<double> getSolution(int eqIdx = 0) "<<eqIdx<<" ";
        for (int c = 0; c<n; c++) {
            sol[c] = x[c][eqIdx];
            //std::cout<<c<<" ";
        }
        //std::cout<<std::endl;
        return sol;
    }
	
    /**
     * Returns the current solution for a single mpi process.
     * Gathering and mapping is done in Python
     */
    virtual std::vector<double> getOldSolution(int eqIdx = 0) {
        int n = checkInitialized();
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
    virtual void setSolution( std::vector<double> sol, int eqIdx = 0) {
        int n = checkInitialized();
        for (int c = 0; c<n; c++) {
            x[c][eqIdx] = sol[c] ;
        }
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
            for (const auto& scvf : scvfs(fvGeometry)) {
                if (scvf.boundary()) {
                    double n = problem->neumann(e, fvGeometry, elemVolVars, scvf)[eqIdx];  // [ kg / (m2 s)]
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
                    auto elemVolVars = Dumux::localView(gridVariables->curGridVolVars());
                    elemVolVars.bindElement(e, fvGeometry, x);
                    f += problem->neumann(e, fvGeometry, elemVolVars, scvf)[eqIdx]; // [kg / (m2 s)]
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
        int n = checkInitialized();
        std::vector<double> fluxes;
        fluxes.resize(n);

        auto elemVolVars = Dumux::localView(gridVariables->curGridVolVars());
        auto fvGeometry = Dumux::localView(*gridGeometry); // soil solution -> volume variable
        auto elemFluxVarsCache = Dumux::localView(gridVariables->gridFluxVarsCache());

        // the upwind term to be used for the volume flux evaluation
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
    virtual int pickCell(VectorType pos) {
        checkInitialized();
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
		//std::cout<<"solverbase::pickCell, !entities.empty(), "<<gIdx<<" "<<entities.size()<<" "<<entities[0]<<std::endl<<std::flush;
        }
		
		//std::cout<<"solverbase::pickCell, gridGeometry, "<<gIdx<<" "<<gridGeometry->gridView().comm().max(gIdx)<<std::endl<<std::flush;
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
     * Checks if the problem was initialized, and returns number of local dof
     * i.e. initializeProblem() was called
     */
    virtual int checkInitialized() {
        if (!gridGeometry) {
            throw std::invalid_argument("SolverBase::checkInitialized: Problem not initialized, call initializeProblem first");
        }
        if (this->isBox) {
            return this->gridGeometry->gridView().size(dim);
        } else {
            return this->gridGeometry->gridView().size(0);
        }
    }

    virtual std::map<int, int> getGlobal2localCellIdx() {
        return localCellIdx;
    }
	
    virtual std::vector<int> getLocal2globalPointIdx() {
        return globalPointIdx;
    }
protected:

    using Grid = typename Problem::Grid;
    using FVGridGeometry = typename Problem::FVGridGeometry;
    using SolutionVector = typename Problem::SolutionVector;
    using GridVariables = typename Problem::GridVariables;
    using FluxVariables = typename Problem::FluxVariables;

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

};

/**
 * pybind11
 */
template<class Problem, class Assembler, class LinearSolver, int dim = 3>
void init_solverbase(py::module &m, std::string name) {
    using Solver = SolverBase<Problem, Assembler, LinearSolver, dim>; // choose your destiny
    py::class_<Solver>(m, name.c_str())
            .def(py::init<>()) // initialization
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
            .def("initializeProblem", &Solver::initializeProblem)
            .def("setInitialCondition", &Solver::setInitialCondition, py::arg("init"), py::arg("eqIdx") = 0)
			.def("reset", &Solver::reset)
            // simulation 
            .def("solve", &Solver::solve, py::arg("dt"), py::arg("maxDt") = -1, py::arg("solverVerbose") = false)
            .def("solveNoMPI", &Solver::solveNoMPI, py::arg("dt"), py::arg("maxDt") = -1, py::arg("solverVerbose") = false,
                py::arg("saveBC") = false)
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
            .def("setSolution", &Solver::setSolution)
            .def("getSolutionAt", &Solver::getSolutionAt, py::arg("gIdx"), py::arg("eqIdx") = 0)
            .def("getNeumann", &Solver::getNeumann, py::arg("gIdx"), py::arg("eqIdx") = 0)
            .def("getAllNeumann", &Solver::getAllNeumann, py::arg("eqIdx") = 0)
            .def("getNetFlux", &Solver::getNetFlux, py::arg("eqIdx") = 0)
            .def("pickCell", &Solver::pickCell)
            .def("pick", &Solver::pick)
            // members
            .def_readonly("dimWorld", &Solver::dimWorld) // read only
            .def_readonly("simTime", &Solver::simTime) // read only
            .def_readwrite("verbose", &Solver::verbose) // initial internal time step
            .def_readwrite("ddt", &Solver::ddt) // initial internal time step
            .def_readonly("rank", &Solver::rank) // read only
            .def_readonly("maxRank", &Solver::maxRank) // read only
            .def_readonly("numberOfCells", &Solver::numberOfCells) // read only
            .def_readonly("periodic", &Solver::periodic) // read only
            // useful
            .def("__str__",&Solver::toString)
            .def("checkInitialized", &Solver::checkInitialized);
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
