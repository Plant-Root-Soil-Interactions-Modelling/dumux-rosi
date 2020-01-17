#ifndef SOLVERBASE_H_
#define SOLVERBASE_H_

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/pybindxi/numpy.h>
namespace py = pybind11;

// initialize
#include <dune/common/parallel/mpihelper.hh> // in dune parallelization is realized with MPI
#include <dumux/common/dumuxmessage.hh> // for fun (a static class)
#include <dumux/common/parameters.hh> // global parameter tree with defaults and parsed from args and .input file

// createGrid
#include <dumux/io/grid/gridmanager.hh>

// simulate
#include <dumux/common/timeloop.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/porousmediumflow/richards/newtonsolver.hh>

// getDofIndices, getPointIndices, getCellIndices
#include <dune/grid/utility/globalindexset.hh>

// pick
#include <dumux/common/geometry/intersectingentities.hh>

#include <ostream>
#include <iostream>

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
 * Examples are given in python/ directory
 */
template<class Problem, class Assembler, class LinearSolver>
class SolverBase {
public:

    std::string gridType = "SPGrid"; // <- for better description and warnings, e.g. YaspGrid, AluGrid, FoamGrid, SPGrid
    static const int dim = 3; // Problem::dimWorld;
    bool isBox = Problem::isBox; // numerical method
    const std::vector<std::string> primNames = { "Matric potential [Pa]" };
    static const bool numberOfEquations = 1;

    double simTime = 0;
    double ddt = -1; // internal time step, minus indicates that its uninitialized
    int maxRank = -1; // max mpi rank
    int rank = -1; // mpi rank

    std::vector<std::array<double, numberOfEquations>> solution;

    SolverBase()
    { }

    virtual ~SolverBase()
    { }

    /**
     * Writes the Dumux welcome message, and creates the global Dumux parameter tree from defaults and the .input file
     *
     * Normally you state an input file, that contains all parameters that are needed for the simulation.
     * SolverBase will optionally set most of them dynamically.
     */
    virtual void initialize(std::vector<std::string> args)
    {
        std::vector<char*> cargs;
        cargs.reserve(args.size());
        for(size_t i = 0; i < args.size(); ++i) {
            cargs.push_back(const_cast<char*>(args[i].c_str()));
        } // its a beautiful language
        int argc = cargs.size();
        char** argv  =  &cargs[0];

        if (isBox) { // add DuMux peculiarities
            setParameter("Grid.Overlap","0");
        } else {
            setParameter("Grid.Overlap","1");
        }

        auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

        maxRank = mpiHelper.size();
        rank = mpiHelper.rank();

        if (rank == 0) { // rank is the process number
            Dumux::DumuxMessage::print(/*firstCall=*/true); // print dumux start message
            std::cout << "\n" << toString() << "\n\n" << std::flush; // add my story
        } else {
            if (gridType.compare("FoamGrid") == 0) {
                throw std::invalid_argument("SolverBase::initialize: FoamGrid does not support parallel computation");
            }
        }
        mpiHelper.getCollectiveCommunication().barrier(); // no one is allowed to mess up the message

        Dumux::Parameters::init(argc, argv); // parse command line arguments and input file
    }

    /**
     * Creates the Grid and gridGeometry from the (global) DuMux parameter tree.
     *
     * Parameters known to me are:
     * Grid.UpperRight
     * Grid.LowerLeft
     * Grid.Cells
     * Grid.Periodic
     * Grid.File
     * Grid.Overlap (should = 0 for box, = 1 for CCTpfa), automatically set in SolverBase::initialize
     */
    virtual void createGrid(std::string modelParamGroup = "")
    {
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
        VectorType numberOfCells, std::string periodic = "false false false")
    {
        auto& p = Dumux::Parameters::paramTree(); // had to modify parameters.hh, its private an no way I can pull it out
        std::ostringstream bmin;
        std::ostringstream bmax;
        std::ostringstream cells;
        bmin << boundsMin[0] << " " << boundsMin[1]<< " " << boundsMin[2];
        bmax << boundsMax[0] << " " << boundsMax[1]<< " " << boundsMax[2];
        cells << numberOfCells[0] << " " << numberOfCells[1]<< " " << numberOfCells[2];
        for (int i=0; i<3; i++) {
            if (boundsMin[i] >= boundsMax[i]) {
                throw std::invalid_argument("SolverBase::createGrid: bounds min >= bounds max");
            }
        }
        p["Grid.LowerLeft"] = bmin.str();
        p["Grid.UpperRight"] = bmax.str();
        p["Grid.Cells"] = cells.str();
        p["Grid.Periodic"] = periodic;
        createGrid();
    }

    /**
     * Creates a grid from a file
     *
     * depending on the Grid you choose at compile time it will accept the file type, or not.
     */
    virtual void readGrid(std::string file)
    {
        auto& p = Dumux::Parameters::paramTree();
        p["Grid.File"] = file;
        createGrid();
    }

    /**
     * Returns a rectangular bounding box around the grid geometry
     *
     * [minx, miny, minz, maxx, maxy, maxz]
     */
    virtual std::array<double, 6> getGridBounds()
    {
        auto bMax = gridGeometry->bBoxMax();
        auto bMin = gridGeometry->bBoxMin();
        return std::array<double,6>({bMin[0], bMin[1], bMin[2], bMax[0], bMax[1], bMax[2]});
    }

    /**
     * Writes a parameter into the global Dumux parameter map
     */
    virtual void setParameter(std::string key, std::string value)
    {
        auto& p = Dumux::Parameters::paramTree();
        p[key] = value;
    }

    /**
     * Reads a parameter from the global Dumux parameter map,
     * returns an empty string if value is not set.
     */
    virtual std::string getParameter(std::string key)
    {
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
    virtual void initializeProblem()
    {
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
    }

    /**
     * Simulates the problem for time span dt, with initial time step ddt.
     *
     * Assembler needs a TimeLoop, so i have to create it in each simulate call.
     * (could be improved, but overhead is likely to be small)
     *
     * todo steady state
     */
    virtual void simulate(double dt, double maxDt = -1)
    {
        checkInitialized();
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

        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop); // dynamic
        auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
        using NonLinearSolver = RichardsNewtonSolver<Assembler, LinearSolver>;
        auto nonLinearSolver = std::make_shared<NonLinearSolver>(assembler, linearSolver);
        nonLinearSolver->setVerbose(false);

        timeLoop->start();
        auto xOld = x;
        do {
            ddt = nonLinearSolver->suggestTimeStepSize(timeLoop->timeStepSize());
            timeLoop->setTimeStepSize(ddt); // set new dt as suggested by the newton solver

            assembler->setPreviousSolution(xOld); // set previous solution for storage evaluations

            nonLinearSolver->solve(x, *timeLoop); // solve the non-linear system with time step control

            xOld = x; // make the new solution the old solution

            gridVariables->advanceTimeStep();
            timeLoop->advanceTimeStep(); // advance to the time loop to the next step
            timeLoop->reportTimeStep(); // report statistics of this time step

            problem->setTime(simTime + timeLoop->time()); // pass current time to the problem ddt?

        } while (!timeLoop->finished());

        simTime += dt;
    }

    /**
     * Returns the Dune vertices (vtk points) of the grid for a single mpi process.
     * Gathering and mapping is done in Python.
     */
    virtual std::vector<VectorType> getPoints()
    {
        checkInitialized();
        std::vector<VectorType> points;
        points.reserve(gridGeometry->gridView().size(dim));
        for (const auto& v : vertices(gridGeometry->gridView())) {
            auto p = v.geometry().center();
            points.push_back(make3d(VectorType({p[0], p[1], p[2]})));
        }
        return points;
    }

    /**
     * return the Dune element (vtk cell) centers of the grid for a single mpi process.
     * Gathering and mapping is done in Python.
     */
    virtual std::vector<VectorType> getCellCenters()
    {
        checkInitialized();
        std::vector<VectorType> cells;
        cells.reserve(gridGeometry->gridView().size(0));
        for (const auto& e : elements(gridGeometry->gridView())) {
            auto p = e.geometry().center();
            cells.push_back(make3d(VectorType({p[0], p[1], p[2]})));
        }
        return cells;
    }

    /**
     * Return the Dune element (vtk cell) of the grid as vertex indices.
     * The number of indices are
     * 2 for line, 4 for tetraeder, 8 for hexaedron (todo other relevant objects? add 2D?)
     *
     * This is done for a single process, gathering and mapping is done in Python. TODO
     */
    virtual std::vector<std::vector<int>> getCells()
    {
        checkInitialized();
        std::vector<std::vector<int>> cells;
        cells.reserve(gridGeometry->gridView().size(0));
        for (const auto& e : elements(gridGeometry->gridView())) {
            std::vector<int> cell;
//            gridGeometry->vertexMapper.
//            auto i0 = vMapper.subIndex(e, 0, 1);
//            auto i1 = vMapper.subIndex(e, 1, 1);

        }
        return cells;
    }

    /**
     * Returns the coordinate, where the DOF sit, in the same order like the solution values.
     *
     * DOF sit either at the vertices (points) for box method or
     * element centers (cell centers) for CCTpfa.
     *
     * For a single mpi process. Gathering and mapping is done in Python
     */
    virtual std::vector<VectorType> getDofCoordinates()
    {
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
    virtual std::vector<int> getDofIndices()
    {
        if (isBox) {
            return getPointIndices();
        } else {
            return getCellIndices();
        }
    }

    /**
     * Returns the current solution for a single mpi process.
     * Gathering and mapping is done in Python
     */
    virtual std::vector<std::array<double, numberOfEquations>> getSolution()
    {
        checkInitialized();
        std::vector<std::array<double, numberOfEquations>> sol;
        int n;
        if (isBox) {
            n = gridGeometry->gridView().size(dim);
        } else {
            n = gridGeometry->gridView().size(0);
        }
        // std::cout << "getSolution(): n " << n << ", " << x.size() << "\n" << std::flush;
        sol.resize(n);
        for (int c = 0; c<n; c++) {
            for (int j=0; j<numberOfEquations; j++) {
                sol[c][j] = x[c]; // [j]
            }
        }
        return sol;
    }

    /**
     * picks and element index (cell index)
     * The lucky rank who found it, maps the local index to a global one,
     * and broadcasts to the others
     */
    virtual int pickCell(VectorType pos) // todo! I have to take care about periodicity myself!
    {
        checkInitialized();
        auto& bBoxTree = gridGeometry->boundingBoxTree();
        Dune::FieldVector<double, 3> p({pos[0], pos[1], pos[2]});
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
     * Quick overview
     */
    virtual std::string toString()
    {
        std::ostringstream msg;
        msg << "DuMux Solver using " << gridType << " in " << dim << "D ";
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
            msg << " process rank "<< rank+1 << "/" << maxRank;
        }
        if (simTime>0) {
            msg << "\nSimulation time is "<< simTime/3600/24 << " days, current internal time step is "<< ddt/3600/24 << " days";
        }
        return msg.str();
    }

    /**
     * Checks if the problem was initialized,
     * i.e. initializeProblem() was called
     */
    virtual void checkInitialized()
    {
        if (!gridGeometry) {
            throw std::invalid_argument("SolverBase::checkInitialized: Problem not initialized, call initializeProblem first");
        }
    }

protected:

    VectorType make3d(VectorType p) { ///<- binding is always 3d, lower dimensional Dumux grids are projected to 3d
        switch(dim) {
        case 1: return VectorType({0., 0., p[0]}); // 1D is depth
        case 2: return VectorType({p[0], p[1], 0.}); // 2D, ignore z
        case 3: return VectorType({p[0], p[1], p[2]}); // 3D
        default:
            throw std::invalid_argument("SolverBase::getPoints: Dimension "+ std::to_string(dim) + "D is not supported");
        }
    }

    using Grid = typename Problem::Grid;
    using FVGridGeometry = typename Problem::FVGridGeometry;
    using SolutionVector = typename Problem::SolutionVector;
    using GridVariables = typename Problem::GridVariables;

    using GridData = Dumux::GridData<Grid>;
    using GridView = typename Grid::Traits::LeafGridView;

    std::shared_ptr<Grid> grid;
    std::shared_ptr<GridData> gridData;
    std::shared_ptr<FVGridGeometry> gridGeometry;
    std::shared_ptr<Problem> problem;
    std::shared_ptr<GridVariables> gridVariables;

    std::shared_ptr<Dune::GlobalIndexSet<GridView>> pointIdx;
    std::shared_ptr<Dune::GlobalIndexSet<GridView>> cellIdx;

    SolutionVector x;

};

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
