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

    std::string gridType = "YaspGrid"; // <- for better description and warnings, e.g. YaspGrid, AluGrid, FoamGrid, SPGrid
    int dim = Problem::dimWorld;
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
     * SolverBase will optionally set some of them dynamically.
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

        // add DuMux peculiarities
        if (isBox) {
            setParameter("Grid.Overlap","0");
        } else {
            setParameter("Grid.Overlap","1");
        }

        auto& mpiHelper = Dune::MPIHelper::instance(argc, argv); // of type MPIHelper, or FakeMPIHelper (in mpihelper.hh)

        maxRank = mpiHelper.size();
        rank = mpiHelper.rank();

        if (rank == 0) { // rank is the process number
            Dumux::DumuxMessage::print(/*firstCall=*/true); // print dumux start message
            std::cout << "\n" << toString() << "\n\n" << std::flush; // add my story
        } else {
            if (gridType.compare("FoamGrid") == 0) {
                throw std::invalid_argument("SolverBase::initialize: FoamGrid does not support parallel computation");
            }
            // std::cout << "SolverBase::initialize: MPI working\n"; // for debugging
        }

        Dumux::Parameters::init(argc, argv); // parse command line arguments and input file
    }

    /**
     * Creates a grid from the (global) Dumux parameter tree
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
     * Returns a rectangular bounding box around the geometry
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
     * creates (a) the GridGeometry, (b) the Problem, (c) applies initial conditions (using input file functions)
     * (d) GridVariables (e) resets simtime, and internal time step ddt
     */
    virtual void initializeProblem()
    {
        gridGeometry = std::make_shared<FVGridGeometry>(grid->leafGridView());
        gridGeometry->update();
        problem = std::make_shared<Problem>(gridGeometry);
        int dof = gridGeometry->numDofs();
        x = SolutionVector(dof);

        std::cout << "SolutionVector of rank "<< rank << " has size "  << dof << "\n";

        problem->applyInitialSolution(x); // Dumux way of saying x = problem->applyInitialSolution()

        gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        gridVariables->init(x); // initialize all variables , updates volume variables to the current solution, and updates the flux variable cache
        simTime = 0; // reset
        ddt = -1;

        copySolution(); // to have access to the initial values
    }

    /**
     * Returns the Dune vertices (vtk points) of the grid
     * for a single mpi process. Gathering is done in Python.
     */
    virtual std::vector<VectorType> getPoints()
    {
        checkInitialized();
        std::vector<VectorType> points;
        points.resize(gridGeometry->gridView().size(dim));
        for (const auto& v : vertices(gridGeometry->gridView())) {
            auto vIdx = gridGeometry->vertexMapper().index(v);
            auto p = v.geometry().center();
            points[vIdx] = make3d(VectorType({p[0], p[1], p[2]}));
        }
        return points;
    }

    /**
     * return the Dune element (vtk cell) centers of the grid
     * for a single mpi process. Gathering is done in Python.
     */
    virtual std::vector<VectorType> getCellCenters()
    {
        checkInitialized();
        std::vector<VectorType> cells;
        cells.resize(gridGeometry->gridView().size(0));
        for (const auto& e : elements(gridGeometry->gridView())) {
            auto eIdx = gridGeometry->elementMapper().index(e);
            auto p = e.geometry().center();
            cells[eIdx] = make3d(VectorType({p[0], p[1], p[2]}));
        }
        return cells;
    }

    /**
     * Returns the coordinate, where the DOF sit, in the same order like the solution values.
     *
     * DOF sit either at the vertices (points) for box method or
     * element centers (cell centers) for CCTpfa.
     *
     * For a single mpi process. Gathering is done in Python
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
     * Return the indices of the grid elements or vertices where the DOF sit.
     * Used to sort the coordinates when gathered from the processes,
     * makes no sense to call directly.
     *
     * For a single mpi process. Gathering is done in Python
     */
    virtual std::vector<int> getDofIndices()
    {
        checkInitialized();
        std::vector<int> indices;
        if (isBox) {
            indices.resize(gridGeometry->gridView().size(dim));
            for (const auto& v : vertices(gridGeometry->gridView())) {
                indices.push_back(gridGeometry->vertexMapper().index(v));
            }
        } else {
            indices.reserve(gridGeometry->gridView().size(0));
            for (const auto& e : elements(gridGeometry->gridView())) {
                indices.push_back(gridGeometry->elementMapper().index(e));
            }
        }
        return indices;
    }

    /**
     * Returns the current solution
     * for a single mpi process. Gathering is done in Python
     */
    virtual std::vector<int> getSolution()
    {
        checkInitialized();
        std::vector<int> indices;
        if (isBox) {
            indices.resize(gridGeometry->gridView().size(dim));
            for (const auto& v : vertices(gridGeometry->gridView())) {
                indices.push_back(gridGeometry->vertexMapper().index(v));
            }
        } else {
            indices.reserve(gridGeometry->gridView().size(0));
            for (const auto& e : elements(gridGeometry->gridView())) {
                indices.push_back(gridGeometry->elementMapper().index(e));
            }
        }
        return indices;
    }

    /**
     * simulates the problem for time span dt, with initial time step ddt.
     *
     * Assembler needs a TimeLoop, so i have to create it in each simulate call.
     * (could be improved, but overhead is likely to be small)
     *
     * todo steady state
     */
    virtual void simulate(double dt, double maxDt = -1) {
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
     * picks and element index (cell index)
     */
    virtual int pickCell(VectorType pos) // todo? do I have to take care about periodicity?
    {
        checkInitialized();
        auto& bBoxTree = gridGeometry->boundingBoxTree();
        Dune::FieldVector<double, 3> p({pos[0], pos[1], pos[2]});
        auto entities = Dumux::intersectingEntities(p, bBoxTree);
        if (entities.empty()) {
            return -1;
        }
        auto element = bBoxTree.entitySet().entity(entities[0]);
        return gridGeometry->elementMapper().index(element);
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
        case 1: return VectorType({0., 0., p[0]}); // 1D
        case 2: return VectorType({p[0], p[1], 0.}); // 2D, ignore z
        case 3: return VectorType({p[0], p[1], p[2]}); // 3D
        default:
            throw std::invalid_argument("SolverBase::getPoints: Dimension "+ std::to_string(dim) + "D is not supported");
        }
    }

//    void copySolution() { ///< copies the solution
//
//        int n = gridGeometry->numDofs();
//        n = gridGeometry->gridView().comm().sum(n);
//        std::cout << "The processor wide DOF are less then: " << n << " \n" << std::flush;
//
//        std::vector<std::array<double, numberOfEquations>> s;
//        s.resize(n);
//        std::array<double, numberOfEquations> ini;
//        std::fill(ini.begin(), ini.end(), -1.e16);
//        std::fill(s.begin(), s.end(), ini);
//
//        if (isBox) { // DOF are located at the vertices
//            int c = 0;
//            for (const auto& v :vertices(gridGeometry->gridView())) {
//                auto vIdx = gridGeometry->vertexMapper().index(v);
//                for (int j=0; j<numberOfEquations; j++) {
//                    s[vIdx][j] = x[j][c]; // rhs, local index ok ?
//                }
//                c++; // increase local index (?)
//            }
//        } else {
//            int c = 0;
//            for (const auto& e :elements(gridGeometry->gridView())) {
//                auto eIdx = gridGeometry->vertexMapper().index(e);
//                for (int j=0; j<numberOfEquations; j++) {
//                    s[eIdx][j] = x[j][c]; // rhs, local index ok ?
//                }
//                c++; // increase local index (?)
//            }
//        }
//        /* collect from mpi (a desperate approach...)*/
//        for (int i = 0; i< n; i++) {
//            // if (i%100==0) std::cout << "i" << i << "\n";
//            for (int j=0; j<numberOfEquations; j++) {
//                s[i][j] = gridGeometry->gridView().comm().max(s[i][j]);
//            }
//        }
//    }

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

//using SolutionVector =  Dune::BlockVector<GetPropType<TypeTag, Properties::PrimaryVariables>>; // in fvproperties


