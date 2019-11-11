#ifndef SOLVERBASE_H_
#define SOLVERBASE_H_

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/pybindxi/numpy.h>
namespace py = pybind11;

#include <dumux/io/grid/gridmanager.hh> // grid types
#include <dune/grid/common/mcmgmapper.hh> // the element and vertex mappers
#include <dumux/common/defaultmappertraits.hh> // nothingness

// initialize
#include <dune/common/parallel/mpihelper.hh> // in dune parallelization is realized with MPI
#include <dumux/common/dumuxmessage.hh> // for fun (a static class)
#include <dumux/common/parameters.hh> // global parameter tree with defaults and parsed from args and .input file

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

    int dim = Problem::dimWorld;
    bool isBox = Problem::isBox;

    SolverBase() { };

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

    //    void createGrid() { }

    //	virtual void createGrid(VectorType boundsMin, VectorType boundsMax, VectorType numberOfCells, std::string periodic = "false false false");
    //    virtual void createGrid(std::string file);
    //
    //    virtual void initializeProblem();
    //
    //    virtual std::vector<VectorType> getPoints(); // dune grid tutorial may provide these two
    //    virtual std::vector<VectorType> getCells();
    //
    ////    virtual void initialConditions(); // TODO
    ////    virtual void boundaryConditions(); // TODO
    ////    virtual void sourceTerm(); // TODO
    //
    //    virtual void simulate(double dt); //
    //
    //    int pickIndex(VectorType pos); // vertex or element index
    //    double solutionAt(VectorType pos);
    //
    //    double simTime = 0;
    //    double ddt = -1; // suggested time step

    static const bool numberOfEquations = 1;
    std::array<std::vector<double>,numberOfEquations> initialValues;
    std::array<std::vector<double>, numberOfEquations> solution;

protected:

    //    using Grid = typename Problem::Grid;
    //    using GridData = Dumux::GridData<Grid>;

    //    using FVGridGeometry = typename Problem::FVGridGeometry;
    //    using SolutionVector = typename Problem::SolutionVector;
    //    using GridVariables = typename Problem::GridVariables;
    //
    //    std::shared_ptr<Grid> grid;
    //    std::shared_ptr<GridData> gridData;
    //    std::shared_ptr<FVGridGeometry> gridGeometry;
    //    std::shared_ptr<Problem> problem;
    //    std::shared_ptr<GridVariables> gridVariables;
    //
    //    SolutionVector x;

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
