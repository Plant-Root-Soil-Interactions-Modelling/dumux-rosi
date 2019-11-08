#ifndef SOLVERBASE_H_
#define SOLVERBASE_H_

#include "pybind11/include/pybind11/pybind11.h"
#include "pybind11/include/pybind11/stl.h"
#include "pybind11/include/pybind11/numpy.h"
namespace py = pybind11;

#include <dumux/io/grid/gridmanager.hh> // grid types

using VectorType = py::array_t<double>; // numpy point and vector type for the binding

using Grid = Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>;
/**
 * lacking polymorphism I have not found a way to make this dynamic,
 * you need to choose the grid at compile time,
 * and I don't know how to pass it via CMakeLists.txt building the Python binding
 */
using GridManager = Dumux::GridManager<Grid>;
using GridData = Dumux::GridData<Grid>;;

/**
 * Derived class will pass ownership
 */

class GridManagerFix :public GridManager {
public:
	using GridManager::GridManager;
	std::shared_ptr<Grid> gridPtr() {
		return gridPtr_;
	}
};


/**
 * Dumux as a solver with a simple Python interface
 */
class SolverBase {

public:

	virtual ~SolverBase() { };

    virtual void initialize(std::vector<std::string> args);

	virtual void createGrid(VectorType boundsMin, VectorType boundsMax, VectorType numberOfCells);
    virtual void createGrid(std::string file);

    virtual std::vector<VectorType> getPoints(); // dune grid tutorial may provide these two
    virtual std::vector<VectorType> getCells();

//    virtual void initialConditions();
//    virtual void boundaryConditions();

    virtual void simulate(); //

    int pickIndex(VectorType pos); // vertex or element index
    double solutionAt(VectorType pos);

    double simTime = 0;
    std::vector<double> initialValues;
    std::vector<double> solution;

private:

    std::shared_ptr<Grid> grid; // i would like to have a grid
    std::shared_ptr<GridData> gridData;

};

/**
 * Python binding of the (Dumux) solver base class
 */
PYBIND11_MODULE(solverbase, m) {
    py::class_<SolverBase>(m, "SolverBase")
    	.def("initialize", &SolverBase::initialize)
        .def("createGrid", (void (SolverBase::*)(VectorType, VectorType, VectorType)) &SolverBase::createGrid)
        .def("createGrid", (void (SolverBase::*)(std::string)) &SolverBase::createGrid)
    	.def("getPoints", &SolverBase::getPoints) // vtk naming
    	.def("getCells", &SolverBase::getCells) // vtk naming
        .def("simulate", &SolverBase::simulate);
}

#endif

