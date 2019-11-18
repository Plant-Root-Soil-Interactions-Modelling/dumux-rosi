#ifndef RICHARDS_YASP_SOLVER_H_
#define RICHARDS_YASP_SOLVER_H_

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/pybindxi/numpy.h>
namespace py = pybind11;


#include <config.h>
//
#include <iostream>
#include <memory>
#include <algorithm>

// Dune
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>

// Dumux
#include <dumux/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>
#include <dumux/common/geometry/boundingboxtree.hh>
#include <dumux/common/geometry/geometricentityset.hh>
#include <dumux/common/geometry/intersectingentities.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/linear/amgbackend.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/io/grid/gridmanager.hh>

#include "solverbase.hh"

#include "../soil_richards/richardsproblem.hh" // the problem class. Defines some TypeTag types and includes its spatialparams.hh class
#include "../soil_richards/propertiesYasp.hh" // the property system related stuff (to pass types, used instead of polymorphism)
#include "../soil_richards/properties_nocoupling.hh" // dummy types for replacing the coupling types

// define the type tag for this problem
using TypeTag = Dumux::Properties::TTag::RichardsBox; // RichardsCC, RichardsBox

// the problem
using Problem = Dumux::RichardsProblem<TypeTag>;

// Pick assembler and solver
using Assembler = Dumux::FVAssembler<TypeTag, Dumux::DiffMethod::numeric>;
using LinearSolver = Dumux::AMGBackend<TypeTag>;

// Name of the choice of Problem, Assembler, LinearSolver
std::string name = "RichardsYaspSolver";

/**
 * Adds solver functionality, that specifically makes sense for Richards equation
 */
class RichardsYaspSolver : public SolverBase<Problem, Assembler, LinearSolver> {
public:

    virtual ~RichardsYaspSolver()
    { }

    /**
     * Total water volume in domain
     */
    virtual double getWaterVolume()
    {
        checkInitialized();
//        double cVol = 0.;
//        for (const auto& element : Dune::elements(gridGeometry->gridView())) { // soil elements
//            auto fvGeometry = Dumux::localView(gridGeometry); // soil solution -> volume variable
//            fvGeometry.bindElement(element);
//            auto elemVolVars = Dumux::localView(gridVariables->curGridVolVars());
//            elemVolVars.bindElement(element, fvGeometry, x);
//            for (const auto& scv : Dumux::scvs(fvGeometry)) {
//                cVol += elemVolVars[scv].saturation(0)*scv.volume();
//            }
//
//        }
//        return cVol;
        return 0;
    }

};

using Solver = RichardsYaspSolver;

/**
 * Python binding of the Dumux solver base class
 */
PYBIND11_MODULE(richards_yasp_solver, m) {

    py::class_<Solver>(m, name.c_str())
        .def(py::init<>())
        .def("initialize", &Solver::initialize)
        .def("createGrid", (void (Solver::*)(std::string)) &Solver::createGrid, py::arg("modelParamGroup") = "") // overloads, defaults
     	.def("createGrid", (void (Solver::*)(VectorType, VectorType, VectorType, std::string)) &Solver::createGrid,
     	   py::arg("boundsMin"), py::arg("boundsMax"), py::arg("numberOfCells"), py::arg("periodic") = "false false false") // overloads, defaults
        .def("readGrid", &Solver::readGrid)
        .def("getGridBounds", &Solver::getGridBounds)
        .def("setParameter", &Solver::setParameter)
        .def("getParameter", &Solver::getParameter)
        .def("initializeProblem", &Solver::initializeProblem)
    	.def("getPoints", &Solver::getPoints) // vtk naming
    	.def("getCellCenters", &Solver::getCellCenters) // vtk naming
    	.def("getDofCoordinates", &Solver::getDofCoordinates)
        .def("getDofIndices", &Solver::getDofIndices)
        .def("simulate", &Solver::simulate)
        .def("getSolution", &Solver::getSolution)
        .def("pickCell", &Solver::pickCell)
        .def_readonly("solution", &Solver::solution) // read only
        .def_readonly("simTime", &Solver::simTime) // read only
        .def_readonly("rank", &Solver::rank) // read only
        .def_readonly("maxRank", &Solver::maxRank) // read only
        .def_readwrite("ddt", &Solver::ddt) // initial internal time step
    	.def("__str__",&Solver::toString)
        .def("checkInitialized", &Solver::checkInitialized)
        // added by class specialization
        .def("getWaterVolume",&Solver::getWaterVolume);
}

#endif

