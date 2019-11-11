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
#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/istl/io.hh> // debug vector/matrix output
#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dune/grid/common/mcmgmapper.hh>

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
#include <dumux/assembly/fvassembler.hh>
#include <dumux/io/grid/gridmanager.hh>

#include "solverbase.hh"

#include "richardsproblem.hh" // the problem class. Defines some TypeTag types and includes its spatialparams.hh class
#include "propertiesYasp.hh" // the property system related stuff (to pass types, used instead of polymorphism)
#include "properties_nocoupling.hh" // dummy types for replacing the coupling types

#include <dumux/linear/amgbackend.hh>
#include <dumux/assembly/fvassembler.hh>

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
 * Python binding of the Dumux solver base class
 */
PYBIND11_MODULE(richards_yasp_solver, m) {

    using Solver = SolverBase<Problem, Assembler, LinearSolver>;

    py::class_<Solver>(m, "RichardsYaspSolver")
        .def(py::init<>())
        .def("initialize", &Solver::initialize)
        .def("createGrid", (void (Solver::*)()) &Solver::createGrid) // (void (Solver::*)())
     	.def("createGrid", (void (Solver::*)(VectorType, VectorType, VectorType, std::string)) &Solver::createGrid)
        .def("createGrid", (void (Solver::*)(std::string)) &Solver::createGrid)
        .def("setParameter", &Solver::setParameter)
        .def("initializeProblem", &Solver::initializeProblem)
    	.def("getPoints", &Solver::getPoints) // vtk naming
    	.def("getCells", &Solver::getCells) // vtk naming
    	.def("getDof", &Solver::getDof) // vtk naming
		.def("simulate", &Solver::simulate)
        .def_readwrite("initialValues", &Solver::initialValues)
        .def_readwrite("solution", &Solver::solution)
    	.def_readwrite("ddt", &Solver::ddt)
    	.def("__str__",&Solver::toString);
}

#endif

