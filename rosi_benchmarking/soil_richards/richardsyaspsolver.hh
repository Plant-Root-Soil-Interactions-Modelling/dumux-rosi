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


/**
 * Python binding of the Dumux solver base class
 */
PYBIND11_MODULE(richards_yasp_solver, m) {

    using RichardsYaspSolver = SolverBase<Problem, Assembler, LinearSolver>;

    py::class_<RichardsYaspSolver>(m, "RichardsYaspSolver")
        .def(py::init<>())
   	   // .def("initialize", &RichardsYaspSolver::initialize)
//        .def("createGrid", &RichardsYaspSolver::createGrid) // (void (RichardsYaspSolver::*)())
//    	.def("createGrid", (void (RichardsYaspSolver::*)(VectorType, VectorType, VectorType, std::string)) &RichardsYaspSolver::createGrid)
//        .def("createGrid", (void (RichardsYaspSolver::*)(std::string)) &RichardsYaspSolver::createGrid)
//    	.def("getPoints", &RichardsYaspSolver::getPoints) // vtk naming
//    	.def("getCells", &RichardsYaspSolver::getCells) // vtk naming
//        .def("simulate", &RichardsYaspSolver::simulate)
        .def_readwrite("initialValues", &RichardsYaspSolver::initialValues)
        .def_readwrite("solution", &RichardsYaspSolver::solution);
}

#endif

