#ifndef PYSOLVER_H_
#define PYSOLVER_H_

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/pybindxi/numpy.h>
#include <dune/pybindxi/functional.h>
namespace py = pybind11;

#include "richardssp.hh" // includes solverbase # TODO seperate richards from problem, assembler, linearsolver defs

PYBIND11_MODULE(dumux_rosi, m) {
	init_solverbase<Problem, Assembler, LinearSolver>(m, "SolverBase"); // types are defined in richardssp
	init_richardssp<Problem, Assembler, LinearSolver>(m, "RichardsSP");
}

#endif

