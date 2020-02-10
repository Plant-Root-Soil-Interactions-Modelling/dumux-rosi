#ifndef PYSOLVER_H_
#define PYSOLVER_H_

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/pybindxi/numpy.h>
namespace py = pybind11;

template<class Problem, class Assembler, class LinearSolver>
void init_solverbase(py::module &);
// void init_richardssp(py::module &);
//void init_rosi(py::module &);

#include "richardssp.hh" // includes solverbase
// #include "RootSoilInteraction.hh"

PYBIND11_MODULE(dumux_rosi, m) {
	init_solverbase<Problem, Assembler, LinearSolver>(m); // types are defined in richardssp
	//init_richardssp(m);
	//init_rosi(m);
}

#endif

