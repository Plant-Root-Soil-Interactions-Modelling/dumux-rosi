#ifndef PYTHON_RICHARDS_CYL_SOLVER_H_
#define PYTHON_RICHARDS_CYL_SOLVER_H_

#include "richards.hh" // most includes are in solverbase

// writeDumuxVTK
#include <dumux/io/vtkoutputmodule.hh>


/**
 * Adds solver functionality, that specifically makes sense for Richards equation
 */
template<class Problem, class Assembler, class LinearSolver, int dim = 1>
class RichardsCyl : public Richards<Problem, Assembler, LinearSolver, dim> {
public:

    std::vector<double> ls; // local sink [kg/s]

    virtual ~RichardsCyl() { }

    /**
     * Calls parent, additionally turns file output off
     */
    void initialize(std::vector<std::string> args) override {
    	Richards<Problem, Assembler, LinearSolver, dim>::initialize(args);
        this->setParameter("Problem.EnableGravity", "false"); // important in 1d axial-symmetric problem
        this->setParameter("Soil.Problem.EnableGravity", "false"); // important in 1d axial-symmetric problem
    }

    // TODO getWaterVolume needs adjusting

};

/**
 * pybind11
 */
template<class Problem, class Assembler, class LinearSolver, int dim = 3>
void init_richards_cyl(py::module &m, std::string name) {
    using RichardsFoam = RichardsCyl<Problem, Assembler, LinearSolver>;
	py::class_<RichardsFoam, SolverBase<Problem, Assembler, LinearSolver, dim>>(m, name.c_str())
   .def(py::init<>())
   .def("setSource", &RichardsFoam::setSource)
   .def("setCriticalPressure", &RichardsFoam::setCriticalPressure)
   .def("getSolutionHead", &RichardsFoam::getSolutionHead, py::arg("eqIdx") = 0)
   .def("getWaterContent",&RichardsFoam::getWaterContent)
   .def("getWaterVolume",&RichardsFoam::getWaterVolume)
   .def("writeDumuxVTK",&RichardsFoam::writeDumuxVTK);
}


#endif
