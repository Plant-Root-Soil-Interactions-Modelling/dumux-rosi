#ifndef PYTHON_RICHARDS10_SOLVER_H_
#define PYTHON_RICHARDS10_SOLVER_H_

// most includes are in solverbase
#include "solverbase.hh"

#include "richards.hh" // most includes are in solverbase
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
// #include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

// writeDumuxVTK
#include <dumux/io/vtkoutputmodule.hh>



/**
 * Adds solver functionality, that specifically makes sense for Richards equation
 */
template<class Problem, class Assembler, class LinearSolver, int dim = 3>
class Richards10 : public Richards<Problem, Assembler, LinearSolver, dim> {
public:

    std::vector<double> getCSS1_out() {
    	return this->problem->getCSS1_(); // 
    }
    std::vector<double> getRF_out() {
    	return this->problem->getRF_(); // 
    }

};

/**
 * pybind11
 */
template<class Problem, class Assembler, class LinearSolver>
void init_richards_10(py::module &m, std::string name) {
    using Richards_ = Richards10<Problem, Assembler, LinearSolver>;
	py::class_<Richards_, SolverBase<Problem, Assembler, LinearSolver>>(m, name.c_str())
   .def(py::init<>())
   .def("initialize", &Richards_::initialize, py::arg("args_") = std::vector<std::string>(0), 
													py::arg("verbose") = true,py::arg("doMPI") = true)
   .def("setSource", &Richards_::setSource, py::arg("sourceMap"), py::arg("eqIdx") = 0)
   .def("applySource", &Richards_::applySource)
   .def("setCriticalPressure", &Richards_::setCriticalPressure)
   .def("setInitialConditionHead", &Richards_::setInitialConditionHead)
   .def("getSolutionHead", &Richards_::getSolutionHead, py::arg("eqIdx") = 0)
   .def("getSolutionHeadAt", &Richards_::getSolutionHeadAt, py::arg("gIdx"), py::arg("eqIdx") = 0)
   .def("getWaterContent",&Richards_::getWaterContent)
   .def("getSaturation",&Richards_::getSaturation)
   .def("getWaterVolume",&Richards_::getWaterVolume)
   .def("getVelocity1D", &Richards_::getVelocity1D)
   .def("writeDumuxVTK",&Richards_::writeDumuxVTK)
   .def("setRegularisation",&Richards_::setRegularisation)
   .def("setTopBC",&Richards_::setTopBC)
   .def("setBotBC",&Richards_::setBotBC)
   .def("setSTopBC",&Richards_::setSTopBC)
   .def("setSBotBC",&Richards_::setSBotBC)
   .def("getAvgDensity",&Richards_::getAvgDensity)
   .def("getCSS1_out",&Richards_::getCSS1_out)
   .def("getRF_out",&Richards_::getRF_out)
   .def("getAvgDensity",&Richards_::getAvgDensity);

}


#endif
