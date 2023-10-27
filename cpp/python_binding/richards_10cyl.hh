#ifndef PYTHON_RICHARDS_10CYL_SOLVER_H_
#define PYTHON_RICHARDS_10CYL_SOLVER_H_

#include "richards.hh" // most includes are in solverbase
#include "richards_cyl.hh" // most includes are in solverbase

// writeDumuxVTK
#include <dumux/io/vtkoutputmodule.hh>


/**
 * Adds solver functionality, that specifically makes sense for Richards equation
 */
template<class Problem, class Assembler, class LinearSolver, int dim = 1>
class Richards10Cyl : public RichardsCyl<Problem, Assembler, LinearSolver, dim> {
public:

    virtual ~Richards10Cyl() { }

    std::vector<double> getCSS1_out() {
    	return this->problem->getCSS1_(); // 
    }
    std::vector<double> getRF_out() {
    	return this->problem->getRF_(); // 
    }
    
    /**
     * set verbose
     */
    virtual void setVerbose(int verbose) {
        this->checkInitialized();
    	this->problem->verbose = verbose;
    }

};

    
/**
 * pybind11
 */
template<class Problem, class Assembler, class LinearSolver, int dim = 3>
void init_richards_10cyl(py::module &m, std::string name) {
    using RichardsFoam = Richards10Cyl<Problem, Assembler, LinearSolver>;
	py::class_<RichardsFoam, SolverBase<Problem, Assembler, LinearSolver, dim>>(m, name.c_str())
   .def(py::init<>())
   .def("initialize", &RichardsFoam::initialize, py::arg("args_") = std::vector<std::string>(0), 
													py::arg("verbose") = true,py::arg("doMPI") = true)
   .def("initializeProblem", &RichardsFoam::initializeProblem)

   .def("setVerbose",&RichardsFoam::setVerbose)
   .def("setSource", &RichardsFoam::setSource, py::arg("sourceMap"), py::arg("eqIdx") = 0)
   .def("setCriticalPressure", &RichardsFoam::setCriticalPressure)
   .def("getSolutionHead", &RichardsFoam::getSolutionHead, py::arg("eqIdx") = 0)
   .def("getSolutionHeadAt", &RichardsFoam::getSolutionHeadAt, py::arg("gIdx"), py::arg("eqIdx") = 0)
   .def("getWaterContent",&RichardsFoam::getWaterContent)
   .def("getSaturation",&RichardsFoam::getSaturation)
   .def("getKrw",&RichardsFoam::getKrw)
   .def("getVelocity1D", &RichardsFoam::getVelocity1D)
   .def("writeDumuxVTK",&RichardsFoam::writeDumuxVTK)
   .def("setRegularisation",&RichardsFoam::setRegularisation)
   .def("setTopBC",&RichardsFoam::setTopBC)
   .def("setBotBC",&RichardsFoam::setBotBC)
   .def("setSTopBC",&RichardsFoam::setSTopBC)
   .def("setSBotBC",&RichardsFoam::setSBotBC)

   .def("getInnerFlux",&RichardsFoam::getInnerFlux, py::arg("eqIdx") = 0)
   .def("getOuterFlux",&RichardsFoam::getOuterFlux, py::arg("eqIdx") = 0)
   .def("getInnerHead",&RichardsFoam::getInnerHead, py::arg("shift") = 0)
   .def("getInnerSolutes",&RichardsFoam::getInnerSolutes, py::arg("shift") = 0, py::arg("compId") = 0)
   .def("setRootSystemBC",&RichardsFoam::setRootSystemBC)

   .def_readonly("outerIdx",&RichardsFoam::outerIdx)
   .def_readonly("rIn",&RichardsFoam::rIn)
   .def_readonly("rOut",&RichardsFoam::rOut)
   .def("getCSS1_out",&RichardsFoam::getCSS1_out)
   .def("getRF_out",&RichardsFoam::getRF_out)
   .def("getAvgDensity",&RichardsFoam::getAvgDensity);
}



#endif
