#ifndef PYTHON_RICHARDS4_SOLVER_H_
#define PYTHON_RICHARDS4_SOLVER_H_

// most includes are in solverbase
//#include "solverbase.hh"

#include "richards.hh" // most includes are in solverbase
// #include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
// // #include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
// #include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

// // writeDumuxVTK
// #include <dumux/io/vtkoutputmodule.hh>



/**
 * Adds solver functionality, that specifically makes sense for Richards equation
 */
template<class Problem, class Assembler, class LinearSolver, int dim = 3>
class Richards4 : public Richards<Problem, Assembler, LinearSolver, dim> {
public:
    using NumEqVector = typename Problem::NumEqVector;
	
	
	void setComputeDtCSS2(const std::function<double(double,double)>& s)
    {
        this->problem->computeDtCSS2 = s;
    }
    
	double computeDtCSS2(double CSW, double CSS2)
    {
        return this->problem->computeDtCSS2(CSW, CSS2);
    }
    
	double computeInitCSS2(double CSS1, double CSW)
    {
        return this->problem->computeInitCSS2_(CSS1, CSW);
    }
	double computeCSS1(double bulkSoilDensity, double C_S_W, int dofIndex)
	{	
		return this->problem->computeCSS1(bulkSoilDensity, C_S_W, dofIndex);
	}
	
	/**
     * set verbose
     */
    virtual void setVerbose(int verbose) {
        this->checkGridInitialized();
    	this->problem->verbose = verbose;
    }
	
	
};

/**
 * pybind11
 */
template<class Problem, class Assembler, class LinearSolver>
void init_richards_4(py::module &m, std::string name) {
    using Richards_ = Richards4<Problem, Assembler, LinearSolver>;
	py::class_<Richards_, SolverBase<Problem, Assembler, LinearSolver>>(m, name.c_str())
   .def(py::init<>())
   .def("initialize", &Richards_::initialize, py::arg("args_") = std::vector<std::string>(0), 
													py::arg("verbose") = true,py::arg("doMPI") = true)
   .def("setSource", &Richards_::setSource, py::arg("sourceMap"), py::arg("eqIdx") = 0)
   .def("setComputeDtCSS2",&Richards_::setComputeDtCSS2)
   .def("computeDtCSS2",&Richards_::computeDtCSS2)
   .def("computeInitCSS2",&Richards_::computeInitCSS2)
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
   .def("computeCSS1",&Richards_::computeCSS1)
   .def_readonly("dimWorld", &Richards_::dimWorld);

}


#endif
