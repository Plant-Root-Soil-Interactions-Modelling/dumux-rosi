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

    virtual ~RichardsCyl() { }

    /**
     * Calls parent, additionally turns gravity off
     */
    void initialize(std::vector<std::string> args_ = std::vector<std::string>(0), bool verbose = true) override {
    	Richards<Problem, Assembler, LinearSolver, dim>::initialize(args_, verbose);
        this->setParameter("Problem.EnableGravity", "false"); // important in 1d axial-symmetric problem
        this->setParameter("Soil.Problem.EnableGravity", "false"); // important in 1d axial-symmetric problem
    }

    /**
     * Adds inner and outer radius and element indices
     * TODO currently only for 1D
     */
	virtual void initializeProblem() {
		SolverBase<Problem, Assembler, LinearSolver, dim>::initializeProblem();
		auto minMax = this->getGridBounds();
		rIn = minMax[0];
		innerIdx = this->pick({ rIn });
		rOut = minMax[1];
    	outerIdx = this->pick({ rOut });
	}

    /**
     * [ kg / (m^2 \cdot s)]
     */
    double getInnerFlux(int eqIdx = 0) {
    	return this->getNeumann(innerIdx, eqIdx); //  [ kg / (m^2 \cdot s)]
    }

    /**
     *
     */
    double getOuterFlux(int eqIdx = 0) {
    	return this->getNeumann(outerIdx, eqIdx); //  [ kg / (m^2 \cdot s)]
    }

    /**
     * Gets the pressure head at the inner boundary [cm]
     */
    double getInnerHead() {
    	return this->getSolutionHeadAt(innerIdx);
    }

    /**
     * Changes the exact root system BC (for coupling) in initialized problem (e.g. within the simulation loop)
     * @parm params 	x0, x1, kr, kx, length
     */
    void setRootSystemBC(std::vector<double> params) {
    	if (this->problem->bcBotType_!=6) {
    		std::cout << "RichardsCyl::setRootSystemBC() warning, wrong bcBotTyp is set (!=6) "<< this->problem->bcBotType_<< "\n";
    	}
    	this->problem->bcBotValues_= params;
    }

    // TODO getWaterVolume needs adjusting

    int innerIdx = -1;
    int outerIdx = -1;
    double rIn = 0.;
    double rOut = 0.;

};

/**
 * pybind11
 */
template<class Problem, class Assembler, class LinearSolver, int dim = 3>
void init_richards_cyl(py::module &m, std::string name) {
    using RichardsFoam = RichardsCyl<Problem, Assembler, LinearSolver>;
	py::class_<RichardsFoam, SolverBase<Problem, Assembler, LinearSolver, dim>>(m, name.c_str())
   .def(py::init<>())
   .def("initialize", &RichardsFoam::initialize, py::arg("args_") = std::vector<std::string>(0), py::arg("verbose") = true)
   .def("initializeProblem", &RichardsFoam::initializeProblem)

   .def("setSource", &RichardsFoam::setSource, py::arg("sourceMap"), py::arg("eqIdx") = 0)
   .def("setCriticalPressure", &RichardsFoam::setCriticalPressure)
   .def("getSolutionHead", &RichardsFoam::getSolutionHead, py::arg("eqIdx") = 0)
   .def("getSolutionHeadAt", &RichardsFoam::getSolutionHeadAt, py::arg("gIdx"), py::arg("eqIdx") = 0)
   .def("getWaterContent",&RichardsFoam::getWaterContent)
   .def("getWaterVolume",&RichardsFoam::getWaterVolume)
   .def("getVelocity1D", &RichardsFoam::getVelocity1D)
   .def("writeDumuxVTK",&RichardsFoam::writeDumuxVTK)
   .def("setRegularisation",&RichardsFoam::setRegularisation)
   .def("setTopBC",&RichardsFoam::setTopBC)
   .def("setBotBC",&RichardsFoam::setBotBC)

   .def("getInnerFlux",&RichardsFoam::getInnerFlux, py::arg("eqIdx") = 0)
   .def("getOuterFlux",&RichardsFoam::getOuterFlux, py::arg("eqIdx") = 0)
   .def("getInnerHead",&RichardsFoam::getInnerHead)
   .def("setRootSystemBC",&RichardsFoam::setRootSystemBC)

   .def_readonly("innerIdx",&RichardsFoam::innerIdx)
   .def_readonly("outerIdx",&RichardsFoam::outerIdx)
   .def_readonly("rIn",&RichardsFoam::rIn)
   .def_readonly("rOut",&RichardsFoam::rOut);
}


#endif
