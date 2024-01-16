#ifndef PYTHON_RICHARDS_10CYL_SOLVER_H_
#define PYTHON_RICHARDS_10CYL_SOLVER_H_

//#include "richards.hh" // most includes are in solverbase
#include "richards_cyl.hh" // most includes are in 
//#include "richards_10.hh" // most includes are in solverbase

// writeDumuxVTK
//#include <dumux/io/vtkoutputmodule.hh>


/**
 * Adds solver functionality, that specifically makes sense for Richards equation
 */
template<class Problem, class Assembler, class LinearSolver, int dim = 1>
class Richards10Cyl : public RichardsCyl<Problem, Assembler, LinearSolver, dim> 
{
	public:

    using NumEqVector = typename Problem::NumEqVector;
    virtual ~Richards10Cyl() { }
	
	
	void setFaceGlobalIndexSet(std::map<int,int>  faceIdx)
	{
		this->problem->setFaceGlobalIndexSet(faceIdx);// = faceIdx;
	}
	
	std::vector<double> getCellVolumes(){
		return this->problem->cellVolumesCyl;
	}
	
	double segLength(){
		return this->problem->segLength;// m
	}
	
    std::vector<NumEqVector> getFluxScvf10c() {
    	return this->problem->getFluxScvf10c_(); // 
    }
    std::vector<int> idxScv4FluxScv_10c() {
    	return this->problem->idxScv4FluxScv_10c_(); // 
    }
    std::vector<double> getCSS1_out() {
    	return this->problem->getCSS1_(); // 
    }
    std::vector<double> getRF_out() {
    	return this->problem->getRF_(); // 
	}
    
    std::vector<double> getSorp() {
    	return this->problem->getSorp_(); // 
    }
    std::vector<double> getReac_CSS2() {
    	return this->problem->getReac_CSS2_(); // 
    }
	
	void resetSetFaceFlux(){this->problem->resetSetFaceFlux();}
	
	std::vector<NumEqVector> getProblemFlux_10c()
	{	
		return this->problem->getFluxScvf10c_();
	}
	std::vector<NumEqVector> getProblemSource_10c()
	{	
		return this->problem->getSource10c_();
	}
	
	//other values to reset?
    virtual void resetInnerVals() {this->problem->setCSS1_(CSS1_saved);}
    virtual void saveInnerVals() { CSS1_saved = this->problem->getCSS1_();}
    virtual void resetInnerValsManual() {this->problem->setCSS1_(CSS1_savedManually);}
    virtual void saveInnerValsManual() { CSS1_savedManually = this->problem->getCSS1_();}
    

	double computeRF(double C_S_W, double theta, double svc_volume)
	{	
		return this->problem->computeRF( C_S_W,  theta,  svc_volume);
	}
	double computeCSS1(double C_S_W, double theta, double svc_volume)
	{	
		return this->problem->computeCSS1( C_S_W,  theta,  svc_volume);
	}
	std::vector<double> computeRFs()
	{	
		std::vector<double> C_S_fr = this->getSolution(1); // [mol / mol] wat
		std::vector<double> theta = this->getWaterContent();// [m3/m3]
		//std::vector<double> svc_volume = this->problem->cellVolumesCyl; // [m3]
		std::vector<double> RFs(C_S_fr.size());// - 
		for(int rfIdx = 0; rfIdx < RFs.size(); rfIdx ++)
		{
			double svc_volume =  this->problem->getCellVolumesCyl(rfIdx);
			RFs.at(rfIdx) = this->problem->computeRF(C_S_fr.at(rfIdx)*this->molarDensityWat_m3, // [mol/m3]
													theta.at(rfIdx), 
													svc_volume);//m3
		}
		return RFs;
	}
	std::vector<double> computeCSS1s()
	{	
		std::vector<double> C_S_fr = this->getSolution(1); // [mol / mol] wat
		std::vector<double> theta = this->getWaterContent();// [m3/m3]
		//std::vector<double> svc_volume = this->problem->cellVolumesCyl; // [m3]
		std::vector<double> CSS1s(C_S_fr.size());// - 
		for(int rfIdx = 0; rfIdx < CSS1s.size(); rfIdx ++)
		{
			double svc_volume =  this->problem->getCellVolumesCyl(rfIdx);
			CSS1s.at(rfIdx) = this->problem->computeCSS1(C_S_fr.at(rfIdx)*this->molarDensityWat_m3, // [mol/m3]
													theta.at(rfIdx), 
													svc_volume);//m3
		}
		return CSS1s;
	}
	std::vector<double> CSS1_saved;
	std::vector<double> CSS1_savedManually;
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
template<class Problem, class Assembler, class LinearSolver, int dim = 1>
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

    .def_readwrite("BC_in_vals", &RichardsFoam::BC_in_vals) 
    .def_readwrite("BC_out_vals", &RichardsFoam::BC_out_vals) 
    .def_readwrite("BC_time", &RichardsFoam::BC_time) 
    .def_readwrite("BC_ddt", &RichardsFoam::BC_ddt) 
   .def_readonly("outerIdx",&RichardsFoam::outerIdx)
   .def_readonly("rIn",&RichardsFoam::rIn)
   .def_readonly("rOut",&RichardsFoam::rOut)
   .def_readonly("dimWorld", &RichardsFoam::dimWorld)
   .def("segLength", &RichardsFoam::segLength)
   .def("getFluxScvf10c",&RichardsFoam::getFluxScvf10c)
   .def("idxScv4FluxScv_10c",&RichardsFoam::idxScv4FluxScv_10c)
   .def("computeRF",&RichardsFoam::computeRF)
   .def("computeCSS1",&RichardsFoam::computeCSS1)
   .def("computeRFs",&RichardsFoam::computeRFs)
   .def("computeCSS1s",&RichardsFoam::computeCSS1s)
   .def("getCSS1_out",&RichardsFoam::getCSS1_out)
   .def("getRF_out",&RichardsFoam::getRF_out)
   .def("getSorp",&RichardsFoam::getSorp)
   .def("getReac_CSS2",&RichardsFoam::getReac_CSS2)
   .def("getProblemFlux_10c",&RichardsFoam::getProblemFlux_10c)
   .def("getAvgDensity",&RichardsFoam::getAvgDensity)
   ;
}



#endif
