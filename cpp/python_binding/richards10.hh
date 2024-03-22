#ifndef PYTHON_RICHARDS10_SOLVER_H_
#define PYTHON_RICHARDS10_SOLVER_H_

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
class Richards10 : public Richards<Problem, Assembler, LinearSolver, dim> {
public:
    using NumEqVector = typename Problem::NumEqVector;
	
	
	void setComputeDtCSS2(const std::function<double(double,double,double)>& s)
    {
        this->problem->computeDtCSS2 = s;
    }
    
	double computeDtCSS2(double CSS1, double CSW, double CSS2)
    {
        return this->problem->computeDtCSS2(CSS1, CSW, CSS2);
    }
    
	double computeInitCSS2(double CSS1, double CSW)
    {
        return this->problem->computeInitCSS2_(CSS1, CSW);
    }
    
	void setFaceGlobalIndexSet(std::map<int,int>  faceIdx)
	{
		this->problem->setFaceGlobalIndexSet(faceIdx);// = faceIdx;
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
	
    std::vector<std::array<double,9>> getFluxScvf() {
		std::vector<NumEqVector> fluxScvf = this->problem->getFluxScvf10c_();
		std::vector<std::array<double,9>> fluxScvfOut(fluxScvf.size());
		for(size_t elem = 0; elem < fluxScvf.size(); elem ++)
		{
			for(size_t elem_ = 0; elem_ < fluxScvf.size(); elem_ ++)
			{
				fluxScvfOut.at(elem).at(elem_) = fluxScvf.at(elem)[elem_];
			}
		}
    	return fluxScvfOut; // 
    }
    std::vector<int> getFluxScvf_idx() {
    	return this->problem->idxScv4FluxScv_10c_(); // 
    }
	
	void resetSetFaceFlux(){this->problem->resetSetFaceFlux();}
	
	std::vector<NumEqVector> getProblemFlux()
	{	
		return this->problem->getFluxScvf10c_();
	}
	
	std::vector<NumEqVector> getProblemSource()
	{	
		return this->problem->getSource10c_();
	}
	double computeRF(double C_S_W, double theta, double svc_volume)
	{	
		return this->problem->computeRF(C_S_W, theta, svc_volume);
	}
	double computeCSS1(double C_S_W, double theta, double svc_volume)
	{	
		return this->problem->computeCSS1(C_S_W, theta, svc_volume);
	}
	std::vector<double> computeRFs()
	{	
		std::vector<double> C_S_fr = this->getSolution(1); // [mol / mol] wat
		std::vector<double> theta = this->getWaterContent();// [m3/m3]
		std::vector<double> svc_volume = this->getCellVolumes(); // [m3]
		std::vector<double> RFs(C_S_fr.size());// - 
		for(int rfIdx = 0; rfIdx < RFs.size(); rfIdx ++)
		{
			RFs.at(rfIdx) = this->problem->computeRF(C_S_fr.at(rfIdx)*this->molarDensityWat_m3, // [mol/m3]
													theta.at(rfIdx), 
													svc_volume.at(rfIdx));//m3
		}
		return RFs;
	}
	std::vector<double> computeCSS1s()
	{	
		std::vector<double> C_S_fr = this->getSolution(1); // [mol / mol] wat
		std::vector<double> theta = this->getWaterContent();// [m3/m3]
		std::vector<double> svc_volume = this->getCellVolumes(); // [m3]
		std::vector<double> CSS1s(C_S_fr.size());// - 
		for(int rfIdx = 0; rfIdx < CSS1s.size(); rfIdx ++)
		{
			CSS1s.at(rfIdx) = this->problem->computeCSS1(C_S_fr.at(rfIdx)*this->molarDensityWat_m3, // [mol/m3]
													theta.at(rfIdx), 
													svc_volume.at(rfIdx));//m3
		}
		return CSS1s;// mol/m3 zone 1
	}
	
	// other values to reset?
    virtual void resetInnerVals() {this->problem->setCSS1_(CSS1_saved);}
    virtual void saveInnerVals() { CSS1_saved = this->problem->getCSS1_();}
    virtual void resetInnerValsManual() {this->problem->setCSS1_(CSS1_savedManually);}
    virtual void saveInnerValsManual() { CSS1_savedManually = this->problem->getCSS1_();}
	
	std::vector<double> CSS1_saved;
	std::vector<double> CSS1_savedManually;
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
void init_richards_10(py::module &m, std::string name) {
    using Richards_ = Richards10<Problem, Assembler, LinearSolver>;
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
   .def("getFluxScvf",&Richards_::getFluxScvf)
   .def("getFluxScvf_idx",&Richards_::getFluxScvf_idx)
   //.def("getAvgDensity",&Richards_::getAvgDensity)
   .def("getCSS1_out",&Richards_::getCSS1_out)
   .def("getRF_out",&Richards_::getRF_out)
   .def("getSorp",&Richards_::getSorp)
   .def("getReac_CSS2",&Richards_::getReac_CSS2)
   //.def("getProblemFlux_10c",&Richards_::getProblemFlux_10c)
   .def("computeRF",&Richards_::computeRF)
   .def("computeCSS1",&Richards_::computeCSS1)
   .def("computeRFs",&Richards_::computeRFs)
   .def("computeCSS1s",&Richards_::computeCSS1s)
   .def_readonly("dimWorld", &Richards_::dimWorld);

}


#endif
