#ifndef PYTHON_RICHARDS_22CYL_SOLVER_H_
#define PYTHON_RICHARDS_22CYL_SOLVER_H_

#include "richards_cyl.hh" // most includes are in solverbase
// writeDumuxVTK
//#include <dumux/io/vtkoutputmodule.hh>


/**
 * Adds solver functionality, that specifically makes sense for Richards equation
 */
template<class Problem, class Assembler, class LinearSolver, int dim = 1>
class Richards22Cyl : public RichardsCyl<Problem, Assembler, LinearSolver, dim> 
{
	public:
	
	//how many reaction terms will be saved, I need a better place to load this
	static constexpr int numReactionTermsgeneral = 8; //general reaction terms, such as how much of the solutes were taken up by the microbes
	static constexpr int numReactionTermsMicrobes = 11; //reaction terms for each microbial species
	static constexpr int numReactionTerms=numReactionTermsgeneral+numReactionTermsMicrobes;

    using NumEqVector = typename Problem::NumEqVector;
    virtual ~Richards22Cyl() { }
	
	
	
	 /**
     * The volume [m3] of each element (vtk cell)
     *
     * This is done for a single process, gathering and mapping is done in Python.
     */
    virtual std::vector<double> getCellVolumesCyl() {
		int numCells = this->checkGridInitialized();
        std::vector<double> volumes;
		volumes.resize(numCells);
        for (int i =0; i<volumes.size(); i++)
            volumes.at(i) = this->problem->getCellVolumesCyl(i);
        return volumes;
    }
    
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
    
	std::vector<double> getCellVolumes(){
		return this->problem->cellVolumesCyl;
	}
	
	double segLength(){
		return this->problem->segLength;// m
	}
	
	

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
	/**
     * set verbose
     */
    virtual void setVerbose(int verbose) {
        this->checkGridInitialized();
    	this->problem->verbose = verbose;
    }
	
    std::vector<int> getSTopBCType() {
    	this->checkGridInitialized();
    	return this->problem->bcSTopType_;
    }
    std::vector<int> getSBotBCType() {
    	this->checkGridInitialized();
    	return this->problem->bcSBotType_;
    }
    std::vector<double> getSTopBCValue() {
    	this->checkGridInitialized();
    	return this->problem->bcSTopValue_;
    }
    std::vector<double> getSBotBCValue() {
    	this->checkGridInitialized();
    	return this->problem->bcSBotValue_;
    }
	
	// virtual void initializeProblem(double maxDt = -1) {
		// // get a pointer to this object as its base version
		// RichardsCyl<Problem, Assembler, LinearSolver, dim>* baseclass;
		// baseclass=this;
		// // call initializeProblem for the baseclass
		// this.RichardsCyl<Problem, Assembler, LinearSolver, dim>::initializeProblem(maxDt);
		// // also assign memory to the saves of the reaction rates
		// this->reactionRates.assign(numReactionTerms,std::vector<double>(this->gridGeometry->numScv(),0.));
	// }
	virtual void assignReactionRates(){
		this->reactionRates.assign(numReactionTerms,std::vector<double>(this->gridGeometry->numScv(),0.));
	}
	
	//overload this function from solverbase to also gather all the reaction terms
	/**
     * For a single mpi process. Gathering is done in Python
	 * [ kg /m^3/s] or [ mol / m^3/s]
     */
    virtual void getScvSourcesAtT(double dt, double outer_dt) {
		this->checkGridInitialized();
		
        auto fvGeometry = Dumux::localView(*this->gridGeometry);
		
        auto elemVolVars = Dumux::localView(this->gridVariables->curGridVolVars());
		

        for (const auto& e : elements(this->gridGeometry->gridView())) { //, Dune::Partitions::interior

            fvGeometry.bind(e);
			
            elemVolVars.bind(e, fvGeometry, this->x);
			
			for (const auto& scv : scvs(fvGeometry))
            {
				double pos0 = 1;
				if(this->dimWorld == 1){
					pos0 = scv.center()[0]; 
				}
				NumEqVector scvfSource_(0.0), scvfSource_dummy(0.0);
				std::vector<double> reactionRates_;
				scvfSource_ = this->problem->source(e, fvGeometry, elemVolVars, scv); // [ kg / (m^3 \cdot s)] or [ mol / (m^3 \cdot s)]
				const auto& volVars = elemVolVars[scv];
				reactionRates_=this->problem->bioChemicalReaction(e, scvfSource_dummy, volVars, pos0, scv);				
				for(int eqIdx = 0; eqIdx < this->nEV.size(); eqIdx ++)
				{
					this->scvSources[eqIdx][scv.dofIndex()] += scvfSource_[eqIdx]/pos0*dt/outer_dt;
				}
				for(int reacIdx = 0; reacIdx < numReactionTerms; reacIdx ++)
				{
					this->reactionRates[reacIdx][scv.dofIndex()] += reactionRates_[reacIdx]/pos0*dt/outer_dt;
				}
			}
		}
		
	}
		
};

	

    
/**
 * pybind11
 */
template<class Problem, class Assembler, class LinearSolver, int dim = 1>
void init_richards_22cyl(py::module &m, std::string name) {
    using RichardsFoam = Richards22Cyl<Problem, Assembler, LinearSolver>;
	py::class_<RichardsFoam, SolverBase<Problem, Assembler, LinearSolver, dim>>(m, name.c_str())
   .def(py::init<>())
   .def("initialize", &RichardsFoam::initialize, py::arg("args_") = std::vector<std::string>(0), 
													py::arg("verbose") = true,py::arg("doMPI") = true)
   .def("initializeProblem", &RichardsFoam::initializeProblem)

   .def("setComputeDtCSS2",&RichardsFoam::setComputeDtCSS2)
   .def("computeDtCSS2",&RichardsFoam::computeDtCSS2)
   .def("computeInitCSS2",&RichardsFoam::computeInitCSS2)
   .def("setVerbose",&RichardsFoam::setVerbose)
   .def("setSource", &RichardsFoam::setSource, py::arg("sourceMap"), py::arg("eqIdx") = 0)
   .def("setCriticalPressure", &RichardsFoam::setCriticalPressure)
   .def("getSolutionHead", &RichardsFoam::getSolutionHead, py::arg("eqIdx") = 0)
   .def("getSolutionHeadAt", &RichardsFoam::getSolutionHeadAt, py::arg("gIdx"), py::arg("eqIdx") = 0)
   .def("getWaterContent",&RichardsFoam::getWaterContent)
   .def("getSaturation",&RichardsFoam::getSaturation)
   .def("getVelocity1D", &RichardsFoam::getVelocity1D)
   .def("writeDumuxVTK",&RichardsFoam::writeDumuxVTK)
   .def("setRegularisation",&RichardsFoam::setRegularisation)

   .def("setTopBC",&RichardsFoam::setTopBC)
   .def("setBotBC",&RichardsFoam::setBotBC)
   .def("setSTopBC",&RichardsFoam::setSTopBC)
   .def("setSBotBC",&RichardsFoam::setSBotBC)
   
   .def("getSTopBCType",&RichardsFoam::getSTopBCType)
   .def("getSBotBCType",&RichardsFoam::getSBotBCType)
   .def("getSTopBCValue",&RichardsFoam::getSTopBCValue)
   .def("getSBotBCValue",&RichardsFoam::getSBotBCValue)
   
   .def("getInnerFlux",&RichardsFoam::getInnerFlux, py::arg("eqIdx") = 0)
   .def("getOuterFlux",&RichardsFoam::getOuterFlux, py::arg("eqIdx") = 0)
   .def("getInnerHead",&RichardsFoam::getInnerHead, py::arg("shift") = 0)
   .def("getInnerSolutes",&RichardsFoam::getInnerSolutes, py::arg("shift") = 0, py::arg("compId") = 1)
   .def("segLength", &RichardsFoam::segLength)
   .def("computeRF",&RichardsFoam::computeRF)
   .def("computeCSS1",&RichardsFoam::computeCSS1)
   .def("computeRFs",&RichardsFoam::computeRFs)
   .def("computeCSS1s",&RichardsFoam::computeCSS1s)
   .def("getAvgDensity",&RichardsFoam::getAvgDensity)
   .def("numComp",&RichardsFoam::numComp)
   .def("getMobility", &RichardsFoam::getMobility)
   .def("getViscosity", &RichardsFoam::getViscosity)
   .def("getPressureHead", &RichardsFoam::getPressureHead)
   .def("getPressure", &RichardsFoam::getPressure)
   .def("getConductivity",&RichardsFoam::getConductivity) 

   .def_readonly("innerIdx",&RichardsFoam::innerIdx)
   .def_readonly("outerIdx",&RichardsFoam::outerIdx)
   .def_readonly("rIn",&RichardsFoam::rIn)
   .def_readonly("rOut",&RichardsFoam::rOut);
}



#endif
