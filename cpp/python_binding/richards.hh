#ifndef PYTHON_RICHARDS_SOLVER_H_
#define PYTHON_RICHARDS_SOLVER_H_

// most includes are in solverbase
#include "solverbase.hh"

#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
// #include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

// writeDumuxVTK
#include <dumux/io/vtkoutputmodule.hh>



/**
 * Adds solver functionality, that specifically makes sense for Richards equation
 */
template<class Problem, class Assembler, class LinearSolver, int dim = 3>
class Richards : public SolverBase<Problem, Assembler, LinearSolver, dim> {
public:

    using MaterialLaw = Dumux::EffToAbsLaw<Dumux::RegularizedVanGenuchten<double>>;
    using MaterialLawParams = typename MaterialLaw::Params;

    virtual ~Richards() { }

    /**
     * Calls parent, additionally turns file output off
     */
    void initialize(std::vector<std::string> args_, bool verbose = true) override {
        SolverBase<Problem, Assembler, LinearSolver, dim>::initialize(args_, verbose);
        this->setParameter("Soil.Output.File", "false");
    }

    /**
     * Sets the source term of the problem.
     *
     * The source is given per cell (Dumux element),
     * as a map with global element index as key, and source as value
     *
     * for simplicity avoiding mpi broadcasting or scattering
     *
     * @param sourceMap 		for each global cell index the source or sink in [kg/s]
     * @param eqIdx				the equation index (default = 0)
     */
    virtual void setSource(const std::map<int, double>& sourceMap, int eqIdx = 0) {
    	this->checkInitialized();
        int n = this->gridGeometry->gridView().size(0);
        std::shared_ptr<std::vector<double>> ls = std::make_shared<std::vector<double>>(n);
        std::fill(ls->begin(), ls->end(), 0.);
        for (const auto& e : elements(this->gridGeometry->gridView())) { // local elements
            int gIdx = this->cellIdx->index(e); // global index
            auto eIdx = this->gridGeometry->elementMapper().index(e);
            if (sourceMap.count(gIdx)>0) {
                //std::cout << "rank: "<< this->rank << " setSource: global index " << gIdx << " local index " << eIdx << "\n" << std::flush;
                ls->at(eIdx) = sourceMap.at(gIdx);
            }
        }
        this->problem->setSource(ls, eqIdx);
    }

    /**
     * Applies source term (operator splitting) to the current solution
     *
     * Limits soil pressure with wilting point from below
     *
     * @param dt 			current time step [s]
     * @param soilFluxes 	fluxes between root and soil [kg/s] (negative fluxes are a sink)
     * @param critP 		critical Pressure or wilting point [Pa]
     */
    virtual void applySource(double dt, const std::map<int, double>& soilFluxes, double critP) {
    	if (this->isBox) {
            throw std::invalid_argument("SolverBase::setInitialCondition: Not implemented yet (sorry)");
        } else {
        	if (cellVolume.size()==0) { // buffer cell volumes [m3]
        		cellVolume = this->getCellVolumes();
        	}
        	const auto& params = this->problem->spatialParams();
        	auto wc = this->getWaterContent();
            for (const auto& e : Dune::elements(this->gridGeometry->gridView())) {  // local elements
            	auto eIdx = this->gridGeometry->elementMapper().index(e);
            	int gIdx = this->cellIdx->index(e); // global index
            	if (soilFluxes.count(gIdx)>0) {
            		const auto& param = params.materialLawParams(e);
                	double sink = soilFluxes.at(gIdx)/1000.; // [kg / s] -> [m3]
                	double new_wc = (wc[eIdx]*cellVolume[eIdx] + dt*sink) / cellVolume[eIdx];
                	// new_wc -> effective saturation
                    double sw = new_wc/params.porosity(e); // water_content -> saturation
                	double p = pRef_-MaterialLaw::pc(param, sw);
                	if (p<critP) {
                    	this->x[eIdx] = critP;
                	} else {
                		this->x[eIdx] = p;
                	}
            	}
            }
        }
    }

    /**
     * Sets critical pressure, used for constantFlux, constantFluxCyl, or atmospheric boundary conditions,
     * to limit maximal the flow.
     *
     * @param critical 		the critical pressure or wilting point [Pa]
     */
    void setCriticalPressure(double critical) {
    	this->checkInitialized();
    	this->problem->criticalPressure(critical); // problem is defined in solverbase.hh
    }

    /**
     * Sets the initial conditions, for a MPI process (TODO for more than 1 equation)
     *
     *  @param init         globally shared initial data, sorted by global index [Pa]
     */
    virtual void setInitialConditionHead(std::vector<double> init) {
        if (this->isBox) {
            throw std::invalid_argument("Richards::setInitialCondition: Not implemented yet (sorry)");
        } else {
            for (const auto& e : Dune::elements(this->gridGeometry->gridView())) {
                int eIdx = this->gridGeometry->elementMapper().index(e);
                int gIdx = this->cellIdx->index(e);
                this->x[eIdx] = toPa(init[gIdx]);
            }
        }
    }

    /**
     * Returns the current solution for a single mpi process in cm pressure head at a global cell index
     * @see SolverBase::getSolutionAt
     */
    virtual double getSolutionHeadAt(int gIdx, int eqIdx = 0) {
        double sol = this->getSolutionAt(gIdx, eqIdx); // Pa
        return toHead(sol);
    }

    /**
     * Returns the current solution for a single mpi process in cm pressure head. @see SolverBase::getSolution
     * Gathering and mapping is done in Python
     */
    virtual std::vector<double> getSolutionHead(int eqIdx = 0) {
        std::vector<double> sol = this->getSolution(eqIdx);
        std::transform(sol.begin(), sol.end(), sol.begin(), std::bind1st(std::plus<double>(), -pRef_));
        std::transform(sol.begin(), sol.end(), sol.begin(), std::bind1st(std::multiplies<double>(), 100. / rho_ / g_));
        return sol;
    }

    /*
     * TODO setLayers(std::map<int, int> l)
     */

    /**
     * Returns the current solution for a single mpi process.
     * Gathering and mapping is done in Python
     */
    virtual std::vector<double> getWaterContent() {
    	int n =  this->checkInitialized();
        std::vector<double> theta;
        theta.reserve(n);
        for (const auto& element : Dune::elements(this->gridGeometry->gridView())) { // soil elements
            double t = 0;
            auto fvGeometry = Dumux::localView(*this->gridGeometry); // soil solution -> volume variable
            fvGeometry.bindElement(element);
            auto elemVolVars = Dumux::localView(this->gridVariables->curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, this->x);
            int c = 0;
            for (const auto& scv : scvs(fvGeometry)) {
                c++;
                t += elemVolVars[scv].waterContent();
            }
            theta.push_back(t/c); // mean value
        }
        return theta;
    }


    /**
     * Returns the current solution for a single mpi process.
     * Gathering and mapping is done in Python
     */
    virtual std::vector<double> getSaturation() {
    	int n =  this->checkInitialized();
        std::vector<double> s;
        s.reserve(n);
        for (const auto& element : Dune::elements(this->gridGeometry->gridView())) { // soil elements
            double t = 0;
            auto fvGeometry = Dumux::localView(*this->gridGeometry); // soil solution -> volume variable
            fvGeometry.bindElement(element);
            auto elemVolVars = Dumux::localView(this->gridVariables->curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, this->x);
            int c = 0;
            for (const auto& scv : scvs(fvGeometry)) {
                c++;
                t += elemVolVars[scv].saturation();
            }
            s.push_back(t/c); // mean value
        }
        return s;
    }

    /**
     * Returns the total water volume [m3] within the domain, TODO wrong because of overlapping cells!
     */
    virtual double getWaterVolume()
    {
        this->checkInitialized();
        double cVol = 0.;
        for (const auto& element : Dune::elements(this->gridGeometry->gridView())) { // soil elements
            auto fvGeometry = Dumux::localView(*this->gridGeometry); // soil solution -> volume variable
            fvGeometry.bindElement(element);
            auto elemVolVars = Dumux::localView(this->gridVariables->curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, this->x);
            for (const auto& scv : scvs(fvGeometry)) {
                cVol += elemVolVars[scv].waterContent()*scv.volume();
            }
        }
        return this->gridGeometry->gridView().comm().sum(cVol);
    }

	/**
	 * Return the darcy (?) velocity in a 1D model TODO not working even in 1D TODO (for nD we would need to multiply with the outer normals)
	 *
	 * For a single mpi process. Gathering is done in Python
	 */
	virtual std::vector<double> getVelocity1D() {
		int n = this->checkInitialized();
		std::vector<double> v;
		v.resize(n);
		for (const auto& e : Dune::elements(this->gridGeometry->gridView())) { // soil elements
			double f = 0.;
			auto fvGeometry = Dumux::localView(*this->gridGeometry); // soil solution -> volume variable
			fvGeometry.bindElement(e);
			for (const auto& scvf : scvfs(fvGeometry)) {
				auto elemVolVars = Dumux::localView(this->gridVariables->curGridVolVars());
				elemVolVars.bindElement(e, fvGeometry, this->x);
				f += this->problem->neumann(e, fvGeometry, elemVolVars, scvf)[0]/1000.; // [kg / (m2 s)] -> [m/s]
			}
			v[this->cellIdx->index(e)] = f;
		}
		return v;
	}

    /**
     * Uses the Dumux VTK Writer
     */
    virtual void writeDumuxVTK(std::string name)
    {
        Dumux::VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*this->gridVariables, this->x, name);
        // using VelocityOutput = PorousMediumFlowVelocityOutput<GridVariables> // <- can't get this type without TTAG :-(
        // vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
        vtkWriter.addVolumeVariable([](const auto& volVars){ return volVars.pressure(); }, "p (Pa)");
        vtkWriter.addVolumeVariable([](const auto& volVars){ return volVars.saturation(); }, "saturation");
        vtkWriter.write(0.0);
    }

    /**
     * Call to change default setting (of 1.e-6 for both)
     *
     * pcEps    capillary pressure regularisation
     * krEps 	relative permeabiltiy regularisation
     */
    virtual void setRegularisation(double pcEps, double krEps) {
    	this->checkInitialized();
    	this->problem->setRegularisation(pcEps, krEps);
    }


    /**
     * Changes boundary condition in initialized problem (e.g. within the simulation loop)
     * Dirichlet types: value [cm]
     * Neumann types: value [cm/day] = [cm3/cm2/day]
     */
    void setTopBC(int type, double value) {
    	this->checkInitialized();
    	this->problem->bcTopType_ = type;
    	this->problem->bcTopValues_[0] = value;
    }

    /**
     * Changes boundary condition in initialized problem (e.g. within the simulation loop)
     * Dirichlet types: value [cm]
     * Neumann types: value [cm/day] = [cm3/cm2/day]
     */
    void setBotBC(int type, double value) {
    	this->checkInitialized();
    	this->problem->bcBotType_ = type;
    	this->problem->bcBotValues_[0] = value;
    }

    /**
     * Changes solute boundary condition in initialized problem (e.g. within the simulation loop)
     * Dirichlet types: value [cm] TODO Units?
     * Neumann types: value [cm/day] = [cm3/cm2/day]
     */
    void setSTopBC(std::vector<int> types, std::vector<double> values) {
    	this->checkInitialized();
    	this->problem->bcSTopType_ = types;
    	this->problem->bcSTopValue_ = values;
    }

    /**
     * Changes solute boundary condition in initialized problem (e.g. within the simulation loop)
     * Dirichlet types: value [cm] TODO Units?
     * Neumann types: value [cm/day] = [cm3/cm2/day]
     */
    void setSBotBC(std::vector<int> types, std::vector<double> values) {
    	this->checkInitialized();
    	this->problem->bcSBotType_ = types;
    	this->problem->bcSBotValue_ = values;
    }

protected:

    std::vector<double> cellVolume;

    using SolutionVector = typename Problem::SolutionVector;
    using GridVariables = typename Problem::GridVariables;

	static constexpr double g_ = 9.81; // cm / s^2 (for type conversions)
	static constexpr double rho_ = 1.e3; // kg / m^3 (for type conversions)
	static constexpr double pRef_ = 1.e5; // Pa

	//! cm pressure head -> Pascal
	static double toPa(double ph) {
		return pRef_ + ph / 100. * rho_ * g_;
	}
	//! Pascal -> cm pressure head
	static double toHead(double p) {
		return (p - pRef_) * 100. / rho_ / g_;
	}

};

/**
 * pybind11
 */
template<class Problem, class Assembler, class LinearSolver>
void init_richards(py::module &m, std::string name) {
    using Richards_ = Richards<Problem, Assembler, LinearSolver>;
	py::class_<Richards_, SolverBase<Problem, Assembler, LinearSolver>>(m, name.c_str())
   .def(py::init<>())
   .def("initialize", &Richards_::initialize, py::arg("args_") = std::vector<std::string>(0), py::arg("verbose") = true)
   .def("setSource", &Richards_::setSource, py::arg("sourceMap"), py::arg("eqIdx") = 0)
   .def("applySource", &Richards_::applySource)
   .def("setCriticalPressure", &Richards_::setCriticalPressure)
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
   .def("setSBotBC",&Richards_::setSBotBC);

}


#endif
