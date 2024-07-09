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

    using GlobalPosition = typename Problem::GlobalPosition; 
	
    virtual ~RichardsCyl() { }

    /**
     * Calls parent, additionally turns gravity off
     */
    virtual void initialize(std::vector<std::string> args_ = std::vector<std::string>(0), bool verbose = true, bool doMPI = true) override {
    	Richards<Problem, Assembler, LinearSolver, dim>::initialize(args_, verbose, doMPI);
        this->setParameter("Problem.EnableGravity", "false"); // important in 1d axial-symmetric problem
    }

    /**
     * Adds inner and outer radius and element indices
     * TODO currently only for 1D
     */
	virtual void initializeProblem(double maxDt = -1) {
		SolverBase<Problem, Assembler, LinearSolver, dim>::initializeProblem(maxDt);
		auto minMax = this->getGridBounds();
		rIn = minMax[0];
		innerIdx = this->pick({ rIn });
		rOut = minMax[1];
    	outerIdx = this->pick({ rOut });
	}
    
	
    /**
     * for 1d3d coupling, useful to get from dumux 
	 * the boundary flows at each time step
	 * more accurate than getNeumann()
	 * currently only implemented for richards_cyl
     */
	 
    void clearSaveBC() {    
        BC_in_vals.clear();
        BC_out_vals.clear();
        BC_time.clear();
        BC_ddt.clear();
        potential_in_vals.clear();
        kr_in_vals.clear();
    }
    void doSaveBC(double ddt_current) {
		int numC = this->numComp();
        std::vector<double> BC_in_vals_i(numC);
        std::vector<double> BC_out_vals_i(numC);
        std::vector<double> potential_in_vals_i(numC);
		
            
        for(int nc = 0.; nc < numC; nc++)
        {
			potential_in_vals_i.at(nc) = this->getInnerPot(nc);// cm or mol/mol
            BC_in_vals_i.at(nc) = this->getInnerFlux(nc)/this->rIn;// [ mol / (m^2 \cdot s)]_axissymmetric / [axyssimetric factor] = [ mol / (m^2 * s)]
            BC_out_vals_i.at(nc) = this->getOuterFlux(nc)/this->rOut;// [ mol / (m^2 \cdot s)]_axissymmetric  / [axyssimetric factor] = [ mol / (m^2 * s)]
        }
        BC_in_vals.push_back(BC_in_vals_i);
        BC_out_vals.push_back(BC_out_vals_i);
		potential_in_vals.push_back(potential_in_vals_i);
        BC_ddt.push_back(ddt_current);// s
		kr_in_vals.push_back(this->getInnerKrw());//  [1/s]
    }
	
	
    /**
     * soil conductance at inner face
     */
    virtual double getInnerKrw() {
        double Krmean = 0.;
		int gIdx = innerIdx;
        if (this->localCellIdx.count(gIdx)>0) {
            int eIdx = this->localCellIdx[gIdx];
            auto e = this->gridGeometry->element(eIdx);
            auto fvGeometry = Dumux::localView(*this->gridGeometry); // soil solution -> volume variable
            fvGeometry.bindElement(e);
            auto elemVolVars = Dumux::localView(this->gridVariables->curGridVolVars());
            elemVolVars.bindElement(e, fvGeometry, this->x);
			
			int c = 0; double t = 0.;
			GlobalPosition ePos = e.geometry().center();
				
			double kc = this->problem->spatialParams().hydraulicConductivity(e); //  [m/s]
			auto materialLaw_ = this->problem->materialLaw(e);
			int h2OIdx = 0;
            for (const auto& scvf : scvfs(fvGeometry)) {
				GlobalPosition pos = scvf.center();
				if ( this->onUpperBoundary_(pos) || this->onLowerBoundary_(pos) ) {
					c++;
					double dz = std::fabs(ePos[this->dimWorld - 1] - pos[this->dimWorld - 1]); // m	
					//auto& volVars = elemVolVars[scvf.insideScvIdx()];
					double s =  elemVolVars[scvf.insideScvIdx()].saturation(h2OIdx);
					double krw = materialLaw_.krw(s);//	The relative permeability for the wetting phase [between 0 and 1]
					
					t += krw * kc/dz; // 1/s
				}

            }
			Krmean = t/c; // mean value
        }
        return Krmean; // so clever
    }
	
    /**
     * [ kg / (m^2 \cdot s)]
     */
    double getInnerPot(int eqIdx = 0) {
		if (eqIdx == 0)
		{
			return this->getSolutionHeadAt(innerIdx); //  [ cm]
		} else{
			return this->getSolutionAt(innerIdx,eqIdx); //  [mol/mol] or [kg/kg]
		}
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
    double getInnerHead(int shift = 0) {
    	return this->getSolutionHeadAt(innerIdx+shift);
    }

    /**
     * Gets the concentration at the inner boundary [cm]
     */
    double getInnerSolutes(int shift = 0, int compId = 1) {
        return this->getSolutionAt(innerIdx+shift,compId);
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
    /**
     * Gets the pressure head at the inner boundary [cm]
     */
    virtual double getWaterVolume() {
        std::vector<double> volumes = this->getCellVolumesCyl();
        std::vector<double> water_content = this->getWaterContent();
        double s = 0.;
        for (int i =0; i<volumes.size(); i++)
            s += volumes.at(i)*water_content.at(i);
        return s;
    }

    int innerIdx = -1;
    int outerIdx = -1;
    double rIn = 0.;
    double rOut = 0.;

    std::vector<std::vector<double>> BC_in_vals;
    std::vector<std::vector<double>> BC_out_vals;
    std::vector<std::vector<double>> potential_in_vals;
    std::vector<double> BC_time;
    std::vector<double> BC_ddt;
    std::vector<double> kr_in_vals;
};

/**
 * pybind11
 */
template<class Problem, class Assembler, class LinearSolver, int dim = 1>
void init_richards_cyl(py::module &m, std::string name) {
    using RichardsFoam = RichardsCyl<Problem, Assembler, LinearSolver>;
	py::class_<RichardsFoam, SolverBase<Problem, Assembler, LinearSolver, dim>>(m, name.c_str())
   .def(py::init<>())
   .def("initialize", &RichardsFoam::initialize, py::arg("args_") = std::vector<std::string>(0),
											py::arg("verbose") = true, py::arg("doMPI") = true)
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
   .def("setSTopBC",&RichardsFoam::setSTopBC)
   .def("setSBotBC",&RichardsFoam::setSBotBC)

   .def("getInnerFlux",&RichardsFoam::getInnerFlux, py::arg("eqIdx") = 0)
   .def("getOuterFlux",&RichardsFoam::getOuterFlux, py::arg("eqIdx") = 0)
   .def("getInnerHead",&RichardsFoam::getInnerHead, py::arg("shift") = 0)
   .def("getInnerSolutes",&RichardsFoam::getInnerSolutes, py::arg("shift") = 0, py::arg("compId") = 1)
   .def("setRootSystemBC",&RichardsFoam::setRootSystemBC)

   .def("numComp",&RichardsFoam::numComp)

    .def_readwrite("BC_in_vals", &RichardsFoam::BC_in_vals) 
    .def_readwrite("BC_out_vals", &RichardsFoam::BC_out_vals) 
    .def_readwrite("BC_time", &RichardsFoam::BC_time) 
    .def_readwrite("BC_ddt", &RichardsFoam::BC_ddt) 
    .def_readwrite("kr_in_vals", &RichardsFoam::kr_in_vals) 
    .def_readwrite("potential_in_vals", &RichardsFoam::potential_in_vals) 

   .def_readonly("innerIdx",&RichardsFoam::innerIdx)
   .def_readonly("outerIdx",&RichardsFoam::outerIdx)
   .def_readonly("rIn",&RichardsFoam::rIn)
   .def_readonly("rOut",&RichardsFoam::rOut);
}


#endif
