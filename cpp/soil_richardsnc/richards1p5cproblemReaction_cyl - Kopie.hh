// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef RICHARDS1P5C_PROBLEM_HH
#define RICHARDS1P5C_PROBLEM_HH
#include <algorithm>
#include <vector>
#include <dumux/porousmediumflow/problem.hh> // base class

#include "../soil_richards/richardsparams.hh"

#include <dune/common/exceptions.hh>
namespace Dumux {

/*!
 * RichardsProblem:
 * Uses Dumux as an easy to use Richards equation solver,
 * where most parameters can be set dynamically
 */
template <class TypeTag>
class Richards1P5CProblem : public PorousMediumFlowProblem<TypeTag>
{
public:
	
	double decayTestShare = 0.;
	
	double getDecay(){return decayTestShare;}
	void setDecay(double input ) const 
	{
		const_cast<double&>(decayTestShare) = input;
	}
	
	double decayTestShareCs = 0.;
	
	double getDecayCs(){return decayTestShareCs;}
	void setDecayCs(double input ) const 
	{
		const_cast<double&>(decayTestShareCs) = input;
	}
	
	double getTheta(){return theta;}
	void setTheta(double input ) const 
	{
		const_cast<double&>(theta) = input;
	}
	
	// exports, used by the binding
	using Grid = GetPropType<TypeTag, Properties::Grid>;
	using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
	using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
	using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
	using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;

	// other
	using GridView = GetPropType<TypeTag, Properties::GridView>;
	using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
	using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
	using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
	using FVElementGeometry = typename FVGridGeometry::LocalView;
	using SubControlVolume = typename FVGridGeometry::SubControlVolume;
	using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
	using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
	using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
	using Scalar = GetPropType<TypeTag, Properties::Scalar>;
	using Element = typename GridView::template Codim<0>::Entity;
	using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
	using MaterialLaw = typename GetPropType<TypeTag, Properties::SpatialParams>::MaterialLaw;
	using MaterialLawParams = typename MaterialLaw::Params;
	using PointSource = GetPropType<TypeTag, Properties::PointSource>;
	using CouplingManager= GetPropType<TypeTag, Properties::CouplingManager>;
	using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
	using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
	using EffectiveDiffusivityModel = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
	
	
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
	
    static constexpr int numFluidComps = FluidSystem::numComponents;
    static constexpr int numSolidComps = SolidSystem::numComponents;
    static constexpr int numInertSolidComps =  SolidSystem::numInertComponents;

	enum {
		pressureIdx = 0, // index of primary variables
		h2OIdx = FluidSystem::liquidPhaseIdx, // fluid index
		soluteIdx = 1, // 1st solute index
		mucilIdx = 2, // mucil index
		CoAIdx = 3, // active oligotrophes
		CoDIdx = 4, // dormant oligo
		CcAIdx = 5, // active copiotrophes
		CcDIdx = 6, // dormant copio
		soilIdx = SolidSystem::mainCompIdx + numFluidComps,
		
		//don t think EqIndx != pvdIdx when we have just 1 phase
		conti0EqIdx = pressureIdx, // indices of the equations
		transportEqIdx = conti0EqIdx + soluteIdx,
		//transportMucilEqIdx = conti0EqIdx + mucilIdx,//is this correct?

		dimWorld = GridView::dimensionworld,
		
		
		numFluidSolutes = numFluidComps -  soluteIdx,
		numSolidSolutes = numSolidComps -  1,
		numSolutes = numFluidSolutes + numSolidSolutes,


		//!!!!!!!!!!!!
		numComponents_ = numFluidComps + numSolidSolutes,//ignore the soil as component
		//!!!!!!!!!!!!
		
		isBox = GetPropType<TypeTag, Properties::FVGridGeometry>::discMethod == DiscretizationMethod::box
	};

	enum BCTypes {
		constantPressure = 1,
		constantConcentration = 1,
		constantFlux = 2,
		constantFluxCyl = 3,
		atmospheric = 4,
		freeDrainage = 5,
		outflow = 6,
		linear = 7,
		michaelisMenten = 8,
		managed = 9,
		outflowCyl = 10
	};

	enum GridParameterIndex {
		materialLayerNumber = 0
	};

	/*!
	 * \brief Constructor: constructed in the main file
	 */
	Richards1P5CProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
	: PorousMediumFlowProblem<TypeTag>(fvGridGeometry) {

		gravityOn_ = false;//Dumux::getParam<bool>("Problem.EnableGravity", true);

		source_.resize(numComponents_); // numComponents_ equations (currently hard coded, where can I get the value?)
		
		 verbose =  getParam<int>("Problem.verbose", 0);
		 toFile =  getParam<bool>("Problem.toFile", false);
		
			
		
		if(verbose>0)
		{
			const auto& myParams = Parameters::paramTree();
			myParams.report();
		}
		// if(toFile)
		// {
			// const auto& myParams = Parameters::paramTree();
			// const auto& paramKeys = myParams.getValueKeys();
			// for (int i = 1; i < paramKeys.size(); ++i)
			// {
				// myfile_ <<paramKeys[i]<<" "<<myParams[paramKeys[i]]<<std::endl;
			// }
		// }
		// Components
		//std::cout<<"all components "<<typeid(FluidSystem).name()<<" "<<numComponents_<<" "<<numFluidComps<<" "<<numSolidComps<<" "<<numInertSolidComps<<std::endl;
		for(int i = 0; i < numComponents_; i++)//all components except h2o
		{
			//std::cout<<"PorousMediumFlowProblem ID"<<i<<" ";
			
			source_[i] = nullptr;
			if(i ==h2OIdx)
			{
				
				// BC			
				bcTopType_ = getParam<int>("Soil.BC.Top.Type", outflow); 
				bcBotType_ = getParam<int>("Soil.BC.Bot.Type", outflow);
				
				//const auto& myParams = Parameters::paramTree();
				//myParams.report();
				bcTopValues_.at(i) =  getParam<double>("Soil.BC.Top.Value", 0.);
				bcBotValues_.at(i) =  getParam<double>("Soil.BC.Bot.Value", 0.);
				//std::cout<<"start h2o "<<bcTopValues_.at(i)<<" "<<bcBotValues_.at(i)<<" ";
				
				//IC
				//std::cout<<"initialSoil_"<<std::endl;
				initialSoil_.at(i) = InputFileFunction("Soil.IC", "P", "Z", 0., this->spatialParams().layerIFF()); // [cm]([m]) pressure head, conversions hard coded
				// Precipitation & Evaporation, and solute input
				if (bcTopType_==atmospheric) {
					componentInput_.at(i) = InputFileFunction("Climate", "Precipitation", "Time", 0.); // cm/day (day)
				}
			}else{
				
				// BC			
				bcSTopType_.at(i - soluteIdx) = getParam<int>("Soil.BC.Top.C"+std::to_string(i)+"Type", outflow); 
				bcSBotType_.at(i - soluteIdx) = getParam<int>("Soil.BC.Bot.C"+std::to_string(i)+"Type", outflow);
				bcSTopValue_.at(i - soluteIdx) = getParam<Scalar>("Soil.BC.Top.C"+std::to_string(i)+"Value", 0.);
				bcSBotValue_.at(i - soluteIdx) = getParam<Scalar>("Soil.BC.Bot.C"+std::to_string(i)+"Value", 0.);
				//std::cout<<"solute start "<<bcSTopType_.at(i - soluteIdx)<<" "<<bcSBotType_.at(i - soluteIdx)
				//<<" "<< bcSTopValue_.at(i - soluteIdx)<<" "<<bcSBotValue_.at(i - soluteIdx)<<" ";
				
				//IC
				initialSoil_.at(i) = InputFileFunction("Soil.IC", "C"+std::to_string(i), "C"+std::to_string(i)+"Z", 
													0., this->spatialParams().layerIFF()); // kg/m2 or mol/m2
				if (bcSTopType_.at(i - soluteIdx)==managed) {
					componentInput_.at(i) = InputFileFunction(std::to_string(i)+".Component.Managed", "Input", "Time", 0.); // cm/day (day)
					
				}
				
				// Buffer power
				b_.at(i - soluteIdx) = getParam<Scalar>(std::to_string(i)+".Component.BufferPower", 0.);
				freundlichN_.at(i - soluteIdx) = getParam<Scalar>(std::to_string(i)+".Component.FreundlichN", 0.);
				freundlichK_.at(i - soluteIdx) = getParam<Scalar>(std::to_string(i)+".Component.FreundlichK", 0.);
				
				// Uptake params (not used)
				vMax_.at(i - soluteIdx)  =  getParam<Scalar>("RootSystem.Uptake.Vmax", 0.)/24./3600.*1e1; //  [g cm-2 day-1] -> [kg m-2 s-1]
				km_.at(i - soluteIdx)    = getParam<Scalar>("RootSystem.Uptake.Km", 0.)*1e3;  // [g cm-3] -> [kg m-3]
				sigma_.at(i - soluteIdx) = getParam<Scalar>("RootSystem.Uptake.ActiveTransport", 0.); // 1 for passive transport, 0 for active transport
			}
			
			
			componentInput_.at(i).setVariableScale(1./(24.*60.*60.)); //day-> s  
			Scalar g2kg = 1/1000 ;
			Scalar m2_2_cm2 = 10000;
			Scalar unitConversion = useMoles ? m2_2_cm2 : m2_2_cm2 * g2kg; //something else needed? 
			componentInput_.at(i).setFunctionScale(unitConversion/(24.*60.*60.)); // g/(cm2 day) or mol/(cm2 day)  -> kg/(m²*s) or mol/(m²*s) 
			//std::cout<<std::endl;
		}

		criticalPressure_ = getParam<double>("Soil.CriticalPressure", -1.e4); // cm
		criticalPressure_ = getParam<double>("Climate.CriticalPressure", criticalPressure_); // cm
		
		double m3_2_cm3 = 1e6;
		 v_maxL = getParam<double>("Soil.v_maxL", v_maxL_)/(24.*60.*60.); //Maximum reaction rate of enzymes targeting large polymers s
		 K_L = getParam<double>("Soil.K_L", K_L_ ) * m3_2_cm3; //[mol cm-3 soil] * [cm3/m3]=> [mol m-3 soil] 
		 
		 m_max = std::vector<double>{getParam<double>("Soil.m_maxO", m_max_[0])/(24.*60.*60.),
									getParam<double>("Soil.m_maxC", m_max_[1])/(24.*60.*60.)};	//Maximum reaction rate 
		 //for troubleshooting , can have m_maxBis != m_max
		 m_maxBis = std::vector<double>{getParam<double>("Soil.m_maxOBis", m_max[0]*(24.*60.*60.))/(24.*60.*60.),
							getParam<double>("Soil.m_maxCBis", m_max[1]*(24.*60.*60.))/(24.*60.*60.)};	//Maximum reaction rate 
		 micro_max = std::vector<double>{getParam<double>("Soil.micro_maxO", micro_max_[0])/(24.*60.*60.),
										getParam<double>("Soil.micro_maxC", micro_max_[1])/(24.*60.*60.)}; //Maximum reaction rate 
		 k_S = std::vector<double>{getParam<double>("Soil.k_SO", k_S_[0])/(24.*60.*60.),
									getParam<double>("Soil.k_SC", k_S_[1])/(24.*60.*60.)}; //[mol soil / mol C soil / s] 
		 k_SBis = std::vector<double>{getParam<double>("Soil.k_SOBis", k_S[0]*(24.*60.*60.))/(24.*60.*60.),
							getParam<double>("Soil.k_SCBis", k_S[1]*(24.*60.*60.))/(24.*60.*60.)};	//[mol soil / mol C soil / s] 
		 k_D = std::vector<double>{getParam<double>("Soil.k_DO",k_D_[0])/(24.*60.*60.),
									getParam<double>("Soil.k_DC",k_D_[1])/(24.*60.*60.)}; // [s-1] 
		 k_R = std::vector<double>{getParam<double>("Soil.k_RO",k_R_[0])/(24.*60.*60.),
									getParam<double>("Soil.k_RC",k_R_[1])/(24.*60.*60.)}; // [s-1] 
		 
		 
		 beta = std::vector<double>{getParam<double>("Soil.betaO", beta_[0]),
									getParam<double>("Soil.betaC", beta_[1])}; //
		 k_growth =std::vector<double>{getParam<double>("Soil.k_growthO", k_growth_[0]),
									getParam<double>("Soil.k_growthC", k_growth_[1])}; //
		 k_growthBis = getParam<double>("Soil.k_growthBis", k_growthBis); //

		 C_S_W_thres = std::vector<double>{getParam<double>("Soil.C_S_W_thresO", C_S_W_thres_[0]) * m3_2_cm3,
										getParam<double>("Soil.C_S_W_thresC", C_S_W_thres_[1]) * m3_2_cm3}; //mol C/m3 soil water
		 	 

		 k_phi = getParam<double>("Soil.k_phi", k_phi_); //
		 k_decay = getParam<double>("Soil.k_decay", k_decay_); // 
		 k_decay2 = getParam<double>("Soil.k_decay2", k_decay2_);//
		 k_decay3 = getParam<double>("Soil.k_decay3", k_decay3);//
		 extra = getParam<double>("Soil.extra", extra);//  
		 extra2 = getParam<double>("Soil.extra2", extra2);// 
		m_maxBisO = getParam<double>("Soil.m_maxBisO2", 0. );
		m_maxBis_Cs = getParam<double>("Soil.m_maxBis_Cs", 0.  );
		 

		// Output
		std::string filestr = this->name() + "_1p5cProblem.txt"; // output file
		myfile_.open(filestr.c_str());
		std::cout << "Richards1P5CProblemR_cyl constructed: bcTopType " << bcTopType_ << ", " << bcTopValues_.at(0) << "; bcBotType "
				<<  bcBotType_ << ", " << bcBotValues_.at(0) 
				<< ", gravitation " << gravityOn_ <<", Critical pressure "
				<< criticalPressure_ //<< "\n" << "Sorption:" << "buffer power "<< b_[0] << ", Freundlich " << freundlichK_[0] << ", " <<
				//freundlichN_[0] << "\n" 
				<< std::flush;
	}

	/**
	 * \brief Eventually, closes output file
	 */
	~Richards1P5CProblem() {
		std::cout << "closing file \n";
		myfile_.close();
	}

	/*!
	 * \brief Temperature [K] within a finite volume. This problem assumes a temperature of 10 degrees Celsius.
	 *
	 * called EnergyVolumeVariablesImplementation::updateTemperature(...) in porousmediumflow/nonisothermal/volumevariables.hh,
	 * included by porousmediumflow/volumevariables.hh,
	 *
	 * todo this makes very little sense for isothermal!
	 *
	 * overwrites PorousMediumFlowProblem::temperature (compiles without, throws exception of base class)
	 */
	Scalar temperature() const {
		return 273.15 + 10; // -> 10°C
	}

	/*!
	 * \brief Reference pressure [Pa] of the non-wetting. This problem assumes a constant reference pressure of 1 bar.
	 *
	 * called by porousmediumflow/richards/volumevariables.hh
	 */
	Scalar nonWettingReferencePressure() const {
		return pRef_;
	}

	/**
	 * The buffer power for a scv for a volVar (linear in this implementation), equals $\rho_b K_d$ in Eqn (4) in phosphate draft
	 *
	 * used by my the modified localresidual.hh (see dumux-rosi/dumux/porousmediumflow/compositional)
	 */
	Scalar bufferPower(const SubControlVolume& scv, const VolumeVariables& volVars, int compIdx = 0) const {
		return 0;
		int compIdx_ = std::max(std::min(compIdx - soluteIdx,static_cast<int>( b_.size())),0);
		//std::cout<<"bufferPower "<<compIdx<<" "<< compIdx_<<std::endl;//1 or 0
		
		if (b_.at(compIdx_)>0.) {
			//std::cout<<"b_.at(compIdx_); "<<b_.at(compIdx_)<<std::endl;
			return b_.at(compIdx_);
		} else {
			if (freundlichK_.at(compIdx_)==0.) {
				//std::cout<<"(freundlichK_.at(compIdx_)==0.)"<<std::endl;
				return 0.;
			}
			if (volVars.massFraction(h2OIdx, compIdx) <= 0) {
				//std::cout<<"volVars.massFraction(h2OIdx, compIdx) "<<volVars.massFraction(h2OIdx, compIdx)<<std::endl;
				return 0.;
			} else {
                Scalar c = volVars.massFraction(h2OIdx, compIdx)*volVars.density(h2OIdx); // mass concentration
				//std::cout<<"bulk_etc "<<(bulkDensity_*freundlichK_.at(compIdx_)*std::pow(c*1e3, freundlichN_.at(compIdx_))*1e-6/c)<<std::endl;
                return bulkDensity_*freundlichK_.at(compIdx_)*std::pow(c*1e3, freundlichN_.at(compIdx_))*1e-6/c; //kg/m3 = 1e3 mg/l ; mg/kg = 1e-6 kg/kg
			}
		}
	}

	/*!
	 * \copydoc FVProblem::initial
	 *
	 * called by FVProblem::applyInitialSolution(...)
	 */
	template<class Entity>
	PrimaryVariables initial(const Entity& entity) const {
		auto eIdx = this->fvGridGeometry().elementMapper().index(entity);
		Scalar z = entity.geometry().center()[dimWorld - 1];
		PrimaryVariables v(0.0);
		v[pressureIdx] = toPa_(initialSoil_[h2OIdx].f(z,eIdx));
		if(verbose>0)
		{
			 std::cout<<"PrimaryVariables initial(1p5cproblem) "<<z<<" "<<v[pressureIdx];
		}
		// if(toFile)
		// {
			 // myfile_<<"PrimaryVariables initial(1p5cproblem) "<<z<<" "<<v[pressureIdx];
		// }
		for(int i = soluteIdx; i<numComponents_;i++)//solutes
		{
			v[i] = initialSoil_.at(i).f(z,eIdx);
			if(verbose>0)
			{
				 std::cout<<" result : "<<v[i];
			}
			// if(toFile)
			// {
				 // myfile_<<" result : "<<v[i];
			// }
		}
		if(verbose>0)
		{
			 std::cout<<std::endl;
		}
		// if(toFile)
		// {
			 // myfile_<<std::endl;
		// }
		
		return v;
	}

	/*!
	 * @copydoc FVProblem::boundaryTypesAtPos
	 *
	 * discretization dependent, e.g. called by BoxElementBoundaryTypes::boundaryTypes(...)
	 * when?
	 */
	BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const {
		BoundaryTypes bcTypes; ///pressureIdx, bcTopType_, bcSTopType_, ...
		bcTypes.setAllNeumann();
        if (onUpperBoundary_(globalPos)) { // top, bot, or outer bc
            if (bcTopType_ == constantPressure) {
                bcTypes.setDirichlet(pressureIdx);
			 
														  
												  
            }
			for(int i = soluteIdx; i<numComponents_;i++)//solutes
			{
				if (bcSTopType_.at(i - soluteIdx) == constantConcentration) 
				{
					bcTypes.setDirichlet(i);
				}
			}
            
        } else if (onLowerBoundary_(globalPos)) { // top, bot, or outer bc
            if (bcBotType_ == constantPressure) {
                bcTypes.setDirichlet(pressureIdx);//,conti0EqIdx
            }
			for(int i = soluteIdx; i<numComponents_;i++)//solutes
			{
				if (bcSBotType_.at(i - soluteIdx) == constantConcentration) 
				{
					bcTypes.setDirichlet(i);
				}
			}
            
        }
		return bcTypes;
	}

    /*!
     * \copydoc FVProblem::dirichletAtPos
     *
     * dirchlet(...) is called by the local assembler, e.g. BoxLocalAssembler::evalDirichletBoundaries
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const {
		DUNE_THROW(Dune::InvalidStateException, "Do not use Dirichlet BC please");
        PrimaryVariables values;
        if (onUpperBoundary_(globalPos)) 
		{ // top bc
            switch (bcTopType_) 
			{
				case constantPressure: values[pressureIdx] = toPa_(bcTopValues_[pressureIdx]); break;
				default: DUNE_THROW(Dune::InvalidStateException, "Top boundary type Dirichlet: unknown boundary type");
            }
			for(int i = soluteIdx; i<numComponents_;i++)//solutes
			{
				switch (bcSTopType_.at(i - soluteIdx)) 
				{
					case constantConcentration: values[i] = bcSTopValue_.at(i - soluteIdx)*rho_; break;
					default: DUNE_THROW(Dune::InvalidStateException, "Top boundary type Dirichlet: unknown boundary type");
				}
            }
        } else if (onLowerBoundary_(globalPos)) 
		{ // bot bc
            switch (bcBotType_) 
			{
				case constantPressure: values[pressureIdx] = toPa_(bcBotValues_[pressureIdx]); break;
				default: DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Dirichlet: unknown boundary type");
            }
			for(int i = soluteIdx; i<numComponents_;i++)//solutes
			{
				switch (bcSBotType_.at(i - soluteIdx)) 
				{
					case constantConcentration: values[i] = bcSBotValue_.at(i - soluteIdx)*rho_; break;
					default: DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Dirichlet: unknown boundary type");
				}
            }
        }
        // values.setState(Indices::bothPhases); /// <--?????
        return values;
    }

	/*!
	 * \copydoc FVProblem::neumann // [kg/(m²*s)] or [mol/(m²*s)]
	 *
	 * called by BoxLocalResidual::evalFlux,
	 * negative = influx, mass flux in \f$ [ (kg or mol) / (m^2 \cdot s)] \f$// effective diffusion coefficient !!!!!!!!!!!!
	 */
	NumEqVector neumann(const Element& element,
			const FVElementGeometry& fvGeometry,
			const ElementVolumeVariables& elemVolVars,
			const SubControlVolumeFace& scvf) const {
		
		NumEqVector flux;
		GlobalPosition pos = scvf.center();
		auto& volVars = elemVolVars[scvf.insideScvIdx()];
		

		/*
		 *  WATER
		 */
		double f = 0.; // return value [kg m-2 s-1)] or [mol m-2 s-1]
		int pos0 = 1;
		if(dimWorld == 1){pos0 =pos[0]; }
		if ( onUpperBoundary_(pos) || onLowerBoundary_(pos) ) {

			Scalar s = volVars.saturation(h2OIdx);
			Scalar kc = this->spatialParams().hydraulicConductivity(element); //  [m/s]
			// = conductivity AT SATURATION. decrease with wat. content, via krw value 
			//(see material/fluidmatrix interaction /vangenuchten.hh)
			MaterialLawParams params = this->spatialParams().materialLawParams(element);
			Scalar p = MaterialLaw::pc(params, s) + pRef_;//water pressure?
			Scalar h = -toHead_(p); // todo why minus -pc?
			GlobalPosition ePos = element.geometry().center();
			
			Scalar dz = 100 * 2 * std::fabs(ePos[dimWorld - 1] - pos[dimWorld - 1]); // m->cm
			Scalar krw = MaterialLaw::krw(params, s);//	The relative permeability for the wetting phase [between 0 and 1]
			
			//std::cout<<"NumEqVector neumann(), dimworld: "<<dimWorld<<", pos: ";
			//for(int ii = 0; ii < dimWorld; ii++){std::cout<<pos[ii]<<" ";}std::cout<<" ePos ";
			//for(int ii = 0; ii < dimWorld; ii++){std::cout<<ePos[ii]<<" ";}
			//std::cout<<"dz "<<dz<<", bcTopType_ "<<bcTopType_<<" "<<bcSTopType_.at(0)<<" "<<bcSTopType_.at(1);
			//std::cout<<" kc "<<kc<<" krw "<<krw<<std::endl;

			//useMole fraction or mass fraction? 
			//[kg/m3] or [mol/m3]
			Scalar rhoW = useMoles ? volVars.molarDensity(h2OIdx) : volVars.density(h2OIdx) ;
			double cm2m = 1./100. ; //[m/cm]
			double unitConversion = cm2m; //something else needed? 
			if (onUpperBoundary_(pos)) { // top bc
				switch (bcTopType_) {
                case constantPressure: {
                    f = rhoW * kc * ((h - bcTopValues_[pressureIdx]) / dz - gravityOn_)*pos[0] *unitConversion; // maximal inflow
                    //std::cout << "!";
                    break;
                }
				case constantFlux: { // with switch for maximum in- or outflow
					f = -bcTopValues_[pressureIdx]*rhoW/(24.*60.*60.) * unitConversion; // cm/day -> kg/(m²*s) or 
					if (f < 0) { // inflow
						Scalar imax = rhoW * kc * ((h - 0.) / dz - gravityOn_); // maximal inflow
						f = std::max(f, imax)*pos0;
					} else { // outflow
						Scalar omax = rhoW * krw * kc * ((h - criticalPressure_) / dz - gravityOn_); // maximal outflow (evaporation)
						f = std::min(f, omax)*pos0;
					}
					break;
				}
				case constantFluxCyl: { // with switch for maximum in- or outflow
					//useMoles:  [cm/d] * [mol/m^3] * [d/s] * [m/cm] = [m/s] *  [mol/m^3] = [mol /(m^2 * s)]
					//useMass:  [cm/d] * [kg/m^3] * [d/s] * [m/cm] = [m/s] *  [kg/m^3] = [kg /(m^2 * s)]
					f = -bcTopValues_[pressureIdx]*rhoW/(24.*60.*60.) * unitConversion * pos[0];
					if (f < 0) { // inflow
					//[mol/m^3] * [m/s] * [cm/cm] = mol/(m2 * s)
					
					//kg/m3 * m2 
						Scalar imax = rhoW * kc * ((h - 0.) / dz - gravityOn_)*pos[0]; // maximal inflow
						f = std::max(f, imax);
					} else { // outflow
						Scalar omax = rhoW * krw * kc * ((h - criticalPressure_) / dz - gravityOn_)* pos[0]; // maximal outflow (evaporation)
						f = std::min(f, omax);
					}
					break;
				}
				case atmospheric: { // atmospheric boundary condition (with surface run-off) // TODO needs testing & improvement
					Scalar prec = -componentInput_[h2OIdx].f(time_);//-precipitation_.f(time_);
					if (prec < 0) { // precipitation
						Scalar imax = rhoW * kc * ((h - 0.) / dz - gravityOn_); // maximal infiltration
						f = std::max(prec, imax)*pos0;
					} else { // evaporation
						Scalar emax = rhoW * krw * kc * ((h - criticalPressure_) / dz - gravityOn_); // maximal evaporation
						f = std::min(prec, emax)*pos0;
					}
					break;
				}
				default: DUNE_THROW(Dune::InvalidStateException, "Top boundary type Neumann (water) unknown type: "+std::to_string(bcTopType_));
				}
			} else if (onLowerBoundary_(pos)) { // bot bc
				switch (bcBotType_) {
                case constantPressure: {
                    f = rhoW * kc * ((h - bcBotValues_[pressureIdx]) / dz - gravityOn_)* pos[0]; // maximal inflow
//                    Scalar omax = rhoW * krw * kc * ((h - criticalPressure_) / dz - gravityOn_); // maximal outflow (evaporation)
//                    f = std::min(f, omax); rho_
                    break;
                }
				case constantFlux: { // with switch for maximum in- or outflow
					f = -bcBotValues_[pressureIdx]*rhoW/(24.*60.*60.) * unitConversion; // cm/day -> kg/(m²*s) or mol/(m2*s)
					if (f < 0) { // inflow
						Scalar imax = rhoW * kc * ((h - 0.) / dz - gravityOn_); // maximal inflow
						f = std::max(f, imax)*pos0;
					} else { // outflow
						Scalar omax = rhoW * krw * kc * ((h - criticalPressure_) / dz - gravityOn_); // maximal outflow (evaporation)
						f = std::min(f, omax)*pos0;
					}
					break;
				}
				case constantFluxCyl: { // with switch for maximum in- or outflow
					f = -bcBotValues_[pressureIdx]*rhoW/(24.*60.*60.) * unitConversion * pos[0];
					
					//std::cout<<" water flow "<<f<<" "<<bcBotValues_[pressureIdx]<<" rhoW "<<rhoW<<" useMoles "<<useMoles<<" molarDensity "<<volVars.molarDensity(h2OIdx)
					//<<" density "<<volVars.density(h2OIdx)<<" unit conversion "<<unitConversion<<" "<<pos[0]<<std::endl;
					if (f < 0) { // inflow
						Scalar imax = rhoW * kc * ((h - 0.) / dz - gravityOn_)* pos[0]; // maximal inflow
						f = std::max(f, imax);
					} else { // outflow
						Scalar omax = rhoW * krw * kc * ((h - criticalPressure_) / dz - gravityOn_)* pos[0]; // maximal outflow (evaporation)
						// std::cout << " f " << f*1.e9  << ", omax "<< omax << ", value " << bcBotValue_.at(0) << ", crit "  << criticalPressure_ << ", " << pos[0] << "\n";
						f = std::min(f, omax);
					}
					break;
				}
				case freeDrainage: { // holds when useMoles?
					f = krw * kc * rhoW *pos0; // * 1 [m]
					break;
				}
				default: DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Neumann (water) unknown: "+std::to_string(bcBotType_));
				}
			}
			
		}
		
		flux[h2OIdx] = f;

		/*
		 * SOLUTES
		 */
		
		Scalar rhoW = useMoles ? volVars.molarDensity(h2OIdx) : volVars.density(h2OIdx) ;
		double g2kg = 1./1000. ;
		double m2_2_cm2 = 10000;
		double unitConversion = useMoles ? m2_2_cm2 : m2_2_cm2 * g2kg; //something else needed? 
		for(int i = soluteIdx;i<numComponents_;i++)
		{
			int i_s = i - soluteIdx;//for vectors which do not have a value for the H2O primary variable
			Scalar massOrMolFraction = useMoles? volVars.moleFraction(0, i) : volVars.massFraction(0, i);
			if (onUpperBoundary_(pos)) { // top bc Solute
				//std::cout<<"neumann solute, upper BC "<<bcSTopType_.at(i_s)<<" ";
				switch (bcSTopType_.at(i_s)) {
				case constantConcentration: {
					GlobalPosition ePos = element.geometry().center();
					Scalar dz = 2 * std::fabs(ePos[dimWorld - 1] - pos[dimWorld - 1]);
					//!!! att component param set
					static const Scalar d = getParam<Scalar>(std::to_string(i)+".Component.LiquidDiffusionCoefficient"); // m2 / s
					Scalar porosity = this->spatialParams().porosity(element);
					Scalar de = EffectiveDiffusivityModel::effectiveDiffusivity(porosity, volVars.saturation(h2OIdx) ,d);
					flux[i] = (de * (massOrMolFraction*rhoW-bcSTopValue_.at(i_s)*rhoW) / dz + f * massOrMolFraction)*pos0;
					break;

				}
				case constantFlux: {
					//flux[i] = -bcSTopValue_.at(i_s)*rhoW/(24.*60.*60.)/100*pos0; // cm/day -> kg/(m²*s)
					//usemoles:
					// mol/(cm2 * d)  * [d/s]  * cm2/m2 =  mol/(m2 * s)
					// use mass:
					// g/(cm2/d) *[d/s] * [kg/g] * cm2/m2 = kg/(m2 * s)
					flux[i] = -bcSTopValue_.at(i_s)/(24.*60.*60.)* unitConversion; // g/cm2/day || mol/cm2/day  -> kg/(m²*s) || mol/(m²*s)
					break;
				}
				case constantFluxCyl: {
					//usemole:
					//[mol/(cm2 * s)]  = [mol/(cm2 * d)] * [cm2/m2] * d/s
					//flux[i] = -bcSTopValue_.at(i_s)*rhoW/(24.*60.*60.)/100*pos[0]; // cm/day -> kg/(m²*s)
					flux[i] = -bcSTopValue_.at(i_s)/(24.*60.*60.)*unitConversion*pos[0]; // g/cm2/day || mol/cm2/day -> kg/(m²*s) || mol/(m²*s)
					//std::cout<<"constantfluxCyl "<<flux[i];
					break;
				}
				case outflow: {
					// std::cout << "f " << f << ", "  << volVars.massFraction(0, i) << "=" << f*volVars.massFraction(0, i) << "\n";
					flux[i] = f * massOrMolFraction;//*pos0;
					break;
				}
				case outflowCyl: {
					// std::cout << "f " << f << ", "  << volVars.massFraction(0, i) << "=" << f*volVars.massFraction(0, i) << "\n";
					flux[i] = f * massOrMolFraction *pos0;
					break;
				}
				case linear: {
					flux[i] = vMax_.at(i_s)  * massOrMolFraction*pos0;
					break;
				}
				case michaelisMenten: {
					flux[i] = (vMax_.at(i_s)  * (std::max(massOrMolFraction,0.)*rhoW)/(km_.at(i_s)  + std::max(massOrMolFraction,0.)*rhoW))*pos0;
					break;
				}
				case managed: {
					Scalar input = componentInput_.at(i).f(time_);
					flux[i] = input*pos0;
					break;
				}
				default:
					DUNE_THROW(Dune::InvalidStateException, "Top boundary type Neumann (solute) unknown: "+std::to_string(bcSTopType_.at(i_s)));
				}
			} else if (onLowerBoundary_(pos)) { // bot bc Solute
				switch (bcSBotType_.at(i_s)) {
				case constantConcentration: {
					GlobalPosition ePos = element.geometry().center();
					Scalar dz = 2 * std::fabs(ePos[dimWorld - 1] - pos[dimWorld - 1]);
					static const Scalar d = getParam<Scalar>(std::to_string(i)+".Component.LiquidDiffusionCoefficient"); // m2 / s
					Scalar porosity = this->spatialParams().porosity(element);
					Scalar de = EffectiveDiffusivityModel::effectiveDiffusivity(porosity, volVars.saturation(h2OIdx) ,d);
					//diffusion + advection (according to wat flux computed above)
					flux[i] =(de * (massOrMolFraction*rhoW-bcSBotValue_.at(i_s)*rhoW) / dz + f * massOrMolFraction)*pos0;
					//[kg_solute/(m²*s)] = [m2 / s] * ([kg_solute/kg_tot] * [kg_tot / m^3_tot] - kg_solute/ m^3_tot)/m + [kg_tot/(m²*s)] * [kg_solute/kg_tot]
					// std::cout << d*1.e9 << ", "<< de*1.e9 << ", " << volVars.massFraction(0, i) << ", " << bcSBotValue_ << ", " << flux[i]*1.e9  << "\n";
					break;
				}
				case constantFlux: {
					//flux[i] = -bcSBotValue_.at(i_s)*rhoW/(24.*60.*60.)/100*pos0; // cm/day -> kg/(m²*s)
					flux[i] = -bcSBotValue_.at(i_s)/(24.*60.*60.)*unitConversion; // g/cm2/day -> kg/(m²*s)
					//[kg_solute/(m²*s)] = [cm_solute/day] * (kg_tot / m^3_tot)* s/day * m/cm ??
					break;
				}
				case constantFluxCyl: {
					//flux[i] = -bcSBotValue_.at(i_s)*rhoW/(24.*60.*60.)/100*pos[0]; // cm/day -> kg/(m²*s)
					flux[i] = -bcSBotValue_.at(i_s)/(24.*60.*60.)*unitConversion*pos[0]; // g/cm2/day -> kg/(m²*s)
					break;
				}
				case outflow: {//?? free outflow??
					// std::cout << "f " << f*1.e6 << ", "  << volVars.massFraction(0, i) << "=" << f*volVars.massFraction(0, i) << "\n";
					flux[i] = f * massOrMolFraction;
					break;
				}
				case outflowCyl: {
					// std::cout << "f " << f << ", "  << volVars.massFraction(0, i) << "=" << f*volVars.massFraction(0, i) << "\n"; rho_
					flux[i] = f * massOrMolFraction*pos0;
					break;
				}
				case linear: {
					flux[i] = vMax_.at(i_s) * massOrMolFraction*pos0;
					break;
				}
				case michaelisMenten: {
					flux[i] = (vMax_.at(i_s) * (std::max(massOrMolFraction,0.)*rhoW)/(km_.at(i_s) + std::max(massOrMolFraction,0.)*rhoW))*pos0;
					break;
				}
				default: DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Neumann (solute): unknown error");
				}
			} else {
				flux[i] = 0.;
			}
		  
							 
		}

		return flux;
	}

	/*!
	 * \copydoc FVProblem::source
	 *
	 * called by FVLocalResidual:computeSource(...)
	 *
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f
     */
	NumEqVector source(const Element &element, const FVElementGeometry& fvGeometry, const ElementVolumeVariables& elemVolVars,
			const SubControlVolume &scv) const {
		NumEqVector source;
		GlobalPosition pos = scv.center();
  
															  
		auto eIdx = this->spatialParams().fvGridGeometry().elementMapper().index(element);
		for(int i = 0;i < numComponents_;i++)
		{			
			if (source_[i] != nullptr) {
				source[i] = source_[i]->at(eIdx)/scv.volume() * pos[0];
			}else{source[i] = 0.;}												 
		}		
        const auto& volVars = elemVolVars[scv];
		bioChemicalReaction(source, volVars, pos[0]);
		
		return source;
	}
	
		/*!
	 *
     * E.g. for the mol balance that would be a mass rate in \f$ [ mol / (m^3 \cdot s)] \f
     */
	void bioChemicalReaction(NumEqVector &q, const VolumeVariables &volVars, double pos0 ) const
	{
		 //DUNE_THROW(Dune::InvalidStateException, "Biochemical stop");
		//[mol solution / m3 solution]
		//55 mol/kg water * 997 kg/m3
		 const auto massOrMoleDensity = [](const auto& volVars, const int compIdx, const bool isFluid)
        {
			double mOMD = isFluid ? (useMoles ? volVars.molarDensity(compIdx) : volVars.density(compIdx) ):
					(useMoles ? volVars.solidComponentMolarDensity(compIdx) : volVars.solidComponentDensity(compIdx) ); 
			return mOMD;
		};

        const auto massOrMoleFraction= [](const auto& volVars, const int phaseIdx, const int compIdx, const bool isFluid)
        {
			double mOMF = isFluid ?( useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx) ): 
					(useMoles ? volVars.solidMoleFraction(compIdx) : volVars.solidMassFraction(compIdx)); 
			return mOMF;
		};
		 // (mol Soil / m3 soil) * (m3 space/m3 space - m3 pores/m3 space)  = (mol Soil / m3 soil) * (m3 soil/m3 space) = mol soil / m3 space
		double bulkSoilDensity = massOrMoleDensity(volVars, soilIdx -  numFluidComps , false) * (1 - volVars.porosity());
		// if(verbose>0)
		// {
			// std::cout<<"bulk soil "<<bulkSoilDensity <<" "<< massOrMoleDensity(volVars, soilIdx -  numFluidComps , false) <<" "<<(1 - volVars.porosity())<<std::endl;
		// }
		double C_LfrW_temp = massOrMoleFraction(volVars,0, mucilIdx, true);
		double C_LfrW =  std::max(C_LfrW_temp, 0.);					//mol C/mol soil water 
		//[mol solution / m3 solution] * [mol solute / mol solution] = [mol solute / m3 solution]
		double C_L_W = massOrMoleDensity(volVars, h2OIdx, true) * C_LfrW;								//mol C/m3 soil water
		//double C_L = C_L_W * volVars.saturation(h2OIdx) * volVars.porosity();							//mol C/m3 bulk soil
		//double C_LfrSoil =  C_L_W* volVars.saturation(h2OIdx) * volVars.porosity()/bulkSoilDensity;		//mol C/mol solid soil 
		
		double C_SfrW = std::max(massOrMoleFraction(volVars,0, soluteIdx, true), 0.);					//mol C/mol soil water
		double C_S_W = massOrMoleDensity(volVars, h2OIdx, true) * C_SfrW;								//mol C/m3 soil water
		//double C_S = C_S_W * volVars.saturation(h2OIdx) * volVars.porosity();							//mol C/m3 bulk soil
		//[mol /mol soil] = [mol / m3 solution] * [m3 solution/m3 pores] * [m3 pores/m3 space] / [mol soil / m3 space]
		//double C_SfrSoil = C_S_W* volVars.saturation(h2OIdx) * volVars.porosity()/bulkSoilDensity;		//mol C/mol solid soil
		
		std::vector<double> C_afrSoil(2), C_a(2), C_dfrSoil(2), C_d(2);
		for(int i = 0; i < 2; i++)// 0 = oligo, 1 = copio
		{
			int CxAIdx = (i == 0) ? CoAIdx : CcAIdx;
			int CxDIdx = (i == 0) ? CoDIdx : CcDIdx;
			C_afrSoil[i] =  std::max(massOrMoleFraction(volVars,0, CxAIdx - numFluidComps, false), 0.);//mol C/mol solid soil
			C_a[i]  = bulkSoilDensity * C_afrSoil[i] ;														//mol C/m3 bulk soil
			C_dfrSoil[i]  =  std::max(massOrMoleFraction(volVars,0, CxDIdx - numFluidComps, false), 0.);//mol C/mol solid soil
			C_d[i] = bulkSoilDensity * C_dfrSoil[i] ;	// mol C / m3 bulk soil													//mol C/m3 bulk soil
		}
		double C_Oafr = std::max(massOrMoleFraction(volVars,0, CoAIdx - numFluidComps, false), 0.) ; // mol C / mol solid	
		
        //	depolymerisation large polymer to small polymers
		//	[s-1] * ([mol C/m3 bulk soil]/([mol C/m3 bulk soil]*[mol C/m3 bulk soil])) * [mol C/m3 bulk solid]
		// double F_depoly = v_maxL * (C_L/(K_L+ C_L)) * C_a[0]  ; //mol C/(m^3 bulk soil *s)
		//	[s-1] * ([mol C/m3 water]/([mol C/m3 water]*[mol C/m3 water])) * [mol C/m3 bulk solid]
		double F_depoly = v_maxL * (C_L_W/(K_L+ C_L_W)) * C_a[0];// * volVars.saturation(h2OIdx) * volVars.porosity() ; //mol C/(m^3 bulk soil *s)
		
		double F_uptake_S = 0.;
		double F_decay = 0.;	//mol C/(m^3 bulk soil *s)
		double F_growth_S = 0.;	
		
		std::vector<double> F_uptake_S_A(2), F_uptake_S_D(2), F_decay_A(2), F_decay_D(2) ;
		std::vector<double> F_growth(2), F_deact(2), F_react(2), phi(2);
		
		for(int i = 0; i < 2; i++)// 0 = oligo, 1 = copio
		{		
			//				Uptake			
			// ([s-1] * [mol C solute / m3 water] * [m3 water / mol C soil / s])/([s-1] + [mol C solute / m3 water] * [m3 water / mol C soil / s]) * [mol C_oX / m3 space] 
			// [s-1] *([-])/([-] + [-]) * [mol C_oX / m3 space] = [mol C_oX / m3 space /s]
			F_uptake_S_A[i] = (m_max[i] * C_S_W * k_SBis[i])/(m_max[i] + C_S_W * k_SBis[i]) * C_a[i] ;			//mol C/(m^3 bulk soil *s)
			F_uptake_S_D[i] = (m_max[i] * C_S_W * k_SBis[i])/(m_max[i] + C_S_W * k_SBis[i]) * beta[i] * C_d[i] ; //mol C/(m^3 bulk soil *s)
			
			//				Decay
			// [mol C microb / m3 bulk soil /s] = [s-1] * [mol C microb / m3 bulk soil] - [mol C microb / m3 bulk soil /s]
			F_decay_A[i] = m_maxBis[i]  * C_a[i]  - F_uptake_S_A[i] ;			//mol C/(m^3 bulk soil *s)
			F_decay_D[i] = m_maxBis[i]  * beta[i]  * C_d[i]  - F_uptake_S_D[i] ;	//mol C/(m^3 bulk soil *s)
			
			//				Other
			
			// ([s-1] * [mol C solute / m3 water] * [m3 water / mol C / s])/([s-1] + [mol C solute / m3 water] * [m3 water / mol C soil / s]) * [mol C_oX / m3 space] 
			// [s-1] *([-])/([-] + [-]) * [mol C_oX / m3 space] = [mol C_oX / m3 space /s]
			F_growth[i] = (micro_max[i] * C_S_W * k_S[i])/(micro_max[i] + C_S_W * k_S[i]) * C_a[i] ;		//mol C/(m^3 bulk soil *s)
			if(verbose==3)//||(massOrMoleFraction(volVars,0, mucilIdx, true)<0.))
			{
				std::cout<<"F_growth["<<i<<"] " << std::scientific<<std::setprecision(20)
				<<micro_max[i] <<" "<< C_S_W <<" "<< k_S[i] <<" "<< C_a[i]<<" "<<F_growth[i]<<std::endl;
				std::cout<<"F_decay_A["<<i<<"] " << std::scientific<<std::setprecision(20)
				<<m_maxBis[i] <<" "<< C_a[i]  <<" "<< F_uptake_S_A[i]  <<" "
				<< F_decay_A[i] <<std::endl;
			}
			phi[i] = 1/(1 + std::exp((C_S_W_thres[i] - C_S_W)/(k_phi * C_S_W_thres[i])));								// - 
			// [-] * [1/s] * [mol C/m3 bulk soil]
			F_deact[i]  = (1 - phi[i] ) * k_D[i]  * C_a[i] ;			//mol C/(m^3 bulk soil *s)
			F_react[i]  = phi[i]  * k_R[i]  * C_d[i] ;				//mol C/(m^3 bulk soil *s)
			
			F_uptake_S += F_uptake_S_A[i] + F_uptake_S_D[i] ;	//mol C/(m^3 bulk soil *s)
			F_decay += F_decay_A[i] + F_decay_D[i];	
			F_growth_S += (1/k_growth[i])*F_growth[i];
		}
		
		
		//Att: using here absolute saturation. should we use the effective? should we multiply by pos?
		//[mol solute / m3 space/s] 
		
		q[soluteIdx] += (  + F_depoly + (1 - k_decay2)*F_decay - F_uptake_S -  F_growth_S)* pos0 ;//* volVars.saturation(h2OIdx) * volVars.porosity();// /FluidSystem::molarMass(soluteIdx)
		q[mucilIdx]  += (-F_depoly +  k_decay2 * F_decay) * pos0;// * volVars.saturation(h2OIdx) * volVars.porosity();// /FluidSystem::molarMass(mucilIdx) 
		
		q[CoAIdx] += (  - extra + F_growth[0] - F_deact[0] + F_react[0] - (1/k_decay)*F_decay_A[0]) * pos0;
		q[CoDIdx] += ( - extra2 + F_deact[0] - F_react[0] - (1/k_decay)*F_decay_D[0]) * pos0;
		
		q[CcAIdx] += (   F_growth[1] - F_deact[1] + F_react[1] - (1/k_decay)*F_decay_A[1]) * pos0;
		q[CcDIdx] += (F_deact[1] - F_react[1] - (1/k_decay)*F_decay_D[1]) * pos0;
		
	}

	/*!
	 * \brief Applies a vector of point sources. The point sources
	 *        are possibly solution dependent.
	 *
	 * \param pointSources A vector of Dumux::PointSource s that contain
              source values for all phases and space positions.
	 *
	 * For this method, the \a values method of the point source
	 * has to return the absolute mass rate in kg/s. Positive values mean
	 * that mass is created, negative ones mean that it vanishes.
	 */
	template<class PointSource>
	void addPointSources(std::vector<PointSource>& pointSources) const {
		if (couplingManager_!=nullptr) {
			pointSources = couplingManager_->bulkPointSources();
		}
	}

	/*!
	 * \brief Evaluate the point sources (added by addPointSources) TODO
	 *        for all phases within a given sub-control-volume.
	 *
	 * This is the method for the case where the point source is
	 * solution dependent and requires some quantities that
	 * are specific to the fully-implicit method.
	 *
	 * For this method, the \a values() method of the point sources returns
	 * the absolute rate mass generated or annihilate in kg/s. Positive values mean
	 * that mass is created, negative ones mean that it vanishes.
	 */
	 //not used??
	template<class ElementVolumeVariables>
	void pointSource(PointSource& source,
			const Element &element,
			const FVElementGeometry& fvGeometry,
			const ElementVolumeVariables& elemVolVars,
			const SubControlVolume &scv) const {

		PrimaryVariables sourceValue(0.);
		std::cout<<"template<class ElementVolumeVariables>"<<std::endl;
		DUNE_THROW(Dune::InvalidStateException, "template<class ElementVolumeVariables>");

	}

	/*!
	 * Sets the current simulation time (within the simulation loop) for atmospheric look up [s]
	 *
	 * eventually, called in the main file (example specific, richards.cc)
	 */
	void setTime(Scalar t, Scalar dt) {
		time_ = t;
		dt_ = dt; // currently unused
	}

	/*!
	 * Source per element index \f$ [ kg / s)] \f$
	 *
	 * eventually, called in the main file (example specific, richards.cc)
	 */
	void setSource(std::shared_ptr<std::vector<double>> s, int eqIdx = 0) {
		source_[eqIdx] = s;
	}

	//! Set the coupling manager
	void setCouplingManager(CouplingManager* cm) {
		couplingManager_ = cm;
	}

	/*!
	 * sets the critical pressure for evaporation [cm] (default = -10000 cm)
	 *
	 *  eventually, called in the main file (example specific, richards.cc)
	 */
	void criticalPressure(Scalar p) {
		criticalPressure_ = p;
	}

	/**
	 * Sets boundary fluxes according to the last solution TODO
	 */
	void postTimeStep(const SolutionVector& sol, const GridVariables& gridVars) {
		bc_flux_upper = 0.;
		bc_flux_lower = 0.;
		int uc = 0;
		int lc = 0;
		for (const auto& e :elements(this->fvGridGeometry().gridView())) {
			auto fvGeometry = localView(this->fvGridGeometry());
			fvGeometry.bindElement(e);
			auto elemVolVars = localView(gridVars.curGridVolVars());
			elemVolVars.bindElement(e, fvGeometry, sol);
			for (const auto& scvf :scvfs(fvGeometry)) { // evaluate root collar sub control faces
				auto p = scvf.center();
				if (onUpperBoundary_(p)) { // top
				//std::cout<<"postTimeStep upper: "<<onUpperBoundary_(p)<<" "<<uc<<std::endl;
					bc_flux_upper += neumann(e, fvGeometry, elemVolVars, scvf);
					uc++;
				} else if (onLowerBoundary_(p)) { // bottom
				//std::cout<<"postTimeStep lower: "<<onLowerBoundary_(p)<<" "<<lc<<std::endl;
					bc_flux_lower += neumann(e, fvGeometry, elemVolVars, scvf);
					lc++;
				}
			}
		}
		bc_flux_upper /= uc;
		bc_flux_lower /= lc;
	}

	/*!
	 * Writes the actual boundary fluxes (top and bottom) into a text file. Call postTimeStep before using it.
	 */
	void writeBoundaryFluxes() {
		myfile_ << time_ << ", " << bc_flux_upper << ", " << bc_flux_lower << "\n";
	}

	/**
	 * debug info TODO make meaningful for 2c
	 */
	void computeSourceIntegral(const SolutionVector& sol, const GridVariables& gridVars) const {
		NumEqVector source(0.0);
		for (const auto& element : elements(this->fvGridGeometry().gridView())) {
			auto fvGeometry = localView(this->fvGridGeometry());
			fvGeometry.bindElement(element);
			auto elemVolVars = localView(gridVars.curGridVolVars());
			elemVolVars.bindElement(element, fvGeometry, sol);
			for (auto&& scv : scvs(fvGeometry)) {
				auto pointSources = this->scvPointSources(element, fvGeometry, elemVolVars, scv);
				pointSources *= scv.volume()*elemVolVars[scv].extrusionFactor();
				source += pointSources;
			}
		}
		std::cout << "Global integrated source (soil): " << source[h2OIdx] << " (kg/s) / " << source[h2OIdx]*3600*24*1000 << " (g/day)" << '\n';
	}

	/**
	 * Forwards to spatialParams
	 *
     * Call to change default setting (of 1.e-6 for both)
     *
     * pcEps    capillary pressure regularisation
     * krEps 	relative permeabiltiy regularisation
     */
    void setRegularisation(double pcEps, double krEps) {
    	this->spatialParams().setRegularisation(pcEps,krEps);
    }


	// BC, direct access for Python binding (setTopBC, setBotBC, in richards.hh)
	int bcTopType_;
	int bcBotType_;
	std::vector<double> bcTopValues_ = std::vector<double>(1);
	std::vector<double> bcBotValues_ = std::vector<double>(1);

	std::vector<int> bcSTopType_ = std::vector<int>(numSolutes); // one solute
	std::vector<int> bcSBotType_ = std::vector<int>(numSolutes);
	std::vector<double> bcSTopValue_ = std::vector<double>(numSolutes);
	std::vector<double> bcSBotValue_= std::vector<double>(numSolutes);
	

private:

	//! cm pressure head -> Pascal
	Scalar toPa_(Scalar ph) const {
		Scalar ph2 = pRef_ + ph / 100. * rho_ * g_;
		//std::cout<<"toPa_ "<<ph<<" "<<ph2<<" "<<rho_<<" "<<g_<<std::endl;
		return ph2;
	}

	//! Pascal -> cm pressure head
	Scalar toHead_(Scalar p) const {
		Scalar p2 = (p - pRef_) * 100. / rho_ / g_;
		//std::cout<<"toHead_ "<<p<<" "<<p2<<" "<<rho_<<" "<<g_<<std::endl;
		return p2;
	}

	//! true if on the point lies on the left boundary
	bool onLeftBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[0] > this->fvGridGeometry().bBoxMin()[0] - eps_;
	}
	//! true if on the point lies on the right boundary
	bool onRightBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_;
	}
	
	//! true if on the point lies on the upper boundary
	bool onUpperBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[dimWorld - 1] > this->fvGridGeometry().bBoxMax()[dimWorld - 1] - eps_;
	}

	//! true if on the point lies on the upper boundary
	bool onLowerBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[dimWorld - 1] < this->fvGridGeometry().bBoxMin()[dimWorld - 1] + eps_;
	}

	// Initial
	std::vector<InputFileFunction> initialSoil_ = std::vector<InputFileFunction>(numComponents_);
								 

	bool gravityOn_;

	// Source
	std::vector<std::shared_ptr<std::vector<double>>> source_; // [kg/s]
	CouplingManager* couplingManager_ = nullptr;

	//InputFileFunction precipitation_;
	std::vector<InputFileFunction> componentInput_ = std::vector<InputFileFunction>(numComponents_);
	Scalar criticalPressure_; // cm
	Scalar time_ = 0.;
	Scalar dt_ = 0.;

	std::ofstream myfile_;
	int verbose;
	bool toFile;
	//not used for anything except post-simulation analysis I think
	NumEqVector bc_flux_upper = NumEqVector(0.);
	NumEqVector bc_flux_lower = NumEqVector(0.);

	static constexpr Scalar eps_ = 1.e-7;
	
	// static constexpr Scalar g_ = 9.81; // cm / s^2 (for type conversions)
	// static constexpr Scalar rho_ = 1.e3; // kg / m^3 (for type conversions)
	// static constexpr Scalar pRef_ = 1.e5; // Pa
	static constexpr Scalar g_ = 9.80665; // cm / s^2 (for type conversions)
	//of pure water and of solution (as low solute content)
	static constexpr Scalar rho_ = 1.e3; // kg / m^3 (for type conversions)
	static constexpr Scalar pRef_ = 101300; // Pa

	// Uptake params
	std::vector<Scalar> vMax_ = std::vector<Scalar>(numSolutes); // Michaelis Menten Parameter [kg m-2 s-1]
	std::vector<Scalar> km_ = std::vector<Scalar>(numSolutes);  // Michaelis Menten Parameter  [kg m-3]
	std::vector<Scalar> sigma_ = std::vector<Scalar>(numSolutes);// 1 for passive transport, 0 for active transport

	std::vector<Scalar> b_ = std::vector<Scalar>(numSolutes); // buffer power
	std::vector<Scalar> freundlichK_ = std::vector<Scalar>(numSolutes); // Freundlich parameters
	std::vector<Scalar> freundlichN_ = std::vector<Scalar>(numSolutes);
	Scalar bulkDensity_ = 1.4;//g/cm3 // TODO check with Mai, buffer power (1+b) or b
	
	// default value in CPB units
	double  v_maxL_ = 5e5; //Maximum reaction rate of enzymes targeting large polymers [d-1]
	double  K_L_ = 10e-3 ; //Half-saturation coefficients of enzymes targeting large polymers [mol cm-3 soil]
	std::vector<double>  m_max_{0.01,0.001}; //[d-1] Maximum maintenance rate coefficient for the corresponding microbial groups
	std::vector<double>  micro_max_{1,10}; //[d-1] Maximum growth rate coefficient for the corresponding microbial group
	std::vector<double>  beta_{0.1,0.3}; //[-] Reduction factor of maintenance requirements in dormant state for the corresponding microbial group
	std::vector<double>  k_S_{ 50000, 30000}; // [mol soil / mol C soil / d] Specific substrate affinity to small polymers for the corresponding microbial group
	std::vector<double>  C_S_W_thres_{ 0.001/(12*1000), 0.01/(12*1000)}; //[mol C/cm3 soil water] Threshold concentration for reactivation and deactivation for the corresponding microbial groups
	double  k_phi_= 0.1; //[-] Sharpness parameter for the switch function from active to dormancy
	std::vector<double>  k_D_{ 1, 5}; //[d-1] Deactivation rate coefficient for the corresponding microbial group
	std::vector<double>  k_R_{ 1, 5}; //[d-1] Reactivation rate coefficient for the corresponding microbial group
	double  k_decay_= 0.75; //[-] Maintenance yield
	double  k_decay2_= 0.5;//[-] Proportion of large polymers formed from dead microbial biomass due to maintenance
	std::vector<double>  k_growth_{ 0.8, 0.6};//[-] Growth yield on small polymers for the corresponding microbial groups
	
	//pagel (2020): the average total initial microbial biomass was 1.67 × 10−4 mg g−1 (C soil−1) * *bulkDensity_ g/cm3
	double  v_maxL ; //Maximum reaction rate of enzymes targeting large polymers [s-1]
	double  K_L  ; //Half-saturation coefficients of enzymes targeting large polymers [kg m-3 soil] or [mol m-3 soil] 
	std::vector<double>  m_max ; //[s-1] Maximum maintenance rate coefficient for the corresponding microbial groups
	std::vector<double>  m_maxBis ; //[s-1] Maximum maintenance rate coefficient for the corresponding microbial groups
	std::vector<double>  micro_max ; //[s-1] Maximum growth rate coefficient for the corresponding microbial group
	std::vector<double>  beta ; //[-] Reduction factor of maintenance requirements in dormant state for the corresponding microbial group
	std::vector<double>  k_S ; // [mol soil / mol C soil / s] Specific substrate affinity to small polymers for the corresponding microbial group
	std::vector<double>  k_SBis ; // [mol soil / mol C soil / s] Specific substrate affinity to small polymers for the corresponding microbial group
	std::vector<double>  C_S_W_thres ; //[mol C/m3 soil water] Threshold concentration for reactivation and deactivation for the corresponding microbial groups
	double  k_phi ; //[-] Sharpness parameter for the switch function from active to dormancy
	std::vector<double>  k_D ; //[s-1] Deactivation rate coefficient for the corresponding microbial groups
	std::vector<double>  k_R ; //[s-1] Reactivation rate coefficient for the corresponding microbial group
	double  k_decay ; //[-] Maintenance yield
	double  k_decay2;//[-] Proportion of large polymers formed from dead microbial biomass due to maintenance
	double  k_decay3 = 1;// for troubleshooting
	std::vector<double>  k_growth;//[-] Growth yield on small polymers for the corresponding microbial groups
	double k_growthBis = 1; // dummy variable for troubleshooting
	double extra = 0;
	double extra2 = 0;
	double m_maxBisO;
	double m_maxBis_Cs;
	double theta;
	//from Magdalena:  have just rechecked all the solute units by looking if the mass of exuded C equals 
	//the mass of C in the soil domain during the simulation and realized that the unit of s.getSolution_(EqIdx)  
	//must be g/cm^3 (you already mentioned that this was not clear).

};

} //end namespace Dumux

#endif