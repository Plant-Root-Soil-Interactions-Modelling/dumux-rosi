// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef RICHARDS1P3C_PROBLEM_HH
#define RICHARDS1P3C_PROBLEM_HH
#include <algorithm>
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
class Richards1P3CProblem : public PorousMediumFlowProblem<TypeTag>
{
public:

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
	using EffectiveDiffusivityModel = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;

	enum {
		pressureIdx = 0, // index of primary variables
		h2OIdx = pressureIdx, // fluid index
		soluteIdx = 1, // 1st solute index
		mucilIdx = 2, // mucil index
		
		//don t think EqIndx != pvdIdx when we have just 1 phase
		conti0EqIdx = pressureIdx, // indices of the equations
		transportEqIdx = conti0EqIdx + soluteIdx,
		//transportMucilEqIdx = conti0EqIdx + mucilIdx,//is this correct?

		dimWorld = GridView::dimensionworld,
		
		//!!!!!!!!!!!!
		numComponents_ = 3 ,//Edit this to have N components instead of 3
		//!!!!!!!!!!!!
		
		numSolutes = numComponents_ -  soluteIdx,


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
		managed = 9
	};

	enum GridParameterIndex {
		materialLayerNumber = 0
	};

	/*!
	 * \brief Constructor: constructed in the main file
	 */
	Richards1P3CProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
	: PorousMediumFlowProblem<TypeTag>(fvGridGeometry) {

		gravityOn_ = Dumux::getParam<bool>("Problem.EnableGravity", true);

		source_.resize(numComponents_); // numComponents_ equations (currently hard coded, where can I get the value?)
		// Components
		for(int i = 0; i < numComponents_; i++)//all components except h2o
		{
			
			source_[i] = nullptr;
			if(i ==h2OIdx)
			{
				
				// BC			
				bcTopType_ = getParam<int>("Soil.BC.Top.Type", outflow); 
				bcBotType_ = getParam<int>("Soil.BC.Bot.Type", outflow);
				
				const auto& myParams = Parameters::paramTree();
				myParams.report();
				bcTopValues_.at(i) =  getParam<double>("Soil.BC.Top.Value", 0.);
				bcBotValues_.at(i) =  getParam<double>("Soil.BC.Bot.Value", 0.);
				
				//IC
				std::cout<<"initialSoil_"<<std::endl;
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
				
				//IC
				initialSoil_.at(i) = InputFileFunction("Soil.IC", "C"+std::to_string(i), "C"+std::to_string(i)+"Z", 
													0., this->spatialParams().layerIFF()); // kg/m2
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
			
			
			componentInput_.at(i).setVariableScale(1./(24.*60.*60.)); // s -> day
			componentInput_.at(i).setFunctionScale(10./(24.*60.*60.)); // g/(cm2 day) -> kg/(m²*s)
			
		}

		criticalPressure_ = getParam<double>("Soil.CriticalPressure", -1.e4); // cm
		criticalPressure_ = getParam<double>("Climate.CriticalPressure", criticalPressure_); // cm
		
		
		 v_maxL = getParam<Scalar>("Soil.v_maxL", 0.1); //Maximum reaction rate of enzymes targeting large polymers [d-1]
		 K_L = getParam<Scalar>("Soil.K_L", 10e-6 ); //Half-saturation coefficients of enzymes targeting large polymers [g cm-3 soil]
		 C_Oa = getParam<Scalar>("Soil.C_oa", 2.338e-06 ); //concentration of active oligotrophic biomass [g (C)cm-3(soil-1)]

												  
		// Output
		std::string filestr = this->name() + "_1p3cProblem.csv"; // output file
		myfile_.open(filestr.c_str());
		std::cout << "Richards1P3CProblem constructed: bcTopType " << bcTopType_ << ", " << bcTopValues_.at(0) << "; bcBotType "
				<<  bcBotType_ << ", " << bcBotValues_.at(0) 
				<< ", gravitation " << gravityOn_ <<", Critical pressure "
				<< criticalPressure_ << "\n" << "Sorption:" << "buffer power "<< b_[0] << ", Freundlich " << freundlichK_[0] << ", " <<
				freundlichN_[0] << "\n" 
				<< std::flush;
		
	}

	/**
	 * \brief Eventually, closes output file
	 */
	~Richards1P3CProblem() {
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
		for(int i = soluteIdx; i<numComponents_;i++)//solutes
		{
			v[i] = initialSoil_.at(i).f(z,eIdx);
		}
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
	 * \copydoc FVProblem::neumann // [kg/(m²*s)]
	 *
	 * called by BoxLocalResidual::evalFlux,
	 * negative = influx, mass flux in \f$ [ kg / (m^2 \cdot s)] \f$// effective diffusion coefficient !!!!!!!!!!!!
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
		double f = 0.; // return value [kg m-2 s-1)]
		if ( onUpperBoundary_(pos) || onLowerBoundary_(pos) ) {

			Scalar s = volVars.saturation(h2OIdx);
			Scalar kc = this->spatialParams().hydraulicConductivity(element); //  [m/s]
			MaterialLawParams params = this->spatialParams().materialLawParams(element);
			Scalar p = MaterialLaw::pc(params, s) + pRef_;//water pressure?
			Scalar h = -toHead_(p); // todo why minus -pc?
			//std::cout<<"NumEqVector neumann() "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<std::endl;
			GlobalPosition ePos = element.geometry().center();
			//std::cout<<"ePos "<<dimWorld<<" "<<ePos[0]<<" "<<ePos[1]<<" "<<ePos[2]<<std::endl;
			Scalar dz = 100 * 2 * std::fabs(ePos[dimWorld - 1] - pos[dimWorld - 1]); // m->cm
			Scalar krw = MaterialLaw::krw(params, s);//	The relative permeability for the wetting phase

			if (onUpperBoundary_(pos)) { // top bc
				switch (bcTopType_) {
                case constantPressure: {
                    f = rho_ * kc * ((h - bcTopValues_[pressureIdx]) / dz - gravityOn_)*pos[0]; // maximal inflow
                    //std::cout << "!";
                    break;
                }
				case constantFlux: { // with switch for maximum in- or outflow
					f = -bcTopValues_[pressureIdx]*rho_/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
					if (f < 0) { // inflow
						Scalar imax = rho_ * kc * ((h - 0.) / dz - gravityOn_); // maximal inflow
						f = std::max(f, imax);
					} else { // outflow
						Scalar omax = rho_ * krw * kc * ((h - criticalPressure_) / dz - gravityOn_); // maximal outflow (evaporation)
						f = std::min(f, omax);
					}
					break;
				}
				case constantFluxCyl: { // with switch for maximum in- or outflow
					f = -bcTopValues_[pressureIdx]*rho_/(24.*60.*60.)/100 * pos[0];
					if (f < 0) { // inflow
						Scalar imax = rho_ * kc * ((h - 0.) / dz - gravityOn_)*pos[0]; // maximal inflow
						f = std::max(f, imax);
					} else { // outflow
						Scalar omax = rho_ * krw * kc * ((h - criticalPressure_) / dz - gravityOn_)* pos[0]; // maximal outflow (evaporation)
						f = std::min(f, omax);
					}
					break;
				}
				case atmospheric: { // atmospheric boundary condition (with surface run-off) // TODO needs testing & improvement
					Scalar prec = -componentInput_[h2OIdx].f(time_);//-precipitation_.f(time_);
					if (prec < 0) { // precipitation
						Scalar imax = rho_ * kc * ((h - 0.) / dz - gravityOn_); // maximal infiltration
						f = std::max(prec, imax);
					} else { // evaporation
						Scalar emax = rho_ * krw * kc * ((h - criticalPressure_) / dz - gravityOn_); // maximal evaporation
						f = std::min(prec, emax);
					}
					break;
				}
				default: DUNE_THROW(Dune::InvalidStateException, "Top boundary type Neumann (water) unknown type: "+std::to_string(bcTopType_));
				}
			} else if (onLowerBoundary_(pos)) { // bot bc
				switch (bcBotType_) {
                case constantPressure: {
                    f = rho_ * kc * ((h - bcBotValues_[pressureIdx]) / dz - gravityOn_)* pos[0]; // maximal inflow
//                    Scalar omax = rho_ * krw * kc * ((h - criticalPressure_) / dz - gravityOn_); // maximal outflow (evaporation)
//                    f = std::min(f, omax);
                    break;
                }
				case constantFlux: { // with switch for maximum in- or outflow
					f = -bcBotValues_[pressureIdx]*rho_/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
					if (f < 0) { // inflow
						Scalar imax = rho_ * kc * ((h - 0.) / dz - gravityOn_); // maximal inflow
						f = std::max(f, imax);
					} else { // outflow
						Scalar omax = rho_ * krw * kc * ((h - criticalPressure_) / dz - gravityOn_); // maximal outflow (evaporation)
						f = std::min(f, omax);
					}
					break;
				}
				case constantFluxCyl: { // with switch for maximum in- or outflow
					f = -bcBotValues_[pressureIdx]*rho_/(24.*60.*60.)/100 * pos[0];
					if (f < 0) { // inflow
						Scalar imax = rho_ * kc * ((h - 0.) / dz - gravityOn_)* pos[0]; // maximal inflow
						f = std::max(f, imax);
					} else { // outflow
						Scalar omax = rho_ * krw * kc * ((h - criticalPressure_) / dz - gravityOn_)* pos[0]; // maximal outflow (evaporation)
						// std::cout << " f " << f*1.e9  << ", omax "<< omax << ", value " << bcBotValue_.at(0) << ", crit "  << criticalPressure_ << ", " << pos[0] << "\n";
						f = std::min(f, omax);
					}
					break;
				}
				case freeDrainage: {
					f = krw * kc * rho_; // * 1 [m]
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
		
		for(int i = soluteIdx;i<numComponents_;i++)
		{
			int i_s = i - soluteIdx;//for vectors which do not have a value for the H2O primary variable
			if (onUpperBoundary_(pos)) { // top bc Solute
				switch (bcSTopType_.at(i_s)) {
				case constantConcentration: {
					GlobalPosition ePos = element.geometry().center();
					Scalar dz = 2 * std::fabs(ePos[dimWorld - 1] - pos[dimWorld - 1]);
					//!!! att component param set
					static const Scalar d = getParam<Scalar>(std::to_string(i)+".Component.LiquidDiffusionCoefficient"); // m2 / s
					Scalar porosity = this->spatialParams().porosity(element);
					Scalar de = EffectiveDiffusivityModel::effectiveDiffusivity(porosity, volVars.saturation(h2OIdx) ,d);
					flux[i] = de * (volVars.massFraction(0, i)*rho_-bcSTopValue_.at(i_s)*rho_) / dz + f * volVars.massFraction(0, i);
					break;

				}
				case constantFlux: {
					flux[i] = -bcSTopValue_.at(i_s)*rho_/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
					break;
				}
				case constantFluxCyl: {
					flux[i] = -bcSTopValue_.at(i_s)*rho_/(24.*60.*60.)/100*pos[0]; // cm/day -> kg/(m²*s)
					break;
				}
				case outflow: {
					// std::cout << "f " << f << ", "  << volVars.massFraction(0, i) << "=" << f*volVars.massFraction(0, i) << "\n";
					flux[i] = f * volVars.massFraction(0, i);
					break;
				}
				case linear: {
					flux[i] = vMax_.at(i_s)  * volVars.massFraction(0, i);
					break;
				}
				case michaelisMenten: {
					flux[i] = vMax_.at(i_s)  * (std::max(volVars.massFraction(0, i),0.)*rho_)/(km_.at(i_s)  + std::max(volVars.massFraction(0, i),0.)*rho_);
					break;
				}
				case managed: {
					Scalar input = componentInput_.at(i).f(time_);
					flux[i] = input;
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
					flux[i] =de * (volVars.massFraction(0, i)*rho_-bcSBotValue_.at(i_s)*rho_) / dz + f * volVars.massFraction(0, i);
					//[kg_solute/(m²*s)] = [m2 / s] * ([kg_solute/kg_tot] * [kg_tot / m^3_tot] - kg_solute/ m^3_tot)/m + [kg_tot/(m²*s)] * [kg_solute/kg_tot]
					// std::cout << d*1.e9 << ", "<< de*1.e9 << ", " << volVars.massFraction(0, i) << ", " << bcSBotValue_ << ", " << flux[i]*1.e9  << "\n";
					break;
				}
				case constantFlux: {
					flux[i] = -bcSBotValue_.at(i_s)*rho_/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
					//[kg_solute/(m²*s)] = [cm_solute/day] * (kg_tot / m^3_tot)* s/day * m/cm ??
					break;
				}
				case constantFluxCyl: {
					flux[i] = -bcSBotValue_.at(i_s)*rho_/(24.*60.*60.)/100*pos[0]; // cm/day -> kg/(m²*s)
					break;
				}
				case outflow: {//?? free outflow??
					// std::cout << "f " << f*1.e6 << ", "  << volVars.massFraction(0, i) << "=" << f*volVars.massFraction(0, i) << "\n";
					flux[i] = f * volVars.massFraction(0, i);
					break;
				}
				case linear: {
					flux[i] = vMax_.at(i_s) * volVars.massFraction(0, i);
					break;
				}
				case michaelisMenten: {
					flux[i] = vMax_.at(i_s) * (std::max(volVars.massFraction(0, i),0.)*rho_)/(km_.at(i_s) + std::max(volVars.massFraction(0, i),0.)*rho_);
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
  
															  
		auto eIdx = this->spatialParams().fvGridGeometry().elementMapper().index(element);
		for(int i = 0;i < numComponents_;i++)
		{			
			if (source_[i] != nullptr) {
				source[i] = source_[i]->at(eIdx)/scv.volume();
			}else{source[i] = 0.;}												 
		}		
        const auto& volVars = elemVolVars[scv];
		bioChemicalReaction(source, volVars);
		
		return source;
	}
	
	
	void bioChemicalReaction(NumEqVector &q, const VolumeVariables &volVars) const
	{
		//depolymerisation large polymer to small polymers
		Scalar C_L = volVars.moleFraction(h2OIdx, mucilIdx);//g/g I think
		Scalar F_depoly = v_maxL * (C_L /(K_L+C_L)) * C_Oa ;
		q[soluteIdx] += F_depoly;
		q[mucilIdx] -= F_depoly;
		
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
		return pRef_ + ph / 100. * rho_ * g_;
	}

	//! Pascal -> cm pressure head
	Scalar toHead_(Scalar p) const {
		return (p - pRef_) * 100. / rho_ / g_;
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
	//not used for anything except post-simulation analysis I think
	NumEqVector bc_flux_upper = NumEqVector(0.);
	NumEqVector bc_flux_lower = NumEqVector(0.);

	static constexpr Scalar eps_ = 1.e-7;
	static constexpr Scalar g_ = 9.81; // cm / s^2 (for type conversions)
	//of pure water and of solution (as low solute content)
	static constexpr Scalar rho_ = 1.e3; // kg / m^3 (for type conversions)
	static constexpr Scalar pRef_ = 1.e5; // Pa

	// Uptake params
	std::vector<Scalar> vMax_ = std::vector<Scalar>(numSolutes); // Michaelis Menten Parameter [kg m-2 s-1]
	std::vector<Scalar> km_ = std::vector<Scalar>(numSolutes);  // Michaelis Menten Parameter  [kg m-3]
	std::vector<Scalar> sigma_ = std::vector<Scalar>(numSolutes);// 1 for passive transport, 0 for active transport

	std::vector<Scalar> b_ = std::vector<Scalar>(numSolutes); // buffer power
	std::vector<Scalar> freundlichK_ = std::vector<Scalar>(numSolutes); // Freundlich parameters
	std::vector<Scalar> freundlichN_ = std::vector<Scalar>(numSolutes);
	Scalar bulkDensity_ = 1.4;//g/cm3 // TODO check with Mai, buffer power (1+b) or b
	
	Scalar v_maxL = 0.1; //Maximum reaction rate of enzymes targeting large polymers [d-1]
	Scalar K_L = 10e-6 ; //Half-saturation coefficients of enzymes targeting large polymers [g cm-3 soil]
	Scalar C_Oa = 2.338e-06 ; //concentration of active oligotrophic biomass [g (C)cm-3(soil-1)]
	//pagel (2020): the average total initial microbial biomass was 1.67 × 10−4 mg g−1 (C soil−1) => 1.67*10-7 *bulkDensity_ g/cm3
	
	//from Magdalena:  have just rechecked all the solute units by looking if the mass of exuded C equals 
	//the mass of C in the soil domain during the simulation and realized that the unit of s.getSolution_(EqIdx)  
	//must be g/cm^3 (you already mentioned that this was not clear).

};

} //end namespace Dumux

#endif
