// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef RICHARDS_PROBLEM_SCHROEDER_HH
#define RICHARDS_PROBLEM_SCHROEDER_HH

#include <dumux/porousmediumflow/problem.hh> // base class

#include "richardsparams.hh"
#include <dumux/external/brent/brent.hpp>                  //T.S.: Brent algorithm to find roots of function



namespace Dumux {

/*!
 * RichardsProblem:
 * Uses Dumux as an easy to use Richards equation solver,
 * where most parameters can be set dynamically
 */
template <class TypeTag>
class RichardsProblem : public PorousMediumFlowProblem<TypeTag>
{
public:

	// exports, used by the binding
	using Grid = GetPropType<TypeTag, Properties::Grid>;
	using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
	using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
	using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

	// other
	using GridView = typename FVGridGeometry::GridView;
	using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
	using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
	using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
	using FVElementGeometry = typename FVGridGeometry::LocalView;
	using SubControlVolume = typename FVGridGeometry::SubControlVolume;
	using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
	using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
	using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
	using Scalar = GetPropType<TypeTag, Properties::Scalar>;
	using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
	using Element = typename GridView::template Codim<0>::Entity;
	using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
	using MaterialLaw = typename GetPropType<TypeTag, Properties::SpatialParams>::MaterialLaw;
	using MaterialLawParams = typename MaterialLaw::Params;

	using PointSource = GetPropType<TypeTag, Properties::PointSource>;
	using CouplingManager= GetPropType<TypeTag, Properties::CouplingManager>;

	enum {
		// copy some indices for convenience
		pressureIdx = Indices::pressureIdx,
		conti0EqIdx = Indices::conti0EqIdx,
		bothPhases = Indices::bothPhases,
		// world dimension
		dimWorld = GridView::dimensionworld,
		// discretization method
		isBox = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethods::Box
	};

	enum BCTypes {
		constantPressure = 1,
		constantFlux = 2,
		constantFluxCyl = 3,
		atmospheric = 4,
		freeDrainage = 5
	};

	enum GridParameterIndex {
		materialLayerNumber = 0
	};

	/*!
	 * \brief Constructor: constructed in the main file
	 */
	RichardsProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
	: PorousMediumFlowProblem<TypeTag>(fvGridGeometry) {

		gravityOn_ = Dumux::getParam<bool>("Problem.EnableGravity", true);

		// BC
		bcTopType_ = getParam<int>("Soil.BC.Top.Type");
		bcBotType_ = getParam<int>("Soil.BC.Bot.Type");
		bcTopValue_ = getParam<Scalar>("Soil.BC.Top.Value",0.);
		bcBotValue_ = getParam<Scalar>("Soil.BC.Bot.Value",0.);

		criticalPressure_ = getParam<double>("Soil.CriticalPressure", -1.e4); // cm
		criticalPressure_ = getParam<double>("Climate.CriticalPressure", criticalPressure_); // cm
		// precipitation
		if (bcTopType_==atmospheric) {
			precipitation_ = InputFileFunction("Climate", "Precipitation", "Time", 0.); // cm/day (day)
			precipitation_.setVariableScale(1./(24.*60.*60.)); // s -> day
			precipitation_.setFunctionScale(rho_/(24.*60.*60.)/100); // cm/day -> kg/(m²*s)
		}
		// IC
		initialSoil_ = InputFileFunction("Soil.IC", "P", "Z", 0., this->spatialParams().layerIFF()); // [cm]([m]) pressure head, conversions hard coded
		// Output
		std::string filestr = this->name() + ".csv"; // output file
		writeFile_ = getParam<bool>("Soil.Output.File", true);
		if (writeFile_) {
			myfile_.open(filestr.c_str());
		}
		std::cout << "RichardsProblem constructed: bcTopType " << bcTopType_ << ", " << bcTopValue_ << "; bcBotType "
				<<  bcBotType_ << ", " << bcBotValue_ << ",  Output File " << writeFile_
				<< ", Critical pressure " << criticalPressure_ << " gravity " << gravityOn_ << "\n" << std::flush;
	}

	/**
	 * \brief Eventually, closes output file
	 */
	~RichardsProblem() {
		if (writeFile_) {
			std::cout << "closing file \n";
			myfile_.close();
		}
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

	/*!
	 * \copydoc FVProblem::initial
	 *
	 * called by FVProblem::applyInitialSolution(...)
	 */
	template<class Entity>
	PrimaryVariables initial(const Entity& entity) const {
		auto eIdx = this->gridGeometry().elementMapper().index(entity);
		Scalar z = entity.geometry().center()[dimWorld - 1];
		// std::cout << "initial " << z << ", " << initialSoil_.f(z,eIdx) << " \n";
		PrimaryVariables v(0.0);
		v[pressureIdx] = toPa_(initialSoil_.f(z,eIdx));
		v.setState(bothPhases);
		return v;
	}

	/*!
	 * @copydoc FVProblem::boundaryTypesAtPos
	 *
	 * discretization dependent, e.g. called by BoxElementBoundaryTypes::boundaryTypes(...)
	 * when?
	 */
	BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const {
		BoundaryTypes bcTypes;
		if (onUpperBoundary_(globalPos)) { // top or outer bc
			switch (bcTopType_) {
			case constantPressure: bcTypes.setAllDirichlet(); break;
			case constantFlux: bcTypes.setAllNeumann(); break;
			case constantFluxCyl: bcTypes.setAllNeumann(); break;
			case atmospheric: bcTypes.setAllNeumann(); break;
			default: DUNE_THROW(Dune::InvalidStateException,"Top boundary type not implemented");
			}
		} else if (onLowerBoundary_(globalPos)) { // bot or inner bc
			switch (bcBotType_) {
			case constantPressure: bcTypes.setAllDirichlet(); break;
			case constantFlux: bcTypes.setAllNeumann(); break;
			case constantFluxCyl: bcTypes.setAllNeumann(); break;
			case freeDrainage: bcTypes.setAllNeumann(); break;
			default: DUNE_THROW(Dune::InvalidStateException,"Bottom boundary type not implemented");
			}
		} else {
			bcTypes.setAllNeumann(); // no top not bottom is no flux
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
		if (onUpperBoundary_(globalPos)) { // top bc
			switch (bcTopType_) {
			case constantPressure: values[Indices::pressureIdx] = toPa_(bcTopValue_); break;
			default: DUNE_THROW(Dune::InvalidStateException, "Top boundary type Dirichlet: unknown boundary type");
			}
		} else if (onLowerBoundary_(globalPos)) { // bot bc
			switch (bcBotType_) {
			case constantPressure: values[Indices::pressureIdx] = toPa_(bcBotValue_); break;
			default: DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Dirichlet: unknown boundary type");
			}
		}
		values.setState(Indices::bothPhases);
		return values;
	}

	/*!
	 * \copydoc FVProblem::neumann // [kg/(m²*s)]
	 *
	 * called by BoxLocalResidual::evalFlux,
	 * mass flux in \f$ [ kg / (m^2 \cdot s)] \f$
	 * Negative values mean influx.
	 */
	NumEqVector neumann(const Element& element,
			const FVElementGeometry& fvGeometry,
			const ElementVolumeVariables& elemVolVars,
			const SubControlVolumeFace& scvf) const {

		NumEqVector flux;
		double f = 0.; // return value
		GlobalPosition pos = scvf.center();

		if ( onUpperBoundary_(pos) || onLowerBoundary_(pos) ) {

			Scalar s = elemVolVars[scvf.insideScvIdx()].saturation(0);
			Scalar kc = this->spatialParams().hydraulicConductivity(element); //  [m/s]
			MaterialLawParams params = this->spatialParams().materialLawParams(element);
			Scalar p = MaterialLaw::pc(params, s) + pRef_; // [Pa]
			Scalar h = -toHead_(p); // cm
			GlobalPosition ePos = element.geometry().center();
			Scalar dz = 100 * 2 * std::fabs(ePos[dimWorld - 1] - pos[dimWorld - 1]); // cm
			Scalar krw = MaterialLaw::krw(params, s);

			if (onUpperBoundary_(pos)) { // top bc
				switch (bcTopType_) {
				case constantFlux: { // with switch for maximum in- or outflow
					f = -bcTopValue_*rho_/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
					if (f < 0) { // inflow
						Scalar imax = rho_ * kc * ((h - 0.) / dz - gravityOn_); // maximal inflow
						//std::cout << "in:" << f <<", " << imax <<"\n";
						f = std::max(f, imax);
					} else { // outflow
						Scalar omax = rho_ * krw * kc * ((h - criticalPressure_) / dz - gravityOn_); // maximal outflow (evaporation)
						f = std::min(f, omax);
					}
					break;
				}
				case constantFluxCyl: { // with switch for maximum in- or outflow
					f = -bcTopValue_*rho_/(24.*60.*60.)/100 * pos[0];  // cm/day -> kg/(m²*s) (Eqns are multiplied by cylindrical radius)
					if (f < 0) { // inflow
						Scalar imax = rho_ * kc * ((h - 0.) / dz - gravityOn_)* pos[0]; // maximal inflow
						f = std::max(f, imax);
					} else { // outflow
						Scalar omax = rho_ * krw * kc * ((h - criticalPressure_) / dz - gravityOn_)* pos[0]; // maximal outflow (evaporation)
						f = std::min(f, omax);
					}
					break;
				}
				case atmospheric: { // atmospheric boundary condition (with surface run-off)
					Scalar prec = -precipitation_.f(time_);
					if (prec < 0) { // precipitation
						// std::cout << "in" << "\n";
						Scalar imax = rho_ * kc * ((h - 0.) / dz - gravityOn_); // maximal infiltration
						f = std::max(prec, imax);
					} else { // evaporation
						// std::cout << "out" << "\n";
						Scalar emax = rho_ * krw * kc * ((h - criticalPressure_) / dz - gravityOn_); // maximal evaporation
						f = std::min(prec, emax);
					}
					break;
				}
				default: DUNE_THROW(Dune::InvalidStateException, "Top boundary type Neumann: unknown error");
				}
			} else if (onLowerBoundary_(pos)) { // bot bc
				switch (bcBotType_) {
				case constantFlux: { // with switch for maximum in- or outflow
					f = -bcBotValue_*rho_/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
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
					f = -bcBotValue_*rho_/(24.*60.*60.)/100 * pos[0];
					if (f < 0) { // inflow
						Scalar imax = rho_ * kc * ((h - 0.) / dz - gravityOn_)* pos[0]; // maximal inflow
						f = std::max(f, imax);
					} else { // outflow
						Scalar omax = rho_ * krw * kc * ((h - criticalPressure_) / dz - gravityOn_)* pos[0]; // maximal outflow (evaporation)
						// std::cout << " f " << f*1.e9 << ", omax "<< omax << ", value " << bcBotValue_ << ", crit "  << criticalPressure_ << ", " << pos[0] << "\n";
						f = std::min(f, omax);
					}
					break;
				}
				case freeDrainage: {
					f = krw * kc * rho_; // * 1 [m]
					break;
				}
				default: DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Neumann: unknown error");
				}
			}
		}

		flux[conti0EqIdx] = f;
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
		if ((source_ != nullptr)) {
			auto eIdx = this->spatialParams().gridGeometry().elementMapper().index(element);
			return source_->at(eIdx)/scv.volume();
		} else {
			return 0.;
		}
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
	 * \brief Evaluate the point sources (added by addPointSources)
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

    // T.S: Function definition: integration by hand (calculate matric-flux-potential based on the currenct absolute pressure)
    const Scalar pc_to_MFP(const Element &element, Scalar lower,  Scalar pressure3D_pc, int n) const
    {
        MaterialLawParams params = this->spatialParams().materialLawParams(element);
        Scalar kc = this->spatialParams().hydraulicConductivity(element); // [m/s]
        std::function<double(double)> f = [=] (double x) { return MaterialLaw::krw(params, MaterialLaw::sw(params, x))*kc*86400; };
        return CPlantBox::Function::quad(f, pressure3D_pc, lower, n);

//        Scalar cumSum =0;
//        for (int i=0; i<n+1; i++)
//        {
//            Scalar xi = pressure3D_pc +i*dx; // pc value for sw call
//            Scalar funValue = MaterialLaw::sw(params, xi);
//            Scalar funValue2 = MaterialLaw::krw(params, funValue);
//            Scalar rectangleArea = funValue2*dx*kc; // height * base length
//            cumSum += rectangleArea*86400;
//            // CHECK UNITS! MFP from cm²/s into cm²/day for comparison with python script, assumption was that MFP should be m²/day, results indicate otherwise
//        }
//        return cumSum;
    }

	template<class ElementVolumeVariables>
	void pointSource(PointSource& source,
			const Element &element,
			const FVElementGeometry& fvGeometry,
			const ElementVolumeVariables& elemVolVars,
			const SubControlVolume &scv) const {
		if (couplingManager_!=nullptr) {
			// compute source at every integration point
			const Scalar pressure3D = couplingManager_->bulkPriVars(source.id())[Indices::pressureIdx];
			const Scalar pressure1D = couplingManager_->lowDimPriVars(source.id())[Indices::pressureIdx];
			const auto& spatialParams = couplingManager_->problem(Dune::index_constant<1>{}).spatialParams();
			const auto lowDimElementIdx = couplingManager_->pointSourceData(source.id()).lowDimElementIdx();
			const Scalar kr = spatialParams.kr(lowDimElementIdx);                                    // T.S: [m/Pa/s]
			const Scalar rootRadius = spatialParams.radius(lowDimElementIdx); // T.S: [m]
			// relative soil permeability
			const auto krel = 1.0;
			// sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
			const auto density = 1000;

            Scalar sourceValue = 2 * M_PI *krel*rootRadius * kr *(pressure1D - pressure3D)*density; //T.S: [kg / s m], deleted const
			source = sourceValue*source.quadratureWeight()*source.integrationElement();

            /// SCHROEDER IMPLEMENTATION
            const Scalar gradients = getParam<Scalar>("Schroeder.gradients");
            // Switch in Input-File (gradients = 1 in [Soil.IC] enables Schroeder, gradients = 0 disables it)

            if (gradients == 1 && sourceValue < 0) {
            // switch to enable/disable Schroeder and check for macroscopic flow from soil to root (Schroeder only used if source value < 0)

                // STEP 1) CALCULATE MFP_SOIL
                // Integration of soil hydraulic conductivity K(h) from -15.000 cm to current pressure head of soil element

                const Scalar pressure3D_pc = -pressure3D + pRef_;
                // upper integration boundary hc, soil point-source pressure [pc]
                const Scalar lowBound = -15000;
                const Scalar lowBound_pc = -toPa_(lowBound) +pRef_ ;
                // lower integration boundary (-15.000 cm) for MFP [pc]
                const int n = getParam<int>("Schroeder.n");
                // integration-steps (10000 gives good results for clay & loam, 40000 needed for sand & still not perfect)
                // hydraulic conductivity of soil voxel
                const Scalar MFP_soil = pc_to_MFP(element, lowBound_pc, pressure3D_pc, n);
                // MFP of source-point soil voxel, call to integration function


                // STEP 2) CALCULATE MFP_ROOT (according to non-stressed equation of Schroeder)

                //calculation of r_out
                const Scalar segment_length = source.integrationElement(); // [m]
                // length of point-source segment in voxel (cut at voxel boundaries)
                const Scalar rootsystem_volume_inElement = couplingManager_->lowDimVolume(element);
                // total volume of rootsystem in voxel
                const Scalar cell_volume = element.geometry().volume();
                // volume of soil voxel
                const Scalar segment_volume = (M_PI * rootRadius * rootRadius) * segment_length;
                // volume of single segment
                const Scalar t = segment_volume / rootsystem_volume_inElement;
                // proportionality factor
                const Scalar targetV = t * cell_volume;
                // targetVolume of segment
                const Scalar r_out = sqrt((targetV + segment_volume) / (M_PI * segment_length)) * 100; // [cm]
                // r_out = radius of bulk soil cylinder assigned to a segment
                const Scalar rho = r_out / (rootRadius*100); // [cm/cm]

                //flux densities at outer and inner boundary
                const Scalar q_out = 0;
                // flux at bulk soil cylinder with radius r_out, for now we assume no-flux
                const Scalar q_root = (-1*sourceValue) / 1000 * 1000000 / 100 *86400 /(2 * M_PI *rootRadius*100);
                // flux at root cylinder with radius r_root. sourceValue Conversion from Dumux units [kg s⁻¹ m⁻1] into Schroeder units [cm³ cm⁻² d⁻¹] <=> [cm/day]
                // this has to be *-1, MFP_soil < MFP_root

                // r, radial coordinate for MFP calculation. For us always r = root, double definition for readibility of MFP_nostress_root equation
                const Scalar r = rootRadius*100; // [cm]
                const Scalar r_root = rootRadius*100; // [cm]
                //MFP at root_surface according to non-stressed equation of Schroeder et al. 2008 (equation 4) [cm²/d]
                const Scalar MFP_nostress_root = MFP_soil + (q_root * r_root - q_out *r_out) * (pow(r,2) / pow(r_root, 2) / (2*(1-pow(rho,2)))
                + pow(rho,2) / (1-pow(rho,2)) * (log(r_out/r) -0.5)) + q_out * r_out * log(r / r_out);

                // STEP 3) TRANSFER MFP at root-surface back to a pressure value
                // parameters for Brent algorithm (finds zero of a function in a bracketing interval)
                double tolerance = brent::r8_epsilon ( );
                // error-tolerance parameter of brent-algorithm
                auto MFP_to_pressure3D = brent::funcLambda([=](double x) { return  pc_to_MFP(element, lowBound_pc, x, n) - MFP_nostress_root; });
                double z = brent::zero (lowBound_pc, 0, tolerance, MFP_to_pressure3D);
                const Scalar pressure3D_pc_new = z;
                const Scalar pressure3D_new = -1*(z-pRef_);

                // STEP 4) PASS NEW PRESSURE3D TO SOURCE-TERM
                MaterialLawParams params = this->spatialParams().materialLawParams(element);
                const Scalar pressure3D_s_new = MaterialLaw::sw(params, pressure3D_pc_new);
                const Scalar krw = MaterialLaw::krw(params, pressure3D_s_new); // pass new pressure3D_s to calculate krw
                const Scalar krw_scaled_to_rootRadius = krw / rootRadius;      // soil hydraulic conductivity scaled to root radius (radius used as proxy)
                const Scalar kmin = std::min(kr, krw_scaled_to_rootRadius);    // minimum of conductivity defines sourceValue
                sourceValue = 2 * M_PI *krel* rootRadius * kmin * (pressure1D - pressure3D_new)*density;    //* kr exchanged with kmin
                source = sourceValue*source.quadratureWeight()*source.integrationElement();

                if(sourceValue > 0) {
                    //discards Schroeder if it leads to inversion of flow (e.g. macroscopic-flow soil => root, schroeder-flow root => soil. Jan thinks this clause may be exluded)
                    source = 0;
                    }

                //Prints (can be enabled / disabled via Schroeder print in coupled input-file)
                const Scalar print = getParam<Scalar>("Schroeder.print");
                if (print == 1) {
                    std::cout << " soilproblem sourceID_soil:" << source.id() << "\n" << " MFP_soil = " << MFP_soil <<  "   MFP_no_stress_root= " << MFP_nostress_root << "   deltaMFP= "
                    << MFP_soil-MFP_nostress_root << "  soil_r_out= " << r_out << "   root_radius= " << rootRadius*100 << "   q_root= " << q_root << "\n"
                    << " pressure3D_head(soil)= " << toHead_(pressure3D) << "   pressure3D_head(root_surface)= " << toHead_(pressure3D_new) << "\n"
                    << " pressure3D_Pa(soil)  = " << pressure3D << "  pressure3D_Pa(root_surface) = " << pressure3D_new << "\n"
                    << " pressure1D_Pa(root)  = " << pressure1D << "  pressure1D_head(root)       = " << toHead_(pressure1D) << "\n" << "\n";
                    }
                }
		     }

// default again
   else {
			source = 0;
		}
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
	 * Source per element index \f$ [ kg / s ] \f$
	 *
	 * eventually, called in the main file (example specific, richards.cc)
	 */
	void setSource(std::shared_ptr<std::vector<double>> s, int eqIdx = 0) {
		source_ = s;
	}

	/*!
	 * sets the critical pressure for evaporation [cm] (default = -10000 cm)
	 *
	 *  eventually, called in the main file (example specific, richards.cc)
	 */
	void criticalPressure(Scalar p) {
		criticalPressure_ = p;
	}

	//! Set the coupling manager
	void setCouplingManager(CouplingManager* cm) {
		couplingManager_ = cm;
	}
	/**
	 * Sets boundary fluxes according to the last solution
	 */
	void postTimeStep(const SolutionVector& sol, const GridVariables& gridVars) {
		bc_flux_upper = 0.;
		int uc = 0;
		bc_flux_lower = 0.;
		int lc = 0;
		for (const auto& e :elements(this->gridGeometry().gridView())) {
			auto fvGeometry = localView(this->gridGeometry());
			fvGeometry.bindElement(e);
			auto elemVolVars = localView(gridVars.curGridVolVars());
			elemVolVars.bindElement(e, fvGeometry, sol);
			for (const auto& scvf :scvfs(fvGeometry)) { // evaluate root collar sub control faces
				try { // will throw an exception, if boundary type is Dirichlet
					auto p = scvf.center();
					if (onUpperBoundary_(p)) { // top
						bc_flux_upper += neumann(e, fvGeometry, elemVolVars, scvf);
						uc++;
					} else if (onLowerBoundary_(p)) { // bottom
						bc_flux_lower += neumann(e, fvGeometry, elemVolVars, scvf);
						lc++;
					}
				} catch (...) {
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
		if (writeFile_) {
			myfile_ << time_ << ", " << bc_flux_upper << ", " << bc_flux_lower << "\n";
		}
	}

	/**
	 * debug info
	 */
	void computeSourceIntegral(const SolutionVector& sol, const GridVariables& gridVars) const {
		NumEqVector source(0.0);
		for (const auto& element : elements(this->gridGeometry().gridView())) {
			auto fvGeometry = localView(this->gridGeometry());
			fvGeometry.bindElement(element);
			auto elemVolVars = localView(gridVars.curGridVolVars());
			elemVolVars.bindElement(element, fvGeometry, sol);
			for (auto&& scv : scvs(fvGeometry)) {
				auto pointSources = this->scvPointSources(element, fvGeometry, elemVolVars, scv);
				pointSources *= scv.volume()*elemVolVars[scv].extrusionFactor();
				source += pointSources;
			}
		}
		std::cout << "Global integrated source (soil): " << source << " (kg/s) / " << source*3600*24*1000 << " (g/day)" << '\n';
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

	// BC, direct access for Python binding (setBcTop, setBcBot)
	int bcTopType_;
	int bcBotType_;
	double bcTopValue_;
	double bcBotValue_;

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
		return globalPos[dimWorld - 1] > this->gridGeometry().bBoxMax()[dimWorld - 1] - eps_;
	}

	//! true if on the point lies on the upper boundary
	bool onLowerBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[dimWorld - 1] < this->gridGeometry().bBoxMin()[dimWorld - 1] + eps_;
	}

	// Initial
	InputFileFunction initialSoil_;
	bool gravityOn_;

	// Source
	std::shared_ptr<std::vector<double>> source_; // [kg/s]
	CouplingManager* couplingManager_ = nullptr;

	InputFileFunction precipitation_;
	Scalar criticalPressure_; // cm
	Scalar time_ = 0.;
	Scalar dt_ = 0.;

	bool writeFile_ = true;
	std::ofstream myfile_;
	Scalar bc_flux_upper = 0.;
	Scalar bc_flux_lower = 0.;

	static constexpr Scalar eps_ = 1.e-7;
	static constexpr Scalar g_ = 9.81; // cm / s^2 (for type conversions)
	static constexpr Scalar rho_ = 1.e3; // kg / m^3 (for type conversions)
	static constexpr Scalar pRef_ = 1.e5; // Pa

};

} //end namespace Dumux

#endif
