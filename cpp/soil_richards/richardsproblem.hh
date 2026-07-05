// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:

#ifndef RICHARDS_PROBLEM_HH
#define RICHARDS_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/method.hh>

#include "../soil_richards/richardsparams.hh"

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/io/inputfilefunction.hh>

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

    using ParentType = PorousMediumFlowProblem<TypeTag>;

	// exports, used by the binding
	using Grid = GetPropType<TypeTag, Properties::Grid>;
	using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
	using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
	using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
	using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;

	// other
	using GridView = typename GridGeometry::GridView;
	using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
	using NumEqVector = typename Dumux::NumEqVector<PrimaryVariables>;
	using FVElementGeometry = typename GridGeometry::LocalView;
	using SubControlVolume = typename GridGeometry::SubControlVolume;
	using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
	using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
	//using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

	using Scalar = GetPropType<TypeTag, Properties::Scalar>;
	using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
	using Element = typename GridView::template Codim<0>::Entity;
	using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

	using PcKrSwCurve = Dumux::FluidMatrix::VanGenuchtenDefault<Scalar>;
	// VanGenuchtenNoReg<double>  // both of type TwoPMaterialLaw // using MaterialLaw = typename GetPropType<TypeTag, Properties::SpatialParams>::MaterialLaw;
	using BasicParams = typename PcKrSwCurve::BasicParams;
    using EffToAbsParams = typename PcKrSwCurve::EffToAbsParams;
    using RegularizationParams = typename PcKrSwCurve::RegularizationParams;

	using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;

	using PointSource = GetPropType<TypeTag, Properties::PointSource>;
	// using CouplingManager= GetPropType<TypeTag, Properties::CouplingManager>;


	enum { isBox = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethods::box };
    static constexpr bool useMoles = true;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static constexpr int numFluidComps = ModelTraits::numFluidComponents();

	enum {
		// copy some indices for convenience
		pressureIdx = Indices::pressureIdx,
		conti0EqIdx = Indices::conti0EqIdx,
		h2OIdx = FluidSystem::liquidPhaseIdx, // fluid index
		// world dimension
		dimWorld = GridView::dimensionworld
		// discretization method
	};
	
	enum BCTypes {
		constantPressure = 1,
		constantFlux = 2,
		constantFluxCyl = 3, // for cylindrical models
		atmospheric = 4,
		freeDrainage = 5,
		rootSystemExact = 6 // for cylindrical models with coupling to root system
	};

	enum GridParameterIndex {
		materialLayerNumber = 0
	};

	/*!
	 * \brief Constructor: constructed in the main file
	 */
	RichardsProblem(std::shared_ptr<const GridGeometry> gridGeometry)
	: ParentType(gridGeometry) {

		gravityOn_ = Dumux::getParam<bool>("Problem.EnableGravity", true);
		temperatureK = Dumux::getParam<double>("Problem.temperatureK", 273.15+10);
		 verbose =  getParam<int>("Problem.verbose", 0);

		// BC
		bcTopType_ = getParam<int>("Soil.BC.Top.Type");
		bcBotType_ = getParam<int>("Soil.BC.Bot.Type");
		bcTopValues_.push_back(getParam<Scalar>("Soil.BC.Top.Value",0.));
		bcBotValues_.push_back(getParam<Scalar>("Soil.BC.Bot.Value",0.));

		criticalPressure_ = getParam<double>("Soil.CriticalPressure", -1.e4); // cm
		criticalPressure_ = getParam<double>("Climate.CriticalPressure", criticalPressure_); // cm
		sourceSlope_ = getParam<double>("Soil.SourceSlope", -1.); // cm, negative value disables regularisation
		// precipitation
		if (bcTopType_==atmospheric) {
			precipitation_ = InputFileFunction("Climate", "Precipitation", "Time", 0.); // cm/day (day)
			precipitation_.setVariableScale(1./(24.*60.*60.)); // s -> day for time
			if(useMoles){
				precipitation_.setFunctionScale(rho_molar/(24.*60.*60.)/100); // cm/day) -> mol/(m²*s)
			} else {
				precipitation_.setFunctionScale(rho_/(24.*60.*60.)/100); // cm/day -> kg/(m²*s)
			}
		}
		// IC
		initialSoil_ = InputFileFunction("Soil.IC", "P", "Z", 0., this->spatialParams().layerIFF()); // [cm]([m]) pressure head, conversions hard coded
		// Output
		std::string filestr = this->name() + ".csv"; // output file
		writeFile_ = getParam<bool>("Soil.Output.File", true);
		if (writeFile_) {
			myfile_.open(filestr.c_str());
		}
		if(verbose){
			std::cout << "RichardsProblem constructed: bcTopType " << bcTopType_ << ", " << bcTopValues_.at(0) << "; bcBotType "
					<<  bcBotType_ << ", " << bcBotValues_.at(0) << ",  Output File " << writeFile_
					<< ", Critical pressure " << criticalPressure_ << " gravity " << gravityOn_ << "\n" << std::flush;
		}
		computedCellVolumesCyl = false;
		if(dimWorld == 1)
		{
			segLength  = Dumux::getParam<double>("Problem.segLength")/100;
			this->setVolumesCyl(computeCellVolumesCyl_());
			computedCellVolumesCyl = true;
		}else{
			segLength  = -1.;
		}
	
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

	double getCellVolumesCyl(int dofIndex) const
	{
		
		if(!computedCellVolumesCyl) // if we want to change the segLEngth without recreating the object?
		{
			DUNE_THROW(Dune::InvalidStateException, "getCellVolumesCyl: cellVolumesCyl not set");
		}
		return cellVolumesCyl.at(dofIndex);
	}

    void setVolumesCyl(std::vector<double> cellvoles) const
	{
		const_cast<std::vector<double>&>(cellVolumesCyl) = cellvoles;
    }
    /**
     * The volume [m3] of each element (vtk cell)
     *
     * This is done for a single process, gathering and mapping is done in Python.
     */
    std::vector<double> computeCellVolumesCyl_() const
	{
        std::vector<double> vols;
		auto points = this->getPoints_();//get the vertices == faces of 1D domain
        for (int i = 0; i < (points.size()-1); i++) {
			double rIn = points.at(i).at(0);
			double rOut = points.at(i + 1).at(0);
            vols.push_back((rOut*rOut - rIn*rIn)* M_PI * segLength);
        }
        return vols;
    }

    /**
     * Returns the Dune vertices (vtk points) of the grid for a single mpi process.
     * Gathering and mapping is done in Python.
     */
    std::vector<std::array<double,1>> getPoints_() const
	{
		//auto eIdx = this->gridGeometry().elementMapper().index(entity);
		//Scalar z = entity.geometry().center()[dimWorld - 1];
		
        std::vector<std::array<double,1>> points;
        points.reserve(this->gridGeometry().gridView().size(dimWorld));
        for (const auto& v : vertices(this->gridGeometry().gridView())) {
            auto p = v.geometry().center();
            std::array<double,1> vp;
            for (int i=0; i<dimWorld; i++) { // found no better way
                vp[i] = p[i];
            }
            points.push_back(vp);
        }
        return points;
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
		return temperatureK; //273.15 + 10; // -> 10°C
	}
	
	/*!
	 * \brief Reference pressure [Pa] of the non-wetting. This problem assumes a constant reference pressure of 1 bar.
	 *
	 * called by porousmediumflow/richards/volumevariables.hh
	 */
	Scalar nonWettingReferencePressure() const {
		return pRef_;
	}
	Scalar nonwettingReferencePressure() const {
		return nonWettingReferencePressure();
	}

	Scalar massOrMoleDensity(const auto& volVars, const int compIdx, const bool isFluid) const
	{
		return isFluid ? (useMoles ? rho_molar : volVars.density(compIdx) ):
				(useMoles ? volVars.solidComponentMolarDensity(compIdx) : volVars.solidComponentDensity(compIdx) ); 
	}

	Scalar massOrMoleFraction(const auto& volVars, const int phaseIdx, const int compIdx, const bool isFluid) const
	{
		return isFluid ?( useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx) ): 
				(useMoles ? volVars.solidMoleFraction(compIdx) : volVars.solidMassFraction(compIdx)); 
	}
	
	/**
	 * The buffer power for a scv for a volVar (linear in this implementation), equals $\rho_b K_d$ in Eqn (4) in phosphate draft
	 *
	 * used by my the modified localresidual.hh (see dumux-rosi/dumux/porousmediumflow/compositional)
	 */
	Scalar bufferPower(const VolumeVariables& volVars, int compIdx, int dofIndex) const {
		
		switch(compIdx)
		{
			case h2OIdx:{
				DUNE_THROW(Dune::InvalidStateException, "bufferPower used for water");
				break;
			}
			default:
			{
				return 0.;
			}
		}
		
	}
	
	double getFreundlichK_(const VolumeVariables& volVars, int dofIndex) const 
	{
		return 0.;
	}
	
	double getB_(const VolumeVariables& volVars, int dofIndex) const 
	{
		return 0.;
	}

	double computeB(const VolumeVariables& volVars, double C_S_W, int dofIndex) const
	{
		return 0.;
	}
	
	double computeCSS1(double bulkSoilDensity, double C_S_W, int dofIndex) const
	{// mol		
	
		return 0.;
	}

    std::function<double(double,double)> computeDtCSS2 =
        std::bind(&RichardsProblem::computeDtCSS2_, this, std::placeholders::_1, std::placeholders::_2); 
        
    double computeInitCSS2_(double CSS1, double CSW) // mol/m3
    {
        return 0.;
    }
    double computeDtCSS2_(double CSW, double CSS2) // mol/m3/s
    {
        return 0.;
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
		// v.setState(pressureIdx); // todo ???
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
			case rootSystemExact: bcTypes.setAllNeumann(); break;
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
			case constantPressure: values[Indices::pressureIdx] = toPa_(bcTopValues_[0]); break;
			default: DUNE_THROW(Dune::InvalidStateException, "Top boundary type Dirichlet: unknown boundary type");
			}
		} else if (onLowerBoundary_(globalPos)) { // bot bc
			switch (bcBotType_) {
			case constantPressure: values[Indices::pressureIdx] = toPa_(bcBotValues_[0]); break;
			default: DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Dirichlet: unknown boundary type");
			}
		}
		// values.setState(pressureIdx); // unused???
		return values;
	}


	PcKrSwCurve materialLaw(const Element& element) const {
	    const BasicParams& basicParams = this->spatialParams().basicParams(element);
	    const EffToAbsParams& effToAbsParams = this->spatialParams().effToAbsParams(element);
	    const RegularizationParams& regularizationParams = this->spatialParams().regularizationParams(element);
	    return PcKrSwCurve(basicParams, effToAbsParams, regularizationParams);
	}

	/*!
	 * \copydoc FVProblem::neumann // [mol/(m²*s)]
	 *
	 * called by BoxLocalResidual::evalFlux,
	 * mass flux in \f$ [ mol / (m^2 \cdot s)] \f$
	 * Negative values mean influx.
	 */
	NumEqVector neumann(const Element& element,
			const FVElementGeometry& fvGeometry,
			const ElementVolumeVariables& elemVolVars,
			const ElementFluxVariablesCache&  fluxCache,
			const SubControlVolumeFace& scvf) const {

		NumEqVector flux;
		double f = 0.; // return value
		GlobalPosition pos = scvf.center();
		double pos0 = 1;
		if((dimWorld == 1)&&(!this->spatialParams().useExtrusion)){
			pos0 = pos[0]; 
		}
		auto& volVars = elemVolVars[scvf.insideScvIdx()];

		if ( onUpperBoundary_(pos) || onLowerBoundary_(pos) ) {

			// kg / m^3  or mol / m3
			Scalar rhoW = useMoles ? rho_molar : rho_ ;
			
			//Scalar s = volVars.saturation(0);
			Scalar kc = this->spatialParams().hydraulicConductivity(element); //  [m/s] // volVars.permeability() * rho_ * g_/volVars.viscosity(h2OIdx);//

			//PcKrSwCurve materialLaw_ = materialLaw(element);
			//Scalar p = materialLaw_.pc(s) + pRef_; // [Pa]
			Scalar h = volVars.pressureHead(0);//-toHead_(p); // cm
			GlobalPosition ePos = element.geometry().center();
			Scalar dz = 100 * std::fabs(ePos[dimWorld - 1] - pos[dimWorld - 1]); // m-> cm (*2 ?)
			Scalar krw = volVars.relativePermeability();//Scalar krw = materialLaw_.krw(s); // [1]

			if (onUpperBoundary_(pos)) { // top bc
				switch (bcTopType_) {
				case constantFlux: { // with switch for maximum in- or outflow
					f = -bcTopValues_[0]*rhoW/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
					if (f < 0.) { // inflow
						Scalar imax = rhoW * kc * ((h - 0.) / dz - gravityOn_); // maximal inflow
//						 std::cout << "in:" << f <<", " << imax <<"\n";
						f = std::max(f, imax);
					} else { // outflow
						Scalar omax = rhoW *  kc * krw * ((h - criticalPressure_) / dz - gravityOn_); // maximal outflow (evaporation)
						// std::cout << "outflow " << f*1.e6 << ", " << omax*1.e6 << " krw " << krw*1.e6 << "\n";
						f = std::min(f, omax);
					}
					break;
				}
				case constantFluxCyl: { // upper = outer, with switch for maximum in- or outflow
					f = -bcTopValues_[0]*rhoW/(24.*60.*60.)/100 * pos0;  // [cm /day] -> [kg/(m²*s)] (Eqns are multiplied by cylindrical radius)
					if (f < 0.) { // inflow
						Scalar imax = rhoW * kc * ((h - 0.) / dz - gravityOn_)* pos0; // maximal inflow
						f = std::max(f, imax);
					} else { // outflow
						Scalar omax = rhoW *  kc * krw *((h - criticalPressure_) / dz - gravityOn_)* pos0; // maximal outflow (evaporation)
						f = std::min(f, omax);
					}
					break;
				}
				case atmospheric: { // atmospheric boundary condition (with surface run-off)
					Scalar prec = -precipitation_.f(time_);
					if (prec < 0.) { // precipitation
						// std::cout << "in" << "\n";
						Scalar imax = rhoW * kc * ((h - 0.) / dz - gravityOn_); // maximal infiltration
						f = std::max(prec, imax);
					} else { // evaporation
						Scalar emax = rhoW * krw * kc * ((h - criticalPressure_) / dz - gravityOn_); // maximal evaporation

						
						// // std::cout << "out" << ", at " << h << " cm \n";
					    // Scalar p2 = toPa_(-10000);
					    // // Scalar h3 = 0.5*(h + criticalPressure_);
					    // // Scalar p3 = toPa_(h);
					    // Scalar s2 = materialLaw_.sw(-(p2- pRef_));
                        // // Scalar s3 = MaterialLaw::sw(params, -(p3- pRef_)) ;
					    // // std::cout << s2 << "\n";
					    // Scalar krw2 = materialLaw_.krw(s2);
					    // // Scalar krw3 = MaterialLaw::krw(params, s3);
                        // Scalar arithmetic = 0.5*(krw2+krw); // arithmetic currently best
					    // // Scalar harmonic = 2*krw2*krw/(krw2+krw);
						// Scalar emax = rhoW * kc * arithmetic *((h - criticalPressure_) / dz + gravityOn_); // maximal evaporation KRW???
						f = std::min(prec, emax);
					}
					break;
				}
				default: DUNE_THROW(Dune::InvalidStateException, "Top boundary type Neumann: unknown error");
				}
			} else if (onLowerBoundary_(pos)) { // bot bc
				switch (bcBotType_) {
				case constantFlux: { // with switch for maximum in- or outflow
					f = -bcBotValues_[0]*rhoW/(24.*60.*60.)/100. * pos0; // [cm /day] * [mol water / m^3] * [d/s] * [m/cm] -> [kg/(m²*s)]
					if (f < 0.) { // inflow
						Scalar imax = rhoW * kc * ((h - 0.) / dz - gravityOn_) * pos0; // maximal inflow
						imax = std::min(imax, 0.); // must stay negative
						f = std::max(f, imax);
					} else { // outflow
						Scalar omax = rhoW * kc * krw *((h - criticalPressure_) / dz - gravityOn_) * pos0; // maximal outflow (evaporation)
						// std::cout << "outflow " << f << ", " << omax << "\n";
						omax = std::max(omax, 0.); // must stay positive
						f = std::min(f, omax);
					}
					break;
				}
				case constantFluxCyl: { // lower = inner, with switch for maximum in- or outflow
					f = -bcBotValues_[0]*rhoW/(24.*60.*60.)/100. * pos0; // [cm /day] -> [kg/(m²*s)]  (Eqns are multiplied by cylindrical radius)
					if (f < 0.) { // inflow
						Scalar imax = rhoW * kc * ((h - (-10.)) / dz - gravityOn_)* pos0; // maximal inflow
						imax = std::min(imax, 0.); // must stay negative
						f = std::max(f, imax);
					} else { // outflow
						Scalar omax = rhoW * kc * krw *((h - criticalPressure_) / dz - gravityOn_)* pos0; // maximal outflow (evaporation)
//						std::cout << " f " << f*1.e9 << ", omax "<< omax*1.e9  << ", value " << bcBotValues_.at(0)
//								<< ", crit "  << criticalPressure_ << ", " << pos[0] << ", krw " << krw <<"\n";
						omax = std::max(omax, 0.); // must stay positive
						f = std::min(f, omax);
					}
					break;
				}
				case rootSystemExact: { // [x0, x1, self.rs.kr_f(age, type_), self.rs.kx_f(age, type_), self.radii[j], l])

					double x0 = bcBotValues_[0]; // node matric potential [cm]
					double x1 = bcBotValues_[1]; // node matric potential [cm]
					double kr = bcBotValues_[2]; // [1 day-1]
					double kx = bcBotValues_[3]; // [cm3 day-1]
					double l = bcBotValues_[4]; // [cm]
					double h0 = bcBotValues_[5]; // matric potential at root soil interface, rsx [cm] (not dynamic, as python)
					// double h0 = h; // dynamic
					double a = 100.*pos[0]; // [cm]

					double f_ = 2 * a * M_PI * kr;
					double tau = std::sqrt(f_ / kx);  // sqrt(c) [cm-1]
					double d = std::exp(-tau * l) - std::exp(tau * l);  //  det
					double fExact = f_ * (1. / (tau * d)) * (x0 - h0 + x1 - h0) * (2. - std::exp(-tau * l) - std::exp(tau * l)); // h or h0?
					f =  fExact / (2 * a * M_PI * l);  // [cm3 / cm2 / day]
					f *= (rhoW*1.e-2) / (24.*3600.); // [cm3 / cm2 / day] -> [kg/(m2*s)]
					f *= -pos[0]; // cylindrical coordinates

					Scalar omax = rhoW * kc * krw *((h - criticalPressure_) / dz - gravityOn_)* pos[0]; // maximal outflow (evaporation)
					f = std::min(f, omax);

					// classical approximation
//					f = kr*(0.5*(x0+x1)-h); // [cm/day]
//					f *= (rhoW*1.e-2) / (24.*3600.); // [cm3 / cm2 / day] -> [kg/(m2*s)]
//					f = std::min(f, 0.);
//					f *= -pos[0]; // cylindrical coordinates

					break; // [kg/(m²*s)]
				}
				case freeDrainage: {
					f = krw * kc * rhoW  *pos0; // [m/s] * [1] * [mol / m3] = [mol / (m2 * s)]  //* 1 [m]
					break;
				}
				default: DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Neumann: unknown error");
				}
			}
		}
		flux[conti0EqIdx] = f;
		return flux;
	}
	struct InvalidElemSol
    {
        template<class Index>
        double operator[] (const Index i) const
        { static_assert(AlwaysFalse<Index>::value, "Solution-dependent material parameters not supported with analytical differentiation"); return 0.0; }
    };
	
	template<class PartialDerivativeMatrices>
	void addRobinFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                             const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& curElemVolVars,
                             const ElementFluxVariablesCache& elemFluxVarsCache,
                             const SubControlVolumeFace& scvf) const
{
    const auto insideScvIdx = scvf.insideScvIdx();
    const auto& insideScv = fvGeometry.scv(insideScvIdx);
    const auto& insideVolVars = curElemVolVars[insideScvIdx];

    GlobalPosition pos = scvf.center();
    double pos0 = 1.;
    if ((dimWorld == 1) && (!this->spatialParams().useExtrusion))
        pos0 = pos[0];

    // only free drainage has a nonzero derivative w.r.t. pw . TODO: add derivatives for the other BCs
    if (onLowerBoundary_(pos) && bcBotType_ == freeDrainage) {

        static const auto rhoW = useMoles
            ? insideVolVars.molarDensity(0)
            : insideVolVars.density(0);

        const auto kc = this->spatialParams().hydraulicConductivity(element); // [m/s]

        // material law derivatives: dkrw/dsw * dsw/dpw
        const auto insideFluidMatrixInteraction =
            this->spatialParams().fluidMatrixInteraction(element, insideScv, InvalidElemSol{});

        const auto sw  = insideVolVars.saturation(0);
        const auto pc  = insideVolVars.capillaryPressure();
        const auto dkrw_dsw = insideFluidMatrixInteraction.dkrw_dsw(sw);
        const auto dsw_dpw  = -insideFluidMatrixInteraction.dsw_dpc(pc); // dsw/dpc * dpc/dpw, dpc/dpw = -1

        // df/dpw = rhoW * kc * (dkrw/dsw * dsw/dpw) * pos0
        derivativeMatrices[insideScvIdx][conti0EqIdx][0] += rhoW * kc * dkrw_dsw * dsw_dpw * pos0;
    }
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
		double pos0 = 1.; 
        double svc_volume = scv.volume();
        int dofIndex = scv.dofIndex();
		if (dimWorld == 1)
		{
            svc_volume = getCellVolumesCyl(dofIndex);//with 1d model, need to evaluate manually the volume of the cell.
                            // for simplicity, we give directly source as [ mol / (m^3 \cdot s)] for now
			if(!this->spatialParams().useExtrusion){
				pos0 = pos[0];
			}
		}
		if ((source_ != nullptr)) {
			if (sourceSlope_<0.) {
				auto eIdx = this->spatialParams().gridGeometry().elementMapper().index(element);
				return source_->at(eIdx)/svc_volume * pos0;
			} else {
			    Scalar s = elemVolVars[scv].saturation();
	            Scalar p = materialLaw(element).pc(s) + pRef_; // [Pa]
	            Scalar h = -toHead_(p); // cm
	            auto eIdx = this->spatialParams().gridGeometry().elementMapper().index(element);
	            if (h<criticalPressure_) {
	                return 0.;
	            } else if (h<=criticalPressure_+sourceSlope_) { //  h in [crit, crit+slope]
	                double theta = (h - criticalPressure_)/sourceSlope_;
	                // std::cout << "source(): " << h << ", "<< theta << "\n" << std::flush;
	                return theta* source_->at(eIdx)/svc_volume * pos0;
	            } else  {
	                return source_->at(eIdx)/svc_volume * pos0;
	            }
			}
		} else {
			return 0.;
		}
	}
	
		/*!
	 *
     * E.g. for the mol balance that would be a mass rate in \f$ [ mol / (m^3 \cdot s)] \f
     */
	void bioChemicalReaction(const Element &element, 
								NumEqVector &q, const VolumeVariables &volVars, double pos0, 
								const SubControlVolume &scv) const
	{}
    
    
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
//		if (couplingManager_!=nullptr) {
//			pointSources = couplingManager_->bulkPointSources();
//		}
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
	template<class ElementVolumeVariables>
	void pointSource(PointSource& source,
			const Element &element,
			const FVElementGeometry& fvGeometry,
			const ElementVolumeVariables& elemVolVars,
			const SubControlVolume &scv) const {
//		if (couplingManager_!=nullptr) {
//			// compute source at every integration point
//			const Scalar pressure3D = couplingManager_->bulkPriVars(source.id())[Indices::pressureIdx];
//			const Scalar pressure1D = couplingManager_->lowDimPriVars(source.id())[Indices::pressureIdx];
//			const auto& spatialParams = couplingManager_->problem(Dune::index_constant<1>{}).spatialParams();
//			const auto lowDimElementIdx = couplingManager_->pointSourceData(source.id()).lowDimElementIdx();
//			const Scalar kr = spatialParams.kr(lowDimElementIdx);
//			const Scalar rootRadius = spatialParams.radius(lowDimElementIdx);
//			// relative soil permeability
//			const auto krel = 1.0;
//			// sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
//			const auto density = 1000;
//			const Scalar sourceValue = 2 * M_PI *krel*rootRadius * kr *(pressure1D - pressure3D)*density;
//			source = sourceValue*source.quadratureWeight()*source.integrationElement();
//			//std::cout << "pointSource " << source.id() << ": " << sourceValue << " -> " << sourceValue*source.quadratureWeight()*source.integrationElement() << "\n";
//		} else {
//			source = 0;
//		}
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
//	void setCouplingManager(CouplingManager* cm) {
//		couplingManager_ = cm;
//	}
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
						//bc_flux_upper += neumann(e, fvGeometry, elemVolVars, scvf);
						uc++;
					} else if (onLowerBoundary_(p)) { // bottom
						//bc_flux_lower += neumann(e, fvGeometry, elemVolVars, scvf);
						lc++;
					}
				} catch (...) {
				}
			}
		}
		bc_flux_upper /= uc;
		bc_flux_lower /= lc;
	}

//    NumEqVector neumann(const Element& element,
//            const FVElementGeometry& fvGeometry,
//            const ElementVolumeVariables& elemVolVars,
//            const ElementFluxVariablesCache&  fluxCache,
//            const SubControlVolumeFace& scvf)


	/*!
	 * Writes the actual boundary fluxes (top and bottom) into a text file. Call postTimeStep before using it.
	 */
	void writeBoundaryFluxes() {
		if (writeFile_) {
			myfile_ << time_ << ", " << bc_flux_upper << ", " << bc_flux_lower << "\n";
		}
	}

	/**
	 * Debug info
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

    /**
     * Forwards to spatialParams
     */
    void addVanGenuchtenDomain(double minx, double miny, double minz, double maxx, double maxy, double maxz, int layerIndex) {
        this->spatialParams().addVanGenuchtenDomain(minx, miny, minz, maxx, maxy, maxz, layerIndex);
    }

    /**
     * Forwards to spatialParams
     */
    void changeVanGenuchtenSet(int vgIndex, double qr, double qs, double alpha, double n, double ks){
        this->spatialParams().changeVanGenuchtenSet(vgIndex, qr, qs, alpha, n, ks);
    }

	void setFaceGlobalIndexSet(std::map<int,int>  faceIdx_)
	{
		faceIdx = faceIdx_;											   
	}
	// BC, direct access for Python binding (setBcTop, setBcBot)
	int bcTopType_;
	int bcBotType_;
	std::vector<double> bcTopValues_ = std::vector<double>(0);
	std::vector<double> bcBotValues_ = std::vector<double>(0);


	// dummy for richardsnc types (needed by richards.hh)
	std::vector<int> bcSTopType_ = std::vector<int>(0); // zero solutes
	std::vector<int> bcSBotType_ = std::vector<int>(0);
	std::vector<double> bcSTopValue_ = std::vector<double>(0);
	std::vector<double> bcSBotValue_= std::vector<double>(0);

	double segLength;//m
	std::vector<double> cellVolumesCyl;
	int nCells_all;
	int nFaces;
	bool computedCellVolumesCyl = false;
	std::map<int,int>  faceIdx;
    int dzScaling;
	int verbose;
    static constexpr Scalar eps_ = 1.e-7;
	static constexpr Scalar g_ = 9.81; // cm / s^2 (for type conversions)
	Scalar rho_ = FluidSystem::H2O::liquidDensity(0,0); // 1.e3; // kg water / m^3 (for type conversions)
	Scalar rho_molar = FluidSystem::H2O::liquidMolarDensity(0,0); // mol water / m^3 (for type conversions)
	static constexpr Scalar pRef_ = 1.e5; // Pa

	//! true if on the point lies on the upper boundary
	bool onLowerBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[dimWorld - 1] < this->gridGeometry().bBoxMin()[dimWorld - 1] + eps_;
	}
	
	bool verbose_assembler = false;
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


	//! true if on the point lies on the left boundary
	bool onLeftBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[0] > this->gridGeometry().bBoxMin()[0] - eps_;
	}
	//! true if on the point lies on the right boundary
	bool onRightBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
	}
	// Initial
	InputFileFunction initialSoil_;
	bool gravityOn_;

	// Source
	std::shared_ptr<std::vector<double>> source_; // [kg/s]
	//CouplingManager* couplingManager_ = nullptr;

	InputFileFunction precipitation_;
	Scalar criticalPressure_; // cm
	Scalar time_ = 0.;
	Scalar dt_ = 0.;
	Scalar sourceSlope_ = -1.; // slope for regularization, negative values disables regularisation

	bool writeFile_ = true;
	std::ofstream myfile_;
	Scalar bc_flux_upper = 0.;
	Scalar bc_flux_lower = 0.;

	double temperatureK;
	
};

} //end namespace Dumux

#endif
