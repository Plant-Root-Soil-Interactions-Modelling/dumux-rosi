// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef RICHARDS1P2C_PROBLEM_HH
#define RICHARDS1P2C_PROBLEM_HH

#include <dumux/porousmediumflow/problem.hh> // base class

#include "../soil_richards/richardsparams.hh"

namespace Dumux {

/*!
 * RichardsProblem:
 * Uses Dumux as an easy to use Richards equation solver,
 * where most parameters can be set dynamically
 */
template <class TypeTag>
class Richards1P2CProblem : public PorousMediumFlowProblem<TypeTag>
{
public:

    // exports, used by the binding
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

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

    enum {
        pressureIdx = 0, // index of primary variables
        h2OIdx = 0, // fluid index
        soluteIdx = 1, // solute index
        conti0EqIdx = 0, // indices of the equations
        transportEqIdx = 1,

        dimWorld = GridView::dimensionworld,

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
		michaelisMenten = 8
    };

    enum GridParameterIndex {
        materialLayerNumber = 0
    };

    /*!
     * \brief Constructor: constructed in the main file
     */
    Richards1P2CProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : PorousMediumFlowProblem<TypeTag>(fvGridGeometry) {

		gravityOn_ = Dumux::getParam<bool>("Problem.EnableGravity", true);

		source_.resize(2); // 2 equations (currently hard coded, where can I get the value?)
		source_[0] = nullptr;
		source_[1] = nullptr;

        // BC
        bcTopType_ = getParam<int>("Soil.BC.Top.Type");
        bcBotType_ = getParam<int>("Soil.BC.Bot.Type");
        bcTopValue_ = getParam<Scalar>("Soil.BC.Top.Value",0.);
        bcBotValue_ = getParam<Scalar>("Soil.BC.Bot.Value",0.);

        // Component
        bcSTopType_ = getParam<int>("Soil.BC.Top.SType", outflow);
        bcSBotType_ = getParam<int>("Soil.BC.Bot.SType", outflow);
        bcSTopValue_ = getParam<Scalar>("Soil.BC.Top.CValue", 0.);
        bcSBotValue_ = getParam<Scalar>("Soil.BC.Bot.CValue", 0.);

		criticalPressure_ = getParam<double>("Soil.CriticalPressure", -1.e4); // cm
		criticalPressure_ = getParam<double>("Climate.CriticalPressure", criticalPressure_); // cm
        // Precipitation & Evaporation
        if (bcTopType_==atmospheric) {
            precipitation_ = InputFileFunction("Climate", "Precipitation", "Time", 0.); // cm/day (day)
            precipitation_.setVariableScale(1./(24.*60.*60.)); // s -> day
            precipitation_.setFunctionScale(1.e3/(24.*60.*60.)/100); // cm/day -> kg/(m²*s)
        }
        // IC
        initialSoilP_ = InputFileFunction("Soil.IC", "P", "Z", 0., this->spatialParams().layerIFF()); // [cm]([m]) pressure head, conversions hard coded
        initialSoilC_ = InputFileFunction("Soil.IC", "C", "CZ", 0., this->spatialParams().layerIFF()); //

        // Uptake params
        vMax_ =  getParam<Scalar>("RootSystem.Uptake.Vmax", 6.2e-11/(10./24./3600.))*10./24./3600.; //  [g cm-2 day-1] -> [kg m-2 s-1]
        km_ = getParam<Scalar>("RootSystem.Uptake.Km", 3.1e-9 /1000. )*1000.;  // [g cm-3] -> [kg m-3]
        sigma_ = getParam<Scalar>("RootSystem.Uptake.ActiveTransport", 0.); // 1 for passive transport, 0 for active transport

        // Buffer power
        b_ = getParam<Scalar>("Component.BufferPower", 0.);

        // Output
        std::string filestr = this->name() + ".csv"; // output file
        myfile_.open(filestr.c_str());
		std::cout << "RichardsProblem constructed: bcTopType " << bcTopType_ << ", " << bcTopValue_ << "; bcBotType "
				<<  bcBotType_ << ", " << bcBotValue_  << ", gravitation " << gravityOn_ <<", Critical pressure "
				<< criticalPressure_ << "\n" << std::flush;
    }

    /**
     * \brief Eventually, closes output file
     */
    ~Richards1P2CProblem() {
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
    Scalar getBufferPower(const SubControlVolume& scv, const VolumeVariables& volVars) const {
        return b_;
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
        v[pressureIdx] = toPa_(initialSoilP_.f(z,eIdx));
        v[soluteIdx] = initialSoilC_.f(z,eIdx);
        // std::cout << v[soluteIdx] << "\n";//////////////////////////////////////////////////////////////
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
        bcTypes.setAllNeumann();
        return bcTypes;
    }

    /*!
     * \copydoc FVProblem::neumann // [kg/(m²*s)]
     *
     * called by BoxLocalResidual::evalFlux,
     * negative = influx, mass flux in \f$ [ kg / (m^2 \cdot s)] \f$
     */
    NumEqVector neumann(const Element& element,
        const FVElementGeometry& fvGeometry,
        const ElementVolumeVariables& elemVolVars,
        const SubControlVolumeFace& scvf) const {

        NumEqVector flux;
        GlobalPosition pos = scvf.center();
        auto& volVars = elemVolVars[scvf.insideScvIdx()];

        // WATER
		double f = 0.; // return value
		if ( onUpperBoundary_(pos) || onLowerBoundary_(pos) ) {

			Scalar s = volVars.saturation(0);
			Scalar kc = this->spatialParams().hydraulicConductivity(element); //  [m/s]
			MaterialLawParams params = this->spatialParams().materialLawParams(element);
			Scalar p = MaterialLaw::pc(params, s) + pRef_;
			Scalar h = -toHead_(p); // todo why minus -pc?
			GlobalPosition ePos = element.geometry().center();
			Scalar dz = 100 * 2 * std::abs(ePos[dimWorld - 1] - pos[dimWorld - 1]); // m->cm
			Scalar krw = MaterialLaw::krw(params, s);

			if (onUpperBoundary_(pos)) { // top bc
				switch (bcTopType_) {
				case constantFlux: { // with switch for maximum in- or outflow
					f = -bcTopValue_*rho_/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
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
					f = -bcTopValue_*rho_/(24.*60.*60.)/100 * pos[0];
					if (f < 0) { // inflow
						Scalar imax = rho_ * kc * ((h - 0.) / dz - gravityOn_)* pos[0]; // maximal inflow
						f = std::max(f, imax);
					} else { // outflow
						Scalar omax = rho_ * krw * kc * ((h - criticalPressure_) / dz - gravityOn_)* pos[0]; // maximal outflow (evaporation)
						f = std::min(f, omax);
					}
					break;
				}
				case atmospheric: { // atmospheric boundary condition (with surface run-off) // TODO needs testing & improvement
					Scalar prec = -precipitation_.f(time_);
					if (prec < 0) { // precipitation
						Scalar imax = rho_ * kc * ((h - 0.) / dz - gravityOn_); // maximal infiltration
						f = std::max(prec, imax);
					} else { // evaporation
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
						// std::cout << " f " << f*1.e9  << ", omax "<< omax << ", value " << bcBotValue_ << ", crit "  << criticalPressure_ << ", " << pos[0] << "\n";
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

        // SOLUTE
        if (onUpperBoundary_(pos)) { // top bc Solute
          switch (bcSTopType_) {
            case constantConcentration: {
                flux[transportEqIdx] = f* (bcSTopValue_ - volVars.massFraction(0, soluteIdx)) ; // TODO ???
                break;
            }
            case constantFlux: {
                flux[transportEqIdx] = -bcSTopValue_*rho_/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
                break;
            }
            case outflow: {
                flux[transportEqIdx] = f * volVars.massFraction(0, soluteIdx);
                break;
            }
            case linear: {
            	flux[transportEqIdx] = vMax_ * volVars.massFraction(0, soluteIdx);
                break;
            }
            case michaelisMenten: {
            	flux[transportEqIdx] = vMax_ * (volVars.massFraction(0, soluteIdx)*rho_)/(km_ + volVars.massFraction(0, soluteIdx)*rho_);
                break;
            }
            default:
                DUNE_THROW(Dune::InvalidStateException, "Top boundary type Neumann: unknown error");
            }
        } else if (onLowerBoundary_(pos)) { // bot bc Solute
            switch (bcSBotType_) {
            case constantConcentration: {
                flux[transportEqIdx] = f * (bcSBotValue_ - volVars.massFraction(0, soluteIdx)) ; // TODO ???
                break;
            }
            case constantFlux: {
                flux[transportEqIdx] = -bcSBotValue_*rho_/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
                break;
            }
            case outflow: {
                flux[transportEqIdx] = f * volVars.massFraction(0, soluteIdx);
                break;
            }
            case linear: { // -2 * M_PI * rootRadius * vMax_ * soilC * density/(km_ + soilC * density);
            	flux[transportEqIdx] = vMax_ * volVars.massFraction(0, soluteIdx);
                break;
            }
            case michaelisMenten: {
            	flux[transportEqIdx] = vMax_ * (volVars.massFraction(0, soluteIdx)*rho_)/(km_ + volVars.massFraction(0, soluteIdx)*rho_);
                break;
            }
            default: DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Neumann: unknown error");
            }
        } else {
            flux[transportEqIdx] = 0.;
        }

        return flux;
    }

    /*!
     * \copydoc FVProblem::source
     *
     * called by FVLocalResidual:computeSource(...)
     */
    NumEqVector source(const Element &element, const FVElementGeometry& fvGeometry, const ElementVolumeVariables& elemVolVars,
        const SubControlVolume &scv) const {

    	auto eIdx = this->spatialParams().fvGridGeometry().elementMapper().index(element);
    	double h = 0.;
    	if (source_[h2OIdx] != nullptr) {
        	h = source_[h2OIdx]->at(eIdx);
    	}
        double s = 0.;
    	if (source_[soluteIdx] != nullptr) {
        	s = source_[soluteIdx]->at(eIdx);
    	}
        return NumEqVector({ h, s });
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
    template<class ElementVolumeVariables>
    void pointSource(PointSource& source,
        const Element &element,
        const FVElementGeometry& fvGeometry,
        const ElementVolumeVariables& elemVolVars,
        const SubControlVolume &scv) const {

        PrimaryVariables sourceValue(0.);

        if (couplingManager_!=nullptr) { // compute source at every integration point

            const Scalar soilP = couplingManager_->bulkPriVars(source.id())[pressureIdx];
            const Scalar tipP = couplingManager_->lowDimPriVars(source.id())[pressureIdx];
            const auto& spatialParams = couplingManager_->problem(Dune::index_constant<1>{}).spatialParams();
            const auto lowDimElementIdx = couplingManager_->pointSourceData(source.id()).lowDimElementIdx();
            const Scalar kr = spatialParams.kr(lowDimElementIdx);
            const Scalar rootRadius = spatialParams.radius(lowDimElementIdx);
            // relative soil permeability
            const auto krel = 1.0;
            // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
            sourceValue[h2OIdx] = 2 * M_PI *krel*rootRadius * kr *(tipP - soilP)*rho_;
            sourceValue[h2OIdx] *= source.quadratureWeight()*source.integrationElement();

            Scalar tipC = couplingManager_ ->lowDimPriVars(source.id())[soluteIdx]; // units [1], fraction
            Scalar soilC = couplingManager_ ->bulkPriVars(source.id())[soluteIdx]; // units [1], fraction
            Scalar passiveUptake;
            if (sourceValue[h2OIdx]>0) {
                passiveUptake = 2 * M_PI * rootRadius * kr * (tipP - soilP) * rho_ * tipC;
            } else {
                passiveUptake = 2 * M_PI * rootRadius * kr * (tipP - soilP) * rho_ * soilC;
            }
            // Active uptake based on Michaelis Menten
            Scalar activeUptake = -2 * M_PI * rootRadius * vMax_ * soilC * rho_/(km_ + soilC * rho_); // todo times root element length

            // choose active or passive
            sourceValue[transportEqIdx] = (sigma_*activeUptake + (1.-sigma_)*passiveUptake) *source.quadratureWeight()*source.integrationElement();

            source = sourceValue;
        }

        source = sourceValue; // return value as reference
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
                    bc_flux_upper += neumann(e, fvGeometry, elemVolVars, scvf);
                    uc++;
                } else if (onLowerBoundary_(p)) { // bottom
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
    InputFileFunction initialSoilP_;
    InputFileFunction initialSoilC_;

    // BC
    int bcTopType_;
    int bcBotType_;
    Scalar bcTopValue_;
    Scalar bcBotValue_;
	bool gravityOn_;

    int bcSTopType_;
    int bcSBotType_;
    Scalar bcSTopValue_;
    Scalar bcSBotValue_;

    // Source
    std::vector<std::shared_ptr<std::vector<double>>> source_;
    CouplingManager* couplingManager_ = nullptr;

    InputFileFunction precipitation_;
    Scalar criticalPressure_; // cm
    Scalar time_ = 0.;
    Scalar dt_ = 0.;

    std::ofstream myfile_;
    NumEqVector bc_flux_upper = NumEqVector(0.);
    NumEqVector bc_flux_lower = NumEqVector(0.);

    static constexpr Scalar eps_ = 1.e-7;
    static constexpr Scalar g_ = 9.81; // cm / s^2 (for type conversions)
    static constexpr Scalar rho_ = 1.e3; // kg / m^3 (for type conversions)
    static constexpr Scalar pRef_ = 1.e5; // Pa

    // Uptake params
    Scalar vMax_; // Michaelis Menten Parameter [kg m-2 s-1]
    Scalar km_;  // Michaelis Menten Parameter  [kg m-3]
    Scalar sigma_;// 1 for passive transport, 0 for active transport

    Scalar b_; // buffer power

};

} //end namespace Dumux

#endif
