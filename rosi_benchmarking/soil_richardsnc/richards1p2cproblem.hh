// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef DUMUX_RICHARDS1P2C_PROBLEM_HH
#define DUMUX_RICHARDS1P2C_PROBLEM_HH

#include <dumux/porousmediumflow/problem.hh> // base class
#include <dumux/porousmediumflow/richardsnc/model.hh>

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
        atmospheric = 4,
        freeDrainage = 5,
        outflow = 5
    };

    enum GridParameterIndex {
        materialLayerNumber = 0
    };

    /*!
     * \brief Constructor: constructed in the main file
     */
    Richards1P2CProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : PorousMediumFlowProblem<TypeTag>(fvGridGeometry) {

        // BC
        bcTopType_ = getParam<int>("Soil.BC.Top.Type"); // todo type as a string might be nicer
        bcBotType_ = getParam<int>("Soil.BC.Bot.Type");
        bcTopValue_ = getParam<Scalar>("Soil.BC.Top.Value",0.);
        bcBotValue_ = getParam<Scalar>("Soil.BC.Bot.Value",0.);
//        if (bcTopType_ == constantPressure ) {
//            bcTopValue_ = toPa_(bcTopValue_);
//        }
//        if (bcBotType_ == constantPressure ) {
//            bcBotValue_ = toPa_(bcBotValue_);
//        }

        // Component
        bcSTopType_ = getParam<int>("Soil.BC.Top.SType", outflow); // todo type as a string might be nicer
        bcSBotType_ = getParam<int>("Soil.BC.Bot.SType", outflow);
        bcSTopValue_ = getParam<Scalar>("Soil.BC.Top.CValue", 0.);
        bcSBotValue_ = getParam<Scalar>("Soil.BC.Bot.CValue", 0.);

        criticalPressure_ = getParam<double>("Climate.CriticalPressure", -1.e4); // cm
        // Precipitation & Evaporation
        if (bcTopType_==atmospheric) {
            precipitation_ = InputFileFunction("Climate", "Precipitation", "Time", 0.); // cm/day (day)
            precipitation_.setVariableScale(1./(24.*60.*60.)); // s -> day
            precipitation_.setFunctionScale(1.e3/(24.*60.*60.)/100); // cm/day -> kg/(m²*s)
        }
        // IC
        initialSoilP_ = InputFileFunction("Soil.IC", "P", "Z", 0., this->spatialParams().layerIFF()); // [cm]([m]) pressure head, conversions hard coded
        initialSoilC_ = InputFileFunction("Soil.IC", "C", "CZ", 0., this->spatialParams().layerIFF()); // [cm]([m]) pressure head, conversions hard coded
        // Output
        std::string filestr = this->name() + ".csv"; // output file
        myfile_.open(filestr.c_str());
        std::cout << "RichardsProblem constructed: bcTopType " << bcTopType_ << ", " << bcTopValue_ << "; bcBotType "
            <<  bcBotType_ << ", " << bcBotValue_ << "\n" << std::flush;
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
    NumEqVector neumann(const Element& e,
        const FVElementGeometry& fvGeometry,
        const ElementVolumeVariables& elemVolVars,
        const SubControlVolumeFace& scvf) const {

        NumEqVector flux;
        GlobalPosition pos = scvf.center();
        auto& volVars = elemVolVars[scvf.insideScvIdx()];

        if (onUpperBoundary_(pos)) { // top bc Water
            switch (bcTopType_) {
            case constantPressure: {
                double dist = (pos - fvGeometry.scv(scvf.insideScvIdx()).center()).two_norm();
                double v = this->spatialParams().permeability(e) * (volVars.pressure() - bcTopValue_) / dist; // [kg/s]
                flux[conti0EqIdx] = v / scvf.area();
                break;
            }
            case constantFlux: {
                flux[conti0EqIdx] = -bcTopValue_*rho_/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
                break;
            }
            case atmospheric: { // atmospheric boundary condition (with surface run-off) // TODO needs testing & improvement
                Scalar Kc = this->spatialParams().hydraulicConductivity(e); //  [m/s]
                Scalar s = volVars.saturation(0);
                Scalar h = toHead_(volVars.pressure());
                GlobalPosition ePos = e.geometry().center();
                Scalar dz = 100 * std::abs(ePos[dimWorld - 1] - pos[dimWorld - 1]); // cm
                Scalar prec = -precipitation_.f(time_);
                if (prec < 0) { // precipitation
                    Scalar imax = rho_ * Kc * ((h - 0.) / dz - 1.); // maximal infiltration
                    Scalar v = std::max(prec, imax);
                    flux[conti0EqIdx] = v;
                } else { // evaporation
                    MaterialLawParams params = this->spatialParams().materialLawParams(e);
                    Scalar krw = MaterialLaw::krw(params, s);
                    Scalar emax = rho_ * krw * Kc * ((h - criticalPressure_) / dz - 1.); // maximal evaporation
                    Scalar v = std::min(prec, emax);
                    // std::cout << prec << ", " << emax << ", " << h << "\n";
                    flux[conti0EqIdx] = v;
                }
                break;
            }
            default:
                DUNE_THROW(Dune::InvalidStateException, "Top boundary type Neumann: unknown error");
            }
        } else if (onLowerBoundary_(pos)) { // bot bc Water
            switch (bcBotType_) {
            case constantPressure: {
                Scalar Kc = this->spatialParams().hydraulicConductivity(e); //  [m/s]
                Scalar s = volVars.saturation(0);
                MaterialLawParams params = this->spatialParams().materialLawParams(e);
                Scalar h = toHead_(volVars.pressure());
                Scalar dz = 100 * std::abs(e.geometry().center()[dimWorld - 1] - pos[dimWorld - 1]); // cm
                flux[conti0EqIdx] = - rho_ * Kc * ((bcBotValue_ - h) / dz - 1.);
                break;
            }
            case constantFlux: {
                flux[conti0EqIdx] = -bcBotValue_*rho_/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
                break;
            }
            case freeDrainage: {
                Scalar Kc = this->spatialParams().hydraulicConductivity(e); // [m/s]
                Scalar s = volVars.saturation(0);
                MaterialLawParams params = this->spatialParams().materialLawParams(e);
                Scalar krw = MaterialLaw::krw(params, s);
                flux[conti0EqIdx] = krw * Kc * rho_; // [kg/(m^2 s)]
                break;
            }
            default:
                DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Neumann: unknown error");
            }
        } else {
            flux[conti0EqIdx] = 0.;
        }

        if (onUpperBoundary_(pos)) { // top bc Solute
          switch (bcSTopType_) {
            case constantConcentration: {
                double dist = (pos - fvGeometry.scv(scvf.insideScvIdx()).center()).two_norm();
                double v = this->spatialParams().permeability(e) * (volVars.massFraction(0, soluteIdx) - bcSTopValue_) / dist; // ??? TODO [kg/s]
                flux[transportEqIdx] = v / scvf.area();
                break;
            }
            case constantFlux: {
                flux[transportEqIdx] = -bcSTopValue_*rho_/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
                break;
            }
            case outflow: {
                flux[transportEqIdx] = flux[conti0EqIdx] * volVars.massFraction(0, soluteIdx);
                break;
            }
            default:
                DUNE_THROW(Dune::InvalidStateException, "Top boundary type Neumann: unknown error");
            }
        } else if (onLowerBoundary_(pos)) { // bot bc Solute
            switch (bcSBotType_) {
            case constantConcentration: {
                double dist = (pos - fvGeometry.scv(scvf.insideScvIdx()).center()).two_norm();
                double v = this->spatialParams().permeability(e) * (volVars.massFraction(0, soluteIdx) - bcSBotValue_) / dist; //  ??? TODO [kg/s]
                flux[transportEqIdx] = v / scvf.area();
                break;
            }
            case constantFlux: {
                flux[transportEqIdx] = -bcSBotValue_*rho_/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
                break;
            }
            case outflow: {
                flux[transportEqIdx] = flux[conti0EqIdx] * volVars.massFraction(0, soluteIdx);
                break;
            }
            default:
                DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Neumann: unknown error");
            }
        } else {
            flux[transportEqIdx] = 0.;
        }
        flux[transportEqIdx] = 0.; // TODO remvove
        return flux;
    }

    /*!
     * \copydoc FVProblem::source
     *
     * called by FVLocalResidual:computeSource(...)
     */
    NumEqVector source(const Element &element, const FVElementGeometry& fvGeometry, const ElementVolumeVariables& elemVolVars,
        const SubControlVolume &scv) const {
        if ((source_ != nullptr)) {
            auto eIdx = this->spatialParams().fvGridGeometry().elementMapper().index(element);
            double v = scv.volume();
            auto s = source_->at(eIdx);
            s[h2OIdx] = s[h2OIdx]/v;
            s[soluteIdx] = s[soluteIdx]/v;
            return s;
        } else {
            return NumEqVector(0.);
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
            const Scalar pressure3D = couplingManager_->bulkPriVars(source.id())[pressureIdx];
            const Scalar pressure1D = couplingManager_->lowDimPriVars(source.id())[pressureIdx];
            const auto& spatialParams = couplingManager_->problem(Dune::index_constant<1>{}).spatialParams();
            const auto lowDimElementIdx = couplingManager_->pointSourceData(source.id()).lowDimElementIdx();
            const Scalar kr = spatialParams.kr(lowDimElementIdx);
            const Scalar rootRadius = spatialParams.radius(lowDimElementIdx);
            // relative soil permeability
            const auto krel = 1.0;
            // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
            const auto density = 1000;
            sourceValue[h2OIdx] = 2 * M_PI *krel*rootRadius * kr *(pressure1D - pressure3D)*density;
            sourceValue[h2OIdx] *= source.quadratureWeight()*source.integrationElement();
            sourceValue[soluteIdx] = 0.; // TODO this might change (on both sides) if we couple richards1c and roots1pnc
            //std::cout << "pointSource " << source.id() << ": " << sourceValue << " -> " << sourceValue*source.quadratureWeight()*source.integrationElement() << "\n";
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
     * todo use smart pointer
     */
    void setSource(std::vector<NumEqVector>* s) {
        source_ = s;
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
     * debug info TODO
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
        std::cout << "Global integrated source (soil): " << source << " (kg/s) / " << source*3600*24*1000 << " (g/day)" << '\n';
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
    int bcSTopType_;
    int bcSBotType_;
    Scalar bcSTopValue_;
    Scalar bcSBotValue_;

    // Source
    std::vector<NumEqVector>* source_ = nullptr;
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

};

} //end namespace Dumux

#endif
