// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef DUMUX_RICHARDS_PROBLEM_HH
#define DUMUX_RICHARDS_PROBLEM_HH

#include <dumux/porousmediumflow/problem.hh> // base class

#include "richardsparams.hh"

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
        isBox = GetPropType<TypeTag, Properties::FVGridGeometry>::discMethod == DiscretizationMethod::box
    };

    enum BCTypes {
        constantPressure = 1,
        constantFlux = 2,
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

        //gravity
        //std::vector<bool> gravity = getParam<std::vector<bool>>("Problem.EnableGravity");
        //gravity  = getParam<bool>("Problem.EnableGravity");
        gravity = Dumux::getParam<bool>("Problem.EnableGravity", 1);
            
        // BC
        bcTopType_ = getParam<int>("Soil.BC.Top.Type"); // todo type as a string might be nicer
        bcBotType_ = getParam<int>("Soil.BC.Bot.Type");
        bcTopValue_ = getParam<Scalar>("Soil.BC.Top.Value",0.);
        bcBotValue_ = getParam<Scalar>("Soil.BC.Bot.Value",0.);

        criticalPressure_ = getParam<double>("Soil.CriticalPressure", -1.e4); // cm
        criticalPressure_ = getParam<double>("Climate.CriticalPressure", criticalPressure_); // cm
        // precipitation
        if (bcTopType_==atmospheric) {
            precipitation_ = InputFileFunction("Climate", "Precipitation", "Time", 0.); // cm/day (day)
            precipitation_.setVariableScale(1./(24.*60.*60.)); // s -> day
            precipitation_.setFunctionScale(1.e3/(24.*60.*60.)/100); // cm/day -> kg/(m²*s)
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
            <<  bcBotType_ << ", " << bcBotValue_ << "\n" << std::flush;
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
        auto eIdx = this->fvGridGeometry().elementMapper().index(entity);
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
            case constantPressure:
                bcTypes.setAllDirichlet();
                break;
            case constantFlux:
                bcTypes.setAllNeumann();
                break;
            case atmospheric:
                bcTypes.setAllNeumann();
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException,"Top boundary type not implemented");
            }
        } else if (onLowerBoundary_(globalPos)) { // bot or inner bc
            switch (bcBotType_) {
            case constantPressure:
                bcTypes.setAllDirichlet();
                break;
            case constantFlux:
                bcTypes.setAllNeumann();
                break;
            case freeDrainage:
                bcTypes.setAllNeumann();
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException,"Bottom boundary type not implemented");
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
            case constantPressure:
                values[Indices::pressureIdx] = toPa_(bcTopValue_);
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException, "Top boundary type Dirichlet: unknown boundary type");
            }
        } else if (onLowerBoundary_(globalPos)) { // bot bc
            switch (bcBotType_) {
            case constantPressure:
                values[Indices::pressureIdx] = toPa_(bcBotValue_);
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Dirichlet: unknown boundary type");
            }
        }
        values.setState(Indices::bothPhases);
        return values;
    }

    /*!
     * \copydoc FVProblem::neumann // [kg/(m²*s)]
     *
     * called by BoxLocalResidual::evalFlux,
     *  mass flux in \f$ [ kg / (m^2 \cdot s)] \f$
     */
    NumEqVector neumann(const Element& element,
        const FVElementGeometry& fvGeometry,
        const ElementVolumeVariables& elemVolVars,
        const SubControlVolumeFace& scvf) const {

        NumEqVector flux;
        GlobalPosition pos = scvf.center();
        if (onUpperBoundary_(pos)) { // top bc
            switch (bcTopType_) {
            case constantFlux: { // with switch for maximum in- or outflow
                Scalar s = elemVolVars[scvf.insideScvIdx()].saturation(0);
                Scalar Kc = this->spatialParams().hydraulicConductivity(element); //  [m/s]
                MaterialLawParams params = this->spatialParams().materialLawParams(element);
                Scalar p = MaterialLaw::pc(params, s) + pRef_;
                Scalar h = -toHead_(p); // todo why minus -pc?
                GlobalPosition ePos = element.geometry().center();
                Scalar dz = 100 * 2 * std::abs(ePos[dimWorld - 1] - pos[dimWorld - 1]); // cm
                Scalar constflux = -bcTopValue_*rho_/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
                if (constflux < 0) { // inflow
                    Scalar imax = rho_ * Kc * ((h - 0.) / dz - gravity); // maximal inflow
                    Scalar v = std::max(constflux, imax);
                    flux[conti0EqIdx] = v;
                } else { // outflow
                    Scalar krw = MaterialLaw::krw(params, s);
                    Scalar omax = rho_ * krw * Kc * ((h - criticalPressure_) / dz - gravity); // maximal outflow (evaporation)
                    Scalar v = std::min(constflux, omax);
                    flux[conti0EqIdx] = v;
                }
                break;
            }
            case atmospheric: { // atmospheric boundary condition (with surface run-off) // TODO needs testing & improvement
                Scalar s = elemVolVars[scvf.insideScvIdx()].saturation(0);
                Scalar Kc = this->spatialParams().hydraulicConductivity(element); //  [m/s]
                MaterialLawParams params = this->spatialParams().materialLawParams(element);
                Scalar p = MaterialLaw::pc(params, s) + pRef_;
                Scalar h = -toHead_(p); // todo why minus -pc?
                GlobalPosition ePos = element.geometry().center();
                Scalar dz = 100 * 2 * std::abs(ePos[dimWorld - 1] - pos[dimWorld - 1]); // cm
                Scalar prec = -precipitation_.f(time_);
                if (prec < 0) { // precipitation
                    Scalar imax = rho_ * Kc * ((h - 0.) / dz - gravity); // maximal infiltration
                    Scalar v = std::max(prec, imax);
                    flux[conti0EqIdx] = v;
                } else { // evaporation
                    Scalar krw = MaterialLaw::krw(params, s);
                    Scalar emax = rho_ * krw * Kc * ((h - criticalPressure_) / dz - gravity); // maximal evaporation
                    Scalar v = std::min(prec, emax);
                    flux[conti0EqIdx] = v;
                }
                break;
            }
            default:
                DUNE_THROW(Dune::InvalidStateException, "Top boundary type Neumann: unknown error");
            }
        } else if (onLowerBoundary_(pos)) { // bot bc
            switch (bcBotType_) {
            case constantFlux: { // with switch for maximum in- or outflow
                Scalar s = elemVolVars[scvf.insideScvIdx()].saturation(0);
                Scalar Kc = this->spatialParams().hydraulicConductivity(element); //  [m/s]
                MaterialLawParams params = this->spatialParams().materialLawParams(element);
                Scalar p = MaterialLaw::pc(params, s) + pRef_;
                Scalar h = -toHead_(p); // todo why minus -pc?
                GlobalPosition ePos = element.geometry().center();
                Scalar dz = 100 * 2 * std::abs(ePos[dimWorld - 1] - pos[dimWorld - 1]); // cm
                Scalar constflux = -bcBotValue_*rho_/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
                if (constflux < 0) { // outflow (drainage)
                    Scalar krw = MaterialLaw::krw(params, s);
                    Scalar omax = rho_ * krw * Kc * ((h - criticalPressure_) / dz - gravity); // maximal outflow
                    Scalar v = std::max(constflux, omax);
                    flux[conti0EqIdx] = v;
                } else { // inflow
                    Scalar imax = rho_ * Kc * ((h - 0.) / dz - gravity); // maximal inflow (capillary rise)
                    Scalar v = std::min(constflux, imax);
                    flux[conti0EqIdx] = v;
                }
                break;
            }
            case freeDrainage: {
                Scalar Kc = this->spatialParams().hydraulicConductivity(element);
                Scalar s = elemVolVars[scvf.insideScvIdx()].saturation(0);
                MaterialLawParams params = this->spatialParams().materialLawParams(element);
                Scalar krw = MaterialLaw::krw(params, s);
                flux[conti0EqIdx] = krw * Kc * rho_; // * 1 [m]
                break;
            }
            default:
                DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Neumann: unknown error");
            }
        } else {
            flux[conti0EqIdx] = 0.;
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
        if ((source_ != nullptr)) {
            auto eIdx = this->spatialParams().fvGridGeometry().elementMapper().index(element);
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
            const Scalar kr = spatialParams.kr(lowDimElementIdx);
            const Scalar rootRadius = spatialParams.radius(lowDimElementIdx);
            // relative soil permeability
            const auto krel = 1.0;
            // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
            const auto density = 1000;
            const Scalar sourceValue = 2 * M_PI *krel*rootRadius * kr *(pressure1D - pressure3D)*density;
            source = sourceValue*source.quadratureWeight()*source.integrationElement();
            //std::cout << "pointSource " << source.id() << ": " << sourceValue << " -> " << sourceValue*source.quadratureWeight()*source.integrationElement() << "\n";
        } else {
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
     * todo use smart pointer
     */
    void setSource(std::vector<double>* s) {
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
        if (writeFile_) {
            myfile_ << time_ << ", " << bc_flux_upper << ", " << bc_flux_lower << "\n";
        }
    }

    /**
     * debug info
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
    InputFileFunction initialSoil_;

    // BC
    int bcTopType_;
    int bcBotType_;
    Scalar bcTopValue_;
    Scalar bcBotValue_;

    // Source
    std::vector<double>* source_ = nullptr;
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
