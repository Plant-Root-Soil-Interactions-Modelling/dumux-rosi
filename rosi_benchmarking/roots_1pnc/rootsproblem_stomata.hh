// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef ROOTS_1P2C_PROBLEM_HH
#define ROOTS_1P2C_PROBLEM_HH

#include <math.h>
#include <map>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/growth/soillookup.hh>

#if DGF
#include "../roots_1p/rootspatialparams_dgf.hh"
#endif
#if ROOTBOX
#include "../roots_1p/rootspatialparams_rb.hh"
#endif

namespace Dumux {

/*!
 * Root system water movement and transpiration driven by hormones.
 *
 * similar Huber et al. 2014
 *
 */
template<class TypeTag>
class RootsOnePTwoCProblem: public PorousMediumFlowProblem<TypeTag> {

    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using CouplingManager= GetPropType<TypeTag, Properties::CouplingManager>;

    enum {
        // indices of primary variables
        pressureIdx = Indices::pressureIdx,

        // component indices
        H2OIdx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::H2OIdx),
        ABAIdx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::ABAIdx),

        // indices of the equations
        contiH2OEqIdx = Indices::conti0EqIdx + H2OIdx,
        transportABAEqIdx = Indices::conti0EqIdx + ABAIdx

    };
    enum {
        isBox = GetPropType<TypeTag, Properties::FVGridGeometry>::discMethod == DiscretizationMethod::box
    };
    enum {
        bcDirichlet = 0,
        bcNeumann = 1
    };

    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static const int dimWorld = GridView::dimensionworld;

public:

    RootsOnePTwoCProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry): ParentType(fvGridGeometry) {

        // initialize fluid system
        FluidSystem::init();

        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions"<<std::endl;
        else
            std::cout<<"problem uses mass fractions"<<std::endl;

        //I.C.
        InputFileFunction sf = InputFileFunction("Soil.IC", "P", "Z", 0.); // [cm]([m])
        sf.setFunctionScale(1.e-2 * rho_ * g_ ); // [cm] -> [Pa], don't forget to add pRef_
        soil_ = new GrowthModule::SoilLookUpTable(sf); // sf is copied by default copy constructor

        if (Dumux::hasParam("RootSystem.Collar.P")) {
            collar_ = InputFileFunction("RootSystem.Collar", "P", "PT", -1.e4); // [cm]([day])
            collar_.setVariableScale(1./(24.*3600)); // [s] -> [day]
            collar_.setFunctionScale(1.e-2 * rho_ * g_); // [cm] -> [Pa], don't forget to add pRef_
            bcType_ = bcDirichlet;
        } else { // ".File " gives transpiration table
            collar_ = InputFileFunction("RootSystem.Collar", "Transpiration", "TranspirationT", 0.02);  // [kg/day]([day])
            collar_.setVariableScale(1./(24.*3600)); // [s] -> [day]
            collar_.setFunctionScale(1./(24.*3600)); // [kg/day] -> [kg/s]
            bcType_ = bcNeumann;
        }
        file_at_.open(this->name() + "_actual_transpiration.txt");
    }

    virtual ~RootsOnePTwoCProblem() {
        delete soil_;
        std::cout << "closing file \n" << std::flush;
        file_at_.close();
    }

    //! evaluates user defined data for vtk fields
    void userData(std::string name, const SolutionVector& sol) {
        const auto& gridView = this->fvGridGeometry().gridView();
        userData_[name] = std::vector<Scalar>(gridView.size(0));
        auto eMapper = this->fvGridGeometry().elementMapper();
        auto vMapper = this->fvGridGeometry().vertexMapper();
        for (const auto& e : elements(gridView)) {
            auto eIdx = eMapper.index(e);
            double d = 0;
            if (name=="kr") {
                d = 1.e4*24.*3600.*this->spatialParams().kr(eIdx); // [m/Pa/s] -> [cm/hPa/day]
            }
            if (name=="kx") {
                d = 1.e10*24.*3600.*this->spatialParams().kx(eIdx); // [m^4/Pa/s] -> [cm^4/hPa/day]
            }
            if (name=="age") {
                d = this->spatialParams().age(eIdx) / 24. / 3600.; // s -> day
            }
            if (name=="order") {
                d = this->spatialParams().order(eIdx);
            }
            if (name=="id") {
                d = this->spatialParams().id(eIdx);
            }
            if (name=="radius") {
                d = this->spatialParams().radius(eIdx);// m
            }
            if (name=="initialPressure") {
                PrimaryVariables d(0);
                d = initialAtPos(e.geometry().center()); // Pa
                double a = 100. * (d[0] - pRef_) / rho_ / g_;  // Pa -> cm
            }
            if (name=="radialFlux") {
                auto geo = e.geometry();
                auto length = geo.volume();
                auto kr = this->spatialParams().kr(eIdx);
                auto a = this->spatialParams().radius(eIdx);
                auto i0 = vMapper.subIndex(e, 0, 1);
                auto i1 = vMapper.subIndex(e, 1, 1);
                auto p = geo.center();
                // kr [m /Pa/s]
                d =  2 * a * M_PI * length* kr * (soil(p) - (sol[i1][0] + sol[i0][0]) / 2); // m^3 / s
                d = 24.*3600*1.e6*d; // [m^3/s] -> [cm^3/day]
            }
            if (name=="axialFlux") {
                auto geo = e.geometry();
                auto length = geo.volume();
                auto kx = this->spatialParams().kx(eIdx);
                auto i0 = vMapper.subIndex(e, 0, 1);
                auto i1 = vMapper.subIndex(e, 1, 1);
                d = kx * ((sol[i1][0] - sol[i0][0]) / length - rho_ * g_); // m^3 / s
                d = 24.*3600*1.e6*d; // [m^3/s] -> [cm^3/day]
            }
            if (name=="p") {
                auto i0 = vMapper.subIndex(e, 0, 1);
                auto i1 = vMapper.subIndex(e, 1, 1);
                d = 0.5 * (sol[i1][0] + sol[i0][0]);
                d = 100. * (d - pRef_) / rho_ / g_;  // Pa -> cm
            }
            userData_[name][eIdx] = d;
        }
    }

    //! vtk fields call back functions
    std::vector<Scalar>& radialFlux() { return userData_["radialFlux"]; }
    std::vector<Scalar>& axialFlux() { return userData_["axialFlux"]; } //
    std::vector<Scalar>& kr() { return userData_["kr"]; }
    std::vector<Scalar>& kx() { return userData_["kx"]; }
    std::vector<Scalar>& age() { return userData_["age"]; }
    std::vector<Scalar>& order() { return userData_["order"]; }
    std::vector<Scalar>& radius() { return userData_["radius"]; }
    std::vector<Scalar>& initialPressure() { return userData_["initialPressure"]; }
    std::vector<Scalar>& id() { return userData_["id"]; }
    std::vector<Scalar>& p() { return userData_["p"]; }

    //! calculates transpiraton, as the sum of radial fluxes (slow but accurate) [cm^3/day]
    Scalar transpiration(const SolutionVector& sol) {
        userData("radialFlux", sol);
        return std::accumulate(userData_["radialFlux"].begin(), userData_["radialFlux"].end(), 0.);
    }

    /*
     * \brief Return the temperature within the domain in [K]. (actually needed? why?)
     */
    Scalar temperature() const {
        return 273.15 + 10; // 10C
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &pos) const {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann(); // default
        if (onUpperBoundary_(pos)) { // root collar
            if (bcType_ == bcDirichlet) {
                bcTypes.setAllDirichlet();
            } else {
                if (!critical_) {
                    bcTypes.setAllNeumann();
                } else {
                    bcTypes.setAllDirichlet();
                }
            }
        } else { // for all other (i.e. root tips)
            bcTypes.setAllNeumann();
        }
        return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &pos) const {
        if (critical_) {
            return PrimaryVariables(criticalCollarPressure_);
        } else {
            return PrimaryVariables(collar_.f(time_)+pRef_);
        }
    }

    /*
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent
     *
     * Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    NumEqVector neumann(const Element& element, const FVElementGeometry& fvGeometry, const ElementVolumeVariables& elemVolVars,
        const SubControlVolumeFace& scvf) const {

        const auto globalPos = scvf.center();
        if (onUpperBoundary_(globalPos)) {
            auto& volVars = elemVolVars[scvf.insideScvIdx()];
            double p = volVars.pressure();
            auto eIdx = this->fvGridGeometry().elementMapper().index(element);
            double kx = this->spatialParams().kx(eIdx);
            auto dist = (globalPos - fvGeometry.scv(scvf.insideScvIdx()).center()).two_norm();
            collarP_ = p;
            potentialTrans_ = collar_.f(time_); // [ kg/s]
            double criticalTranspiration;
            // chemical concentration
            double cL = densityABA/MolarMass * (useMoles ? volVars.moleFraction(0, ABAIdx) :
                                                                    volVars.massFraction(0, ABAIdx)	); // (mol/m3)

            // stomatal conductance definition
            if (collarP_ < p_crit)
            {
                alpha = alphaR + (1 - alphaR)*exp(-(1-cD)*sC*cL - cD)*exp(-sH*(collarP_ - p_crit));
                criticalTranspiration = alpha * potentialTrans_; // Tact = alpha * Tpot
            }
            else
            {
                criticalTranspiration = volVars.density(0)/volVars.viscosity(0) * kx * (p - criticalCollarPressure_) / dist;
            }

            double v = std::min(potentialTrans_, criticalTranspiration);
            actualTrans_ = v;
            neumannTime_ = time_;
            maxTrans_ = criticalTranspiration;
            v /= volVars.extrusionFactor(); // [kg/s] -> [kg/(s*m^2)]
            return NumEqVector(v);
        } else {
            return NumEqVector(0.); // no flux at root tips
        }
    }

    /*!
     * For this method, the return parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    NumEqVector source(const Element &element, const FVElementGeometry& fvGeometry, const ElementVolumeVariables& elemVolVars,
        const SubControlVolume &scv) const {

        NumEqVector values;
        if (couplingManager_==nullptr) {
            auto params = this->spatialParams();
            const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
            Scalar a = params.radius(eIdx); // root radius (m)
            Scalar kr = params.kr(eIdx); //  radial conductivity (m^2 s / kg)
            Scalar phx;
            if (isBox) { // dumux
                phx = elemVolVars[scv.localDofIndex()].pressure(); // kg/m/s^2
            } else {
                phx = elemVolVars[scv.dofIndex()].pressure(); // kg/m/s^2
            }
            Scalar phs = soil(scv.center()); // kg/m/s^2
            values[contiH2OEqIdx] = kr * 2 * a * M_PI * (phs - phx); // m^2/s
            values[contiH2OEqIdx] /= (a * a * M_PI); // 1/s
            values[contiH2OEqIdx] *= rho_; // (kg/s/m^3)
        } else {
            values[contiH2OEqIdx] = 0;
        }
        return values;
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& p) const {
        return PrimaryVariables(soil(p)); // soil(p)
    }

    /**
     * deletes old soil, sets new soil, takes ownership
     */
    void setSoil(CRootBox::SoilLookUp* s) {
        delete soil_;
        soil_ = s;
        std::cout << "setSoil(...): manually changed soil to " << s->toString() << "\n";
    }

    //! soil pressure (called by initial, and source term)
    Scalar soil(const GlobalPosition& p) const {
        auto p2 = CRootBox::Vector3d(p[0] * 100, p[1] * 100, p[2] * 100); // m -> cm
        double d = soil_->getValue(p2);
        return pRef_+d;
    }

    //! sets the current simulation time [s] (within the simulation loop) for collar boundary look up
    void setTime(double t, double dt) {
        // std::cout << "Time " << t << " time step " << dt << "\n";
        this->spatialParams().setTime(t, dt);
        time_ = t;
        dt_ = dt;
    }

    /*!
     * writes the actual transpiration into a text file:
     * 0 time [s], 1 actual transpiration [kg/s], 2 potential transpiration [kg/s], 3 maximal transpiration [kg/s],
     * 4 collar pressure [Pa], 5 calculated actual transpiration [cm^3/day]
     *
     * 0 - 4 work only for neuman bc
     */
    void writeTranspirationRate(const SolutionVector& sol) {
        Scalar trans = this->transpiration(sol); // [cm3/day]
        file_at_ << neumannTime_ << ", " << actualTrans_ << ", " << potentialTrans_ << ", " << maxTrans_ << ", " << collarP_ << ", "
            << trans << ", "<< time_ << " , " << cL << "\n";
    }

    //! if true, sets bc to Dirichlet at criticalCollarPressure (false per default)
    void setCritical(bool b) {
        critical_ = b;
    }

    //! sets the criticalCollarPressure [Pa]
    void criticalCollarPressure(Scalar p) {
        criticalCollarPressure_ = p;
    }

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     *
     * The extrusion factor here makes extrudes the 1d line to a circular tube with
     * cross-section area pi*r^2.
     *
     * called by volumevariables (why there?), no compilation error if you remove it, just wrong results
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element &element, const SubControlVolume &scv, const ElementSolution& elemSol) const {
        const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
        const auto radius = this->spatialParams().radius(eIdx);
        return M_PI*radius*radius;
    }

    /*!
     * \brief Applies a vector of point sources. The point sources
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of PointSource s that contain
             source values for all phases and space positions.
     *
     * For this method, the \a values method of the point source
     * has to return the absolute mass rate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const {
        pointSources = couplingManager_->lowDimPointSources();
    }

    /*!
     * \brief Evaluate the point sources (added by addPointSources)
     *        for all phases within a given sub-control-volume.
     *
     * This is the method for the case where the point source is
     * solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param pointSource A single point source
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub-control volume within the element
     *
     * For this method, the \a values() method of the point sources returns
     * the absolute rate mass generated or annihilate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    template<class ElementVolumeVariables>
    void pointSource(PointSource& source, const Element &element, const FVElementGeometry& fvGeometry,
        const ElementVolumeVariables& elemVolVars, const SubControlVolume &scv) const
    {
        if (couplingManager_!=nullptr) {
            // compute source at every integration point
            const Scalar pressure3D = couplingManager_->bulkPriVars(source.id())[Indices::pressureIdx];
            const Scalar pressure1D = couplingManager_->lowDimPriVars(source.id())[Indices::pressureIdx];
            const auto lowDimElementIdx = couplingManager_->pointSourceData(source.id()).lowDimElementIdx();
            const Scalar kr = this->spatialParams().kr(lowDimElementIdx);
            const Scalar rootRadius = this->spatialParams().radius(lowDimElementIdx);
            // relative soil permeability
            const auto krel = 1.0;//this->couplingManager().relPermSoil(pressure3D);
            // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
            const auto density = 1000;
            PrimaryVariables sourceValue;
            sourceValue[contiH2OEqIdx] = 2* M_PI *krel*rootRadius * kr *(pressure3D - pressure1D)*density; // (kg/m/s)
            sourceValue[contiH2OEqIdx] *= source.quadratureWeight()*source.integrationElement();
            const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
            Scalar RootAge = this->spatialParams().age(eIdx); // RootAge in s
            RootAge /= (24. * 3600.); // days
            if (RootAge <= 1)  // if root segment age <= 1 day
            {
                Scalar tipP_ = pressure1D; // local tip pressure in kg/m/s^2
                if (abs(tipP_) >= abs(p0))
                {
                    if (useMoles) {
                    Msignal = 3.26e-16*(abs(tipP_) - abs(p0))*mi;     //(mol/s) 3.2523e-16 is production rate per dry mass in mol kg-1 Pa-1 s-1
                    }
                    else {
                    Msignal =  Msignal = 3.26e-16*(abs(tipP_) - abs(p0))*mi * MolarMass; // (kg/s)
                    }
                    sourceValue[transportABAEqIdx] = Msignal*source.quadratureWeight()*source.integrationElement(); // Msignal in (mol/s) TODO: ask for units
                }
                else
                {
                    sourceValue[transportABAEqIdx] = 0.;
                }
            }
            source = sourceValue;
        }

            else {
            source = 0;
        }
    }

    /**
     * for debugging
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
        std::cout << "Global integrated source (root): " << source[0] << " (kg/s) / "
            <<                           source[0]*3600*24*1000 << " (g/day)" << '\n';
    }

    //! Set the coupling manager
    void setCouplingManager(CouplingManager* cm) {
        couplingManager_ = cm;
    }

private:

    //! cm pressure head -> Pascal
    Scalar toPa_(Scalar ph) const {
        return pRef_ + ph / 100. * rho_ * g_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const {  // on root collar
        return globalPos[dimWorld - 1] > this->fvGridGeometry().bBoxMax()[dimWorld - 1] - eps_;
    }

    CouplingManager* couplingManager_ = nullptr;

    CRootBox::SoilLookUp* soil_;
    InputFileFunction collar_;
    size_t bcType_;
    double time_ = 0.;
    double dt_ = 0.;
    double criticalCollarPressure_ = -1.4e6;
    bool critical_ = false; // imposes dirichlet strong

    static constexpr Scalar g_ = 9.81; // cm / s^2
    static constexpr Scalar rho_ = 1.e3; // kg / m^3
    static constexpr Scalar pRef_ = 1.e5; // Pa
    static constexpr Scalar eps_ = 1e-6;

    std::ofstream file_at_; // file for actual transpiration
    mutable double neumannTime_ = 0;
    mutable double actualTrans_ = 0;
    mutable double potentialTrans_ = 0;
    mutable double maxTrans_ = 0.;
    mutable double collarP_ = 0.;

    std::map<std::string, std::vector<Scalar>> userData_;

    // chemical signalling variables
    Scalar alphaR = 0; // residual stomatal conductance, taken as 0
    bool cD = getParam<bool>("Control.cD"); // boolean variable: cD = 0 -> interaction between pressure and chemical regulation
    Scalar MolarMass = getParam<Scalar>("Component.MolarMass"); // (kg/mol)
    Scalar densityABA = getParam<Scalar>("Component.Density"); // (kg/m3)

    const Scalar p0 = toPa_(-4500); // cm -> Pa
    const Scalar p_crit = toPa_(-5500); // cm -> Pa
    const Scalar mi= 1.76e-7;  //dry mass = 140 kg_DM/m3, calculated using root tip = 1 cm length, and 0.02 cm radius
    const Scalar sH = 1.02e-6; // (Pa-1) from Huber et. al [2014]
    const Scalar sC = 5e+4; // (m^3/mol) from Huber et. al [2014]

    mutable Scalar cL = 0.;
    mutable Scalar Msignal = 0; // initial production rate of chemicals
    mutable Scalar alpha = 1; // initial stomatal conductance (=1) is stomata is fully open

};

/*
    //! compute the actual transpiration rate
    Scalar computeActualTranspirationRate(const SolutionVector& sol, const GridVariables& gridVars, bool verbose = true) const
    {
        NumEqVector transpirationRate(0.0);
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVars.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, sol);

            for (const auto& scvf : scvfs(fvGeometry))
                if (scvf.boundary())
                    transpirationRate += this->neumann(element, fvGeometry, elemVolVars, scvf)
 *scvf.area()*elemVolVars[scvf.insideScvIdx()].extrusionFactor();
        }
        if (verbose)
        {
            std::cout << "Actual transpiration rate:       " << transpirationRate << " (kg/s) / "
                << transpirationRate[0]*86400*1000 << " (g/day) / "
                << transpirationRate[0]/domainSize_[0]/domainSize_[1]*86400 << " (mm/day)\n"
                << "Potential transpiration rate:    " << potentialTranspirationRate() << " (kg/s) / "
                << potentialTranspirationRate()*86400*1000 << " (g/day) / "
                << potentialTranspirationRate()/domainSize_[0]/domainSize_[1]*86400 << " (mm/day)\n";
        }
        return transpirationRate[0];
    }
 */


} //end namespace Dumux

#endif
