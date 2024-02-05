// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef ROOTS_PROBLEM_1P2C_HH
#define ROOTS_PROBLEM_1P2C_HH

#include <math.h>
#include <cmath>        // std::abs
#include <map>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/growth/soillookup.hh>

//// maybe we will need it for advective flux approx
//#include <dumux/discretization/cctpfa.hh>
//#include <dumux/discretization/ccmpfa.hh>
//#include <dumux/discretization/box.hh>
//#include <dumux/discretization/evalsolution.hh>
//#include <dumux/discretization/evalgradients.hh>

#if DGF
#include "../roots_1p/rootspatialparams_dgf.hh"
#endif
#if ROOTBOX
#include "../roots_1p/rootspatialparams_rb.hh" // TODO
#endif

namespace Dumux {

/*!
 * Root system water movement and transpiration driven by hormones.
 *
 * similar Huber et al. 2014
 *
 */
template<class TypeTag>
class Roots1P2CProblem: public PorousMediumFlowProblem<TypeTag> {

    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename FVGridGeometry::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using CouplingManager= GetPropType<TypeTag, Properties::CouplingManager>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    enum {

        pressureIdx = 0, // indices of primary variables
        h2OIdx = 0, // component indices
        soluteIdx = 1,
        conti0EqIdx = 0,  // indices of the equations
        transportEqIdx = 1, //

        isBox = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethods::Box,

        bcDirichlet = 0,
        bcNeumann = 1
    };

    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static const int dimWorld = GridView::dimensionworld;

public:

    //! Constructor
    Roots1P2CProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry): ParentType(fvGridGeometry) {

        // initialize fluid system
        FluidSystem::init();

        // stating in the console whether mole or mass fractions are used
        if(useMoles) {
            std::cout<<"Roots1P2CProblem uses mole fractions, CURRENTLY WRONG, CARE!"<<std::endl;
        } else {
            std::cout<<"Roots1P2CProblem problem uses mass fractions"<<std::endl;
        }
        std::cout<<"contiH2OEqIdx "<< conti0EqIdx <<"transportEqIdx "<< transportEqIdx << std::endl;

        // for the uncoupled case, a static soil is created
        InputFileFunction sf = InputFileFunction("Soil.IC", "P", "Z", 0.); // [cm]([m])
        sf.setFunctionScale(1.e-2 * rho_ * g_ ); // [cm] -> [Pa], don't forget to add pRef_
        soil_ = new GrowthModule::SoilLookUpTable(sf); // sf is copied by default copy constructor

        if (Dumux::hasParam("RootSystem.Collar.P")) {
            collar_ = InputFileFunction("RootSystem.Collar", "P", "PT", -1.e4); // [cm]([day])
            collar_.setVariableScale(1./(24.*3600)); // [s] -> [day]
            collar_.setFunctionScale(1.e-2 * rho_ * g_); // [cm] -> [Pa], don't forget to add pRef_
            bcType_ = bcDirichlet;
            std::cout << "RootsOnePTwoCProblem(): Warning: transpiration should be predetermined, hormone model is not implemented for Dirichlet BC yet. \n";
        } else { // ".File " gives transpiration table
            collar_ = InputFileFunction("RootSystem.Collar", "Transpiration", "TranspirationT", 0.02);  // [kg/day]([day])
            collar_.setVariableScale(1./(24.*3600)); // [s] -> [day]
            collar_.setFunctionScale(1./(24.*3600)); // [kg/day] -> [kg/s]
            bcType_ = bcNeumann;
        }
        file_at_.open(this->name() + "_actual_transpiration.txt");

        grow_ = getParam<bool>("RootSystem.Grid.Grow", false); // for mimicing growth using root age

        leafVolume_ = InputFileFunction("RootSystem.Leaf", "Volume", "VolumeT", 1.); // [cm^3]([day])
        leafVolume_.setVariableScale(1./(24.*3600)); // [s] -> [day]
        leafVolume_.setFunctionScale(1.e-6); // [cm^3] -> [m^3]

        molarMass = getParam<double>("Component.MolarMass", 0.26432); // (kg/mol)

        // Optionally, give crit values
        critPCollarDirichlet_ = toPa_(getParam<double>("Control.CritCollarP", -1.5e4));  // cm -> Pa

        // Uptake params
        vMax_ =  getParam<Scalar>("RootSystem.Uptake.Vmax", 6.2e-11/(10.*24.*3600.))*10*24*3600; // Michaelis Menten Parameter [g/cm2/day] -> [kg m-2 s-1]
        km_ = getParam<Scalar>("RootSystem.Uptake.Km", 3.1e-9/1000.)*1000.;  // Michaelis Menten Parameter  [g/cm3] - > [kg m-3]
        sigma_ = getParam<Scalar>("RootSystem.Uptake.ActiveTransport", 0.); // 1 for active transport, 0 for passive

        std::cout << "***********************************************\n";
        std::cout << "leafVolume "<< leafVolume_.f(0.) << ", grow " << grow_ << "\n";
        std::cout << "critPCollarDirichlet " << critPCollarDirichlet_ << "\n";
        std::cout << "***********************************************\n";
    }

    //! Destructor - close transpiration file
    virtual ~Roots1P2CProblem() {
        delete soil_;
        std::cout << "closing file \n" << std::flush;
        file_at_.close();
    }

    //! evaluates user defined data for vtk fields
    void userData(std::string name, const SolutionVector& sol) { // , const GridVariables& gridVars
        const auto& gridView = this->gridGeometry().gridView();
        userData_[name] = std::vector<Scalar>(gridView.size(0));
        auto eMapper = this->gridGeometry().elementMapper();
        auto vMapper = this->gridGeometry().vertexMapper();
        for (const auto& e :elements(this->gridGeometry().gridView())) {
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
                d = initialAtPos(e.geometry().center())[pressureIdx]; // Pa
                d = 100. * (d - pRef_) / rho_ / g_;  // Pa -> cm
            }
            if (name=="radialFlux") { // two point approximation
                auto geo = e.geometry();
                auto length = geo.volume();
                auto kr = this->spatialParams().kr(eIdx);
                auto a = this->spatialParams().radius(eIdx);
                auto i0 = vMapper.subIndex(e, 0, 1);
                auto i1 = vMapper.subIndex(e, 1, 1);
                auto p = geo.center();
                // kr [m /Pa/s]
                d =  2 * a * M_PI * length* kr * (soil(p) - (sol[i1][pressureIdx] + sol[i0][pressureIdx]) / 2); // m^3 / s
                d = 24.*3600*1.e6*d; // [m^3/s] -> [cm^3/day]
            }
            if (name=="axialFlux") {
                auto geo = e.geometry();
                auto length = geo.volume();
                auto kx = this->spatialParams().kx(eIdx);
                auto i0 = vMapper.subIndex(e, 0, 1);
                auto i1 = vMapper.subIndex(e, 1, 1);
                d = kx * ((sol[i1][pressureIdx] - sol[i0][pressureIdx]) / length - rho_ * g_); // m^3 / s
                d = 24.*3600*1.e6*d; // [m^3/s] -> [cm^3/day]
            }
            if (name=="pSoil") {
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
    std::vector<Scalar>& p() { return userData_["pSoil"]; }

    //! calculates transpiration, as the sum of radial fluxes [cm^3/day]
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

    /**
     * The buffer power for a scv for a volVar (linear in this implementation), equals $\rho_b K_d$ in Eqn (4) in phosphate draft
     *
     * used by my the modified localresidual.hh (see dumux-rosi/dumux/porousmediumflow/compositional)
     */
    Scalar getBufferPower(const SubControlVolume& scv, const VolumeVariables& volVars) const {
        return 0.;
    }
    
    /**
	 * The buffer power for a scv for a volVar (linear in this implementation), equals $\rho_b K_d$ in Eqn (4) in phosphate draft
	 *
	 * used by my the modified localresidual.hh (see dumux-rosi/dumux/porousmediumflow/compositional)
	 */
	Scalar bufferPower(const SubControlVolume& scv, const VolumeVariables& volVars, int compIdx = 0) const {
        return 0.;
	}

    /*!
     * \brief Evaluate the initial value for a control volume.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& p) const {
        PrimaryVariables priVars;
        priVars[pressureIdx] = soil(p); //
        priVars[soluteIdx] = 0.0;  // todo (if we want some initial hormone state)
        return priVars;
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &pos) const {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann(); // default
        return bcTypes;
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

        NumEqVector flux;
        const auto globalPos = scvf.center();
        if (onUpperBoundary_(globalPos)) { // root collar

            auto eIdx = this->gridGeometry().elementMapper().index(element);
            auto& volVars = elemVolVars[scvf.insideScvIdx()];
            double p = volVars.pressure(); // pressure at the root collar
            double kx = this->spatialParams().kx(eIdx);
            auto dist = (globalPos - fvGeometry.scv(scvf.insideScvIdx()).center()).two_norm();
            double criticalTranspiration = volVars.density(0) * kx * (p - critPCollarDirichlet_) / dist; // [kg/s]
            double potentialTrans = collar_.f(time_); // [kg/s]
            double actTrans = std::min(potentialTrans, criticalTranspiration);// actual transpiration rate [kg/s]
            flux[conti0EqIdx] = actTrans/volVars.extrusionFactor(); // [kg/s] -> [kg/(s*m^2)];

            double fraction = useMoles ? volVars.moleFraction(0, soluteIdx) : volVars.massFraction(0, soluteIdx);
            flux[transportEqIdx] = flux[conti0EqIdx] * fraction; // [kg_aba/(s*m^2)],  convective outflow BC

        } else { // root tip

            flux[conti0EqIdx] = 0.;
            flux[transportEqIdx] = 0.; // solutes cannot leave from the tip

        }
        return flux;
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
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            Scalar a = params.radius(eIdx); // root radius (m)
            Scalar kr = params.kr(eIdx); //  radial conductivity (m^2 s / kg)
            Scalar tipP_;
            if (isBox) { // would elemVolVars[scv] work? not sure...
                tipP_ = elemVolVars[scv.localDofIndex()].pressure(); // kg/m/s^2
            } else {
                tipP_ = elemVolVars[scv.dofIndex()].pressure(); // kg/m/s^2
            }
            Scalar phs = soil(scv.center()); // kg/m/s^2
            values[conti0EqIdx] = kr * 2 * a * M_PI * (phs - tipP_); // m^3/s
            values[conti0EqIdx] /= (a * a * M_PI); // 1/s
            values[conti0EqIdx] *= rho_; // (kg/s/m^3)

            Scalar rootAge = this->spatialParams().age(eIdx) / (24. * 3600.); // days
            if (!grow_) { // for static root system, static root tips should not age
                rootAge -= time_ / (24. * 3600.);
            }
            values[transportEqIdx] = 0.; // TODO

        } else { // couplingManager_ is set, pointSources is used instead of source(...)
            values[conti0EqIdx] = 0;
            values[transportEqIdx] = 0.;
        }
        return values;
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
     * For this method, the \a values() method of the point sources returns
     * the absolute rate mass generated or annihilate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    template<class ElementVolumeVariables>
    void pointSource(PointSource& source, const Element &element, const FVElementGeometry& fvGeometry,
        const ElementVolumeVariables& elemVolVars, const SubControlVolume &scv) const
    {
        PrimaryVariables sourceValue(0.);

        if (couplingManager_!=nullptr) { // compute source at every integration point

            Scalar soilP = couplingManager_->bulkPriVars(source.id())[pressureIdx];
            Scalar tipP = couplingManager_->lowDimPriVars(source.id())[pressureIdx];

            auto eIdx = this->gridGeometry().elementMapper().index(element);
            Scalar kr = this->spatialParams().kr(eIdx); //  [m^2 s/kg]
            Scalar rootRadius = this->spatialParams().radius(eIdx); // [m]
            double density = 1000.; // [kg /m^3]
            sourceValue[conti0EqIdx] = 2* M_PI *rootRadius * kr *(soilP - tipP)*density; // [kg/m/s] TODO check signs, must equal the richardsnc::pointSource
            sourceValue[conti0EqIdx] *= source.quadratureWeight()*source.integrationElement(); // [kg /s]

            Scalar tipC = couplingManager_ ->lowDimPriVars(source.id())[soluteIdx]; // units [1], fraction
            Scalar soilC = couplingManager_ ->bulkPriVars(source.id())[soluteIdx]; // units [1], fraction

            Scalar passiveUptake;
            if (sourceValue[conti0EqIdx]>0) {      // flow from root to soil TODO check signs, must equal the richardsnc::pointSource
                passiveUptake = 2 * M_PI * rootRadius * kr * (tipP - soilP) * density * tipC;
            } else {
                passiveUptake = 2 * M_PI * rootRadius * kr * (tipP - soilP) * density * soilC;
            }
            // Active uptake based on Michaelis Menten
            Scalar activeUptake = -2 * M_PI * rootRadius * vMax_ * soilC * density/(km_ + soilC * density);

            // choose active or passive
            sourceValue[transportEqIdx] = 0.; // (sigma_*activeUptake + (1.-sigma_)*passiveUptake) *source.quadratureWeight()*source.integrationElement();

            source = sourceValue;

        } else { // should not happen...

            std::cout << "RootsOnePTwoCProblem::pointSource(): Coupling manager must be set in main file \n";
            source = sourceValue;

        }
    }

    /**
     * deletes old soil, sets new soil, takes ownership
     * TODO use smart pointers
     */
    void setSoil(CPlantBox::SoilLookUp* s) {
        delete soil_;
        soil_ = s;
        std::cout << "setSoil(...): manually changed soil to " << s->toString() << "\n";
    }

    //! soil pressure (called by initial, and source term)
    Scalar soil(const GlobalPosition& p) const {
        return pRef_+soil_->getValue(CPlantBox::Vector3d(p[0] * 100, p[1] * 100, p[2] * 100)); // m -> cm;
    }

    //! sets the current simulation time [s] (within the simulation loop) for collar boundary look up
    void setTime(double t, double dt) {
        this->spatialParams().setTime(t, dt);
        time_ = t;
        dt_ = dt;
    }

    //! sets the criticalCollarPressure [Pa]
    void criticalCollarPressure(Scalar p) {
        critPCollarDirichlet_ = p;
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
        auto eIdx = this->gridGeometry().elementMapper().index(element);
        double radius = this->spatialParams().radius(eIdx);
        return M_PI*radius*radius;
    }

    /**
     * Sets the cumulative outflow according to the last solution
     */
    void postTimeStep(const SolutionVector& sol, const GridVariables& gridVars) {

        NumEqVector source(0.0);
        mLRate_ = 0.;
        for (const auto& e :elements(this->gridGeometry().gridView())) {

            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(e);
            auto elemVolVars = localView(gridVars.curGridVolVars());
            elemVolVars.bindElement(e, fvGeometry, sol);

            for (const auto& scv :scvs(fvGeometry)) {
                if (couplingManager_!=nullptr) {
                    auto pointSources = this->scvPointSources(e, fvGeometry, elemVolVars, scv); // kg s-1
                    pointSources *= scv.volume()*elemVolVars[scv].extrusionFactor();
                    source += pointSources;
                } else {
                    auto sources = this->source(e, fvGeometry, elemVolVars, scv); // kg m-3 s-1
                    source += sources;
                }
            }

            for (const auto& scvf :scvfs(fvGeometry)) { // evaluate root collar sub control faces

                if (onUpperBoundary_(scvf.center())) { // root collar

                    auto& volVars = elemVolVars[scvf.insideScvIdx()];
                    double p = volVars.pressure();
                    auto eIdx = this->gridGeometry().elementMapper().index(e);
                    double kx = this->spatialParams().kx(eIdx);
                    auto dist = (scvf.center() - fvGeometry.scv(scvf.insideScvIdx()).center()).two_norm();
                    double criticalTranspiration = volVars.density(0) * kx * (p - critPCollarDirichlet_) / dist; // [kg/s]
                    potentialTrans_ = collar_.f(time_); // [kg/s]
                    double cL = mL_ / leafVolume_.f(time_); // mL from last time step [kg], leaf volume at simulation time [m^3]
                    double fraction = useMoles ? volVars.moleFraction(0, soluteIdx) : volVars.massFraction(0, soluteIdx); // [-]
                    neumannTime_ = time_; // [s]
                    actualTrans_ = std::min(potentialTrans_, criticalTranspiration);// actual transpiration rate [kg/s]
                    maxTrans_ = criticalTranspiration; // [kg/s]
                    collarP_ = p; // [Pa]
                    mLRate_ = actualTrans_*fraction; // [kg/s]
                    std::cout << " 1.e6* { cL "<< 1.e6*cL << " mL_ "<< 1.e6*mL_ << " leafVolume " << 1.e6*leafVolume_.f(time_) <<
                        " fraction " << 1.e6*fraction << " }\n\n";
                }
            }
        }

        mL_ += mLRate_*dt_; // integrate rate with old time step, we might need additional decay rate
        mRootRate_ = source[transportEqIdx]; // kg/s
        mRoot_ +=  mRootRate_*dt_; //kg
    }

    /*!
     * Writes the actual transpiration into a text file. Call postTimeStep before using it.
     *
     * 0 time [s], 1 actual transpiration [kg/s], 2 potential transpiration [kg/s], 3 maximal transpiration [kg/s],
     * 4 collar pressure [Pa], 5 - (0.), 6 simtime [s], 7 hormone leaf mass [kg],
     * 8 hormone collar flow rate [kg/s], 9 hormone root system mass [kg] , 10 hormone source rate [kg/s]
     */
    void writeTranspirationRate() {
        file_at_ << neumannTime_ << ", " << actualTrans_ << ", " << potentialTrans_ << ", " << maxTrans_ << ", " << collarP_ << ", "
            << 0. << ", "<< time_ << " , " << mL_ << ", "<< mLRate_  << ", " << mRoot_ << ", " << mRootRate_ << "\n";
    }

    /**
     * for debugging
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
        std::cout << "Source: water " << source[conti0EqIdx]*3600*24*1000 << " (g/day),\n"
            << " solutes " << source[transportEqIdx]*3600*24*1000 << " (g/day)" << '\n';
    }

    //! Set the coupling manager
    void setCouplingManager(CouplingManager* cm) {
        couplingManager_ = cm;
    }

protected:

    //! cm pressure head -> Pascal
    Scalar toPa_(Scalar ph) const {
        return pRef_ + ph / 100. * rho_ * g_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const {  // on root collar
        return globalPos[dimWorld - 1] > this->gridGeometry().bBoxMax()[dimWorld - 1] - eps_;
    }

    CouplingManager* couplingManager_ = nullptr;

    CPlantBox::SoilLookUp* soil_;
    InputFileFunction collar_;
    size_t bcType_;
    double time_ = 0.; // s
    double dt_ = 0.; // s

    static constexpr Scalar g_ = 9.81; // cm / s^2
    static constexpr Scalar rho_ = 1.e3; //1.e3; // kg / m^3
    static constexpr Scalar pRef_ = 1.e5; // Pa
    static constexpr Scalar eps_ = 1e-6;

    std::ofstream file_at_; // file for actual transpiration
    double neumannTime_ = 0;
    double actualTrans_ = 0;
    double potentialTrans_ = 0;
    double maxTrans_ = 0.;
    double collarP_ = 0.;
    double mRootRate_ = 0.; // [kg/s]

    std::map<std::string, std::vector<Scalar>> userData_; // for vtk fields

    bool grow_; // indicates if root segments age, or not

    double critPCollarDirichlet_ = -1.4e6; // -1.4e6;
    double mL_ = 0.; // (kg) mass of hormones in the leaf
    double mRoot_ = 0.; // (kg) mass of hormones in the root system
    double mLRate_ = 0.; // (kg / s) production rate of hormones flowing into the leaf volume
    InputFileFunction leafVolume_; // (m^3)

    Scalar molarMass; // (kg/mol)

    Scalar vMax_; // Michaelis Menten Parameter [kg m-2 s-1]
    Scalar km_; // Michaelis Menten Parameter  [kg m-3]
    Scalar sigma_; // 1 for passive transport, 0 for active transport

};

} //end namespace Dumux

#endif
