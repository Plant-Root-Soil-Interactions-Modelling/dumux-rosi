// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef ROOTS_PROBLEM_HH
#define ROOTS_PROBLEM_HH

#include <math.h>
#include <map>

#include <dumux/porousmediumflow/problem.hh>

#include <dumux/growth/soillookup.hh>

#if DGF
#include "rootspatialparams_cavitation_dgf.hh"
#endif
#if ROOTBOX
#include "rootspatialparams_rb.hh"
#endif

namespace Dumux {

/*!
 * Root Doussan Model
 *
 * with optional coupling to a soil model
 */
template<class TypeTag>
class RootsProblem: public PorousMediumFlowProblem<TypeTag> {

    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename FVGridGeometry::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
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

    enum {
        conti0EqIdx = Indices::conti0EqIdx, // indices of the primary variables
        pressureIdx = Indices::pressureIdx
    };
    enum {
        isBox = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethods::Box
    };
    enum {
        bcDirichlet = 0,
        bcNeumann = 1
    };

    static const int dimWorld = GridView::dimensionworld;

public:

    //! Constructor
    RootsProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry): ParentType(fvGridGeometry) {

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
        criticalCollarPressure_ = toPa_(getParam<double>("RootSystem.Collar.CritCollarP", -1.5e4));  // cm -> Pa
        
        // Optionally give cavitation parameter
        double b = getParam<double>("Control.b", 1.e16); // cm pressure head
        double c = getParam<double>("Control.c", 1);
        this->spatialParams().setParameter(b,c);

    }

    //! Destructor - close transpiration file
    virtual ~RootsProblem() {
        delete soil_;
        std::cout << "closing file \n" << std::flush;
        file_at_.close();
    }

    //! evaluates user defined data for vtk fields
    void userData(std::string name, const SolutionVector& sol) {
        const auto& gridView = this->gridGeometry().gridView();
        userData_[name] = std::vector<Scalar>(gridView.size(0));
        auto eMapper = this->gridGeometry().elementMapper();
        auto vMapper = this->gridGeometry().vertexMapper();
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
                d = initialAtPos(e.geometry().center()); // Pa
                d = 100. * (d - pRef_) / rho_ / g_;  // Pa -> cm
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
                d =  2 * a * M_PI * length* kr * (soil(p) - (sol[i1] + sol[i0]) / 2); // m^3 / s
                d = 24.*3600*1.e6*d; // [m^3/s] -> [cm^3/day]
            }
            if (name=="axialFlux") {
                auto geo = e.geometry();
                auto length = geo.volume();
                auto kx = this->spatialParams().kx(eIdx);
                auto i0 = vMapper.subIndex(e, 0, 1);
                auto i1 = vMapper.subIndex(e, 1, 1);
                d = kx * ((sol[i1] - sol[i0]) / length - rho_ * g_); // m^3 / s
                d = 24.*3600*1.e6*d; // [m^3/s] -> [cm^3/day]
            }
            if (name=="pXylem") { // this is the segment xylem pressure as 0.5*(node1 + node2)
                auto i0 = vMapper.subIndex(e, 0, 1);
                auto i1 = vMapper.subIndex(e, 1, 1);
                d = 0.5 * (sol[i1][0] + sol[i0][0]);
                d = 100. * (d - pRef_) / rho_ / g_;  // Pa -> cm
            }
            // pSoil is set in postTimeStep
            userData_[name][eIdx] = d;
        }
    }

    //! vtk fields call back functions (initialize with userData(name) )
    std::vector<Scalar>& radialFlux() { return userData_["radialFlux"]; }
    std::vector<Scalar>& axialFlux() { return userData_["axialFlux"]; }
    std::vector<Scalar>& kr() { return userData_["kr"]; }
    std::vector<Scalar>& kx() { return userData_["kx"]; }
    std::vector<Scalar>& age() { return userData_["age"]; }
    std::vector<Scalar>& order() { return userData_["order"]; }
    std::vector<Scalar>& radius() { return userData_["radius"]; }
    std::vector<Scalar>& initialPressure() { return userData_["initialPressure"]; }
    std::vector<Scalar>& id() { return userData_["id"]; }
    std::vector<Scalar>& p() { return userData_["pSoil"]; }
    std::vector<Scalar>& pXylem() { return userData_["pXylem"]; }


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
            return criticalCollarPressure_;
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
            auto eIdx = this->gridGeometry().elementMapper().index(element);
            double kx = this->spatialParams().kx(eIdx);
            auto dist = (globalPos - fvGeometry.scv(scvf.insideScvIdx()).center()).two_norm();
            double criticalTranspiration = volVars.density(0) * kx * (p - criticalCollarPressure_) / dist; // [kg/s]
            double actTrans = std::min(collar_.f(time_), criticalTranspiration);
            actTrans /= volVars.extrusionFactor(); // [kg/s] -> [kg/(s*m^2)]
            return NumEqVector(actTrans);
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
        values[conti0EqIdx] = 0.;
        if (couplingManager_==nullptr) {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            auto params = this->spatialParams();
            Scalar age = params.age(eIdx); // root age (s)

            Scalar a = params.radius(eIdx); // root radius (m)
            Scalar kr = params.kr(eIdx); //  radial conductivity (m^2 s / kg)
            Scalar phx;
            if (isBox) { // dumux
                phx = elemVolVars[scv.localDofIndex()].pressure(); // kg/m/s^2
            } else {
                phx = elemVolVars[scv.dofIndex()].pressure(); // kg/m/s^2
            }
            Scalar phs = soil(scv.center()); // kg/m/s^2
            values[conti0EqIdx] = kr * 2 * a * M_PI * (phs - phx); // m^3/s
            values[conti0EqIdx] /= (a * a * M_PI); // 1/s
            values[conti0EqIdx] *= rho_; // (kg/s/m^3)
        }
        return values;
    }

    /*!
     * \brief Applies a vector of point sources. The point sources
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of PointSource s that contain
     * source values for all phases and space positions.
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
        source = 0;
        if (couplingManager_!=nullptr) {
            auto eIdx = couplingManager_->pointSourceData(source.id()).lowDimElementIdx();
            Scalar age = this->spatialParams().age(eIdx);
            if (age>0) {
                // compute source at every integration point
                Scalar pressure3D = couplingManager_->bulkPriVars(source.id())[Indices::pressureIdx];
                Scalar pressure1D = couplingManager_->lowDimPriVars(source.id())[Indices::pressureIdx];
                Scalar kr = this->spatialParams().kr(eIdx);
                Scalar rootRadius = this->spatialParams().radius(eIdx);
                // relative soil permeability
                auto krel = 1.0;//this->couplingManager().relPermSoil(pressure3D);
                // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
                auto density = 1000;
                Scalar sourceValue = 2* M_PI *krel*rootRadius * kr *(pressure3D - pressure1D)*density;
                source = sourceValue*source.quadratureWeight()*source.integrationElement();
            }
        }
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& p) const {
        return PrimaryVariables(soil(p)); // soil(p)
    }

    /**
     * deletes old soil, sets new soil, takes ownership
     * TODO smart pointers
     */
    void setSoil(CPlantBox::SoilLookUp* s) {
        delete soil_;
        soil_ = s;
        std::cout << "setSoil(...): manually changed soil to " << s->toString() << "\n";
    }

    //! soil pressure (called by initial, and source term)
    Scalar soil(const GlobalPosition& p) const {
        auto p2 = CPlantBox::Vector3d(p[0] * 100, p[1] * 100, p[2] * 100); // m -> cm
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
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        const auto radius = this->spatialParams().radius(eIdx);
        return M_PI*radius*radius;
    }

    /**
     * Sets transpiration according to the last solution
     */
    void postTimeStep(const SolutionVector& sol, const GridVariables& gridVars) {

        NumEqVector source(0.0);
        for (const auto& e :elements(this->gridGeometry().gridView())) {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(e);
            auto elemVolVars = localView(gridVars.curGridVolVars());
            elemVolVars.bindElement(e, fvGeometry, sol);
            for (const auto& scvf :scvfs(fvGeometry)) { // evaluate root collar sub control faces
                if (onUpperBoundary_(scvf.center())) { // root collar
                    auto& volVars = elemVolVars[scvf.insideScvIdx()];
                    double p = volVars.pressure();
                    auto eIdx = this->gridGeometry().elementMapper().index(e);
                    double kx = this->spatialParams().kx(eIdx);
                    auto dist = (scvf.center() - fvGeometry.scv(scvf.insideScvIdx()).center()).two_norm();
                    double criticalTranspiration = volVars.density(0) * kx * (p - criticalCollarPressure_) / dist; // [kg/s]
                    potentialTrans_ = collar_.f(time_); // [kg/s]
                    neumannTime_ = time_; // [s]
                    actualTrans_ = std::min(potentialTrans_, criticalTranspiration);// actual transpiration rate [kg/s]
                    maxTrans_ = criticalTranspiration; // [kg/s]
                    collarP_ = p; // [Pa]
                }
            }
            if (couplingManager_!=nullptr) {
                const auto& gridView = this->gridGeometry().gridView();
                if (userData_.count("pSoil")) {
                    userData_["pSoil"].resize(gridView.size(0));
                } else {
                    userData_["pSoil"] = std::vector<Scalar>(gridView.size(0));
                }
                for (const auto& scv : scvs(fvGeometry)) {
                    scvSoilMatricPotential_(e, fvGeometry, elemVolVars, scv); //
                }
            }
        }
    }

    /*!
     * Writes the actual transpiration into a text file. Call postTimeStep before using it.
     *
     * 0 time [s], 1 actual transpiration [kg/s], 2 potential transpiration [kg/s], 3 maximal transpiration [kg/s],
     * 4 collar pressure [Pa], 5 - (0.) 6 simtime [s]
     */
    void writeTranspirationRate() {
        file_at_ << neumannTime_ << ", " << actualTrans_ << ", " << potentialTrans_ << ", " << maxTrans_ << ", " << collarP_ << ", "
            << 0. << ", "<< time_ << "\n";
    }

    /*!
     * much copy paste from FVProblem::scvPointSources
     * stores the matric potential of the soil cell in userData_
     */
    void scvSoilMatricPotential_(const Element &element, const FVElementGeometry& fvGeometry, const ElementVolumeVariables& elemVolVars, const SubControlVolume &scv) {
        auto scvIdx = scv.indexInElement();
        auto key = std::make_pair(this->gridGeometry().elementMapper().index(element), scvIdx);
        if (this->pointSourceMap().count(key)) {
            auto pointSources = this->pointSourceMap().at(key);
            for (auto&& pointSource : pointSources) {
                pointSource.update(this->asImp_(), element, fvGeometry, elemVolVars, scv);
                    auto eIdx = couplingManager_->pointSourceData(pointSource.id()).lowDimElementIdx();
                    Scalar age = this->spatialParams().age(eIdx);
                    if (age>0) {
                        // compute source at every integration point
                        Scalar pressure3D = couplingManager_->bulkPriVars(pointSource.id())[Indices::pressureIdx];
                        userData_["pSoil"].at(eIdx) =  100. * (pressure3D - pRef_) / rho_ / g_;  // Pa -> cm
                    }


            }
        }
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
        std::cout << "Global integrated source (root): " << source << " (kg/s) / " <<  source*3600*24*1000 << " (g/day)" << '\n';
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
        return globalPos[dimWorld - 1] > this->gridGeometry().bBoxMax()[dimWorld - 1] - eps_;
    }

    CouplingManager* couplingManager_ = nullptr;

    CPlantBox::SoilLookUp* soil_;
    InputFileFunction collar_;
    size_t bcType_;
    double time_ = 0.;
    double dt_ = 0.;
    double criticalCollarPressure_ = -1.4e6; // -15290 cm ??
    bool critical_ = false; // imposes dirichlet strong

    static constexpr Scalar g_ = 9.81; // cm / s^2
    static constexpr Scalar rho_ = 1.e3; // kg / m^3
    static constexpr Scalar pRef_ = 1.e5; // Pa
    static constexpr Scalar eps_ = 1e-6;

    std::ofstream file_at_; // file for actual transpiration
    double neumannTime_ = 0;
    double actualTrans_ = 0;
    double potentialTrans_ = 0;
    double maxTrans_ = 0.;
    double collarP_ = 0.;

    std::map<std::string, std::vector<Scalar>> userData_;

};

} //end namespace Dumux

#endif
