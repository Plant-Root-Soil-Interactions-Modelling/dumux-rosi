// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef ROOTS_PROBLEM_HH
#define ROOTS_PROBLEM_HH

#include <math.h>
#include <map>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh> // import for MaterialLaw Schroeder
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>             // import for MaterialLaw Schroeder
#include <dumux/external/brent/brent.hpp>                       //T.S.: Brent algorithm to find roots of function


#include <dumux/porousmediumflow/problem.hh>

#include <dumux/growth/soillookup.hh>

#if DGF
#include "rootspatialparams_dgf.hh"
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
    using MaterialLaw = EffToAbsLaw<RegularizedVanGenuchten<Scalar>>;
    using MaterialLawParams = typename MaterialLaw::Params;


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
            if (name=="pSoil") {
                auto i0 = vMapper.subIndex(e, 0, 1);
                auto i1 = vMapper.subIndex(e, 1, 1);
                d = 0.5 * (sol[i1][0] + sol[i0][0]);
                d = 100. * (d - pRef_) / rho_ / g_;  // Pa -> cm
            }
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

    // Templates for brent-algorithm taken from https://stackoverflow.com/questions/51931479/conversion-between-stdfunctiondoubledouble-to-double-double
    // note: function-builder-base and function builder need to be adapteded with 2x const each

    // T.S: Function definition: integration by hand (calculate matric-flux-potential based on the currenct absolute pressure)
//    const Scalar pc_to_MFP(const auto& bulkElement, const Scalar pressure3D_pc, int n, const Scalar dx, const Scalar kc) const
    const Scalar pc_to_MFP(const auto& element, Scalar lower,  Scalar pressure3D_pc, int n) const
    {
        const auto& soilSpatialParams = couplingManager_->problem(Dune::index_constant<0>{}).spatialParams();
        const auto& params = soilSpatialParams.materialLawParams(element);
        Scalar kc = soilSpatialParams.hydraulicConductivity(element); // [m/s]
        std::function<double(double)> f = [=] (double x) { return MaterialLaw::krw(params, MaterialLaw::sw(params, x))*kc*86400; };
        return CPlantBox::Function::quad(f, pressure3D_pc, lower, n);

//        Scalar cumSum =0;
//        for (int i=0; i<n+1; i++)
//        {
//            const auto& soilSpatialParams = couplingManager_->problem(Dune::index_constant<0>{}).spatialParams();
//            MaterialLawParams params = soilSpatialParams.materialLawParams(bulkElement);
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
    void pointSource(PointSource& source, const Element &element, const FVElementGeometry& fvGeometry,
        const ElementVolumeVariables& elemVolVars, const SubControlVolume &scv) const
    {
        source = 0;
        if (couplingManager_!=nullptr) {
            auto eIdx = couplingManager_->pointSourceData(source.id()).lowDimElementIdx(); // const eingefuegt
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

                // SCHROEDER IMPLEMENTATION
                const Scalar gradients = getParam<Scalar>("Schroeder.gradients");
                // Switch in Input-File (gradients = 1 in [Soil.IC] enables Schroeder, gradients = 0 disables it)

                if (gradients == 1 && sourceValue > 0) {
                // switch to enable/disable Schroeder and check for macroscopic flow from soil to root (Schroeder only used if source value > 0, meaning uptake by root)

                    // STEP 0) GRAB NEEDED VARIABLES FROM RICHARDSPROBLEM + ID-DEFINITIONS
                    const auto& soilSpatialParams = couplingManager_->problem(Dune::index_constant<0>{}).spatialParams();
                    //get soil params into rootsproblem
                    const auto lowDimElementIdx = couplingManager_->pointSourceData(source.id()).lowDimElementIdx();
                    const auto bulkElementIdx = couplingManager_->pointSourceData(source.id()).bulkElementIdx();
                    // get bulkElement in rootsproblem (equivalent to pointSource element in richardsproblem)
                    const auto& soilProblem = couplingManager_->problem(Dune::index_constant<0>{});
                    const auto& gridGeometry = soilProblem.gridGeometry();
                    const auto& bulkElement = gridGeometry.element(bulkElementIdx);
                    MaterialLawParams params = soilSpatialParams.materialLawParams(bulkElement);

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
                    const Scalar MFP_soil = pc_to_MFP(bulkElement, lowBound_pc, pressure3D_pc, n);
                    // MFP of source-point soil voxel, call to integration function
                    //std::cout << " pointSource_root " << source.id() << "MFP_soil (rootsproblem)= " << MFP_soil << "\n";

                    // STEP 2) CALCULATE MFP_ROOT (according to non-stressed equation of Schroeder)

                    //calculation of r_out:
                    const Scalar segment_length = source.integrationElement(); // [m]
                    // length of point-source segment in voxel (cut at voxel boundaries)
                    const Scalar rootsystem_volume_inElement = couplingManager_->lowDimVolume(bulkElement);
                    // total volume of rootsystem in voxel
                    const Scalar cell_volume = bulkElement.geometry().volume();
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
                    const Scalar q_root = sourceValue / 1000 * 1000000 / 100 *86400 /(2 * M_PI *rootRadius*100); // here: -1*sourceValue exchanged with sourceValue
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
                    auto MFP_to_pressure3D = brent::funcLambda([=](double x) { return  pc_to_MFP(bulkElement, lowBound_pc, x, n) - MFP_nostress_root; });
                    double z = brent::zero(lowBound_pc, 0, tolerance, MFP_to_pressure3D);
                    const Scalar pressure3D_pc_new = z;
                    const Scalar pressure3D_new = -1*(z-pRef_);


                    // STEP 4) PASS NEW PRESSURE3D TO SOURCE-TERM
                    const Scalar pressure3D_s_new = MaterialLaw::sw(params, pressure3D_pc_new);
                    const Scalar krw = MaterialLaw::krw(params, pressure3D_s_new); // pass new pressure3D_s to calculate krw
                    const Scalar krw_scaled_to_rootRadius = krw / rootRadius;      // soil hydraulic conductivity scaled to root radius (radius used as proxy)
                    const Scalar kmin = std::min(kr, krw_scaled_to_rootRadius);    // minimum of conductivity defines sourceValue
                    sourceValue = 2 * M_PI *krel* rootRadius * kmin * (pressure3D_new - pressure1D)*density;    //* kr exchanged with kmin
                    source = sourceValue*source.quadratureWeight()*source.integrationElement();

                    if(sourceValue < 0) {
                        //discards Schroeder if it leads to inversion of flow (e.g. macroscopic-flow soil => root, schroeder-flow root => soil. Jan thinks this clause may be exluded)
                        source = 0;
                        }

                    //Prints (can be enabled / disabled via Schroeder print in coupled input-file)
                    const Scalar print = getParam<Scalar>("Schroeder.print");
                    if (print == 1) {
                                    std::cout << "rootsproblem sourceID_root:" << source.id() << "\n" << " MFP_soil = " << MFP_soil <<  "   MFP_no_stress_root= " << MFP_nostress_root << "   deltaMFP= "
                                    << MFP_soil-MFP_nostress_root << "  soil_r_out= " << r_out << "   root_radius= " << rootRadius*100 << "   q_root= " << q_root << "\n"
                                    << " pressure3D_head(soil)= " << toHead_(pressure3D) << "   pressure3D_head(root_surface)= " << toHead_(pressure3D_new) << "\n"
                                    << " pressure3D_Pa(soil)  = " << pressure3D << "  pressure3D_Pa(root_surface) = " << pressure3D_new << "\n" 
                                    << " pressure1D_Pa(root)  = " << pressure1D << "  pressure1D_head(root)       = " << toHead_(pressure1D) << "\n" << "\n";
                        }
                    }
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
//Tobias: inserted toPa_ and toHead here
    //! cm pressure head -> Pascal
	Scalar toPa_(Scalar ph) const {
		return pRef_ + ph / 100. * rho_ * g_;
	}

	//! Pascal -> cm pressure head
	Scalar toHead_(Scalar p) const {
		return (p - pRef_) * 100. / rho_ / g_;
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
