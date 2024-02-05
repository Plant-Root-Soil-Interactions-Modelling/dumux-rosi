// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef DUMUX_ROOT_SPATIALPARAMS_DGF_HH
#define DUMUX_ROOT_SPATIALPARAMS_DGF_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/material/spatialparams/fv1p.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/io/inputfilefunction.hh> // in dumux-rosi

namespace Dumux {

/*!
 * \brief Root spatial parameters class for grid files.
 *
 * use initParameters to initialize the class with data from the grid file
 *
 */
template<class FVGridGeometry, class Scalar>
class RootSpatialParamsDGF
    : public FVSpatialParamsOneP<FVGridGeometry, Scalar, RootSpatialParamsDGF<FVGridGeometry, Scalar>> {
    using ThisType = RootSpatialParamsDGF<FVGridGeometry, Scalar>;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar, ThisType>;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Water = Components::SimpleH2O<Scalar>;

public:

    using PermeabilityType = Scalar; // export permeability type

    RootSpatialParamsDGF(std::shared_ptr<const FVGridGeometry> fvGridGeometry):
        ParentType(fvGridGeometry) {
        // DGF specific,
        // new IBG-3 default: order, brnID, surf [cm2], length [cm], radius [cm], kz [cm4 hPa-1 d-1], kr [cm hPa-1 d-1], emergence time [d], subType, organType
        orderIdx_ = Dumux::getParam<int>("RootSystem.Grid.orderIdx", 0);
        radiusIdx_ = Dumux::getParam<int>("RootSystem.Grid.radiusIdx", 4);
        ctIdx_ = Dumux::getParam<int>("RootSystem.Grid.ctIdx", 7);
        krIdx_ = Dumux::getParam<int>("RootSystem.Grid.krIdx", 6);
        kxIdx_ = Dumux::getParam<int>("RootSystem.Grid.kxIdx", 5);
        idIdx_ = Dumux::getParam<int>("RootSystem.Grid.idIdx", 1);

        // DGF specific
        time0_ = Dumux::getParam<double>("RootSystem.Grid.InitialT", 14)*24*3600; // days -> s
        kr_ = InputFileFunction("RootSystem.Conductivity", "Kr", "KrAge", krIdx_, orderIdx_); // [cm/hPa/day] ([day])
        kr_.setVariableScale(1./(24.*3600.)); // [s] -> [day]
        kr_.setFunctionScale(1.0/(24.*3600.*1000*9.81)); // [cm/hPa/day] -> [m/Pa/s]
        kx_ = InputFileFunction("RootSystem.Conductivity","Kx", "KxAge", kxIdx_, orderIdx_); // [cm^4/hPa/day] ([day])
        kx_.setVariableScale(1./(24.*3600.)); // [s] -> [day]
        kx_.setFunctionScale(1e-6/(24.*3600.*1000*9.81)); // [cm^4/hPa/day] -> [m^4/Pa/s]
        radius_ = InputFileFunction("RootSystem.Grid", "Radius", "RadiusAge", radiusIdx_, orderIdx_); // [cm] ([day])
        radius_.setVariableScale(1./(24.*3600.)); // [s] -> [day]
        radius_.setFunctionScale(1.e-2); // [cm] -> [m]
        ct_ = InputFileFunction("RootSystem", "CreationTime", ctIdx_, orderIdx_); // segment creation time [day], optional grid data, no variable
        ct_.setFunctionScale(24.*3600.); // [day] -> [s]
        order_ = InputFileFunction("RootSystem", "Order", orderIdx_, orderIdx_); // [1], optional grid data, no variable
        id_ = InputFileFunction("RootSystem", "Id", idIdx_, idIdx_);

        // Conductivities and radius of (artificial) shoot
        kx0_ = InputFileFunction("RootSystem.Conductivity" , "ShootKx", "ShootKxAge", 1.); // default = high axial flux
        kx0_.setVariableScale(1./(24.*3600.)); // [s] -> [day]
        kx0_.setFunctionScale(1.0/(24.*3600.*1000*9.81)); // [cm^4/hPa/day] -> [m^4/Pa/s]
        kr0_ = InputFileFunction("RootSystem.Conductivity", "ShootKr", "ShootKrAge", 0.); // default = no radial flux
        kr0_.setVariableScale(1./(24.*3600.)); // [s] -> [day]
        kr0_.setFunctionScale(1e-6/(24.*3600.*1000*9.81)); // [cm/hPa/day] -> [m/Pa/s]
        radius0_ = InputFileFunction("RootSystem.Conductivity", "ShootRadius", "ShootRadiusAge", 0.5); // default = large radius
        radius0_.setVariableScale(1./(24.*3600.)); // [s] -> [day]
        radius0_.setFunctionScale(1.e-2); // [cm] -> [m]
    }

    /*!
     * \brief Return the intrinsic permeability for the current sub-control volume in [m^2].
     *
     * \note Kx has units [m^4/(Pa*s)] so we have to divide by the cross-section area
     *       and multiply with a characteristic viscosity
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
        const SubControlVolume& scv, const ElementSolution& elemSol) const {
        return permeability(element);
    }

    //! simpler interface
    PermeabilityType permeability(const Element& element) const {
        Scalar mu = Water::liquidViscosity(285.15, 1e5); // temperature, pressure
        auto eIdx = this->gridGeometry().elementMapper().index(element);
        Scalar a = this->radius(eIdx);
        Scalar kx = this->kx(eIdx);
        return kx * mu / (M_PI * a * a);
    }

    //! \brief returns the porosity \f$[-]\f$
    Scalar porosityAtPos(const GlobalPosition& globalPos) const {
        return 1.;
    }

    // id of the polyline the segment belongs to [1]
    Scalar id(std::size_t eIdx) const {
        return id_.f(0.,eIdx);
    }

    //! segment root order, or root type [1], starts at zero
    int order(std::size_t eIdx) const {
        int o = (int) order_.f(eIdx);
        if (o<0) {
            std::cout << "RootSpatialParams::order: warning root order is negative for element index "
                << eIdx <<": " << o << ", resuming with root order 0\n" << std::flush;
            o = 0;
        }
        return o;
    }

    //! segment radius [m]
    Scalar radius(std::size_t eIdx) const {
        int o = (int) order_.f(eIdx);
        if (o<0) {
            return radius0_.f(this->age(eIdx), eIdx);
        }
        return radius_.f(this->age(eIdx), eIdx);
    }

    //! segment age [s]
    Scalar age(std::size_t eIdx) const {
        return time0_ - ct_.f(eIdx) + time_;
    }

    //! radial conductivity [m/Pa/s] == [m^2 s/kg]
    Scalar kr(std::size_t eIdx) const {
        int o = (int) order_.f(eIdx);
        if (o<0) {
            return kr0_.f(this->age(eIdx), eIdx);
        }
        return kr_.f(this->age(eIdx), eIdx);
    }

    //! axial conductivity [m^4/Pa/s]
    Scalar kx(std::size_t eIdx) const {
        int o = (int) order_.f(eIdx);
        if (o<0) {
            return kx0_.f(this->age(eIdx), eIdx);
        }
        return kx_.f(this->age(eIdx), eIdx);
    }

    //! set current simulation time, age is time dependent (so sad), kx and kr can be age dependent
    void setTime(double t, double dt) {
        time_ = t;
        dt_ = dt;
    }

    //! Read initial parameters from grid data object
    template<class GridData>
    void initParameters(const GridData& gridData) {
        const auto& fvGridGeometry = this->gridGeometry();
        kr_.setGridData(gridData, fvGridGeometry);
        kx_.setGridData(gridData, fvGridGeometry);
        order_.setGridData(gridData, fvGridGeometry);
        id_.setGridData(gridData, fvGridGeometry);
        radius_.setGridData(gridData, fvGridGeometry);
        ct_.setGridData(gridData, fvGridGeometry);
    }

    //! ignore (because there is no common base class with rootspatialparams_rb.hh)
    template<class RootSystem>
    void updateParameters(const RootSystem& rs) {
        DUNE_THROW(Dune::InvalidStateException, "updateParameters is called for DGF");
    }

private:

    int orderIdx_;
    int radiusIdx_;
    int ctIdx_;
    int krIdx_;
    int kxIdx_;
    int idIdx_;

    InputFileFunction kr_;
    InputFileFunction kx_;
    InputFileFunction order_;
    InputFileFunction radius_;
    InputFileFunction ct_;
    InputFileFunction id_;

    InputFileFunction kx0_;
    InputFileFunction kr0_;
    InputFileFunction radius0_;

    double time_ = 0.;
    double dt_ = 0.;
    double time0_ = 0.; // for calculating age from ct
};

} // end namespace Dumux

#endif
