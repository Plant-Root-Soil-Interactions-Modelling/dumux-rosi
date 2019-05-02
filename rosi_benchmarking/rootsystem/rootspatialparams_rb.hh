// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief The spatial parameters class
 */
#ifndef DUMUX_ROOT_SPATIALPARAMS_RB_HH
#define DUMUX_ROOT_SPATIALPARAMS_RB_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/material/spatialparams/fv1p.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/growth/growthinterface.hh>
#include <dumux/growth/gridgrowth.hh>

#include <dumux/io/inputfilefunction.hh>
#include <dumux/growth/soillookup.hh>



namespace Dumux {

/*!
 * \brief Root spatial parameters class for static CRootBox root systems
 *
 * use initParameters to initialize the class with data from the root system model
 *
 */
template<class FVGridGeometry, class Scalar>
class RootSpatialParamsRB
    : public FVSpatialParamsOneP<FVGridGeometry, Scalar, RootSpatialParamsRB<FVGridGeometry, Scalar>> {
    using ThisType = RootSpatialParamsRB<FVGridGeometry, Scalar>;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar, ThisType>;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Water = Components::SimpleH2O<Scalar>;
    using RootSystem = typename GrowthModule::GrowthInterface<GlobalPosition>; // decoupled from CRootBox

public:
    // export permeability type
    using PermeabilityType = Scalar;

    RootSpatialParamsRB(std::shared_ptr<const FVGridGeometry> fvGridGeometry) :
        ParentType(fvGridGeometry) {
        kr_ = InputFileFunction("RootSystem.Conductivity.Kr", "RootSystem.Conductivity.KrAge", -1, -1);
        kx_ = InputFileFunction("RootSystem.Conductivity.Kx", "RootSystem.Conductivity.KxAge", -1, -1);
        assert(kr_.type() != InputFileFunction::data && "RootSpatialParamsRB: no grid data available for CRootBox root systems"); // todo use kr_.setData(...) by hand
        assert(kx_.type() != InputFileFunction::data && "RootSpatialParamsRB: no grid data available for CRootBox root systems");
        time0_ = getParam<double>("RootSystem.Grid.InitialT")*3600.*24.; // root system initial time
        // shoot
        radii_ = { 1.17/100. };
        orders_ = { 0 };
        ctimes_ = { 0. };
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
        auto eIdx = this->fvGridGeometry().elementMapper().index(element);
        Scalar a = this->radius(eIdx);
        Scalar kx = this->kx(eIdx);
        return kx * mu / (M_PI * a * a);
    }

    //! \brief returns the porosity \f$[-]\f$
    Scalar porosityAtPos(const GlobalPosition& globalPos) const {
        return 1.;
    }

    Scalar order(std::size_t eIdx) const {
        return orders_[eIdx];
    }

    // m
    Scalar radius(std::size_t eIdx) const {
        return radii_[eIdx]; // m
    }

    // s
    Scalar age(std::size_t eIdx) const {
        double a = (time0_+time_) - ctimes_[eIdx];
        return a; // todo there is a crootbox bug, where time0_ << ctimes_ for certain tip segments (???)
    }

    //! radial conductivity [m /Pa/s]
    Scalar kr(std::size_t eIdx) const {
        if (eIdx==0) { // no radial flow at the shoot element
            return 0.;
        }
        return kr_.f(this->age(eIdx)/(24.*3600.), eIdx)*1e-4/(24*3600); // cm / hPa / day -> m / Pa /s
    }

    //! axial conductivity [m^4/Pa/s]
    Scalar kx(std::size_t eIdx) const {
//        if (eIdx==0) { // high axial flow at the shoot element
//            return 1.;
//        }
        return kx_.f(this->age(eIdx)/(24.*3600.), eIdx)*1e-10/(24*3600); // cm^4 / hPa / day -> m^4 / Pa /s
    }

    void setTime(double t) {
        time_ = t;
    }

    //! Update the Root System Parameters (root system must implement GrowthModule::GrowthInterface)
    void updateParameters(const RootSystem& rs) {

        const auto& gridView = this->fvGridGeometry().gridView();
        radii_.resize(gridView.size(0));
        orders_.resize(gridView.size(0));
        ctimes_.resize(gridView.size(0));

        auto segs = rs.newSegments();
        auto segCT = rs.segmentCreationTimes();
        auto segO = rs.segmentOrders();
        auto segRadii = rs.segmentRadii();

        std::cout << "UPDATING " << gridView.size(0) << ": " << segs.size() << " segments " << "\n"<< std::flush;

        for (size_t i = 0; i < segs.size(); i++) {
            size_t rIdx = segs[i][1] - 1; // rootbox segment index = second node index - 1 ==
            size_t eIdx = rs.map2dune(rIdx);
            // std::cout << "updateParameters: age at root index " << rIdx << " element index " << eIdx << " = " <<  segCT[i] << " s = " << segCT[i]/24/3600 << " d \n";
            orders_.at(eIdx) = segO[i];
            radii_.at(eIdx) = segRadii[i];
            ctimes_.at(eIdx) = segCT[i];
        }

        if ((kr_.type() == InputFileFunction::perType) || (kr_.type() == InputFileFunction::tablePerType)) {
            kr_.setData(orders_);
        }
        if ((kx_.type() == InputFileFunction::perType) || (kx_.type() == InputFileFunction::tablePerType)) {
            kx_.setData(orders_);
        }
    }

    //! ignore (because there is no common base class with rootspatialparams_dgf.hh)
    template<class GridData>
    void initParameters(const GridData& gridData) {
        DUNE_THROW(Dune::InvalidStateException, "initParameters is called for a root growth model");
    }

private:
    InputFileFunction kr_;
    InputFileFunction kx_;

    std::vector<double> orders_; // root order, or root type
    std::vector<double> radii_; // [m]
    std::vector<double> ctimes_; // [s]

    double time_ = 0.; // [s]
    double time0_ = 0; // initial time [s]
};

} // end namespace Dumux

#endif
