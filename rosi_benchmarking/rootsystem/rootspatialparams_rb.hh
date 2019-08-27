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
    using GrowingSystem = typename GrowthModule::GrowthInterface<GlobalPosition>; // decoupled from CRootBox

public:

    using PermeabilityType = Scalar; // export permeability type

    RootSpatialParamsRB(std::shared_ptr<const FVGridGeometry> fvGridGeometry) :
        ParentType(fvGridGeometry) {

        kr_ = InputFileFunction("RootSystem.Conductivity", "Kr", "KrAge", 0, 0); // [cm/hPa/day] ([day]), indices are not used, data are set manually in updateParameters
        kr_.setVariableScale(1./(24.*3600.)); // [s] -> [day]
        kr_.setFunctionScale(1.e-4/(24.*3600.)); // [cm/hPa/day] -> [m/Pa/s]

        kx_ = InputFileFunction("RootSystem.Conductivity", "Kx", "KxAge", 0, 0); // [cm^4/hPa/day] ([day]), indices are not used, data are set manually in updateParameters
        kx_.setVariableScale(1./(24.*3600.)); // [s] -> [day]
        kx_.setFunctionScale(1.e-10/(24.*3600.)); // [cm^4/hPa/day] -> [m^4/Pa/s]

        kx0_ = InputFileFunction("RootSystem.Conductivity" , "ShootKx", "ShootKxAge", 1, 0);
        kx0_.setVariableScale(1./(24.*3600.)); // [s] -> [day]
        kx0_.setFunctionScale(1.e-10/(24.*3600.)); // [cm^4/hPa/day] -> [m^4/Pa/s]

        kr0_ = InputFileFunction("RootSystem.Conductivity", "ShootKr", "ShootKxAge", 0, 0);
        kr0_.setVariableScale(1./(24.*3600.)); // [s] -> [day]
        kr0_.setFunctionScale(1.e-4/(24.*3600.)); // [cm/hPa/day] -> [m/Pa/s]

        time0_ = getParam<double>("RootSystem.Grid.InitialT")*3600.*24.; // root system initial time
        double radius0 = getParam<double>("RootSystem.Grid.ShootRadius", 1.5)/100; // [cm] -> [m]
        radii_ = { radius0 }; // vectors will incrementially grow in updateParameters

        /*InputFileFunction radius0 = InputFileFunction("RootSystem.Grid", "ShootRadius", 1.5);
        radius0.setVariableScale(1./(24.*3600.)); // [s] -> [day]
        radius0.setFunctionScale(1.e-2); // [cm] -> [m]
        radii_ = { radius0 };*/

        orders_ = { 0 };
        ctimes_ = { 0. };
        ids_ = { 0 };
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

    // [1] the id of the polyline the segment belongs to
    Scalar id(std::size_t eIdx) const {
        return ids_[eIdx];
    }

    // [1] the order of the root the segment belongs to, starts at zero
    Scalar order(std::size_t eIdx) const {
        return orders_[eIdx];
    }

    // [m]
    Scalar radius(std::size_t eIdx) const {
        return radii_[eIdx]; // m
    }

    // [s]
    Scalar age(std::size_t eIdx) const {
        double a = (time0_+time_) - ctimes_[eIdx];
        return a;
    }

    //! radial conductivity [m /Pa/s]
    Scalar kr(std::size_t eIdx) const {
        if (eIdx==0) { // no radial flow at the shoot element
            return kr0_.f(this->age(eIdx), eIdx);
        }
        return kr_.f(this->age(eIdx), eIdx);
    }

    //! axial conductivity [m^4/Pa/s]
    Scalar kx(std::size_t eIdx) const {
        if (eIdx==0) { // high axial flow at the shoot element
            return kx0_.f(this->age(eIdx), eIdx);
        }
        return kx_.f(this->age(eIdx), eIdx);
    }

    //! sets the simulation time @param t [s]
    void setTime(double t, double dt) {
        time_ = t;
        dt_ = dt;
    }

    //! Update the Root System Parameters (root system must implement GrowthModule::GrowthInterface)
    void updateParameters(const GrowingSystem& rs) {

        const auto& gridView = this->fvGridGeometry().gridView();
        radii_.resize(gridView.size(0));
        orders_.resize(gridView.size(0));
        ctimes_.resize(gridView.size(0));
        ids_.resize(gridView.size(0));

        auto segs = rs.newSegments();
        auto segCT = rs.segmentCreationTimes();
        auto segO = rs.segmentParameter("order");
        auto segRadii = rs.segmentRadii();
        auto segId = rs.segmentParameter("id");

        if (segs.size()!=segCT.size() || segs.size()!=segO.size() || segs.size()!=segRadii.size() || segs.size()!=segId.size()) { // sanity check
            throw Dumux::ParameterException("updateParameters: new segments sizes differ");
        }

        std::cout << "updateParameters: " << gridView.size(0) << ": " << segs.size() << " new segments " << "\n"<< std::flush;

        for (size_t i = 0; i < segs.size(); i++) {
            size_t rIdx = segs[i][1] - 1; // rootbox segment index = second node index - 1
            size_t eIdx = rs.map2dune(rIdx);
            // std::cout << "updateParameters: age at root index " << rIdx << " element index " << eIdx << " = " <<  segCT[i] << " s = " << segCT[i]/24/3600 << " d \n";
            orders_.at(eIdx) = segO[i];
            ids_.at(eIdx) = segId[i];
            radii_.at(eIdx) = segRadii[i];
            ctimes_.at(eIdx) = segCT[i];
            if (segCT[i]<0) { // sanity checks
                throw Dumux::ParameterException("updateParameters: creation time cannot be negative");
            }
            if (segCT[i]>time_+time0_+dt_+1) {// sanity checks
                throw Dumux::ParameterException("updateParameters: creation time cannot be larger than simulation time, "+
                    std::to_string(segCT[i])+">"+std::to_string(time_+time0_));
            }
        }

        // update ctimes: when tips move, there segmentCTs need to be updated
        auto uni = rs.updatedNodeIndices();
        auto cts = rs.updatedNodeCTs();
        for (int i : uni) {
            size_t rIdx = i - 1; // rootbox segment index = node index - 1
            size_t eIdx = rs.map2dune(rIdx);
            ctimes_.at(eIdx) = segCT[i]; // replace time
        }

        if ((kr_.type() == InputFileFunction::perType) || (kr_.type() == InputFileFunction::tablePerType)) {
            kr_.setData(orders_);
        }
        if ((kx_.type() == InputFileFunction::perType) || (kx_.type() == InputFileFunction::tablePerType)) {
            kx_.setData(orders_);
        }
        // std::cout << "updateParameters done\n" << std::flush;
    }

    //! ignore (used in rootspatialparams_dgf.hh)
    template<class GridData>
    void initParameters(const GridData& gridData) {
        DUNE_THROW(Dune::InvalidStateException, "initParameters is called for a root growth model");
    }

private:
    InputFileFunction kr_;
    InputFileFunction kx_;
    InputFileFunction kx0_;
    InputFileFunction kr0_;

    std::vector<double> radii_; // [m]
    std::vector<double> ctimes_; // [s]
    std::vector<double> ids_; // [1]
    std::vector<double> orders_; // root order, or root type [1]


    double time_ = 0.; // [s]
    double dt_ = 0.;
    double time0_ = 0; // initial time [s]

};

} // end namespace Dumux

#endif
