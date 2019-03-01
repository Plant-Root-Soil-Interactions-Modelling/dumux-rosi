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
#ifndef DUMUX_ROOT_SPATIALPARAMS_DGF_HH
#define DUMUX_ROOT_SPATIALPARAMS_DGF_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/material/spatialparams/fv1p.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/io/inputfilefunction.hh>



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

    // parameter indices, where the parameter (if used) is located the dgf file
    enum {
        orderIdx = 0, radiusIdx = 1, ageIdx = 2, krIdx = 3, kxIdx = 4
    };

public:
    // export permeability type
    using PermeabilityType = Scalar;

    RootSpatialParamsDGF(std::shared_ptr<const FVGridGeometry> fvGridGeometry) :
        ParentType(fvGridGeometry) {
        kr_ = InputFileFunction("RootSystem.Conductivity.Kr", "RootSystem.Conductivity.KrAge", krIdx, orderIdx);
        kx_ = InputFileFunction("RootSystem.Conductivity.Kx", "RootSystem.Conductivity.KxAge", kxIdx, orderIdx);
        radius_ = InputFileFunction("RootSystem.Radius", "RootSystem.RadiusAge", radiusIdx, orderIdx);
        age_ = InputFileFunction("RootSystem.Age", ageIdx, orderIdx);
        order_ = InputFileFunction("RootSystem.Order", orderIdx, orderIdx); //
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
        // std::cout << "params " << kx * 1e13 << ", " << a << ", " << mu << "\n";
        return kx * mu / (M_PI * a * a);
    }

    //! \brief returns the porosity \f$[-]\f$
    Scalar porosityAtPos(const GlobalPosition& globalPos) const {
        return 1.;
    }

    Scalar order(std::size_t eIdx) const {
        return order_.f(eIdx);
    }

    Scalar radius(std::size_t eIdx) const {
        return radius_.f(this->age(eIdx), eIdx);
    }

    Scalar age(std::size_t eIdx) const {
        return age_.f(eIdx)+time_;
    }

    Scalar kr(std::size_t eIdx) const {
        return kr_.f(this->age(eIdx), eIdx);
    }

    Scalar kx(std::size_t eIdx) const {
        return kx_.f(this->age(eIdx), eIdx);
    }

    void setTime(double t) {
        time_ = t;
    }

    //! Read initial parameters from grid data object
    template<class GridData>
    void initParameters(const GridData& gridData)
    {
        const auto& fvGridGeometry = this->fvGridGeometry();
        kr_.setGridData(gridData, fvGridGeometry);
        kx_.setGridData(gridData, fvGridGeometry);
        order_.setGridData(gridData, fvGridGeometry);
        radius_.setGridData(gridData, fvGridGeometry);
        age_.setGridData(gridData, fvGridGeometry);
    }

private:
    InputFileFunction kr_;
    InputFileFunction kx_;
    InputFileFunction order_;
    InputFileFunction radius_;
    InputFileFunction age_;

    double time_ = 0.;
};

} // end namespace Dumux

#endif
