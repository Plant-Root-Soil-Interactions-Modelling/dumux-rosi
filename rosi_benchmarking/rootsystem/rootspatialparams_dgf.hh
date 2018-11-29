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

#include <dumux/growth/rootparameters.hh>

namespace Dumux {

/*!
 * \brief Root spatial params
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

    using RootParams = GrowthModule::RootParams<Scalar>;

    enum {
        orderIdx = 0, radiusIdx = 1, ageIdx = 2, krIdx = 3, kxIdx = 4
    };

public:
    // export permeability type
    using PermeabilityType = Scalar;

    RootSpatialParamsDGF(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
    }

    /*!
     * \brief Return the intrinsic permeability for the current sub-control volume in [m^2].
     *
     * \param ipGlobal The integration point
     * \note Kx has units [m^4/(Pa*s)] so we have to divide by the cross-section area
     *       and multiply with a characteristic viscosity
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
        const SubControlVolume& scv, const ElementSolution& elemSol) const
    { return permeability(element); }

    //! simpler interface
    PermeabilityType permeability(const Element& element) const
    {
        auto eIdx = this->fvGridGeometry().elementMapper().index(element);
        Scalar r = radii_[eIdx];
        return this->Kx(eIdx) / (M_PI*r*r) * Components::SimpleH2O<Scalar>::liquidViscosity(285.15, 1e5);
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos the scv center
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return 1.;
    }

    Scalar order(std::size_t eIdx) const
    {
        return orders_[eIdx];
    }

    Scalar radius(std::size_t eIdx) const
    {
        return radii_[eIdx];
    }

    Scalar age(std::size_t eIdx) const {
        return ages_[eIdx];
    }

    Scalar Kr(std::size_t eIdx) const
    {
        return Kr_[eIdx];
    }

    Scalar Kx (std::size_t eIdx) const
    {
        return Kx_[eIdx];
    }

    //! Read initial parameters from grid data object
    template<class GridData>
    void initParameters(const GridData& gridData)
    {
        // extrac the segment params
        const auto& gridView = this->fvGridGeometry().gridView();
        radii_.resize(gridView.size(0));
        orders_.resize(gridView.size(0));
        ages_.resize(gridView.size(0));
        Kx_.resize(gridView.size(0));
        Kr_.resize(gridView.size(0));
        const auto& elementMapper = this->fvGridGeometry().elementMapper();
        for (const auto& element : elements(gridView))
        {
            const auto eIdx = elementMapper.index(element);
            const auto& elemParams = gridData.parameters(element);

            // initialize the previous length to the current length
            ages_[eIdx] = elemParams[ageIdx] / 86400;
            radii_[eIdx] = elemParams[radiusIdx];
            orders_[eIdx] = elemParams[orderIdx];
            const auto& order = orders_[eIdx];
            if (order < 0)
            {
                Kx_[eIdx] = 1.0e-10; // Kx: very high value for shoots to not contribute to resistance
                Kr_[eIdx] = 0.0; // Kr: 0.0 for shoots so that shoot element do not exchange mass with the soil
                radii_[eIdx] = 0.01;
            }
            else
            {
                Kx_[eIdx] = elemParams[kxIdx];
                Kr_[eIdx] = elemParams[krIdx];
            }
        }
    }

    //! Output an analysis of the root system
    void analyseRootSystem() const
    {
        // compute some statistics of the initial root system
        Scalar totalLength = 0.0, totalLengthTop = 0.0, totalLengthBottom = 0.0;
        Scalar totalLengthPrimary = 0.0, totalLengthSecondary = 0.0;
        Scalar rootVolume = 0.0;
        Scalar totalAge = 0.0;
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
            if (orders_[eIdx] >= 0) // exclude shoot
            {
                const auto geo = element.geometry();
                const auto length = geo.volume();
                const auto r = radii_[eIdx];
                totalLength += length;
                rootVolume += length*M_PI*r*r;
                totalAge = std::max(totalAge, ages_[eIdx]);

                if (geo.center()[2] > -0.42)
                    totalLengthTop += length;
                else
                    totalLengthBottom += length;

                if (orders_[eIdx] == 0)
                    totalLengthPrimary += length;
                else
                    totalLengthSecondary += length;
            }
        }

        std::cout << ".........................................................\n"
            << "-- Root system age:            " << totalAge << " days\n"
            << "-- Total length:               " << totalLength << " m\n" << "-- Total length (top 42 cm):   " << totalLengthTop << " m\n"
            << "-- Total length (below 42 cm): " << totalLengthBottom << " m\n" << "-- Total length (primary):     " << totalLengthPrimary << " m\n"
            << "-- Total length (secondary):   " << totalLengthSecondary << " m\n" << "-- Total volume:               " << rootVolume << " m³\n"
            << ".........................................................\n";
    }

private:

    //! Tabularized root conductivity functions (pairs of x, y) for each order
    std::vector<std::pair<std::vector<double>, std::vector<double>>> axialRootConductivity_;
    std::vector<std::pair<std::vector<double>, std::vector<double>>> radialRootConductivity_;

    //! Segment parameters
    std::vector<Scalar> orders_;  //! root orders for output
    std::vector<Scalar> radii_; //! radii for output
    std::vector<Scalar> ages_; //! the current age of the element
    std::vector<Scalar> Kx_;  //! axial conductivities
    std::vector<Scalar> Kr_;  //! radial conductivities

};

} // end namespace Dumux

#endif
