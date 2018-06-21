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
 * \ingroup OneTests
 * \brief The spatial parameters class blood flow problem
 */
#ifndef DUMUX_ROOT_SPATIALPARAMS_HH
#define DUMUX_ROOT_SPATIALPARAMS_HH

#include <dumux/common/parameters.hh>
#include <dumux/material/spatialparams/fv1p.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/growth/rootparameters.hh>

#include <RootSystem.h>

namespace Dumux {

/*!
 * \brief Root spatial params
 */
template<class FVGridGeometry, class Scalar>
class RootSpatialParams
: public FVSpatialParamsOneP<FVGridGeometry, Scalar, RootSpatialParams<FVGridGeometry, Scalar>>
{
    using ThisType = RootSpatialParams<FVGridGeometry, Scalar>;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar, ThisType>;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using RootParams = GrowthModule::RootParams<Scalar>;
    using DGFParamIndices = GrowthModule::DGFParamIndices;

public:

    // export permeability type
    using PermeabilityType = Scalar;

    RootSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        porosity_ = getParam<Scalar>("Root.SpatialParams.Porosity", 0.4);
        constantKx_ = getParam<Scalar>("SpatialParams.Kx", 5.0968e-17);
        constantKr_ = getParam<Scalar>("SpatialParams.Kr", 2.04e-13);
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
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        const Scalar r = rootParams(element).radius;
        return constantKx_ / (M_PI*r*r) * Components::SimpleH2O<Scalar>::liquidViscosity(285.15, 1e5);
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos the scv center
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return porosity_;
    }

    //! Return the root parameters
    const RootParams& rootParams(const Element &element) const
    {
        const auto& elementMapper = this->fvGridGeometry().elementMapper();
        const auto eIdx = elementMapper.index(element);
        int rootId = rootId_[eIdx];
        if (rootId < 0)
            return shootParams_;
        else
            return rootParams_[rootId];
    }

    Scalar previousLength(std::size_t eIdx) const
    {
        return previousLength_[eIdx];
    }

    Scalar radius(std::size_t eIdx) const
    {
        int rootId = rootId_[eIdx];
        if (rootId < 0)
            return shootParams_.radius;
        else
            return rootParams_[rootId].radius;
    }

    Scalar Kr(std::size_t eIdx) const
    {
        return constantKr_;
    }

    //! Access to the radii vector for output
    const std::vector<Scalar>& radii() const
    { return radii_; }

    //! Access to the orders vector for output
    const std::vector<Scalar>& orders() const
    { return orders_; }

    //! Read parameters from the a crootbox rootsystem
    void initParameters(const ::RootSystem& rs)
    {
        // be careful: the total number of roots and the grown number of roots is not the same
        rootParams_.resize(rs.getNumberOfRoots(/*allRoots=*/true));
        // this gets the root that already have at least one segment
        const auto roots = rs.getRoots();
        for (auto&& r : roots)
        {
            auto& rp = rootParams_[r->id];

            // compute the root order
            rp.order = 0;
            auto* r_ = r;
            while (r_->parent != nullptr)
            {
                ++rp.order;
                r_=r_->parent;
            }

            rp.radius = r->param.a * 0.01; // conversion from cm to m
            rp.rootId = r->id;
            rp.axialPerm = constantKx_;
            rp.radialPerm = constantKr_;
            rp.plantId = 0;
        }

        // set the shootParams
        shootParams_.order = 0;
        shootParams_.radius = 0.01;
        shootParams_.rootId = -1;
        shootParams_.axialPerm = constantKx_*1e20;rootParams_
        shootParams_.radialPerm = 0.0;
        shootParams_.plantId = 0;

        const auto shootSegments = rs.getShootSegments();
        const auto segmentRoots = rs.getSegmentsOrigin();
        const auto& gridView = this->fvGridGeometry().gridView();
        rootId_.resize(gridView.size(0));
        radii_.resize(gridView.size(0));
        orders_.resize(gridView.size(0));
        previousLength_.resize(gridView.size(0));
        for (const auto& element : elements(gridView))
        {
            const auto& elementMapper = this->fvGridGeometry().elementMapper();
            const auto eIdx = elementMapper.index(element);

            // initialize the previous length to the current length
            previousLength_[eIdx] = element.geometry().volume();

            // subtract the number of shoot segments to get the normal segment index
            const int segmentIdx = eIdx - shootSegments.size();
            // shoot segments
            if (segmentIdx < 0)
            {
                rootId_[eIdx] = -1;
                const auto& rp = rootParams(element);
                radii_[eIdx] = rp.radius;
                orders_[eIdx] = rp.order;
            }
            // regular segements
            else
            {
                // we hope that the crootbox indices and the dune indices on initialization are the same!
                rootId_[eIdx] = segmentRoots[segmentIdx]->id;
                const auto& rp = rootParams_[rootId_[eIdx]];
                radii_[eIdx] = rp.radius;
                orders_[eIdx] = rp.order;
            }
        }
    }

    //! update the root parameters that aren't set by growth directly
    void updateParameters(const ::RootSystem& rs)
    {
        const auto roots = rs.getRoots();
        for (auto&& r : roots)
        {
            auto& rp = rootParams_[r->id];
            if (rp.isUpdated)
                continue;

            // TODO: compute the permeabilities dependent on age/radius
            rp.axialPerm = constantKx_;
            rp.radialPerm = constantKr_;
            rp.isUpdated = true;
        }

        // update radii for output
        const auto& gridView = this->fvGridGeometry().gridView();
        radii_.resize(gridView.size(0));
        orders_.resize(gridView.size(0));
        for (const auto& element : elements(gridView))
        {
            const auto& elementMapper = this->fvGridGeometry().elementMapper();
            const auto eIdx = elementMapper.index(element);
            const auto& rp = rootParams(element);
            radii_[eIdx] = rp.radius;
            orders_[eIdx] = rp.order;
        }
    }

    //! resize the segment parameters
    void resizeSegmentParams(std::size_t size)
    {
        rootId_.resize(size);
        previousLength_.resize(size);
    }

    //! set root index (needed by growth if eIdx of an element changed after growth)
    void setRootId(unsigned int eIdx, int rootId)
    { rootId_[eIdx] = rootId; }

    //! set the previous length
    void setPreviousLength(unsigned int eIdx, Scalar length)
    { previousLength_[eIdx] = length; }

    //! Resize root params to add more roots
    void rootParamsResize(std::size_t size)
    { rootParams_.resize(size); }

    //! Return a reference to the root parameters of a root id to change them (called from growth module)
    RootParams& getRootParams(int rootId)
    { return rootParams_[rootId]; }

private:
    //! Parameters for every root
    std::vector<RootParams> rootParams_;
    RootParams shootParams_;

    //! Segment paramters
    std::vector<int> rootId_; //! the root id for each segment
    std::vector<Scalar> previousLength_; //! the length of the element at the last time step (new elements return 0 length)

    std::vector<Scalar> radii_; //! radii for output
    std::vector<Scalar> orders_;  //! root orders for output

    Scalar porosity_, constantKx_, constantKr_;
};

} // end namespace Dumux

#endif
