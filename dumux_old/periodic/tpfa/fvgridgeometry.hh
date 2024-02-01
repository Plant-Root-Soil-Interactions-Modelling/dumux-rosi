// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \author Timo Koch <timo.koch@iws.uni-stuttgart.de>
 * \author Timo Koch
 * \ingroup CCTpfaDiscretization
 * \brief The finite volume geometry (scvs and scvfs) for cell-centered TPFA models on a grid view
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element of the grid partition.
 */
#ifndef DUMUX_DISCRETIZATION_CCTPFA_FV_GRID_GEOMETRY_PERIODIC_HH
#define DUMUX_DISCRETIZATION_CCTPFA_FV_GRID_GEOMETRY_PERIODIC_HH

#include <algorithm>

#include <dumux/common/indextraits.hh>
#include <dumux/common/defaultmappertraits.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/basefvgridgeometry.hh>
#include <dumux/discretization/checkoverlapsize.hh>
#include <dumux/discretization/cellcentered/subcontrolvolume.hh>
#include <dumux/discretization/cellcentered/connectivitymap.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/fvelementgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/subcontrolvolumeface.hh>

namespace Dumux {

/*!
 * \ingroup CCTpfaDiscretization
 * \brief The finite volume geometry (scvs and scvfs) for cell-centered TPFA models on a grid view
 *        This builds up the sub control volumes and sub control volume faces
 * \note This class is specialized for versions with and without caching the fv geometries on the grid view
 */
template<class GridView,
         bool enableFVGridGeometryCache = false,
         class Traits = CCTpfaDefaultGridGeometryTraits<GridView> >
class PeriodicCCTpfaFVGridGeometry;

/*!
 * \ingroup CCTpfaDiscretization
 * \brief The finite volume geometry (scvs and scvfs) for cell-centered TPFA models on a grid view
 *        This builds up the sub control volumes and sub control volume faces
 * \note For caching enabled we store the fv geometries for the whole grid view which is memory intensive but faster
 */
template<class GV, class Traits>
class PeriodicCCTpfaFVGridGeometry<GV, true, Traits>
: public BaseFVGridGeometry<PeriodicCCTpfaFVGridGeometry<GV, true, Traits>, GV, Traits>
{
    using ThisType = PeriodicCCTpfaFVGridGeometry<GV, true, Traits>;
    using ParentType = BaseFVGridGeometry<ThisType, GV, Traits>;
    using ConnectivityMap = typename Traits::template ConnectivityMap<ThisType>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using Element = typename GV::template Codim<0>::Entity;

    static const int dim = GV::dimension;
    static const int dimWorld = GV::dimensionworld;

public:
    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, true>;
    //! export the type of sub control volume
    using SubControlVolume = typename Traits::SubControlVolume;
    //! export the type of sub control volume
    using SubControlVolumeFace = typename Traits::SubControlVolumeFace;
    //! export dof mapper type
    using DofMapper = typename Traits::ElementMapper;

    //! export the discretization method this geometry belongs to
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;

    //! The maximum admissible stencil size (used for static memory allocation during assembly)
    static constexpr int maxElementStencilSize = LocalView::maxNumElementScvfs*Traits::maxNumScvfNeighbors + 1;

    //! export the grid view type
    using GridView = GV;

    //! Constructor
    PeriodicCCTpfaFVGridGeometry(const GridView& gridView)
    : ParentType(gridView)
    {
        // Check if the overlap size is what we expect
        if (!CheckOverlapSize<DiscretizationMethod::cctpfa>::isValid(gridView))
            DUNE_THROW(Dune::InvalidStateException, "The cctpfa discretization method needs at least an overlap of 1 for parallel computations. "
                                                     << " Set the parameter \"Grid.Overlap\" in the input file.");
    }

    //! the element mapper is the dofMapper
    //! this is convenience to have better chance to have the same main files for box/tpfa/mpfa...
    const DofMapper& dofMapper() const
    { return this->elementMapper(); }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return scvs_.size();
    }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    {
        return scvfs_.size();
    }

    //! The total number of boundary sub control volume faces
    std::size_t numBoundaryScvf() const
    {
        return numBoundaryScvf_;
    }

    //! The total number of degrees of freedom
    std::size_t numDofs() const
    { return this->gridView().size(0); }

    //! to be called before update to connect elements that are not connected by grid intersections
    //! maps from boundary vertices to all connected elements
    void setExtraConnectivity(const std::unordered_map<GridIndexType, std::vector<GridIndexType>>& connectivity)
    { extraConnectivity_ = connectivity; }

    //! update all fvElementGeometries (do this again after grid adaption)
    void update()
    {
        static_assert(dim == 1, "Only implemented for one-dimensional network grids");

        ParentType::update();

        // clear containers (necessary after grid refinement)
        scvs_.clear();
        scvfs_.clear();
        scvfIndicesOfScv_.clear();
        flipScvfIndices_.clear();

        // determine size of containers
        std::size_t numScvs = numDofs();
        std::size_t numScvf = 0;
        for (const auto& element : elements(this->gridView()))
            numScvf += element.subEntities(1);

        // reserve memory
        scvs_.resize(numScvs);
        scvfs_.reserve(numScvf);
        scvfIndicesOfScv_.resize(numScvs);
        hasBoundaryScvf_.resize(numScvs, false);

        // Build the scvs and scv faces
        GridIndexType scvfIdx = 0;
        numBoundaryScvf_ = 0;
        for (const auto& element : elements(this->gridView()))
        {
            const auto eIdx = this->elementMapper().index(element);
            scvs_[eIdx] = SubControlVolume(element.geometry(), eIdx);

            // the element-wise index sets for finite volume geometry
            std::vector<GridIndexType> scvfsIndexSet;
            scvfsIndexSet.reserve(element.subEntities(1));

            // for network grids there might be multiple intersection with the same geometryInInside
            // we indentify those by the indexInInside for now (assumes conforming grids at branching facets)
            using ScvfGridIndexStorage = typename SubControlVolumeFace::Traits::GridIndexStorage;
            std::vector<ScvfGridIndexStorage> neighborIndices;

            //! first, push inside index in all neighbor sets
            neighborIndices.resize(element.subEntities(1));
            std::for_each(neighborIndices.begin(), neighborIndices.end(), [eIdx] (auto& nIndices) { nIndices.reserve(4); nIndices.push_back(eIdx); });

            // second, insert neighbors
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (intersection.neighbor())
                {
                    const auto nIdx = this->elementMapper().index( intersection.outside() );
                    neighborIndices[intersection.indexInInside()].push_back(nIdx);
                }

                else if (intersection.boundary())
                {
                    const auto vIdx = this->vertexMapper().subIndex(element, intersection.indexInInside(), dim);
                    if (extraConnectivity_.count(vIdx)) // if it is a periodic vertex
                    {
                        const auto& neighbors = extraConnectivity_.at(vIdx);
                        for (const auto nIdx : neighbors)
                            if (nIdx != eIdx) // we already put the inside element at the beginning
                                neighborIndices[intersection.indexInInside()].push_back(nIdx);
                    }
                }
            }

            for (const auto& intersection : intersections(this->gridView(), element))
            {
                // inner sub control volume faces
                if (intersection.neighbor())
                {
                    auto indexInInside = intersection.indexInInside();
                    // check if we already handled this facet
                    if (neighborIndices[indexInInside].empty())
                        continue;
                    else
                    {
                        scvfs_.emplace_back(intersection,
                                            intersection.geometry(),
                                            scvfIdx,
                                            neighborIndices[indexInInside],
                                            false);
                        scvfsIndexSet.push_back(scvfIdx++);
                        neighborIndices[indexInInside].clear();
                    }
                }
                // boundary sub control volume faces
                else if (intersection.boundary())
                {
                    const auto indexInInside = intersection.indexInInside();
                    const auto vIdx = this->vertexMapper().subIndex(element, indexInInside, dim);
                    if (extraConnectivity_.count(vIdx)) // if it is a periodic vertex
                    {
                        scvfs_.emplace_back(intersection,
                                            intersection.geometry(),
                                            scvfIdx,
                                            neighborIndices[indexInInside],
                                            false);
                        scvfsIndexSet.push_back(scvfIdx++);
                        neighborIndices[indexInInside].clear();
                    }

                    // regular boundary intersection
                    else
                    {
                        scvfs_.emplace_back(intersection,
                                            intersection.geometry(),
                                            scvfIdx,
                                            ScvfGridIndexStorage({eIdx, static_cast<GridIndexType>(this->gridView().size(0) + numBoundaryScvf_++)}),
                                            true);
                        scvfsIndexSet.push_back(scvfIdx++);

                        hasBoundaryScvf_[eIdx] = true;
                    }
                }
            }

            // Save the scvf indices belonging to this scv to build up fv element geometries fast
            scvfIndicesOfScv_[eIdx] = scvfsIndexSet;
        }

        // Make the flip index set for network and surface grids
        flipScvfIndices_.resize(scvfs_.size());
        for (auto&& scvf : scvfs_)
        {
            if (scvf.boundary())
                continue;

            flipScvfIndices_[scvf.index()].resize(scvf.numOutsideScvs());
            const auto insideScvIdx = scvf.insideScvIdx();
            // check which outside scvf has the insideScvIdx index in its outsideScvIndices
            for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
                flipScvfIndices_[scvf.index()][i] = findFlippedScvfIndex_(insideScvIdx, scvf.outsideScvIdx(i));
        }

        // build the connectivity map for an effecient assembly
        connectivityMap_.update(*this);
    }

    //! Get a sub control volume with a global scv index
    const SubControlVolume& scv(GridIndexType scvIdx) const
    {
        return scvs_[scvIdx];
    }

    //! Get a sub control volume face with a global scvf index
    const SubControlVolumeFace& scvf(GridIndexType scvfIdx) const
    {
        return scvfs_[scvfIdx];
    }

    //! Get the scvf on the same face but from the other side
    //! Note that e.g. the normals might be different in the case of surface grids
    const SubControlVolumeFace& flipScvf(GridIndexType scvfIdx, std::size_t outsideScvfIdx = 0) const
    {
        return scvfs_[flipScvfIndices_[scvfIdx][outsideScvfIdx]];
    }

    //! Get the sub control volume face indices of an scv by global index
    const std::vector<GridIndexType>& scvfIndicesOfScv(GridIndexType scvIdx) const
    {
        return scvfIndicesOfScv_[scvIdx];
    }

    /*!
     * \brief Returns the connectivity map of which dofs have derivatives with respect
     *        to a given dof.
     */
    const ConnectivityMap &connectivityMap() const
    { return connectivityMap_; }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf(GridIndexType eIdx) const
    { return hasBoundaryScvf_[eIdx]; }

private:
    // find the scvf that has insideScvIdx in its outsideScvIdx list and outsideScvIdx as its insideScvIdx
    GridIndexType findFlippedScvfIndex_(GridIndexType insideScvIdx, GridIndexType outsideScvIdx)
    {
        // go over all potential scvfs of the outside scv
        for (auto outsideScvfIndex : scvfIndicesOfScv_[outsideScvIdx])
        {
            const auto& outsideScvf = this->scvf(outsideScvfIndex);
            for (int j = 0; j < outsideScvf.numOutsideScvs(); ++j)
                if (outsideScvf.outsideScvIdx(j) == insideScvIdx)
                    return outsideScvf.index();
        }

        DUNE_THROW(Dune::InvalidStateException, "No flipped version of this scvf found!");
    }

    //! connectivity map for efficient assembly
    ConnectivityMap connectivityMap_;
    std::unordered_map<GridIndexType, std::vector<GridIndexType>> extraConnectivity_;

    //! containers storing the global data
    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;
    std::vector<std::vector<GridIndexType>> scvfIndicesOfScv_;
    std::size_t numBoundaryScvf_;
    std::vector<bool> hasBoundaryScvf_;

    //! needed for embedded surface and network grids (dim < dimWorld)
    std::vector<std::vector<GridIndexType>> flipScvfIndices_;
};

} // end namespace Dumux

#endif
