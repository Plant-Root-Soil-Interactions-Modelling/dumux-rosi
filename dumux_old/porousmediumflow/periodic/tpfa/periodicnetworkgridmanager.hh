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
#ifndef DUMUX_PERIODIC_NETWORK_GRID_MANAGER_HH
#define DUMUX_PERIODIC_NETWORK_GRID_MANAGER_HH

#include <dune/common/fvector.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/grid/griddata.hh>
#include <dumux/periodic/periodicnetworktransform.hh>

namespace Dumux {

#if HAVE_DUNE_FOAMGRID

/*!
 * \brief Wrapper around the 1d foam grid grid manager to make foamgrids periodic
 */
template<int dimWorld>
class PeriodicNetworkGridManager
{
    static constexpr int dim = 1;
    using GridType = Dune::FoamGrid<dim, dimWorld>;
    using IndexType = typename GridType::LeafGridView::IndexSet::IndexType;
    using Element = typename GridType::template Codim<0>::Entity;
    using GlobalCoordinate = typename Element::Geometry::GlobalCoordinate;
    using VertexMarker = Dune::FieldVector<int, dimWorld>;

    class GridDataType
    {
    public:
        GridDataType(std::unordered_map<IndexType, std::vector<IndexType>>&& connectivity,
                     std::vector<std::vector<double>>&& elementParams,
                     std::shared_ptr<const GridType> grid,
                     Dune::GridFactory<GridType>&& gridFactory)
        : periodicConnectivity_(std::move(connectivity))
        , elementParams_(std::move(elementParams))
        , grid_(grid)
        , gridFactory_(std::move(gridFactory))
        {}

        GridDataType(typename std::shared_ptr<GridData<GridType>> hostGridData)
        : hostGridData_(hostGridData)
        {}

        //! get element parameters from host grid dgf
        const std::vector<double>& parameters(const Element& element) const
        {
            if (!hostGridData_)
                return elementParams_[gridFactory_.insertionIndex(element)];
            else
                return hostGridData_->parameters(element);
        }

        //! create the periodic vertex set given mappers
        template<class ElementMapper, class VertexMapper>
        std::unordered_map<IndexType, std::vector<IndexType>>
        createPeriodicConnectivity(const ElementMapper& elementMapper, const VertexMapper& vertexMapper) const
        {
            if (periodicConnectivity_.empty())
                return periodicConnectivity_;

            std::vector<IndexType> leafElementIndex(grid_->leafGridView().size(0));
            for (const auto& element : elements(grid_->leafGridView()))
                leafElementIndex[gridFactory_.insertionIndex(element)] = elementMapper.index(element);

            std::vector<IndexType> leafVertexIndex(grid_->leafGridView().size(1));
            for (const auto& vertex : vertices(grid_->leafGridView()))
                leafVertexIndex[gridFactory_.insertionIndex(vertex)] = vertexMapper.index(vertex);

            std::unordered_map<IndexType, std::vector<IndexType>> periodicConnectivity;
            for (const auto& entry : periodicConnectivity_)
            {
                auto mappedElementIndices = entry.second;
                std::transform(mappedElementIndices.begin(), mappedElementIndices.end(), mappedElementIndices.begin(),
                               [&](const auto& insertionElementIndex) { return leafElementIndex[insertionElementIndex]; });
                periodicConnectivity[leafVertexIndex[entry.first]] = mappedElementIndices;
            }

            return periodicConnectivity;
        }

    private:
        const std::unordered_map<IndexType, std::vector<IndexType>> periodicConnectivity_;
        const std::vector<std::vector<double>> elementParams_;
        std::shared_ptr<const GridType> grid_;
        Dune::GridFactory<GridType> gridFactory_;

        std::shared_ptr<GridData<GridType>> hostGridData_;
    };

public:
    using Grid = GridType;
    using GridData = GridDataType;

    PeriodicNetworkGridManager(const GlobalCoordinate& lowerLeft,
                               const GlobalCoordinate& upperRight,
                               const std::bitset<dimWorld>& periodic)
    : transformation_(lowerLeft, upperRight, periodic)
    {}

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     * \param lowerLeft the lower left corner of the bounding box the grid should be restricted to
     * \param upperRight the upper right corner of the bounding box the grid should be restricted to
     * \param periodic if the grid is periodic in x and/or y and/or z
     * \param paramGroup the parameter group to preferably look in
     */
    void init(const std::string& paramGroup = "")
    {
        // first create the host grid
        hostGridManager_.init(paramGroup);

        // for non-periodic grid just forward to host grid
        if (transformation_.periodic().none())
        {
            gridData_ = std::make_shared<GridData>(hostGridManager_.getGridData());
            return;
        }

        // get the data, we also have to transfer that to the new grid
        auto gridData = hostGridManager_.getGridData();

        const auto& grid = hostGridManager_.grid();
        const auto& gridView = grid.leafGridView();

        std::cout << "Read non-periodic grid with " << gridView.size(dim) << " vertices and "
                  << gridView.size(0) << " elements." << std::endl;

        std::vector<GlobalCoordinate> verts(gridView.size(dim));
        std::vector<VertexMarker> vertexMarker(gridView.size(dim));
        std::vector<IndexType> faceIndices(gridView.size(dim));
        std::vector<bool> isPeriodicFace(gridView.size(dim), false);
        verts.reserve(verts.size() + verts.size()/10 + 5);
        vertexMarker.reserve(verts.size() + verts.size()/10 + 5);
        faceIndices.reserve(verts.size() + verts.size()/10 + 5);
        isPeriodicFace.reserve(verts.size() + verts.size()/10 + 5);

        Dune::MultipleCodimMultipleGeomTypeMapper<typename Grid::LeafGridView>
            vertexMapper(gridView, Dune::mcmgVertexLayout());

        // add all vertices and compute the marker
        // the marker signifies uniquely in which periodic cell the vertex belongs
        // the marker is the periodic shift, i.e. 1 0 0 is a shift of 1-times the periodic cell length in x direction
        for (const auto& vertex : vertices(gridView))
        {
            const auto vIdx = vertexMapper.index(vertex);
            verts[vIdx] = vertex.geometry().corner(0);
            vertexMarker[vIdx] = transformation_.getVertexMarker(verts[vIdx]);
            faceIndices[vIdx] = vIdx;
        }

        IndexType faceIndexCounter = verts.size();

        using std::abs;
        std::vector<std::vector<unsigned int>> elems;
        std::vector<std::vector<double>> elemData;
        elems.reserve(gridView.size(0) + gridView.size(0)/10 + 5);
        elemData.reserve(gridView.size(0) + gridView.size(0)/10 + 5);
        for (const auto& element : elements(gridView))
        {
            const unsigned int corner0Idx = vertexMapper.subIndex(element, 0, dim);
            const unsigned int corner1Idx = vertexMapper.subIndex(element, 1, dim);

            // the element doesn't intersect any periodic boundary -> insert unchanged
            if (vertexMarker[corner0Idx] == vertexMarker[corner1Idx])
            {
                elems.emplace_back(std::vector<unsigned int>{corner0Idx, corner1Idx});
                elemData.emplace_back(gridData->parameters(element));
            }

            // the element is intersecting a periodic boundary,
            // split at intersection point and insert new vertex
            else
            {
                // find out if one of the vertices is exactly on the boundary,
                // if yes, add only one vertex with the marker of the other vertex
                bool v0OnBoundary = transformation_.onBoundary(verts[corner0Idx] - transformation_.getShift(vertexMarker[corner0Idx]));
                bool v1OnBoundary = transformation_.onBoundary(verts[corner1Idx] - transformation_.getShift(vertexMarker[corner1Idx]));

                if (v0OnBoundary && v1OnBoundary)
                    DUNE_THROW(Dune::InvalidStateException, "Element is too large");
                else if (v0OnBoundary)
                {
                    verts.emplace_back(verts[corner0Idx]);
                    vertexMarker.emplace_back(vertexMarker[corner1Idx]);
                    faceIndices.push_back(corner0Idx);
                    isPeriodicFace[corner0Idx] = true;
                    isPeriodicFace.emplace_back(true);
                    elems.emplace_back(std::vector<unsigned int>{static_cast<unsigned int>(verts.size()-1), corner1Idx});
                    elemData.emplace_back(gridData->parameters(element));
                }
                else if (v1OnBoundary)
                {
                    verts.emplace_back(verts[corner1Idx]);
                    vertexMarker.emplace_back(vertexMarker[corner0Idx]);
                    faceIndices.push_back(corner1Idx);
                    isPeriodicFace[corner1Idx] = true;
                    isPeriodicFace.emplace_back(true);
                    elems.emplace_back(std::vector<unsigned int>{corner0Idx, static_cast<unsigned int>(verts.size()-1)});
                    elemData.emplace_back(gridData->parameters(element));
                }

                // otherwise we compute the intersection point and
                // add two vertices, one for each of the markers
                else
                {
                    // get periodic box of corner 0
                    const auto box = transformation_.getPeriodicBBox(vertexMarker[corner0Idx]);
                    // find intersection point
                    const auto iPos = transformation_.intersectRayBox(box.first, box.second, verts[corner0Idx], verts[corner1Idx]-verts[corner0Idx]);

                    verts.emplace_back(iPos);
                    vertexMarker.emplace_back(vertexMarker[corner0Idx]);
                    faceIndices.emplace_back(faceIndexCounter);
                    isPeriodicFace.emplace_back(true);
                    elems.emplace_back(std::vector<unsigned int>{corner0Idx, static_cast<unsigned int>(verts.size()-1)});
                    elemData.emplace_back(gridData->parameters(element));

                    // some data that is volume based has to be transformed too
                    elemData.back()[2] *= (verts[corner0Idx]-iPos).two_norm()/element.geometry().volume();

                    verts.emplace_back(iPos);
                    vertexMarker.emplace_back(vertexMarker[corner1Idx]);
                    faceIndices.emplace_back(faceIndexCounter);
                    isPeriodicFace.emplace_back(true);
                    elems.emplace_back(std::vector<unsigned int>{static_cast<unsigned int>(verts.size()-1), corner1Idx});
                    elemData.emplace_back(gridData->parameters(element));

                    // some data that is volume based has to be transformed too
                    elemData.back()[2] *= (verts[corner1Idx]-iPos).two_norm()/element.geometry().volume();

                    ++faceIndexCounter;
                }
            }
        }

        // compute connectivities
        std::vector<std::vector<IndexType>> faceElements(verts.size());
        for (std::size_t eIdx = 0; eIdx < elems.size(); ++eIdx)
        {
            faceElements[faceIndices[elems[eIdx][0]]].push_back(eIdx);
            faceElements[faceIndices[elems[eIdx][1]]].push_back(eIdx);
        }

        std::unordered_map<IndexType, std::vector<IndexType>> periodicConnectivity;
        for (std::size_t vIdx = 0; vIdx < verts.size(); ++vIdx)
        {
            if (isPeriodicFace[faceIndices[vIdx]])
            {
                auto& vertexConnectivity = periodicConnectivity[vIdx];
                for (const auto eIdx : faceElements[faceIndices[vIdx]])
                    vertexConnectivity.push_back(eIdx);
            }
        }

        // create grid
        Dune::GridFactory<Grid> factory;

        // move all vertices by periodic shift
        for (unsigned int vIdx = 0; vIdx < verts.size(); ++vIdx)
            factory.insertVertex(verts[vIdx] - transformation_.getShift(vertexMarker[vIdx]));

        for (const auto& element : elems)
            factory.insertElement(Dune::GeometryTypes::line, element);

        grid_ = std::shared_ptr<Grid>(factory.createGrid());
        gridData_ = std::make_shared<GridData>(std::move(periodicConnectivity), std::move(elemData), grid_, std::move(factory));

        std::cout << "Created periodic grid with " << grid_->leafGridView().size(dim) << " vertices and "
                  << grid_->leafGridView().size(0) << " elements." << std::endl;
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    Grid& grid()
    {
        if (grid_)
            return *grid_;
        else
            return hostGridManager_.grid();
    }

    /*!
     * \brief get a shared_ptr to the grid data
     */
    std::shared_ptr<GridData> getGridData() const
    { return gridData_; }

private:
    PeriodicNetworkTransform<GlobalCoordinate> transformation_;
    std::shared_ptr<Grid> grid_;
    GridManager<Grid> hostGridManager_;

    std::shared_ptr<GridData> gridData_;
};

#endif // HAVE_DUNE_FOAMGRID

} // end namespace Dumux

#endif
