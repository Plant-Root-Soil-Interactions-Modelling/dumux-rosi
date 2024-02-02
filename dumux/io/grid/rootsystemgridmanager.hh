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
 * \brief A grid factory using a CRootBox root system
 */
#ifndef DUMUX_ROOTSYSTEM_GRIDMANAGER_HH
#define DUMUX_ROOTSYSTEM_GRIDMANAGER_HH

#include <unordered_map>

#include <dune/common/fvector.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dumux/growth/crootboxinterface.hh>
#include <dumux/growth/rootparameters.hh>
#include <dumux/periodic/periodicnetworktransform.hh>

#include <RootSystem.h>

#if HAVE_DUNE_FOAMGRID

namespace Dumux {

template<class Grid>
class RootSystemGridData
{
    using IndexType = typename Grid::LeafGridView::IndexSet::IndexType;
    using Element = typename Grid::template Codim<0>::Entity;
public:
    //! export type of root parameters
    using RootParameter = GrowthModule::RootParams<double>;

    //! export type of element parameters
    struct ElementParameter
    {
        int id = {};
        double age = {};
        double radius = {};
        double order = {};
        bool isShootElement = false;
    };

    //! constructor (assumes all the data is moved in here)
    RootSystemGridData(std::unordered_map<IndexType, std::vector<IndexType>>&& connectivity,
                       std::vector<ElementParameter>&& elementParams,
                       std::vector<RootParameter>&& rootParams,
                       RootParameter&& shootParams,
                       std::shared_ptr<const Grid> grid,
                       Dune::GridFactory<Grid>&& gridFactory)
    : periodicConnectivity_(std::move(connectivity))
    , elementParams_(std::move(elementParams))
    , rootParams_(std::move(rootParams))
    , shootParams_(std::move(shootParams))
    , grid_(grid)
    , gridFactory_(std::move(gridFactory))
    {}

    //! get element parameters from CRootBox::RootSystem
    const ElementParameter& elementParameters(const Element& element) const
    { return elementParams_[gridFactory_.insertionIndex(element)]; }

    //! get the root parameters from CRootBox::RootSystem
    const std::vector<RootParameter>& rootParameters() const
    { return rootParams_; }

    //! get the root parameters from CRootBox::RootSystem
    const RootParameter& shootParameters() const
    { return shootParams_; }

    //! create the periodic vertex set given mappers
    template<class ElementMapper, class VertexMapper>
    std::unordered_map<IndexType, std::vector<IndexType>>
    createPeriodicConnectivity(const ElementMapper& elementMapper, const VertexMapper& vertexMapper) const
    {
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

    std::vector<ElementParameter> elementParams_; //!< parameter given for a single element
    std::vector<RootParameter> rootParams_; //!< root data for whole roots (possibly several elements until branch)
    RootParameter shootParams_; //!< root data for whole shoot roots (possibly several elements until branch)

    std::shared_ptr<const Grid> grid_;
    Dune::GridFactory<Grid> gridFactory_;
};

/*!
 * \brief Builds a Dune::FoamGrid<1, 3> from a CRootBox::RootSystem
 */
class RootSystemGridManager
{
    static constexpr int dim = 1;
    static constexpr int dimWorld = 3;
    using GridType = Dune::FoamGrid<1, 3>;
    using Element = GridType::template Codim<0>::Entity;
    using GlobalCoordinate = Element::Geometry::GlobalCoordinate;
    using IndexType = GridType::LeafGridView::IndexSet::IndexType;

    using GridDataType = RootSystemGridData<GridType>;
    using RootParams = GridDataType::RootParameter;
    using ElementParams = GridDataType::ElementParameter;
    using VertexMarker = Dune::FieldVector<int, dimWorld>;

public:
    //! export grid type
    using Grid = GridType;
    //! export grid data type
    using GridData = GridDataType;

    /*!
     * \brief Initialize the grid manager (builds the grid and extracts the data)
     * \param rs the CRootBox::RootSystem to extract the grid from
     * \param verbose whether the output is verbose or not
     */
    void init(const CRootBox::RootSystem& rs, bool verbose = false)
    {
        // the grid factory creates the grid
        if (verbose) std::cout << "RootSystemGridFactory: " << std::endl;
        Dune::GridFactory<Grid> factory;

        const auto nodes = rs.getNodes();
        int counter = 0;
        for (const auto& n : nodes)
        {
            if (verbose) std::cout << "-- add vertex " << counter++ <<  " at " << n.toString() << std::endl;
            factory.insertVertex(CRootBoxInterface::convert(n));
        }

        // get shoot and root segments
        const auto shootSegments = rs.getShootSegments();
        const auto segments = rs.getSegments();
        const auto numSegments = shootSegments.size() + segments.size();

        // get some root data
        const auto segmentRoots = rs.getSegmentOrigins(); // the root origin for each segment
        const auto segmentAges = rs.getRootTips(); // the segment ages
        auto rootParams = getRootParams_(rs);
        auto shootParams = getShootParams_(rs);
        std::vector<ElementParams> elemData; elemData.reserve(numSegments);

        IndexType segmentIdx = 0;
        for(const auto& s : shootSegments)
        {
            if (verbose) std::cout << "-- add element with vertices " << s.toString() << std::endl;
            factory.insertElement(Dune::GeometryTypes::line, CRootBoxInterface::convert(s));
            elemData.emplace_back(getElementParams_(segmentIdx++, rs,
                                                    rootParams, segmentRoots, segmentAges,
                                                    shootSegments.size(), shootParams));
        }

        for(const auto& s : segments)
        {
            if (verbose) std::cout << "-- add element with vertices " << s.toString() << std::endl;
            factory.insertElement(Dune::GeometryTypes::line, CRootBoxInterface::convert(s));
            elemData.emplace_back(getElementParams_(segmentIdx++, rs,
                                                    rootParams, segmentRoots, segmentAges,
                                                    shootSegments.size(), shootParams));
        }

        // create the grid
        grid_ = std::shared_ptr<Grid>(factory.createGrid());

        // create the grid data object
        std::unordered_map<IndexType, std::vector<IndexType>> periodicConnectivity; // empty
        gridData_ = std::make_shared<GridData>(std::move(periodicConnectivity),
                                               std::move(elemData), std::move(rootParams), std::move(shootParams),
                                               grid_, std::move(factory));

        std::cout << "Created non-periodic root grid with " << grid_->leafGridView().size(dim) << " vertices and "
                  << grid_->leafGridView().size(0) << " elements." << std::endl;
    }

    /*!
     * \brief Initialize the grid manager (builds the grid and extracts the data)
     * \note this is the overload for constructing a periodic grid (given a periodic cell)
     * \param rs the CRootBox::RootSystem to extract the grid from
     * \param lowerLeft the lower left corner of the periodic cell
     * \param upperRight the upper right corner of the periodic cell
     * \param periodic in which direction the cells is periodic
     * \param verbose whether the output is verbose or not
     */
    void init(const CRootBox::RootSystem& rs,
              const GlobalCoordinate& lowerLeft,
              const GlobalCoordinate& upperRight,
              const std::bitset<dimWorld>& periodic = std::bitset<dimWorld>{},
              bool verbose = false)
    {
        // if there is no periodicity simply forward to the non-periodic member function overload
        if (periodic.none())
            return init(rs, verbose);

        // create a helper class for mapping the grid onto a periodic cell
        PeriodicNetworkTransform<GlobalCoordinate> transformation(lowerLeft, upperRight, periodic);

        // get nodes from CRootBox::RootSystem
        const auto nodes = rs.getNodes();

        // create vertex data structures
        // -- verts, the vector of vertex positions
        // -- vertexMarker, the markers denoting in which periodic cell the vertex is (used to compute the shift)
        // -- faceIndices, unique face indices identifying periodic vertices as one face
        // -- isPeriodicFace, face markers that indentify faces as periodic
        // -- faceIndexCounter, index counter if new faces are created during the mapping
        std::vector<GlobalCoordinate> verts(nodes.size());
        std::vector<VertexMarker> vertexMarker(nodes.size());
        std::vector<IndexType> faceIndices(nodes.size());
        std::vector<bool> isPeriodicFace(nodes.size(), false);
        IndexType faceIndexCounter = verts.size();
        // allocate enough memory
        verts.reserve(verts.size() + verts.size()/10 + 5);
        vertexMarker.reserve(verts.size() + verts.size()/10 + 5);
        faceIndices.reserve(verts.size() + verts.size()/10 + 5);
        isPeriodicFace.reserve(verts.size() + verts.size()/10 + 5);

        // convert the node list to Dune data structure
        for (IndexType nodeIdx = 0; nodeIdx < nodes.size(); ++nodeIdx)
        {
            verts[nodeIdx] = CRootBoxInterface::convert(nodes[nodeIdx]);
            vertexMarker[nodeIdx] = transformation.getVertexMarker(verts[nodeIdx]);
            faceIndices[nodeIdx] = nodeIdx;
        }

        // get shoot segments and regular segments from CRootBox::RootSystem
        const auto shootSegments = rs.getShootSegments();
        const auto segments = rs.getSegments();
        const auto numSegments = shootSegments.size() + segments.size();

        // get some root data
        const auto segmentRoots = rs.getSegmentOrigins(); // the root origin for each segment
        const auto segmentAges = rs.getRootTips(); // the segment ages
        auto rootParams = getRootParams_(rs);
        auto shootParams = getShootParams_(rs);

        // create element data structures
        // -- elems, the elements as vertex lists
        // -- elemData, data associated with the elements
        std::vector<std::vector<unsigned int>> elems;
        std::vector<ElementParams> elemData;
        // allocate enough memory
        elems.reserve(numSegments + numSegments/10 + 5);
        elemData.reserve(numSegments + numSegments/10 + 5);

        IndexType segmentIdx = 0;
        auto handleSegment = [&](const auto& s)
        {
            const auto corners = CRootBoxInterface::convert(s);
            // the element doesn't intersect any periodic boundary -> insert unchanged
            if (vertexMarker[corners[0]] == vertexMarker[corners[1]])
            {
                elems.emplace_back(std::move(corners));
                elemData.emplace_back(getElementParams_(segmentIdx, rs,
                                                        rootParams, segmentRoots, segmentAges,
                                                        shootSegments.size(), shootParams));
            }

            // the element is intersecting a periodic boundary,
            // split at intersection point and insert new vertex
            else
            {
                // find out if one of the vertices is exactly on the boundary,
                // if yes, add only one vertex with the marker of the other vertex
                bool v0OnBoundary = transformation.onBoundary(verts[corners[0]] - transformation.getShift(vertexMarker[corners[0]]));
                bool v1OnBoundary = transformation.onBoundary(verts[corners[1]] - transformation.getShift(vertexMarker[corners[1]]));

                if (v0OnBoundary && v1OnBoundary)
                    DUNE_THROW(Dune::InvalidStateException, "Element is too large");
                else if (v0OnBoundary)
                {
                    verts.emplace_back(verts[corners[0]]);
                    vertexMarker.emplace_back(vertexMarker[corners[1]]);
                    faceIndices.push_back(corners[0]);
                    isPeriodicFace[corners[0]] = true;
                    isPeriodicFace.emplace_back(true);
                    elems.emplace_back(std::vector<unsigned int>{static_cast<unsigned int>(verts.size()-1), corners[1]});
                    elemData.emplace_back(getElementParams_(segmentIdx, rs,
                                                            rootParams, segmentRoots, segmentAges,
                                                            shootSegments.size(), shootParams));
                }
                else if (v1OnBoundary)
                {
                    verts.emplace_back(verts[corners[1]]);
                    vertexMarker.emplace_back(vertexMarker[corners[0]]);
                    faceIndices.push_back(corners[1]);
                    isPeriodicFace[corners[1]] = true;
                    isPeriodicFace.emplace_back(true);
                    elems.emplace_back(std::vector<unsigned int>{corners[0], static_cast<unsigned int>(verts.size()-1)});
                    elemData.emplace_back(getElementParams_(segmentIdx, rs,
                                                            rootParams, segmentRoots, segmentAges,
                                                            shootSegments.size(), shootParams));
                }

                // otherwise we compute the intersection point and
                // add two vertices, one for each of the markers
                else
                {
                    // get periodic box of corner 0
                    const auto box = transformation.getPeriodicBBox(vertexMarker[corners[0]]);
                    // compute the intersection point using directed ray box intersection
                    const auto iPos = transformation.intersectRayBox(box.first, box.second, verts[corners[0]], verts[corners[1]]-verts[corners[0]]);

                    verts.emplace_back(iPos);
                    vertexMarker.emplace_back(vertexMarker[corners[0]]);
                    faceIndices.emplace_back(faceIndexCounter);
                    isPeriodicFace.emplace_back(true);
                    elems.emplace_back(std::vector<unsigned int>{corners[0], static_cast<unsigned int>(verts.size()-1)});
                    elemData.emplace_back(getElementParams_(segmentIdx, rs,
                                                            rootParams, segmentRoots, segmentAges,
                                                            shootSegments.size(), shootParams));

                    verts.emplace_back(iPos);
                    vertexMarker.emplace_back(vertexMarker[corners[1]]);
                    faceIndices.emplace_back(faceIndexCounter);
                    isPeriodicFace.emplace_back(true);
                    elems.emplace_back(std::vector<unsigned int>{static_cast<unsigned int>(verts.size()-1), corners[1]});
                    elemData.emplace_back(getElementParams_(segmentIdx, rs,
                                                            rootParams, segmentRoots, segmentAges,
                                                            shootSegments.size(), shootParams));

                    ++faceIndexCounter;
                }
            }

            ++segmentIdx;
        };

        // first handle the shoot segments
        for (const auto& s : shootSegments)
            handleSegment(s);

        // then the same for the normal segments
        for (const auto& s : segments)
            handleSegment(s);

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

        // the grid factory creates the grid
        Dune::GridFactory<Grid> factory;

        // move all vertices by periodic shift
        for (unsigned int vIdx = 0; vIdx < verts.size(); ++vIdx)
            factory.insertVertex(verts[vIdx] - transformation.getShift(vertexMarker[vIdx]));

        for (const auto& element : elems)
            factory.insertElement(Dune::GeometryTypes::line, element);

        // create the grid
        grid_ = std::shared_ptr<Grid>(factory.createGrid());

        // create the grid data object
        gridData_ = std::make_shared<GridData>(std::move(periodicConnectivity),
                                               std::move(elemData), std::move(rootParams), std::move(shootParams),
                                               grid_, std::move(factory));

        std::cout << "Created periodic grid with " << grid_->leafGridView().size(dim) << " vertices and "
                  << grid_->leafGridView().size(0) << " elements." << std::endl;
    }

    /*!
     * \brief a reference to the grid
     */
    Grid& grid()
    { return *grid_; }

    /*!
     * \brief a shared_ptr to the grid
     */
    std::shared_ptr<Grid> gridPtr()
    { return grid_; }

    /*!
     * \brief create a shared pointer to the grid data
     */
    std::shared_ptr<GridData> getGridData() const
    { return gridData_; }

    /*!
     * \brief Distributes the grid over all processes for a parallel computation.
     */
    void loadBalance()
    {
        DUNE_THROW(Dune::NotImplemented, "Parallel (periodic) network grids");
    }

private:
    //! extract the root params from the root system
    std::vector<RootParams> getRootParams_(const CRootBox::RootSystem& rs) const
    {
        // be careful: the total number of roots and the grown number of roots is not the same
        std::vector<RootParams> rootParams(rs.getNumberOfRoots(/*allRoots=*/true));
        // this gets the root that already have at least one segment
        const auto roots = rs.getRoots();
        for (auto&& root : roots)
        {
            auto& rp = rootParams[root->getId()];

            // compute the root order
            rp.order = 0;
            auto* parentRoot = root;
            while (parentRoot->parent != nullptr)
            {
                ++rp.order;
                parentRoot = parentRoot->parent;
            }

            rp.radius = root->param.a * 0.01; // conversion from cm to m
            rp.rootId = root->getId();
            rp.plantId = 0;
        }

        return rootParams;
    }

    //! extract the shoot params from the root system
    RootParams getShootParams_(const CRootBox::RootSystem& rs) const
    {
        RootParams shootParams;
        shootParams.order = 0;
        shootParams.radius = 0.01;
        shootParams.rootId = -1;
        shootParams.plantId = 0;
        return shootParams;
    }

    //! extract the element params of a CRootBox segment
    ElementParams getElementParams_(unsigned int segmentIdx,
                                    const CRootBox::RootSystem& rs,
                                    const std::vector<RootParams>& rootParams,
                                    const std::vector<CRootBox::Root*>& segmentRoots,
                                    const std::vector<double>& segmentAges,
                                    std::size_t numShootSegments,
                                    const RootParams& shootParams) const
    {
        ElementParams params;

        // shoot segment come first and have all the same parameters
        if (segmentIdx < numShootSegments)
        {
            params.id = -1;
            params.age = rs.getSimTime(); // they where always there
            params.radius = shootParams.radius;
            params.order = shootParams.order;
            params.isShootElement = true;
        }

        // regular segments
        else
        {
            const auto sIdx = segmentIdx-numShootSegments;
            params.id = segmentRoots[sIdx]->getId();
            params.age = std::max(0.0, rs.getSimTime() - segmentAges[sIdx]);
            const auto& rp = rootParams[params.id];
            params.radius = rp.radius;
            params.order = rp.order;
            params.isShootElement = false;
        }

        return params;
    }

    std::shared_ptr<Grid> grid_; //!< a pointer to the grid
    std::shared_ptr<GridData> gridData_; //!< a grid data pointer
};

} // end namespace Dumux

#endif // HAVE_DUNE_FOAMGRID
#endif // DUMUX_ROOTSYSTEM_GRIDMANAGER_HH
