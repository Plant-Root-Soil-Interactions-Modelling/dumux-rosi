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
 * \brief This file contains a growth algorithm using crootbox as a backend
 */
#ifndef DUMUX_CROOTBOX_GROWTH_ALGORITHM_HH
#define DUMUX_CROOTBOX_GROWTH_ALGORITHM_HH

#include <memory>

#include <dune/common/version.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/exceptions.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/entitymap.hh>
#include <dumux/growth/algorithm.hh>
#include <dumux/growth/crootboxinterface.hh>
#include <RootSystem.h>

namespace Dumux {

namespace GrowthModule {

template<class TypeTag>
class CRootBoxGrowthAlgorithm : public GrowthAlgorithm
{
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

    struct Data
    {
        PrimaryVariables priVars;
        unsigned int rootId;
        unsigned int rootOrder;
        double age;
        double previousLength;
    };

    using PersistentContainer = Dune::PersistentContainer<Grid, Data>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using CRB = Dumux::GrowthModule::CRootBoxInterface;

public:
    //! contruct the rootbox interface growth algorithm with a root system
    CRootBoxGrowthAlgorithm(std::shared_ptr<Grid> grid,
                            std::shared_ptr<FVGridGeometry> fvGridGeometry,
                            std::shared_ptr<CRootBox::RootSystem> rs,
                            Problem& problem,
                            SolutionVector& sol)
    : grid_(grid)
    , fvGridGeometry_(fvGridGeometry)
    , rootSystem_(rs)
    , indexToVertex_(*grid, fvGridGeometry->vertexMapper())
    , data_(*grid, 0)
    , problem_(problem)
    , sol_(sol)
    {
        const auto& gv = grid->leafGridView();
        indexMap_.resize(gv.size(Grid::dimension));
        //! TODO: we assume that in the beginning the node indices are the same
        std::iota(indexMap_.begin(), indexMap_.end(), 0);
    }

    //! do a growth step with step size dt
    //! returns how many element were inserted, how many removed, and how many vertices moved.
    //! \param dt the time step size in seconds
    StatusReport grow(double dt) final
    {
        //! keep track of how much we inserted
        unsigned int insertedSegments = 0;
        unsigned int movedNodes = 0;

        // we work on the leaf grid view
        const auto& gv = grid_->leafGridView();

        // remember the old segment and vertex amount
        const auto oldNumSegments = gv.size(0);
        const auto oldNumVertices = gv.size(Grid::dimension);

        // store the old grid data
        storeData_();

        //! grow the crootbox root system
        rootSystem_->simulate(dt*CRB::sToDays);

        //! get nodes that changed their position
        const auto updatedNodeIndices = rootSystem_->getUpdatedNodeIndices();
        const auto updatedNodes = rootSystem_->getUpdatedNodes();
        assert(updatedNodeIndices.size() == updatedNodes.size());

        //! change their position in dune-foamgrid too
        if (!updatedNodes.empty())
        {
            for (unsigned int nIdx = 0; nIdx < updatedNodeIndices.size(); ++nIdx)
                grid_->setPosition(indexToVertex_[getDuneIndex_(updatedNodeIndices[nIdx])],
                                   CRB::convert(updatedNodes[nIdx]));
            // log how many have been moved
            movedNodes = updatedNodes.size();
        }

        //! nodes and segments that have be added
        const auto newNodes = rootSystem_->getNewNodes();
        const auto newNodeIndices = rootSystem_->getNewNodeIndices();
        const auto newSegments = rootSystem_->getNewSegments();
        // for each segment we have a point to the root they belong to
        const auto newSegmentRoots = rootSystem_->getNewSegmentsOrigin();
        const auto newSegmentAges = rootSystem_->getNewNETimes();

        //! add them to the dune-foamgrid too
        if (!newSegments.empty())
        {
            // register new vertices
            std::vector<unsigned int> newNodeIndicesDune(newNodes.size());
            for (unsigned int nIdx = 0; nIdx < newNodes.size(); ++nIdx)
            {
                newNodeIndicesDune[nIdx] = grid_->insertVertex(CRB::convert(newNodes[nIdx]));

                //! Do a check that newly inserted vertices are numbered consecutively in both dune and crootbox
                if(newNodeIndicesDune[nIdx] != newNodeIndices[nIdx])
                    DUNE_THROW(Dune::GridError, "Nodes are not consecutively numbered!");
            }

            // register the new segments and their parameters and their parent root params if non-existent
            problem_.spatialParams().rootParamsResize(rootSystem_->getNumberOfRoots(/*allRoots=*/true));
            newSegmentRootIds_.resize(newSegments.size());
            newSegmentAges_.resize(newSegments.size());
            for (unsigned int sIdx = 0; sIdx < newSegments.size(); ++sIdx)
            {
                const auto* root = newSegmentRoots[sIdx];
                newSegmentRootIds_[sIdx] = root->id;
                newSegmentAges_[sIdx] = std::max(0.0, rootSystem_->getSimTime()-newSegmentAges[sIdx]);
                const auto& s = newSegments[sIdx];
                grid_->insertElement(Dune::GeometryTypes::line, {getDuneIndex_(s.x), getDuneIndex_(s.y)});

                // if the root parameters are not set yet update them here
                auto& rp = problem_.spatialParams().getRootParams(newSegmentRootIds_[sIdx]);
                if (!rp.isInitialized)
                {
                    // compute root order
                    rp.order = 0;
                    auto* r_ = root;
                    while (r_->parent != nullptr)
                    {
                        ++rp.order;
                        r_=r_->parent;
                    }
                    rp.rootId = root->id;
                    rp.radius = root->param.a * 0.01; // conversion from cm to m
                    rp.plantId = 0;
                    rp.isInitialized = true;
                }
            }

            // log how many we inserted
            insertedSegments = newSegments.size();

            // this acutally inserts the new vertices and segments
            grid_->preGrow();
            grid_->grow();

            // update grid geometry (updates the mappers)
            fvGridGeometry_->update();

            // Update the index map
            indexMap_.resize(gv.size(Grid::dimension));
            // update the actual index the new vertices got in the grid in the crootbox index to dune index map
            for (const auto& vertex : vertices(gv))
                if (fvGridGeometry_->vertexMapper().index(vertex) >= oldNumVertices)
                    indexMap_[grid_->growthInsertionIndex(vertex)] = fvGridGeometry_->vertexMapper().index(vertex);

            // also update the index to vertex map
            indexToVertex_.update(fvGridGeometry_->vertexMapper());

            // clear the markers in the problem
            problem_.clearIsNewMarkers();

            // restore the old grid data and initialize new data
            reconstructData_(dt);

            // cleanup grid after growth
            grid_->postGrow();
        }

        // check if all segments could be inserted
        if (oldNumSegments + newSegments.size() != gv.size(0))
            DUNE_THROW(Dune::GridError, "Not all segments could be inserted!");

        return StatusReport{insertedSegments, /*removedSegments=*/0, movedNodes};
    }

private:

    unsigned int getDuneIndex_(int cRootBoxIdx)
    {
        if (cRootBoxIdx >= indexMap_.size())
            return static_cast<unsigned int>(cRootBoxIdx); //! new vertices have to same numbering during insertion
        else
            return indexMap_[cRootBoxIdx];
    }

    /*!
     * Store the data of the old grid
     */
    void storeData_()
    {
        data_.resize();
        const auto& gv = grid_->leafGridView();
        for (const auto& element : elements(gv))
        {
            // put your value in the map
            const auto eIdx = fvGridGeometry_->elementMapper().index(element);
            const auto& rp = problem_.spatialParams().rootParams(element);
            auto& data = data_[element];
            data.priVars = sol_[eIdx];
            data.age = problem_.spatialParams().age(eIdx);
            data.rootId = rp.rootId;
            data.rootOrder = rp.order;
            data.previousLength = problem_.spatialParams().previousLength(eIdx);
        }
    }

    /*!
     * Reconstruct missing primary variables (where elements are created/deleted)
     * \param problem The current problem
     */
    void reconstructData_(double dt)
    {
        data_.resize();
        const auto& gv = grid_->leafGridView();
        sol_.resize(gv.size(0));
        problem_.spatialParams().resizeSegmentParams(gv.size(0));

        for (const auto& element : elements(gv))
        {
            // old elements get their old variables assigned
            if(!element.isNew())
            {
                // get your primary variables from the map
                const auto eIdx = fvGridGeometry_->elementMapper().index(element);
                const auto& data = data_[element];
                sol_[eIdx] = data.priVars;
                // set the root id of the old element
                problem_.spatialParams().setRootId(eIdx, data.rootId);
                // set the previous length of the old element
                problem_.spatialParams().setPreviousLength(eIdx, data.previousLength);
                // set the age
                problem_.spatialParams().setAge(eIdx, data.age + dt*CRB::sToDays);
            }

            // initialize new elements with the pressure of the closest old element
            else
            {
                problem_.setIsNewMarker(element, true);
                const auto newEIdx = fvGridGeometry_->elementMapper().index(element);
                const auto insertionIdx = grid_->growthInsertionIndex(element);
                const auto newRootId = newSegmentRootIds_[insertionIdx];
                const auto newRootAge = newSegmentAges_[insertionIdx];
                // set the root id of the new element
                problem_.spatialParams().setRootId(newEIdx, newRootId);
                // set the previous length to zero
                problem_.spatialParams().setPreviousLength(newEIdx, 0.0);
                // set the age
                problem_.spatialParams().setAge(newEIdx, newRootAge);

                // initialize the oldElement to be found with the new element and walk up the
                // root hierarchy until an old element is found
                auto oldRootOrder = problem_.spatialParams().getRootParams(newRootId).order;
                auto oldElement = element;
                while (oldElement.isNew())
                {
                    for (const auto& intersection : intersections(gv, oldElement))
                    {
                        if (intersection.neighbor())
                        {
                            // potential up-branch element
                            auto nextElement = intersection.outside();

                            if (nextElement.isNew())
                            {
                                // skip new neighbor elements on the same branch that were inserted later
                                // i.e. further towards the branch tip
                                const auto nextIdx = grid_->growthInsertionIndex(nextElement);
                                if (nextIdx > insertionIdx)
                                    continue;

                                // or skip if it's an element of a daughter branch
                                const auto& nextRootParams = problem_.spatialParams().getRootParams(newSegmentRootIds_[nextIdx]);
                                if (nextRootParams.order > oldRootOrder)
                                    continue;

                                // we found the next up branch element
                                oldRootOrder = nextRootParams.order;
                                oldElement = std::move(nextElement);
                                break;
                            }

                            // skip old neighbor elements of daughter branches
                            else
                            {
                                // we get the order of an old element from the persistent container
                                const auto nextRootOrder = data_[nextElement].rootOrder;
                                if (nextRootOrder > oldRootOrder)
                                    continue;

                                // we found the final element we can extract the privars from
                                oldRootOrder = nextRootOrder;
                                oldElement = std::move(nextElement);
                                break;
                            }
                        }
                    }
                }

                if (fvGridGeometry_->elementMapper().index(oldElement) == newEIdx)
                    DUNE_THROW(Dune::InvalidStateException, "No old element was found to extract the primary variables for element: " << newEIdx);

                // extract the privars from the closest old element we found above
                sol_[newEIdx] = data_[oldElement].priVars;
            }
        }

        // reset entries in restrictionmap
        data_.resize(typename PersistentContainer::Value());
        data_.shrinkToFit();
        data_.fill(typename PersistentContainer::Value());

        // clear the root id container for the newt growth step
        newSegmentRootIds_.clear();
        newSegmentAges_.clear();
    }

    std::shared_ptr<Grid> grid_; //! the dune-foamgrid
    std::shared_ptr<FVGridGeometry> fvGridGeometry_; //! fv grid geometry
    std::shared_ptr<CRootBox::RootSystem> rootSystem_; //! the croot box root system

    // list of root ids for new roots
    std::vector<unsigned int> newSegmentRootIds_;
    std::vector<double> newSegmentAges_;

    // index mapping stuff
    EntityMap<Grid, Grid::dimension> indexToVertex_; //! a map from dune vertex indices to vertex entitites
    std::vector<int> indexMap_; //! a map from crootbox node indices to foamgrid indices

    PersistentContainer data_; //! data container with persistent ids as key to transfer element data from old to new grid

    // the data (non-const) to transfer from the old to the new grid
    Problem& problem_;
    SolutionVector& sol_;
};

} // end namespace GridGrowth

} // end namespace Dumux

#endif
