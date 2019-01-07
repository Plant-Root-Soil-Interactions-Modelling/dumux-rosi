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

#include "growthinterface.hh" // holds base class, and crootbox interface

namespace Dumux {

namespace GrowthModule {

/**
 * A class managing the growth of the grid (Dune::FoamGrid<1,3>) according to a some grow model GrowInterface
 */
template<class TypeTag>
class GridGrowth {

    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using GlobalPosition = typename GET_PROP_TYPE(TypeTag, GlobalPosition);
    using PersistentContainer = Dune::PersistentContainer<Grid, PrimaryVariables>; // todo remove macros

public:

    //! contruct the rootbox interface growth algorithm with a root system
    GridGrowth(std::shared_ptr<Grid> grid, std::shared_ptr<FVGridGeometry> fvGridGeometry, std::shared_ptr<GrowInterface> grow, SolutionVector& sol) :
        grid_(grid), fvGridGeometry_(fvGridGeometry), growInterface_(grow), indexToVertex_(*grid, fvGridGeometry->vertexMapper()), data_(*grid, 0), sol_(sol) {
        //const auto& gv = grid->leafGridView();
        // indexMap_.resize(gv.size(Grid::dimension));
        //! TODO: we assume that in the beginning the node indices are the same
        // std::iota(indexMap_.begin(), indexMap_.end(), 0);
    }

    //! \param dt the time step size in seconds
    void grow(double dt) {

        const auto& gv = grid_->leafGridView(); // we work on the leaf grid view

        //        // remember the old segment and vertex amount
        //        const auto oldNumSegments = gv.size(0);
        //        const auto oldNumVertices = gv.size(Grid::dimension);

        // store the old grid data
        storeData_();

        //! get nodes that changed their position
        auto updatedNodeIndices = growInterface_->updatedNodeIndices();
        auto updatedNodes = growInterface_->updatedNodes();
        assert((updatedNodeIndices.size() == updatedNodes.size()) && "GridGrowth: number of updated nodes <> updated indices"); // a little trust...

        //! change their position in dune-foamgrid too
        for (unsigned int i = 0; i < updatedNodeIndices.size(); i++) {
            grid_->setPosition(indexToVertex_[getDuneIndex_(updatedNodeIndices[i])], updatedNodes[nIdx]);
        }
        // i dont quite get what indexToVertex_ is for

        //! nodes and segments that have be added
        auto newNodes = growInterface_->newNodes();
        auto newNodeIndices = growInterface_->newNodeIndices();
        auto newSegments = growInterface_->newSegments();

        //! add them to the dune-foamgrid too
        if (!newSegments.empty()) {
            // register new vertices
            std::vector<size_t> newNodeIndicesDune(newNodes.size()); // todo actually, we do not need to store these
            for (size_t i = 0; i < newNodes.size(); i++) {
                newNodeIndicesDune[i] = grid_->insertVertex(newNodes[i]);
                //! Do a check that newly inserted vertices are numbered consecutively in both dune and crootbox
                if (newNodeIndicesDune[i] != newNodeIndices[i]) {
                    DUNE_THROW(Dune::GridError, "Nodes are not consecutively numbered!");
                }
            }

            for (size_t i = 0; i < newSegments.size(); i++) {
                const auto& s = newSegments[i];
                grid_->insertElement(Dune::GeometryTypes::line, { getDuneIndex_(s.x), getDuneIndex_(s.y) });
            }

            // this acutally inserts the new vertices and segments
            grid_->preGrow();
            grid_->grow();

            // update grid geometry (updates the mappers)
            fvGridGeometry_->update();

            // Update the index map
            indexMap_.resize(gv.size(Grid::dimension));
            // update the actual index the new vertices got in the grid in the crootbox index to dune index map
            for (const auto& vertex : vertices(gv)) {
                if (fvGridGeometry_->vertexMapper().index(vertex) >= oldNumVertices)
                    indexMap_[grid_->growthInsertionIndex(vertex)] = fvGridGeometry_->vertexMapper().index(vertex);
            }

            // also update the index to vertex map
            indexToVertex_.update(fvGridGeometry_->vertexMapper());

            // restore the old grid data and initialize new data
            reconstructData_(dt);

            // cleanup grid after growth
            grid_->postGrow();
        }

    }

private:

    /*!
     * dune vertex index from a crootbox node index
     */
    unsigned int getDuneIndex_(int cRootBoxIdx) {
        if (cRootBoxIdx >= indexMap_.size()) {
            return static_cast<unsigned int>(cRootBoxIdx); //! new vertices have to same numbering during insertion
        } else {
            return indexMap_[cRootBoxIdx];
        }
    }

    /*!
     * Store the primary variables of the old grid
     */
    void storeData_() {
        data_.resize();
        const auto& gv = grid_->leafGridView();
        for (const auto& element : elements(gv)) {
            const auto eIdx = fvGridGeometry_->elementMapper().index(element);
            data_[element] = sol_[eIdx];
        }
    }

    /*!
     * Reconstruct missing primary variables (where elements are created/deleted)
     * \param problem The current problem
     */
    void reconstructData_(double dt) {
        data_.resize();
        const auto& gv = grid_->leafGridView();
        sol_.resize(gv.size(0));

        for (const auto& element : elements(gv)) {
            // old elements get their old variables assigned
            if (!element.isNew()) { // get your primary variables from the map
                const auto eIdx = fvGridGeometry_->elementMapper().index(element);
                sol_[eIdx] = data_[element];
            } else { // initialize new elements with the pressure of the surrounding soil (todo)
                const auto newEIdx = fvGridGeometry_->elementMapper().index(element);
                const auto insertionIdx = grid_->growthInsertionIndex(element);
                sol_[newEIdx] = sol_[0];
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
    std::shared_ptr<GrowInterface<GlobalPosition>> growInterface_;

    //    // list of root ids for new roots
    //    std::vector<unsigned int> newSegmentRootIds_;
    //    std::vector<double> newSegmentAges_;

    // index mapping stuff
    EntityMap<Grid, Grid::dimension> indexToVertex_; //! a map from dune vertex indices to vertex entitites (?)
    std::vector<int> indexMap_; //! a map from crootbox node indices to foamgrid indices

    PersistentContainer data_; //! data container with persistent ids as key to transfer element data from old to new grid

    // the data (non-const) to transfer from the old to the new grid
    SolutionVector& sol_;
};

} // end namespace GridGrowth

} // end namespace Dumux

#endif
