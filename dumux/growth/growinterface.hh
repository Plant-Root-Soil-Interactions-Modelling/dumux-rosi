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
 * \brief This file contains a generic grid growth class
 */
#ifndef DUMUX_GROW_INTERFACE_HH
#define DUMUX_GROW_INTERFACE_HH

#include <RootSystem.h>

namespace Dumux {

namespace GrowthModule {

/**
 * Implement this abstract class to let the FoamGrid<1,3> grow
 */
template<class GlobalPosition>
class GrowInterface {
public:
    virtual ~GrowInterface() { };
    virtual std::vector<size_t> updatedNodeIndices() = 0; ///< Indices of nodes that were updated in the previous time step
    virtual std::vector<GlobalPosition> updatedNodes() = 0; ///< Values of the updated nodes
    virtual std::vector<int> newNodeIndices() = 0; ///< Node indices that were created in the previous time step
    virtual std::vector<GlobalPosition> newNodes() = 0; ///< Nodes created in the previous time step
    virtual std::vector<std::vector<size_t>> newSegments() = 0; ///< Segments created in the previous time step
    virtual std::vector<int> newSegmentOrders() = 0;
    virtual std::vector<double> newSegmentCreationTimes() = 0; ///< Copies a pointer to the root containing the new segments
    virtual std::vector<double> newSegmentRadii() = 0; ///< node ermergence times of segments created in the previous time step
};

/**
 * Implementation for CRootBox
 */
template<class GlobalPosition>
class CRootBoxInterface :public GrowInterface<GlobalPosition> {
public:
    CRootBoxInterface(CRootBox::RootSystem& rs) :rootsystem_(rs) { };
    virtual ~CRootBoxInterface() { }; // nothing to do, rootsystem_ is someone elses problem
    virtual std::vector<size_t> updatedNodeIndices() {
        auto nodes = rootsystem->getUpdatedNodeIndices();
        // todo write converter
        return nodes; // todo
    };
    // for all


private:
    CRootBox::RootSystem& rootsystem_;
};



} // end namespace GridGrowth

} // end namespace Dumux

#endif
