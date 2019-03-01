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
#ifndef DUMUX_GROWTH_INTERFACE_HH
#define DUMUX_GROWTH_INTERFACE_HH

namespace Dumux {

namespace GrowthModule {

/**
 * Implement this abstract class to let the FoamGrid<1,3> grow via GridGrowth::grow
 */
template<class GlobalPosition>
class GrowthInterface {

public:
    virtual ~GrowthInterface() { };

    // run simulation
    virtual void simulate(double dt) = 0; //< simulate the next dt seconds [s]

    // nodes that moved
    virtual std::vector<size_t> updatedNodeIndices() const = 0; //< Indices of nodes that were updated in the previous time step
    virtual std::vector<GlobalPosition> updatedNodes() const = 0; //< Values of the updated nodes [m]

    // new nodes
    virtual std::vector<size_t> newNodeIndices() const = 0; //< Node indices that were created in the previous time step
    virtual std::vector<GlobalPosition> newNodes() const = 0; //< Nodes created in the previous time step [m]

    // new segments
    virtual std::vector<std::array<size_t, 2>> newSegments() const = 0; //< Segments created in the previous time step

    // virtual double getModelData(size_t rIdx, size_t type) = 0;

    virtual std::vector<double> segmentCreationTimes() const = 0; //< Segment emergence time of segment sIdx created in the previous time step [s]
    virtual std::vector<int> segmentOrders() const = 0; //< Segment orders of segment created in the previous time step
    virtual std::vector<double> segmentRadii() const = 0; //< Segment radius of segment sIdx created in the previous time step [m]

    // Mapper
    std::vector<size_t> root2dune; ///< Maps a growth model index to the dune element index

    //! returns the dune element index of the root model index
    size_t map2dune(size_t rIdx) const {
        if (rIdx>root2dune.size()) {
            std::cout << "mapping problmes";
        }
        return root2dune.at(rIdx);
    }

};

} // end namespace GridGrowth

} // end namespace Dumux

#endif
