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
#ifndef DUMUX_GRID_GROWTH_HH
#define DUMUX_GRID_GROWTH_HH

#include <memory>
#include <dune/common/timer.hh>
#include <dumux/common/properties.hh>
#include <dumux/growth/algorithm.hh>

namespace Dumux {

namespace GrowthModule {

class GridGrowth
{
public:
    //! constructor with a grid and a growth algorithm
    GridGrowth(std::shared_ptr<GrowthAlgorithm> growthAlgorithm)
    : growthAlgorithm_(growthAlgorithm)
    {}

    //! do a single time step
    void grow(double dt)
    {
        // reset internal status
        hasGrown_ = false;

        // grow using the chosen algorithm
        Dune::Timer timer;
        std::cout << std::endl << "Beginning growth step with time step size: " << dt << " seconds." << std::endl;
        const auto status = growthAlgorithm_->grow(dt);
        std::cout << "Growth step created " << status.inserted << " new elements"
                  << ", removed " << status.removed << " elements"
                  << ", and moved " << status.moved << " vertices in "
                  << timer.stop() << " seconds." << std::endl << std::endl;

        // set the internal status
        if (status.inserted != 0 || status.removed != 0 || status.moved != 0)
            hasGrown_ = true;
    }

    //! query if the grid has grown in the last growth step
    bool hasGrown() const
    { return hasGrown_; }

private:
    //! the growth algorithm implementation
    std::shared_ptr<GrowthAlgorithm> growthAlgorithm_;

    //! if the grid has grown in the last growth step
    bool hasGrown_;
};

} // end namespace GridGrowth

} // end namespace Dumux

#endif
