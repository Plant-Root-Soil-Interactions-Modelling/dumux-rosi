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
#ifndef DUMUX_GRID_GROWTH_ALGORITHM_HH
#define DUMUX_GRID_GROWTH_ALGORITHM_HH

#include <dumux/common/properties.hh>

namespace Dumux {

namespace GrowthModule {

class GrowthAlgorithm
{
public:
    //! a data structure for a status report from the growth step
    struct StatusReport
    {
        unsigned int inserted;
        unsigned int removed;
        unsigned int moved;
    };

    //! pure abstract base class has virtual destructor
    virtual ~GrowthAlgorithm() = default;

    //! do a growth step with step size dt
    //! returns how many element were inserted, how many removed, and how many vertices moved.
    virtual StatusReport grow(double dt) = 0;
};

class NoGrowthAlgorithm : public GrowthAlgorithm
{
public:
    //! do a growth step with step size dt
    //! returns how many element were inserted, how many removed, and how many vertices moved.
    StatusReport grow(double dt) final
    { return StatusReport{0,0,0}; }
};

} // end namespace GridGrowth

} // end namespace Dumux

#endif
