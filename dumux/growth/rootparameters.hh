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
 * \brief This file contains parameter classes needed for growth
 */
#ifndef DUMUX_ROOT_PARAMETERS_HH
#define DUMUX_ROOT_PARAMETERS_HH

namespace Dumux {

namespace GrowthModule {

//! Parameters given for every root
template <class Scalar>
struct RootParams
{
    Scalar radius; //! root base radius
    int order; //! the order if the root
    int rootId; //! which branch this segment is belonging to
    int plantId; //! which plant this segment is belonging to
    bool isInitialized = false; //! if at least some of the parameters have been set
};

//! Indices to access the parameters in the dgf file
enum DGFParamIndices {
    orderIdx = 0,
    rootIdIdx = 1,
    surfaceIdx = 2,
    massIdx = 3,
    plantIdx = 5
};

} // end namespace GridGrowth

} // end namespace Dumux

#endif
