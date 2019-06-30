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
 * \ingroup ImplicitGridAdaptProperties
 * \ingroup ImplicitGridAdapt
 * \file
 *
 * \brief Defines a type tag and some fundamental properties for
 *        linear solvers
 */
#ifndef DUMUX_IMPLICIT_GRIDGROWTH_PROPERTIES_HH
#define DUMUX_IMPLICIT_GRIDGROWTH_PROPERTIES_HH

#include <dumux/common/basicproperties.hh>

namespace Dumux
{
namespace Properties
{
//! Grid growth type tag
NEW_TYPE_TAG(GridGrowth);

//! Defines if the grid has growth ability
NEW_PROP_TAG(GrowingGrid);

//! Class defining the indicator
NEW_PROP_TAG(GrowthIndicator);

//! Helper handling user data transfer
NEW_PROP_TAG(GrowthHelper);

//! Time step interval for growth
NEW_PROP_TAG(GrowthInterval);
} // namespace Properties
} // namespace Dumux

#endif
