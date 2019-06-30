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
 * \ingroup ImplicitGridGrowthPropertyDefaults
 * \ingroup ImplicitGridGrowth
 * \file
 *
 * \brief Defines a type tag and some fundamental properties for
 *        grid growth
 */
#ifndef DUMUX_IMPLICIT_GRIDGROWTH_PROPERTY_DEFAULTS_HH
#define DUMUX_IMPLICIT_GRIDGROWTH_PROPERTY_DEFAULTS_HH

#include <dumux/common/basicproperties.hh>
#include "gridgrowthproperties.hh"
#include "gridgrowthindicatordefault.hh"
#include "gridgrowthhelperdefault.hh"

namespace Dumux
{
namespace Properties
{

// no growing grid
SET_BOOL_PROP(GridGrowth, GrowingGrid, false);

//standard setting
SET_INT_PROP(GridGrowth, GrowthInterval, 1);

//! Set the default indicator class models for growth
SET_TYPE_PROP(GridGrowth, GrowthIndicator, ImplicitGridGrowthIndicatorDefault<TypeTag>);

//! Set the default indicator class models for growth
SET_TYPE_PROP(GridGrowth, GrowthHelper, ImplicitGridGrowthHelperDefault<TypeTag>);
} // namespace Properties
} // namespace Dumux

#endif
