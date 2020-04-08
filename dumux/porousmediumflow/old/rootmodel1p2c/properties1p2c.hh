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
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \ingroup RootsystemModel
 * \file
 *
 * \brief Defines the properties required for the one-phase fully implicit model.
 */
#ifndef DUMUX_ROOTSYSTEM_PROPERTIES_1P2C_HH
#define DUMUX_ROOTSYSTEM_PROPERTIES_1P2C_HH

//#include <dumux/implicit/box/properties.hh>
#include <dumux/porousmediumflow/1p2c/implicit/properties.hh>
//#include <dumux/implicit/cellcentered/properties.hh>
#include <dumux/implicit/growth/gridgrowthproperties.hh>

namespace Dumux
{
// \{
///////////////////////////////////////////////////////////////////////////
// properties for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit single-phase problems
NEW_TYPE_TAG(RootsystemOnePTwoC, INHERITS_FROM(GridGrowth));
NEW_TYPE_TAG(BoxRootsystemOnePTwoC, INHERITS_FROM(BoxModel, RootsystemOnePTwoC));
NEW_TYPE_TAG(CCRootsystemOnePTwoC, INHERITS_FROM(CCModel, RootsystemOnePTwoC));

}

} // end namespace

#endif
