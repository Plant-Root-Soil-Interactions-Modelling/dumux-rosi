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
#ifndef DUMUX_ROOTSYSTEM_PROPERTIES_DATA_HH
#define DUMUX_ROOTSYSTEM_PROPERTIES_DATA_HH

#include <dumux/implicit/box/properties.hh>
#include <dumux/implicit/cellcentered/properties.hh>
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
NEW_TYPE_TAG(Rootsystem, INHERITS_FROM(GridGrowth));
NEW_TYPE_TAG(BoxRootsystem, INHERITS_FROM(BoxModel, Rootsystem));
NEW_TYPE_TAG(CCRootsystem, INHERITS_FROM(CCModel, Rootsystem));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(Indices); //!< Enumerations for the model
NEW_PROP_TAG(SpatialParams); //!< The type of the spatial parameters object
NEW_PROP_TAG(FluidSystem); //!< The type of the fluid system to use
NEW_PROP_TAG(Fluid); //!< The fluid used for the default fluid system
NEW_PROP_TAG(ProblemEnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(ImplicitMassUpwindWeight); //!< Returns weight of the upwind cell when calculating fluxes
NEW_PROP_TAG(ImplicitMobilityUpwindWeight); //!< Weight for the upwind mobility in the velocity calculation
NEW_PROP_TAG(SpatialParamsForchCoeff); //!< Property for the forchheimer coefficient
// \}
}

} // end namespace

#endif
