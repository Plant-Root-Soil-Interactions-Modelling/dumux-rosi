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
 * \ingroup RichardsTwoCModel
 * \file
 *
 * \brief Contains the property declarations for the RichardsTwoC
 *        model.
 */
#ifndef DUMUX_RICHARDS_2C_BUFFER_RADIALLY_SYMMETRIC_PROPERTIES_HH
#define DUMUX_RICHARDS_2C_BUFFER_RADIALLY_SYMMETRIC_PROPERTIES_HH

#include <dumux/implicit/box/properties.hh>
#include <dumux/implicit/cellcentered/properties.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/properties.hh>

namespace Dumux
{
// \{
///////////////////////////////////////////////////////////////////////////
// properties for the isothermal richards model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit isothermal one-phase two-component problems
NEW_TYPE_TAG(RichardsTwoCBufferRadiallySymmetric);
NEW_TYPE_TAG(BoxRichardsTwoCBufferRadiallySymmetric, INHERITS_FROM(BoxModel, RichardsTwoCBufferRadiallySymmetric));
NEW_TYPE_TAG(CCRichardsTwoCBufferRadiallySymmetric, INHERITS_FROM(CCModel, RichardsTwoCBufferRadiallySymmetric));

//! The type tags for the corresponding non-isothermal problems
NEW_TYPE_TAG(RichardsTwoCBufferRadiallySymmetricNI, INHERITS_FROM(RichardsTwoCBufferRadiallySymmetric, NonIsothermal));
NEW_TYPE_TAG(BoxRichardsTwoCBufferRadiallySymmetricNI, INHERITS_FROM(BoxModel, RichardsTwoCBufferRadiallySymmetricNI));
NEW_TYPE_TAG(CCRichardsTwoCBufferRadiallySymmetricNI, INHERITS_FROM(CCModel, RichardsTwoCBufferRadiallySymmetricNI));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(PhaseIdx); //!< A phase index in to allow that a two-phase fluidsystem is used
NEW_PROP_TAG(NumComponents);   //!< Number of fluid components in the system
NEW_PROP_TAG(Indices); //!< Enumerations for the model
NEW_PROP_TAG(SpatialParams); //!< The type of the spatial parameters
NEW_PROP_TAG(EffectiveDiffusivityModel); //!< The employed model for the computation of the effective diffusivity
NEW_PROP_TAG(FluidSystem); //!< Type of the multi-component relations
NEW_PROP_TAG(FluidState); //!< Type of the fluid state to be used
NEW_PROP_TAG(ImplicitMassUpwindWeight);   //!< The default value of the upwind weight
NEW_PROP_TAG(ImplicitMobilityUpwindWeight); //!< Weight for the upwind mobility in the velocity calculation
NEW_PROP_TAG(ProblemEnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(UseMoles); //!< Defines whether mole (true) or mass (false) fractions are used
NEW_PROP_TAG(Scaling); //!< Defines Scaling of the model
NEW_PROP_TAG(SpatialParamsForchCoeff); //!< Property for the forchheimer coefficient
NEW_PROP_TAG(VtkAddVelocity); //!< Returns whether velocity vectors are written into the vtk output

NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (by default extracted from the spatial parameters)
NEW_PROP_TAG(MaterialLawParams); //!< The type of the parameter object for the material law (by default extracted from the spatial parameters)
NEW_PROP_TAG(UsePH); //!< Defines whether pressure [Pa] (false) or pressure head [cm] (ture) is used

// \}
}

} // end namespace

#endif
