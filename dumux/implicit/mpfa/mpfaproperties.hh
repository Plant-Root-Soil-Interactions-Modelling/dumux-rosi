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
#ifndef DUMUX_MPFA_PROPERTIES_HH
#define DUMUX_MPFA_PROPERTIES_HH

#include <dumux/implicit/properties.hh>

/*!
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \ingroup MpfaModel
 */
namespace Dumux
{
////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the cell-centered mpfa scheme
NEW_TYPE_TAG(MpfaModel, INHERITS_FROM(ImplicitBase));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(MpfaMethod);                           //!< specifies the mpfa method to be used
NEW_PROP_TAG(MpfaMethods);                          //!< indices of the different methods
NEW_PROP_TAG(MpfaInteractionVolume);                //!< the interaction volume class
NEW_PROP_TAG(MpfaBoundaryInteractionVolume);        //!< the interaction volume class of the O method (for boundary treatment within the other methods)
NEW_PROP_TAG(MpfaInteractionVolumeManager);         //!< class to fill interactionvolumes with data
NEW_PROP_TAG(MpfaBoundaryInteractionVolumeManager); //!< class to fill interactionvolumes of the O method (for boundary treatment within the other methods)
NEW_PROP_TAG(MpfaInteractionVolumeContainer);       //!< class to store the interactionvolumes data
NEW_PROP_TAG(MpfaFluxCalculator);                   //!< class to calculate fluxes over sub control volume faces
NEW_PROP_TAG(MpfaBoundaryFluxCalculator);           //!< class to calculate fluxes over sub control volume faces of the O method (for boundary treatment within the other methods)
NEW_PROP_TAG(MpfaContinuityPoint);                  //!< sets the position on which continuity is forced
NEW_PROP_TAG(MpfaFaceTypes);                        //!< a set of indices determining the type of a sub control volume face
NEW_PROP_TAG(FullMpfaOBoundary); 	                //!< determines if the full face is treated with the o method (e.g. for the l method)
NEW_PROP_TAG(BoundaryLayers);                       //!< Property to find interaction volumes on different layers of the boundary
NEW_PROP_TAG(ImplicitIsMpfa);                       //!< Property to determine whether Mpfa is used or not
NEW_PROP_TAG(FacetCoupling);                        //!< specifies whether there is a coupled model on the element facets or not
NEW_PROP_TAG(CouplingStrategies);                   //!< set of indices indicating which coupling strategy is used on the facets
NEW_PROP_TAG(CouplingStrategy);                     //!< defines the actual strategy used in the coupled model
NEW_PROP_TAG(XiFactor);                             //!< coefficient for the coupling conditions in the coupled models
}
}

#include "mpfapropertydefaults.hh"

#endif
