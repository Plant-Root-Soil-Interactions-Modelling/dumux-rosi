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
 * \brief Index names for the RichardsTwoC model.
 */
#ifndef DUMUX_RICHARDS_2C_BUFFER_RADIALLY_SYMMETRIC_INDICES_HH
#define DUMUX_RICHARDS_2C_BUFFER_RADIALLY_SYMMETRIC_INDICES_HH

#include "richards2cbufferproperties.hh"

namespace Dumux
{
// \{

/*!
 * \ingroup RichardsTwoCModel
 * \ingroup ImplicitIndices
 * \brief Index names for the RichardsTwoC model.
 */

template <class TypeTag, int PVOffset = 0>
struct RichardsTwoCBufferRadiallySymmetricIndices
{

    //! Set the default phase used by the fluid system to the first one
    static const int phaseIdx = 0;

    //! Component indices
    static const int phaseCompIdx = phaseIdx;//!< The index of the main component of the considered phase
    //! The index of the transported (minor) component; ASSUMES phase indices of 0 and 1
    static const int transportCompIdx = (unsigned int)(1-phaseIdx);

    // Equation indices
   static const int conti0EqIdx = PVOffset + 0; //!< continuity equation index
   static const int transportEqIdx = PVOffset + 1; //!< transport equation index

    // primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< pressure
    static const int massOrMoleFracIdx = PVOffset + 1; //!< mole fraction of the second component
};

// \}
}

#endif

