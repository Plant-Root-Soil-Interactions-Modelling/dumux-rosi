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
#ifndef DUMUX_RICHARDS_2C_BUFFER_INDICES_HH
#define DUMUX_RICHARDS_2C_BUFFER_INDICES_HH

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
struct RichardsTwoCBufferIndices
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    //////////
    // primary variable indices
    //////////

    //! Primary variable index for the wetting phase pressure
    static const int pwIdx = PVOffset + 0;
    //! Primary variable index for the wetting phase pressure head (used for pressure head formulation
    static const int hIdx = PVOffset + 0;
    static const int massOrMoleFracIdx = PVOffset + 1; //!< mole fraction of the second component

    //////////
    // equation indices
    //////////
    static const int contiEqIdx = PVOffset + 0;//!< continuity equation index
    static const int transportEqIdx = PVOffset + 1; //!< transport equation index

    //////////
    // phase indices
    //////////
    static const int wPhaseIdx = FluidSystem::wPhaseIdx; //!< Index of the wetting phase;
    static const int nPhaseIdx = FluidSystem::nPhaseIdx; //!< Index of the non-wetting phase;

    //! Component indices
    static const int phaseCompIdx = wPhaseIdx;//!< The index of the main component of the considered phase
    //! The index of the transported (minor) component; ASSUMES phase indices of 0 and 1
    static const int transportCompIdx = (unsigned int)(1-wPhaseIdx);

};
// \}

} // end namespace

#endif
