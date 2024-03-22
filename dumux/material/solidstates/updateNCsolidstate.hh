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
 * \ingroup SolidStates
 * \brief Update the solid volume fractions (inert and reacitve) and set them in the solidstate
 */
#ifndef DUMUX_UPDATE_NCSOLID_VOLUME_FRACTION_HH
#define DUMUX_UPDATE_NCSOLID_VOLUME_FRACTION_HH

namespace Dumux {
	

/*!
 * \ingroup SolidStates
 * \brief update the solid volume fractions (inert and reacitve) and set them in the solidstate
 * \note updates the inert components (TODO: these are assumed to come last right now in the solid system!)
 * \note gets the non-inert components from the primary variables
 */
template<class ElemSol, class Problem, class Element, class Scv, class SolidState>
void updateNCSolidstate(const ElemSol& elemSol,
                                const Problem& problem,
                                const Element& element,
                                const Scv& scv,
                                SolidState& solidState,
                                const int solidVolFracOffset)
{    
	//assert(useMoles, "updateNCSolidstate only works if useMoles == true");
	const auto& sp = problem.spatialParams();
    auto&& priVars = elemSol[scv.localDofIndex()];
	solidState.setVolumeFraction(solidState.mainCompIdx, 1 - sp.porosity(element));//main component, takes all the space
	double sumSecondaryFractions = 0.0;
	for (int sCompIdx = 0; sCompIdx < solidState.numComponents - 1; ++sCompIdx)
    {
        solidState.setVolumeFraction(sCompIdx, 0);//assume all the others do not take up space
		int compIdx = sCompIdx + solidVolFracOffset;
		solidState.setMoleFraction(sCompIdx, std::max(0.,priVars[compIdx]) ); //refuse values bellow 0.
        sumSecondaryFractions += std::max(0.,priVars[compIdx]) ; //refuse values bellow 0.
    }
	solidState.setMoleFraction(solidState.mainCompIdx, 1.0 - sumSecondaryFractions);
}

} // end namespace Dumux

#endif
