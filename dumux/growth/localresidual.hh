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
 * \ingroup GrowthModule
 * \brief Local residual considering volume change
 */
#ifndef DUMUX_GROWTH_LOCAL_RESIDUAL_ADAPTER_HH
#define DUMUX_GROWTH_LOCAL_RESIDUAL_ADAPTER_HH

namespace Dumux {
namespace GrowthModule {

/*!
 * \ingroup GrowthModule
 * \brief Local residual considering volume change
 */
template<class Base>
class GrowthLocalResidualAdapter : public Base
{
public:
    using Base::Base;
    using ElementResidualVector = typename Base::ElementResidualVector;

    /*!
     * \brief Compute the storage local residual, i.e. the deviation of the
     *        storage term from zero for instationary problems.
     *
     * \param residual The residual vector to fill
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevElemVolVars The volume averaged variables for all
     *                        sub-control volumes of the element at the previous time level
     * \param curElemVolVars The volume averaged variables for all
     *                       sub-control volumes of the element at the current  time level
     * \param scv The sub control volume the storage term is integrated over
     */
    using Base::evalStorage;
    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables>
    void evalStorage(ElementResidualVector& residual,
                     const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& prevElemVolVars,
                     const ElementVolumeVariables& curElemVolVars,
                     const typename FVElementGeometry::SubControlVolume& scv) const
    {
        const auto& curVolVars = curElemVolVars[scv];
        const auto& prevVolVars = prevElemVolVars[scv];

        // mass balance within the element. this is the
        // \f$\frac{m}{\partial t}\f$ term if using implicit or explicit
        // euler as time discretization.
        //
        // TODO: We might need a more explicit way for
        // doing the time discretization...

        //! Compute storage with the model specific storage residual
        auto prevStorage = this->computeStorage(problem, scv, prevVolVars);
        auto storage = this->computeStorage(problem, scv, curVolVars);

        prevStorage *= problem.preGrowthVolume(element, scv, prevVolVars);
        storage *= problem.postGrowthVolume(element, scv, prevVolVars);

        storage -= prevStorage;
        storage /= this->timeLoop().timeStepSize();

        residual[scv.indexInElement()] += storage;
    }
};

} // end namespace GrowthModule
} // end namespace Dumux

#endif
