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
 * \brief A newton solver specific to the Rootsystem problem.
 */
#ifndef DUMUX_ROOTSYSTEM_NEWTON_CONTROLLER_HH
#define DUMUX_ROOTSYSTEM_NEWTON_CONTROLLER_HH

#include "properties.hh"
#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux {
/*!
 * \ingroup Newton
 * \brief A Rootsystem model specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * and can thus do update smarter than the plain Newton controller.
 */
template <class TypeTag>
class RootsystemNewtonController : public NewtonController<TypeTag>
{
    typedef NewtonController<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

public:
    /*!
     * \brief Constructor
     */
    RootsystemNewtonController(const Problem &problem)
        : ParentType(problem)
    {}

    /* TODO: remove or overwrite if  water stress is not used */

    /*
     * \brief Indicates that we're done solving the non-linear system
     *        of equations.*/
    void newtonEnd()
    {
        this->problem_().setWaterStress();
        bool switchBC = this->problem_().getSwitchBC();

        if (switchBC){
            DUNE_THROW(NumericalProblem,
                           "switch collar BC");
        }
        else{
            if (GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, WriteConvergence))
                this->convergenceWriter_.endTimestep();
        }
    }

    /*
     * \brief Called if the Newton method broke down.
     *
     * This method is called _after_ newtonEnd()*/
    void newtonFail()
    {
        bool switchBC = this->problem_().getSwitchBC();

        if (switchBC){
            this->model_().jacobianAssembler().reassembleAll();
            this->numSteps_ =  this->targetSteps_;
        }
        else {
            this->model_().jacobianAssembler().reassembleAll();
             this->numSteps_ =  this->targetSteps_*2;
        }
        this->problem_().resetSwitchBC();
    }

}; //class
} //namespace

#endif
