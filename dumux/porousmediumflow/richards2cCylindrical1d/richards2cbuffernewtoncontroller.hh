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
 * \brief A newton solver specific to the RichardsTwoC problem.
 */
#ifndef DUMUX_RICHARDS_2C_BUFFER_RADIALLY_SYMMETRIC_NEWTON_CONTROLLER_HH
#define DUMUX_RICHARDS_2C_BUFFER_RADIALLY_SYMMETRIC_NEWTON_CONTROLLER_HH

#include <dune/common/version.hh>

#include "richards2cbufferproperties.hh"

#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux {
/*!
 * \ingroup Newton
 * \brief A RichardsTwoC model specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * and can thus do update smarter than the plain Newton controller.
 */
template <class TypeTag>
class RichardsTwoCBufferRadiallySymmetricNewtonController : public NewtonController<TypeTag>
{
    typedef NewtonController<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { pwIdx = Indices::pressureIdx };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    enum { dim = GridView::dimension };
    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };
    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
    const Dune::ParameterTree &tree = ParameterTree::tree();

public:
    /*!
     * \brief Constructor
     */
    RichardsTwoCBufferRadiallySymmetricNewtonController(const Problem &problem)
        : ParentType(problem)
	{
	   avoidNegativeValue_ = tree.template get<bool>("Newton.AvoidNegativeValue", true);
	}
     /*!
     * \brief Returns true if another iteration should be done.
     *
     * \param uCurrentIter The solution of the current Newton iteration
     */
    bool newtonProceed(const SolutionVector &uCurrentIter)
    {
   		haveNegValue_ = false;
   		for (const auto sol : uCurrentIter)
           {
   		    if (sol[1] < 0.0)
   		        {
   		            haveNegValue_ = true;
                       std::cout << "\nNegative concentrations in solution "<< sol[1]<< " \n";
   		            if (ParentType::numSteps_ >= ParentType::maxSteps_)
   		                return false;
   		            else return true;
   		        }
               //break; // just consider the first element at root surface
           }

        if (ParentType::numSteps_ < 2)
            return true; // we always do at least two iterations
        else if (ParentType::asImp_().newtonConverged()) {
            return false; // we are below the desired tolerance
        }
        else if (ParentType::numSteps_ >= ParentType::maxSteps_) {
            // We have exceeded the allowed number of steps. If the
            // maximum relative shift was reduced by a factor of at least 4,
            // we proceed even if we are above the maximum number of steps.
            if (ParentType::enableShiftCriterion_)
                return ParentType::shift_*4.0 < ParentType::lastShift_;
            else
                return ParentType::reduction_*4.0 < ParentType::lastReduction_;
        }

        return true;
    }

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    bool newtonConverged() const
    {
        if (haveNegValue_)
            if ((avoidNegativeValue_) or (ParentType::numSteps_ <= ParentType::maxSteps_))
            {   std::cout << "\nNegative concentrations in solution, set Newton return to be not converged !!! \n";
                return false;
            }

        if (ParentType::enableShiftCriterion_ && !ParentType::enableResidualCriterion_)
        {
            return ParentType::shift_ <= ParentType::shiftTolerance_;
        }
        else if (!ParentType::enableShiftCriterion_ && ParentType::enableResidualCriterion_)
        {
            if(ParentType::enableAbsoluteResidualCriterion_)
                return ParentType::residual_ <= ParentType::residualTolerance_;
            else
                return ParentType::reduction_ <= ParentType::reductionTolerance_;
        }
        else if (ParentType::satisfyResidualAndShiftCriterion_)
        {
            if(ParentType::enableAbsoluteResidualCriterion_)
                return ParentType::shift_ <= ParentType::shiftTolerance_
                        && ParentType::residual_ <= ParentType::residualTolerance_;
            else
                return ParentType::shift_ <= ParentType::shiftTolerance_
                        && ParentType::reduction_ <= ParentType::reductionTolerance_;
        }
        else
        {
            return ParentType::shift_ <= ParentType::shiftTolerance_
                    || ParentType::reduction_ <= ParentType::reductionTolerance_
                    || ParentType::residual_ <= ParentType::residualTolerance_;
        }

        return false;
    }

    /*!
     * \brief Update the current solution of the newton method
     *
     * This is basically the step
     * \f[ u^{k+1} = u^k - \Delta u^k \f]
     *
     * \param uCurrentIter The solution after the current Newton iteration \f$ u^{k+1} \f$
     * \param uLastIter The solution after the last Newton iteration \f$ u^k \f$
     * \param deltaU The vector of differences between the last
     *               iterative solution and the next one \f$ \Delta u^k \f$
     */
    /* void newtonUpdate(SolutionVector &uCurrentIter,
                      const SolutionVector &uLastIter,
                      const SolutionVector &deltaU)
    {
        ParentType::newtonUpdate(uCurrentIter, uLastIter, deltaU);

        if (!GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, UseLineSearch))
        {
            // do not clamp anything after 5 iterations
            if (this->numSteps_ > 4)
                return;

            // clamp saturation change to at most 20% per iteration
            FVElementGeometry fvGeometry;
            const GridView &gridView = this->problem_().gridView();
            ElementIterator eIt = gridView.template begin<0>();
            const ElementIterator eEndIt = gridView.template end<0>();
            for (; eIt != eEndIt; ++eIt) {
                fvGeometry.update(gridView, *eIt);
                for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
                    int dofIdxGlobal = this->model_().dofMapper().subIndex(*eIt, scvIdx, dofCodim);
#else
                    int dofIdxGlobal = this->model_().dofMapper().map(*eIt, scvIdx, dofCodim);
#endif

                    // calculate the old wetting phase saturation
                    const SpatialParams &spatialParams = this->problem_().spatialParams();
                    const MaterialLawParams &mp = spatialParams.materialLawParams(*eIt, fvGeometry, scvIdx);
                    Scalar pcMin = MaterialLaw::pc(mp, 1.0);
                    Scalar pw = uLastIter[dofIdxGlobal][pwIdx];
                    Scalar pn = std::max(this->problem_().referencePressure(*eIt, fvGeometry, scvIdx),
                                         pw + pcMin);
                    Scalar pcOld = pn - pw;
                    Scalar SwOld = std::max<Scalar>(0.0, MaterialLaw::sw(mp, pcOld));

                    // convert into minimum and maximum wetting phase
                    // pressures
                    Scalar pwMin = pn - MaterialLaw::pc(mp, SwOld - 0.2);
                    Scalar pwMax = pn - MaterialLaw::pc(mp, SwOld + 0.2);

                    // clamp the result
                    pw = uCurrentIter[dofIdxGlobal][pwIdx];
                    pw = std::max(pwMin, std::min(pw, pwMax));
                    uCurrentIter[dofIdxGlobal][pwIdx] = pw;

                }
            }
        }
    }*/
private:
	bool haveNegValue_, avoidNegativeValue_;
};
}

#endif
