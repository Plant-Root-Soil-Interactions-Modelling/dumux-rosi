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
 * \ingroup RichardsModel
 * \brief A Richards model Newton solver.
 */

#ifndef DUMUX_RICHARDS_NEWTON_SOLVER_HH
#define DUMUX_RICHARDS_NEWTON_SOLVER_HH

#include <dumux/common/properties.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/discretization/elementsolution.hh>

namespace Dumux {
/*!
 * \ingroup RichardsModel
 * \brief A Richards model specific Newton solver.
 *
 * This solver 'knows' what a 'physically meaningful' solution is
 * and can thus do update smarter than the plain Newton solver.
 *
 * \todo make this typetag independent by extracting anything model specific from assembler
 *       or from possible ModelTraits.
 */
template <class Assembler, class LinearSolver,
          class Reassembler = PartialReassembler<Assembler>,
          class Comm = Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> >
class RichardsNewtonSolver : public NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>
{
    using Scalar = typename Assembler::Scalar;
    using ParentType = NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>;
    using SolutionVector = typename Assembler::ResidualType;
    //using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    //using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;

    using MaterialLaw = typename Assembler::Problem::SpatialParams::MaterialLaw;
    using Indices = typename Assembler::GridVariables::VolumeVariables::Indices;
    using VolumeVariables = typename Assembler::GridVariables::VolumeVariables;
    static constexpr int numFComponents =  VolumeVariables::FluidSystem::numComponents;
    static constexpr int numSComponents =  VolumeVariables::SolidSystem::numComponents - VolumeVariables::SolidSystem::numInertComponents; 
    
    enum { 
    numComponents = numFComponents + numSComponents,
        //numComponents = VolumeVariables::ModelTraits::numFluidComponents(),
        //numComponents = GetPropType<TypeTag, Properties::ModelTraits>::numComponents(),
        //nElem =  VolumeVariables::numFluidComponents() + VolumeVariables::numSolidComponents() -1 ,
        pressureIdx = Indices::pressureIdx };
    //enum { numComponents    = ModelTraits::numFluidComponents() };

public:
    using ParentType::ParentType;

private:

    /*!
     * \brief Update the current solution of the Newton method
     *
     * \param uCurrentIter The solution after the current Newton iteration \f$ u^{k+1} \f$
     * \param uLastIter The solution after the last Newton iteration \f$ u^k \f$
     * \param deltaU The vector of differences between the last
     *               iterative solution and the next one \f$ \Delta u^k \f$
     */
    void choppedUpdate_(SolutionVector &uCurrentIter,
                        const SolutionVector &uLastIter,
                        const SolutionVector &deltaU) final
    {
        uCurrentIter = uLastIter;
        uCurrentIter -= deltaU;

        // do not clamp anything after 5 iterations
        
            // clamp saturation change to at most 20% per iteration
            const auto& fvGridGeometry = this->assembler().fvGridGeometry();
            for (const auto& element : elements(fvGridGeometry.gridView()))
            {
                auto fvGeometry = localView(fvGridGeometry);
                fvGeometry.bindElement(element);

                for (auto&& scv : scvs(fvGeometry))
                {
                    auto dofIdxGlobal = scv.dofIndex();

                    if (this->numSteps_ <= 4)
                    {
                    
                        // calculate the old wetting phase saturation
                        const auto& spatialParams = this->assembler().problem().spatialParams();
                        const auto elemSol = elementSolution(element, uCurrentIter, fvGridGeometry);
                        const auto& materialLawParams = spatialParams.materialLawParams(element, scv, elemSol);
                        const Scalar pcMin = MaterialLaw::pc(materialLawParams, 1.0);
                        const Scalar pw = uLastIter[dofIdxGlobal][pressureIdx];
                        const Scalar pn = std::max(this->assembler().problem().nonWettingReferencePressure(), pw + pcMin);
                        const Scalar pcOld = pn - pw;
                        const Scalar SwOld = std::max(0.0, MaterialLaw::sw(materialLawParams, pcOld));

                        // convert into minimum and maximum wetting phase
                        // pressures
                        const Scalar pwMin = pn - MaterialLaw::pc(materialLawParams, SwOld - 0.2);
                        const Scalar pwMax = pn - MaterialLaw::pc(materialLawParams, SwOld + 0.2);

                        // clamp the result
                        //using std::min; using std::max;
                        uCurrentIter[dofIdxGlobal][pressureIdx] = std::max(pwMin, std::min(uCurrentIter[dofIdxGlobal][pressureIdx], pwMax));
                    }//
                    
                    for(int eqIdx = pressureIdx +1; eqIdx < numComponents ;  eqIdx++)
                    {
                        // std::cout<<"eqIdx "<<eqIdx<<" numComponents "<<numComponents;
                        // std::cout<<"uCurrentIter[dofIdxGlobal]["<<eqIdx<<"] "<<uCurrentIter[dofIdxGlobal][eqIdx]<<", ";
                        uCurrentIter[dofIdxGlobal][eqIdx] = std::max(0., uCurrentIter[dofIdxGlobal][eqIdx]);//all other elemnents >=0.
                    }
                }
            }
        

        if (this->enableResidualCriterion())
            this->computeResidualReduction_(uCurrentIter);

        else
        {
            // If we get here, the convergence criterion does not require
            // additional residual evalutions. Thus, the grid variables have
            // not yet been updated to the new uCurrentIter.
            this->assembler().updateGridVariables(uCurrentIter);
        }
    }
};

} // end namespace Dumux

#endif
