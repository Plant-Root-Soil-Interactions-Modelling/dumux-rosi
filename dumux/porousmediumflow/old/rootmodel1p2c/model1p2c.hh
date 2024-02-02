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
 *
 * \brief Base class for all models which use the one-phase,
 *        fully implicit model.
 *        Adaption of the fully implicit scheme to the one-phase flow model.
 */

#ifndef DUMUX_ROOTSYSTEM_MODEL_1P2C_HH
#define DUMUX_ROOTSYSTEM_MODEL_1P2C_HH

#include <dumux/porousmediumflow/implicit/velocityoutput.hh>
#include <dumux/implicit/growth/gridgrowthproperties.hh>
//#include "properties.hh"
#include "properties1p2c.hh"

namespace Dumux
{
/*!
 * \ingroup RootsystemBoxModel
 * \brief A single-phase, isothermal flow model using the fully implicit scheme.
 *
 * Single-phase, isothermal flow model, which uses a standard Darcy approach as the
 * equation for the conservation of momentum:
 * \f[
 v = - \frac{\textbf K}{\mu}
 \left(\textbf{grad}\, p - \varrho {\textbf g} \right)
 * \f]
 *
 * and solves the mass continuity equation:
 * \f[
 \phi \frac{\partial \varrho}{\partial t} + \text{div} \left\lbrace - \varrho \frac{\textbf K}{\mu} \left( \textbf{grad}\, p -\varrho {\textbf g} \right) \right\rbrace = q,
 * \f]
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 * The model supports compressible as well as incompressible fluids.
 */
template<class TypeTag >
class RootsystemOnePTwoCModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseModel) ParentType;
    //define FluidSystem
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
    //added
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { phaseIdx = Indices::phaseIdx };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using ScalarField = Dune::BlockVector<Dune::FieldVector<Scalar, 1> >;
    using VectorField = Dune::BlockVector<GlobalPosition>;

public:

    /*!
     * \brief Called by the update() method before it tries to
     *        apply the newton method. This is primarily a hook
     *        which the actual model can overload.
     */
    void updateBegin()
    {
        // call of parenttyp fuction
        ParentType::updateBegin();

        if(GET_PROP_VALUE(TypeTag, GrowingGrid) && this->problem_().gridGrowth().wasGrown())
        {
            this->prevSol() =  this->curSol();
            this->updateBoundaryIndices_();
            this->jacobianAssembler().init(this->problem_());

       }
    }

    /*!
     * \brief \copybrief Dumux::ImplicitModel::addOutputVtkFields
     *
     * Specialization for the RootsystemModel, adding the pressure and
     * the process rank to the VTK writer.
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<double, dimWorld> > VectorField;

        // create the required scalar fields
        unsigned numDofs = this->numDofs();
        auto numElements = this->gridView_().size(0);

        ScalarField *p = writer.allocateManagedBuffer(numDofs);
        ScalarField *pint = writer.allocateManagedBuffer(numDofs);

        ScalarField *Kx = writer.allocateManagedBuffer(numDofs);
        ScalarField *Kr = writer.allocateManagedBuffer(numDofs);
        ScalarField *rootOrder = writer.allocateManagedBuffer(numDofs);
        ScalarField *rootBranch = writer.allocateManagedBuffer(numDofs);
        ScalarField *rootMass = writer.allocateManagedBuffer(numDofs);
        ScalarField *rootRadius = writer.allocateManagedBuffer(numDofs);
        ScalarField *rootSurface = writer.allocateManagedBuffer(numDofs);
        ScalarField *rootLength = writer.allocateManagedBuffer(numDofs);
        // Define fraction of components
        ScalarField *moleFraction0 = writer.allocateManagedBuffer(numDofs);
        ScalarField *moleFraction1 = writer.allocateManagedBuffer(numDofs);
        ScalarField *massFraction0 = writer.allocateManagedBuffer(numDofs);
        ScalarField *massFraction1 = writer.allocateManagedBuffer(numDofs);

        //ScalarField *ph = writer.allocateManagedBuffer(numDofs);
        ScalarField *volumeFlux = writer.allocateManagedBuffer(numElements);
        VectorField *velocity = writer.template allocateManagedBuffer<Scalar, dimWorld>(numElements);
        ImplicitVelocityOutput<TypeTag> velocityOutput(this->problem_());

        if (velocityOutput.enableOutput())
        {
            // initialize velocity field
            for (unsigned int i = 0; i < numDofs; ++i)
            {
                (*velocity)[i] = Scalar(0.0);
            }
        }

        //signed numElements = this->gridView_().size(0);
        ScalarField *rank = writer.allocateManagedBuffer(numElements);

        for (const auto& element : elements(this->gridView_()))
        {
            int eIdx = this->problem_().model().elementMapper().index(element);
            (*rank)[eIdx] = this->gridView_().comm().rank();

            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView_(), element);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(this->problem_(),
                               element,
                               fvGeometry,
                               false /* oldSol? */);

            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
            {
                //int globalIdx = this->dofMapper().subIndex(*eIt, scvIdx, dofCodim);
                int globalIdx = eIdx;

                const SpatialParams &spatialParams = this->problem_().spatialParams();

                (*p)[globalIdx] = elemVolVars[scvIdx].pressure();

                //(*ph)[globalIdx] = elemVolVars[scvIdx].pressureHead(0);

                (*rootRadius)[globalIdx] = spatialParams.rootRadius(element,
                                                                    fvGeometry,
                                                                    scvIdx);
                (*rootSurface)[globalIdx] = spatialParams.rootSurface(element,
                                                                      fvGeometry,
                                                                      scvIdx);
                (*rootBranch)[globalIdx] = spatialParams.rootBranch(element,
                                                                      fvGeometry,
                                                                      scvIdx);
                (*rootOrder)[globalIdx] = spatialParams.rootOrder(element,
                                                                  fvGeometry,
                                                                  scvIdx);
                (*rootMass)[globalIdx] = spatialParams.rootMass(element,
                                                                fvGeometry,
                                                                scvIdx);
                (*Kx)[globalIdx] = spatialParams.Kx(element,
                                                    fvGeometry,
                                                    scvIdx);
                (*Kr)[globalIdx] = spatialParams.Kr(element,
                                                    fvGeometry,
                                                    scvIdx);

                (*rootLength)[globalIdx] = element.geometry().volume();
                // fraction of components
                (*moleFraction0)[globalIdx] = elemVolVars[scvIdx].moleFraction(0);
                (*moleFraction1)[globalIdx] = elemVolVars[scvIdx].moleFraction(1);
                (*massFraction0)[globalIdx] = elemVolVars[scvIdx].massFraction(0);
                (*massFraction1)[globalIdx] = elemVolVars[scvIdx].massFraction(1);
            }


            // velocity output
            if(velocityOutput.enableOutput())
            {
                calculateVolumeFluxVector(elemVolVars, fvGeometry, element, velocityOutput, *velocity);
                (*volumeFlux)[eIdx] = (*velocity)[eIdx].two_norm();
                (*velocity)[eIdx] /= M_PI*(*rootRadius)[eIdx]*(*rootRadius)[eIdx];
            }
        }

        writer.attachDofData(*p, "p", isBox);
        writer.attachDofData(*pint, "pint", isBox);
        //writer.attachDofData(*ph, "pressure head", isBox);
        writer.attachDofData(*rootRadius, "rootRadius", isBox);
        writer.attachDofData(*rootSurface, "rootSurface", isBox);
        writer.attachDofData(*rootLength, "rootLength", isBox);
        writer.attachDofData(*rootOrder, "rootOrder", isBox);
        writer.attachDofData(*rootBranch, "rootBranch", isBox);
        writer.attachDofData(*rootMass, "rootMass", isBox);

        writer.attachDofData(*moleFraction0, "x_" + FluidSystem::componentName(0), isBox);
        writer.attachDofData(*moleFraction1, "x_" + FluidSystem::componentName(1), isBox);

        writer.attachDofData(*massFraction0, "w_" + FluidSystem::componentName(0), isBox);
        writer.attachDofData(*massFraction1, "w_" + FluidSystem::componentName(1), isBox);

        if (velocityOutput.enableOutput())
        {
        //    writer.attachCellData(*velocity,  "velocity (m/s)", dimWorld);
            writer.attachCellData(*volumeFlux,  "volumeFlux (m^3/s)");
        }
        writer.attachCellData(*rank, "process rank");
    }

    // helper function for one-dimensional network grids
    // we always want element velocities. For cell-centered models the
    // velocity outputs computes actually the velocity * cross-section area
    // i.e. a "volume flux vector". In order to get the vectorial velocity
    // we need to divide by the area. The volume flux is the vector's magnitude.
    void calculateVolumeFluxVector(const ElementVolumeVariables& elemVolVars,
                                   const FVElementGeometry& fvGeometry,
                                   const Element& element,
                                   ImplicitVelocityOutput<TypeTag>& velocityOutput,
                                   VectorField& volumeFluxVector)
    {
        auto eIdx = this->problem_().elementMapper().index(element);
        volumeFluxVector[eIdx] = GlobalPosition(0.0);
        if(isBox)
        {
            FluxVariables fluxVars;
            fluxVars.update(this->problem_(),
                            element,
                            fvGeometry,
                            0,
                            elemVolVars);
            auto flux = fluxVars.volumeFlux(0);
            volumeFluxVector[eIdx] = (element.geometry().corner(1) - element.geometry().corner(0)).two_norm();
            volumeFluxVector[eIdx] *= flux;
        }
        else // cell-centered models use the standard output
            velocityOutput.calculateVelocity(volumeFluxVector, elemVolVars, fvGeometry, element, /*phaseIdx=*/0);
    }

    void sourceValues(SolutionVector& sourcevalues)
    {
        int gridsize = this->gridView_().size(0);
        sourcevalues.resize(gridsize);

        ElementIterator eIt = this->gridView_().template begin<0>();
        ElementIterator eEndIt = this->gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt)
        {
            int eIdx = this->problem_().model().elementMapper().index(*eIt);

            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView_(), *eIt);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(this->problem_(),
                               *eIt,
                               fvGeometry,
                               false /* oldSol? */);
            PrimaryVariables values;

            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
            {
                    this->problem_().solDependentSource(values,
                                                        *eIt,
                                                        fvGeometry,
                                                        scvIdx,
                                                        elemVolVars);
            }
            Scalar vol = eIt->geometry().volume();
            sourcevalues[eIdx][0] = values*vol;
        }
    }

};
}

#include "propertydefaults1p2c.hh"

#endif
