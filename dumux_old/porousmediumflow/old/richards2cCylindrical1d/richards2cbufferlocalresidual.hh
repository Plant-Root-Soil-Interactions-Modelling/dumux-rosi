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
 * \brief Element-wise calculation of the residual for the RichardsTwoC fully implicit model.
 */
#ifndef DUMUX_RICHARDS_2C_BUFFER_RADIALLY_SYMMETRIC_LOCAL_RESIDUAL_HH
#define DUMUX_RICHARDS_2C_BUFFER_RADIALLY_SYMMETRIC_LOCAL_RESIDUAL_HH

#include "richards2cbufferproperties.hh"

namespace Dumux
{
/*!
 * \ingroup RichardsTwoCModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the residual for the RichardsTwoC fully implicit model.
 */
template<class TypeTag>
class RichardsTwoCBufferRadiallySymmetricLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        phaseIdx = Indices::phaseIdx,
        transportCompIdx = Indices::transportCompIdx
    };
    // indices of the equations
    enum {
        contiEqIdx = Indices::conti0EqIdx,
        transportEqIdx = Indices::transportEqIdx,
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum {dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    static const bool usePH = GET_PROP_VALUE(TypeTag, UsePH);
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    RichardsTwoCBufferRadiallySymmetricLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
    };

    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub control
     *        volume of a finite volume element for the RichardsTwoC
     *        model.
     *
     * This function should not include the source and sink terms.
     *
     * \param storage Stores the average mass per unit volume for each phase \f$\mathrm{[kg/m^3]}\f$
     * \param scvIdx The sub control volume index of the current element
     * \param usePrevSol Calculate the storage term of the previous solution
     *                   instead of the model's current solution.
     */
    void computeStorage(PrimaryVariables &storage, const int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const VolumeVariables &volVars =
            usePrevSol ?
            this->prevVolVars_(scvIdx) :
            this->curVolVars_(scvIdx);

        storage = 0;

        // partial time derivative of the wetting phase mass
        // pressure head formulation
        storage[contiEqIdx] =
                volVars.saturation(phaseIdx)
                *volVars.porosity()*volVars.coordinatesCenter()[0];;

        if(!useMoles) //mass-fraction formulation
        {
            // storage term of continuity equation - massfractions
            if (!usePH)
                storage[contiEqIdx] *= volVars.density();

            //storage term of the transport equation - massfractions
            storage[transportEqIdx] +=
                volVars.density() * volVars.massFraction(transportCompIdx) *
                (volVars.saturation(phaseIdx)*volVars.porosity()+volVars.buffer())
                *volVars.coordinatesCenter()[0];
        }
        else //mole-fraction formulation
        {
            // storage term of continuity equation- molefractions
            //careful: molarDensity changes with moleFrac!
            if (!usePH)
                storage[contiEqIdx] *= volVars.molarDensity();

            // storage term of the transport equation - molefractions
            storage[transportEqIdx] +=
                volVars.molarDensity()*volVars.moleFraction(transportCompIdx) *
                (volVars.saturation(phaseIdx)*volVars.porosity()+volVars.buffer())
				* volVars.coordinatesCenter()[0];
        }

    }


    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     *
     * \param flux Stores the total mass fluxes over a sub-control volume face
     *             of the current element \f$\mathrm{[kg/s]}\f$
     * \param fIdx The sub control volume face index inside the current
     *                element
     * \param onBoundary A boolean variable to specify whether the flux variables
     *        are calculated for interior SCV faces or boundary faces, default=false
     */
    void computeFlux(PrimaryVariables &flux, const int fIdx, bool onBoundary=false) const
    {
        FluxVariables fluxVars;
        fluxVars.update(this->problem_(),
                               this->element_(),
                               this->fvGeometry_(),
                               fIdx,
                               this->curVolVars_(),
                               onBoundary);

        flux = 0;
        asImp_()->computeAdvectiveFlux(flux, fluxVars);
        asImp_()->computeDiffusiveFlux(flux, fluxVars);
    }

    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a sub-control volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV
     */

    void computeAdvectiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        // data attached to upstream and the downstream vertices
        // of the current phase
        const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx(phaseIdx));
        const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx(phaseIdx));

        //pressure head formulation
        flux[contiEqIdx] =
            fluxVars.volumeFlux(phaseIdx);

        if(!useMoles) //mass-fraction formulation
        {
            // total mass flux - massfraction
            //KmvpNormal is the Darcy velocity multiplied with the normal vector, calculated in 1p2cfluxvariables.hh
            if (!usePH)
                flux[contiEqIdx] *=fluxVars.coordinatesCenter()[0]*
                    ((     massUpwindWeight_)*up.density()
                     +
                     ((1 - massUpwindWeight_)*dn.density()));

            // advective flux of the second component - massfraction

            flux[transportEqIdx] +=
                fluxVars.volumeFlux(phaseIdx) *
                ((    massUpwindWeight_)*up.density()*up.massFraction(transportCompIdx)
                 +
                 (1 - massUpwindWeight_)*dn.density()*dn.massFraction(transportCompIdx))
                 *fluxVars.coordinatesCenter()[0];
        }
        else //mole-fraction formulation
        {
            // total mass flux - molefraction
            //KmvpNormal is the Darcy velocity multiplied with the normal vector, calculated in 1p2cfluxvariables.hh
            if (!usePH)
                flux[contiEqIdx] *=
                    ((      massUpwindWeight_)*up.molarDensity()
                     +
                     ((1 - massUpwindWeight_)*dn.molarDensity()))
                    *fluxVars.coordinatesCenter()[0];

            // advective flux of the second component -molefraction
            flux[transportEqIdx] +=
                fluxVars.volumeFlux(phaseIdx)
                *((    massUpwindWeight_)*up.molarDensity() * up.moleFraction(transportCompIdx)
                 +
                 (1 - massUpwindWeight_)*dn.molarDensity() * dn.moleFraction(transportCompIdx))
                * fluxVars.coordinatesCenter()[0];
        }

    }

    /*!
     * \brief Adds the diffusive flux to the flux vector over
     *        the face of a sub-control volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each phase
     * \param fluxVars The flux variables at the current SCV
     *
     * This function doesn't do anything but may be used by the
     * non-isothermal three-phase models to calculate diffusive heat
     * fluxes
     */
    void computeDiffusiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        // diffusive fluxes
        Scalar tmp(0);
//
        tmp = -(fluxVars.radialConcentrationGrad(transportCompIdx)*fluxVars.face().normal);
        tmp *= fluxVars.porousDiffCoeff();;//* fluxVars.coordinatesCenter()[0];
//
    //    //dispersive flux of second component - massfraction
    //    GlobalPosition normalDisp;
    //    fluxVars.dispersionTensor().mv(fluxVars.face().normal, normalDisp);
    //    tmp -= (normalDisp * fluxVars.concentrationGrad(transportCompIdx));
//
    //    flux[transportEqIdx] += tmp;


    //    //flux[transportEqIdx]=0;
    //    // diffusive flux of second component
    //    if(!useMoles) //mass-fraction formulation
    //    {
//
            // diffusive flux of the second component - massfraction
        //    tmp = -(fluxVars.massFractionGrad(transportCompIdx)*fluxVars.face().normal);
        //    tmp *= fluxVars.porousDiffCoeff()*fluxVars.density();// * fluxVars.molarDensity();

            // dispersive flux of second component - massfraction
            //GlobalPosition normalDisp;
            //fluxVars.dispersionTensor().mv(fluxVars.face().normal, normalDisp);
            //tmp -= (normalDisp * fluxVars.massFractionGrad(transportCompIdx))*fluxVars.density();;
//
    //        // convert it to a mass flux and add it
            flux[transportEqIdx] += tmp;
             //* FluidSystem::molarMass(transportCompIdx);

//
    //        //if ((fluxVars.face().normal[2]!=0) &&
    //        //        (fluxVars.moleFractionGrad(transportCompIdx)*fluxVars.face().normal !=0))
//
    //        //       std::cout << "  computeDiffusiveFlux " << fluxVars.moleFractionGrad(transportCompIdx)
    //        //             <<" "<< fluxVars.face().normal <<" "<< fluxVars.porousDiffCoeff() <<" "
    //        //             << fluxVars.molarDensity() <<" "<<FluidSystem::molarMass(transportCompIdx)<<std::endl;
    //        //        std::cout << "          flux[transportEqIdx] " << flux[transportEqIdx]<< " "
    //        //             << (fluxVars.moleFractionGrad(transportCompIdx)*fluxVars.face().normal) << std::endl;;
    //        //std::cout << "                normalDisp "<< normalDisp << std::endl;
//
    //    }
    //    else //mole-fraction formulation
    //    {
    //        // diffusive flux of the second component - molefraction
    //        tmp = -(fluxVars.moleFractionGrad(transportCompIdx)*fluxVars.face().normal);
    //        tmp *= fluxVars.porousDiffCoeff() * fluxVars.molarDensity();
//
    //        // dispersive flux of second component - molefraction
    //        GlobalPosition normalDisp;
    //        fluxVars.dispersionTensor().mv(fluxVars.face().normal, normalDisp);
    //        tmp -= fluxVars.molarDensity()*
    //            (normalDisp * fluxVars.moleFractionGrad(transportCompIdx));
//
    //        flux[transportEqIdx] += tmp*FluidSystem::molarMass(transportCompIdx);
    //    }
    }
//    /*!
//     * \brief Calculate the source term of all equations.
//     *        The pressure gradient at the center of a SCV is computed
//     *        and the gravity term evaluated.
//     *
//     * \param source The source/sink in the sub control volume for each component
//     * \param scvIdx The local index of the sub-control volume
//     */
//    void computeSource(PrimaryVariables &source, const int scvIdx)
//    {
//        const VolumeVariables &volVars = this->curVolVars_(scvIdx);
//
//        this->problem_().solDependentSource(source,
//                                     this->element_(),
//                                     this->fvGeometry_(),
//                                     scvIdx,
//                                     this->curVolVars_());
//
//        // add contribution from possible point sources
//        this->problem_().scvPointSources(source,
//                                     this->element_(),
//                                     this->fvGeometry_(),
//                                     scvIdx,
//                                     this->curVolVars_());
//
//        Scalar tmp1;
//        tmp1 = volVars.porousDiffCoeff()
//                *volVars.concentrationGrad(transportCompIdx)
//                /volVars.coordinatesCenter()[0];
//        source[transportEqIdx] += tmp1;
//        }

protected:
    Implementation *asImp_()
    {
        return static_cast<Implementation *> (this);
    }

    const Implementation *asImp_() const
    {
        return static_cast<const Implementation *> (this);
    }

private:
    Scalar massUpwindWeight_;


};

}

#endif
