

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
 * \brief Element-wise calculation the local Jacobian for the single-phase,
 *        two-component one-dimensional stokes model in the fully implicit scheme.
 */

#ifndef DUMUX_ROOTSYSTEM_LOCAL_RESIDUAL_1P2C_HH
#define DUMUX_ROOTSYSTEM_LOCAL_RESIDUAL_1P2C_HH

#include <cmath>
#include "properties1p2c.hh"

namespace Dumux
{
/*!
 *
 * \ingroup OneDStokesTwoCModel
 * \ingroup ImplicitLocalResidual
 * \brief Calculate the local Jacobian for the single-phase,
 *        two-component model in the fully implicit scheme.
 *
 *  This class is used to fill the gaps in BaseLocalResidual for the 1dstokes2c flow and transport.
 */
template<class TypeTag>
class RootSystemOnePTwoCLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        //phase index
        phaseIdx = Indices::phaseIdx,
        transportCompIdx = Indices::transportCompIdx
    };
    // indices of the equations
    enum {
        conti0EqIdx = Indices::conti0EqIdx,
        transportEqIdx = Indices::transportEqIdx
    };

    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    RootSystemOnePTwoCLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        upwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
    };

    /*!
     * \brief Evaluate the amount of all conservation quantities
     *        (e.g. phase mass) within a finite volume.
     *
     *        \param storage The mass of the component within the sub-control volume
     *        \param scvIdx The index of the considered face of the sub-control volume
     *        \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &storage, const int scvIdx, const bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];


        Scalar radius = this->problem_().spatialParams().rootRadius(this->element_(),this->fvGeometry_(),scvIdx);
        storage[conti0EqIdx] += M_PI*radius*radius*volVars.density()*volVars.porosity();
        storage = 0;
        if(!useMoles) //mass-fraction formulation
        {
            //storage term of the transport equation - massfractions
            storage[transportEqIdx] += M_PI*radius*radius*volVars.density()*volVars.massFraction(transportCompIdx)*volVars.porosity();
        }
        else //mole-fraction formulation
        {
            // storage term of the transport equation - molefractions
            storage[transportEqIdx] += M_PI*radius*radius*volVars.molarDensity()*volVars.moleFraction(transportCompIdx)*volVars.porosity();
        }

    }

    /*!
     * \brief Evaluate the mass flux over a face of a sub-control
     *        volume.
     *
     *        \param flux The flux over the SCV (sub-control-volume) face for each component
     *        \param faceIdx The index of the considered face of the sub control volume
     *        \param onBoundary A boolean variable to specify whether the flux variables
     *               are calculated for interior SCV faces or boundary faces, default=false
     */
    void computeFlux(PrimaryVariables &flux, const int faceIdx, const bool onBoundary=false) const
    {
        flux = 0;
        FluxVariables fluxVars;
        fluxVars.update(this->problem_(),
                        this->element_(),
                        this->fvGeometry_(),
                        faceIdx,
                        this->curVolVars_(),
                        onBoundary);

        asImp_()->computeAdvectiveFlux(flux, fluxVars, faceIdx);
        asImp_()->computeDiffusiveFlux(flux, fluxVars);
    }

    /*!
     * \brief Evaluate the advective mass flux of all components over
     *        a face of a sub-control volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV
     */
    void computeAdvectiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars, const int faceIdx) const
    {
        ///////////////////////////////////////////////////
        // advective fluxes of all components in all phases
        ///////////////////////////////////////////////////

        const VolumeVariables &up =
            this->curVolVars_(fluxVars.upstreamIdx(phaseIdx));
        const VolumeVariables &dn =
            this->curVolVars_(fluxVars.downstreamIdx(phaseIdx));

        // total mass flux
        flux[conti0EqIdx] += fluxVars.volumeFlux(phaseIdx)*
            ((upwindWeight_)*up.density()
            + (1-upwindWeight_)*dn.density());

        if(!useMoles) //mass-fraction formulation
        {
            // advective flux of the second component - massfraction
            flux[transportEqIdx] += fluxVars.volumeFlux(phaseIdx)*
                ((upwindWeight_)*up.massFraction(transportCompIdx)*up.density()
                + (1-upwindWeight_)*dn.massFraction(transportCompIdx)*dn.density());
            Valgrind::CheckDefined(flux[transportEqIdx]);

        }
        else //mole-fraction formulation
        {

            // advective flux of the second component - molefraction
            flux[transportEqIdx] += fluxVars.volumeFlux(phaseIdx) *
                ((upwindWeight_)*up.moleFraction(transportCompIdx)*up.molarDensity()
                + (1-upwindWeight_)*dn.moleFraction(transportCompIdx)*dn.molarDensity());
            Valgrind::CheckDefined(flux[transportEqIdx]);
        }

        Valgrind::CheckDefined(flux[conti0EqIdx]);
    }

    /*!
     * \brief Adds the diffusive mass flux of all components over
     *        a face of a sub-control volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV
     */
    void computeDiffusiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        Scalar tmp(0);

        // diffusive flux of second component
        if(!useMoles) //mass-fraction formulation
        {
            tmp = fluxVars.diffusiveFlux(transportCompIdx);
            // convert it to a mass flux and add it
            flux[transportEqIdx] += tmp * FluidSystem::molarMass(transportCompIdx);
        }
        else //mole-fraction formulation
        {
            tmp = fluxVars.diffusiveFlux(transportCompIdx);
            flux[transportEqIdx] += tmp;
        }
        Valgrind::CheckDefined(flux[transportEqIdx]);
    }

    /*!
     * \brief Return the temperature given the solution vector of a
     *        finite volume.
     */
    template <class PrimaryVariables>
    Scalar temperature(const PrimaryVariables &priVars)
    { return this->problem_.temperature(); /* constant temperature */ }

    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }

private:
    Scalar upwindWeight_;
};

} // end namespace Dumux

#endif