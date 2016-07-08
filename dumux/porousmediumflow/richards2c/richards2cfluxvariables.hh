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
 * \brief Overwrites the volumeFlux() function from the ImplicitDarcyFluxVariables.
 *
 */
#ifndef DUMUX_RICHARDS_2C_FLUX_VARIABLES_HH
#define DUMUX_RICHARDS_2C_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
//#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include "richards2cproperties.hh"

namespace Dumux
{


/*!
 * \ingroup RichardsTwoCModel
 * \ingroup ImplicitFluxVariables
 * \brief Overwrites the volumeFlux() function from the ImplicitDarcyFluxVariables.
 */
template <class TypeTag>
class RichardsTwoCFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel) EffectiveDiffusivityModel;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;

    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases)} ;
    enum { nPhaseIdx = Indices::nPhaseIdx} ;
    enum { transportCompIdx = Indices::transportCompIdx };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;

public:

    /*!
     * \brief Compute / update the flux variables
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param fIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     * \todo The fvGeometry should be better initialized, passed and stored as an std::shared_ptr
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int fIdx,
                const ElementVolumeVariables &elemVolVars,
                const bool onBoundary = false)
    {
        fvGeometryPtr_ = &fvGeometry;
        onBoundary_ = onBoundary;
        faceIdx_ = fIdx;

        mobilityUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MobilityUpwindWeight);

        asImp_().calculateGradients_(problem, element, elemVolVars);
        asImp_().calculateNormalVelocity_(problem, element, elemVolVars);
        asImp_().calculatePorousDiffCoeff_(problem, element, elemVolVars);
        asImp_().calculateDispersionTensor_(problem, element, elemVolVars);
    }

   /*!
    * \brief Return the pressure potential multiplied with the
    *        intrinsic permeability  and the face normal which
    *        goes from vertex i to vertex j.
    *
    * Note that the length of the face's normal is the area of the
    * phase, so this is not the actual velocity but the integral of
    * the velocity over the face's area. Also note that the phase
    * mobility is not yet included here since this would require a
    * decision on the upwinding approach (which is done in the
    * actual model).
    */
   Scalar KmvpNormal() const
   { return KmvpNormal_; }

   /*!
    * \brief Return the pressure potential multiplied with the
    *        intrinsic permeability as vector (for velocity output).
    */
   GlobalPosition Kmvp() const
   { return Kmvp_; }

   /*!
    * \brief The face of the current sub-control volume. This may be either
    *        an inner sub-control-volume face or a SCV face on the boundary.
    */
   const SCVFace &face() const
   {
       if (onBoundary_)
           return fvGeometry_().boundaryFace[faceIdx_];
       else
           return fvGeometry_().subContVolFace[faceIdx_];
   }

    /*!
     * \brief Return the intrinsic permeability tensor \f$\mathrm{[m^2]}\f$.
     */
    const DimWorldMatrix &intrinsicPermeability() const
    { return K_; }

    /*!
     * \brief Return the dispersion tensor \f$\mathrm{[m^2/s]}\f$.
     */
    const DimWorldMatrix &dispersionTensor() const
    { return dispersionTensor_; }

    /*!
     * \brief Return the pressure potential gradient \f$\mathrm{[Pa/m]}\f$.
     */
    const GlobalPosition &potentialGrad(int phaseIdx) const
    { return potentialGrad_[phaseIdx]; }


    /*!
     * \brief Return the mole-fraction gradient of a component in a phase \f$\mathrm{[mol/mol/m)]}\f$.
     *
     * \param compIdx The index of the considered component
     */
    const GlobalPosition &moleFractionGrad(int compIdx) const
    {
       if (compIdx != 1)
       { DUNE_THROW(Dune::InvalidStateException,
                "The 1p2c model is supposed to need "
                "only the concentration gradient of "
                "the second component!"); }
       return moleFractionGrad_;
    };

    /*!
    * \brief The binary diffusion coefficient for each fluid phase in the porous medium \f$\mathrm{[m^2/s]}\f$.
    */
    Scalar porousDiffCoeff() const
    {
        // TODO: tensorial porousDiffCoeff_usion coefficients
        return porousDiffCoeff_;
    };

    /*!
    * \brief Return viscosity \f$\mathrm{[Pa s]}\f$ of a phase at the integration
    *        point.
    */
    Scalar viscosity() const
    { return viscosity_;}

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ of a phase at the integration
     *        point.
     */
    Scalar molarDensity() const
    { return molarDensity_; }

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ of a phase at the integration
     *        point.
     */
    Scalar density() const
    { return density_; }

    /*!
     * \brief Given the intrinsic permeability times the pressure
     *        potential gradient and SCV face normal for a phase,
     *        return the local index of the upstream control volume
     *        for a given phase.
     *
     *        \param normalFlux The flux over a face of the sub-control volume
     */
    int upstreamIdx(Scalar normalFlux) const
    { return (normalFlux >= 0)?face().i:face().j; }

    /*!
     * \brief Given the intrinsic permeability times the pressure
     *        potential gradient and SCV face normal for a phase,
     *        return the local index of the downstream control volume
     *        for a given phase.
     *
     *        \param normalFlux The flux over a face of the sub-control volume
     */
    int downstreamIdx(Scalar normalFlux) const
    { return (normalFlux > 0)?face().j:face().i; }

    /*!
    * \brief Return the local index of the upstream control volume
    *        for a given phase.
    */
    int upstreamIdx(int phaseIdx) const
    { return upstreamIdx_[phaseIdx]; }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase.
     */
    int downstreamIdx(int phaseIdx) const
    { return downstreamIdx_[phaseIdx]; }

    /*!
     * \brief Return the volumetric flux over a face of a given phase.
     *
     *        This is the calculated velocity multiplied by the unit normal
     *        and the area of the face.
     *        face().normal
     *        has already the magnitude of the area.
     *
     * \param phaseIdx index of the phase
     */
    Scalar volumeFlux(const unsigned int phaseIdx) const
    {
        return volumeFlux_[phaseIdx];
    }

protected:
    //! Returns the implementation of the flux variables (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    /*!
     * \brief Calculation of the potential gradients
     *
     * \param problem The problem
     * \param element The finite element
     * \param elemVolVars The volume variables of the current element
     */
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars)
    {
        // loop over all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
        {
            potentialGrad_[phaseIdx]= 0.0;

            for (unsigned int idx = 0;
                 idx < face().numFap;
                 idx++) // loop over adjacent vertices
            {
                // FE gradient at vertex idx
                const GlobalPosition &feGrad = face().grad[idx];

                // index for the element volume variables
                int volVarsIdx = face().fapIndices[idx];

                // the pressure gradient
                GlobalPosition tmp(feGrad);
                tmp *= elemVolVars[volVarsIdx].fluidState().pressure(phaseIdx);
                potentialGrad_[phaseIdx] += tmp;
            }

            // correct the pressure gradient by the gravitational acceleration
            if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            {
                // ask for the gravitational acceleration at the given SCV face
                GlobalPosition g(problem.gravityAtPos(face().ipGlobal));

                // calculate the phase density at the integration point. we
                // only do this if the wetting phase is present in both cells
                Scalar SI = elemVolVars[face().i].fluidState().saturation(phaseIdx);
                Scalar SJ = elemVolVars[face().j].fluidState().saturation(phaseIdx);
                Scalar rhoI = elemVolVars[face().i].fluidState().density(phaseIdx);
                Scalar rhoJ = elemVolVars[face().j].fluidState().density(phaseIdx);
                Scalar fI = std::max(0.0, std::min(SI/1e-5, 0.5));
                Scalar fJ = std::max(0.0, std::min(SJ/1e-5, 0.5));
                if (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(fI + fJ, 0.0, 1.0e-30))
                    // doesn't matter because no wetting phase is present in
                    // both cells!
                    fI = fJ = 0.5;
                const Scalar density = (fI*rhoI + fJ*rhoJ)/(fI + fJ);

                // make gravity acceleration a force
                GlobalPosition f(g);
                f *= density;

                // calculate the final potential gradient
                potentialGrad_[phaseIdx] -= f;
            } // gravity
        } // loop over all phases
     }

    /*!
     * \brief Actual calculation of the normal Darcy velocities.
     *
     * \param problem The problem
     * \param element The finite element
     * \param elemVolVars The volume variables of the current element
     */
    void calculateNormalVelocity_(const Problem &problem,
                                  const Element &element,
                                  const ElementVolumeVariables &elemVolVars)
    {
        // calculate the mean intrinsic permeability
        const SpatialParams &spatialParams = problem.spatialParams();
        DimWorldMatrix K;
        if (GET_PROP_VALUE(TypeTag, ImplicitIsBox))
        {
            spatialParams.meanK(K,
                                spatialParams.intrinsicPermeability(element,
                                                                    fvGeometry_(),
                                                                    face().i),
                                spatialParams.intrinsicPermeability(element,
                                                                    fvGeometry_(),
                                                                    face().j));
        }
        else
        {
            const Element& elementI = fvGeometry_().neighbors[face().i];
            FVElementGeometry fvGeometryI;
            fvGeometryI.subContVol[0].global = elementI.geometry().center();

            const Element& elementJ = fvGeometry_().neighbors[face().j];
            FVElementGeometry fvGeometryJ;
            fvGeometryJ.subContVol[0].global = elementJ.geometry().center();

            spatialParams.meanK(K,
                                spatialParams.intrinsicPermeability(elementI, fvGeometryI, 0),
                                spatialParams.intrinsicPermeability(elementJ, fvGeometryJ, 0));
        }

        // loop over all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
        {
            // calculate the flux in the normal direction of the
            // current sub control volume face:
            //
            // v = - (K_f grad phi) * n
            // with K_f = rho g / mu K
            //
            // Mind, that the normal has the length of it's area.
            // This means that we are actually calculating
            //  Q = - (K grad phi) dot n /|n| * A


            K.mv(potentialGrad_[phaseIdx], kGradP_[phaseIdx]);
            kGradPNormal_[phaseIdx] = kGradP_[phaseIdx]*face().normal;

            // determine the upwind direction
            if (kGradPNormal_[phaseIdx] < 0)
            {
                upstreamIdx_[phaseIdx] = face().i;
                downstreamIdx_[phaseIdx] = face().j;
            }
            else
            {
                upstreamIdx_[phaseIdx] = face().j;
                downstreamIdx_[phaseIdx] = face().i;
            }

            // obtain the upwind volume variables
            const VolumeVariables& upVolVars = elemVolVars[ upstreamIdx(phaseIdx) ];
            const VolumeVariables& downVolVars = elemVolVars[ downstreamIdx(phaseIdx) ];

            // the minus comes from the Darcy relation which states that
            // the flux is from high to low potentials.
            // set the velocity
            velocity_[phaseIdx] = kGradP_[phaseIdx];
            velocity_[phaseIdx] *= - ( mobilityUpwindWeight_*upVolVars.mobility(phaseIdx)
                    + (1.0 - mobilityUpwindWeight_)*downVolVars.mobility(phaseIdx)) ;

            // set the volume flux
            volumeFlux_[phaseIdx] = velocity_[phaseIdx] * face().normal;
        }// loop all phases
    }

    /*!
    * \brief Calculation of the effective diffusion coefficient.
    *
    *        \param problem The considered problem file
    *        \param element The considered element of the grid
    *        \param elemVolVars The parameters stored in the considered element
    */
    void calculatePorousDiffCoeff_(const Problem &problem,
                                   const Element &element,
                                   const ElementVolumeVariables &elemVolVars)
    {
        const VolumeVariables &volVarsI = elemVolVars[face().i];
        const VolumeVariables &volVarsJ = elemVolVars[face().j];

        const Scalar diffCoeffI = EffectiveDiffusivityModel::effectiveDiffusivity(volVarsI.porosity(),
                                                                     /*sat=*/1.0,
                                                                     volVarsI.diffCoeff());

        const Scalar diffCoeffJ = EffectiveDiffusivityModel::effectiveDiffusivity(volVarsJ.porosity(),
                                                                     /*sat=*/1.0,
                                                                     volVarsJ.diffCoeff());

        // -> harmonic mean
        porousDiffCoeff_ = harmonicMean(diffCoeffI, diffCoeffJ);
    }

    /*!
    * \brief Calculation of the dispersion.
    *
    *        \param problem The considered problem file
    *        \param element The considered element of the grid
    *        \param elemVolVars The parameters stored in the considered element
    */
    void calculateDispersionTensor_(const Problem &problem,
                                    const Element &element,
                                    const ElementVolumeVariables &elemVolVars)
    {
        const VolumeVariables &volVarsI = elemVolVars[face().i];
        const VolumeVariables &volVarsJ = elemVolVars[face().j];

        //calculate dispersivity at the interface: [0]: alphaL = longitudinal disp. [m], [1] alphaT = transverse disp. [m]
        Scalar dispersivity[2];
        dispersivity[0] = 0.5 * (volVarsI.dispersivity()[0] +  volVarsJ.dispersivity()[0]);
        dispersivity[1] = 0.5 * (volVarsI.dispersivity()[1] +  volVarsJ.dispersivity()[1]);

        //calculate velocity at interface: v = -1/mu * vDarcy = -1/mu * K * grad(p)
        GlobalPosition velocity;
        Valgrind::CheckDefined(potentialGrad(0));
        Valgrind::CheckDefined(K_);
        K_.mv(potentialGrad(0), velocity);
        velocity /= - 0.5 * (volVarsI.viscosity(0) + volVarsJ.viscosity(0));

        //matrix multiplication of the velocity at the interface: vv^T
        dispersionTensor_ = 0;
        for (int i=0; i<dim; i++)
            for (int j = 0; j<dim; j++)
                dispersionTensor_[i][j] = velocity[i]*velocity[j];

        //normalize velocity product --> vv^T/||v||, [m/s]
        Scalar vNorm = velocity.two_norm();

        dispersionTensor_ /= vNorm;
        if (vNorm < 1e-20)
            dispersionTensor_ = 0;

        //multiply with dispersivity difference: vv^T/||v||*(alphaL - alphaT), [m^2/s] --> alphaL = longitudinal disp., alphaT = transverse disp.
        dispersionTensor_ *= (dispersivity[0] - dispersivity[1]);

        //add ||v||*alphaT to the main diagonal:vv^T/||v||*(alphaL - alphaT) + ||v||*alphaT, [m^2/s]
        for (int i = 0; i<dim; i++)
            dispersionTensor_[i][i] += vNorm*dispersivity[1];
    }

    //const FVElementGeometry &fvGeometry_;
    // return const reference to fvGeometry
    const FVElementGeometry& fvGeometry_() const
    { return *fvGeometryPtr_; }

    int faceIdx_;
    bool onBoundary_;

    //! pressure potential gradient
    GlobalPosition potentialGrad_[numPhases] ;
    //! mole-fraction gradient
    GlobalPosition moleFractionGrad_;
    //! the effective diffusion coefficent in the porous medium
    Scalar porousDiffCoeff_;

    //! the dispersion tensor in the porous medium
    DimWorldMatrix dispersionTensor_;

    //! the intrinsic permeability tensor
    DimWorldMatrix K_;
    // intrinsic permeability times pressure potential gradient
    GlobalPosition Kmvp_;
    // projected on the face normal
    Scalar KmvpNormal_;
    GlobalPosition  kGradP_[numPhases] ; //!< Permeability multiplied with gradient in potential
    Scalar          kGradPNormal_[numPhases] ;  //!< Permeability multiplied with gradient in potential, multiplied with normal (magnitude=area)
    GlobalPosition  velocity_[numPhases] ;      //!< The velocity as determined by Darcy's law or by the Forchheimer relation
    unsigned int    upstreamIdx_[numPhases] , downstreamIdx_[numPhases]; //!< local index of the upstream / downstream vertex

    //! viscosity of the fluid at the integration point
    Scalar viscosity_;

    //! molar densities of the fluid at the integration point
    Scalar molarDensity_, density_;

    Scalar volumeFlux_[numPhases]; //!< Velocity multiplied with normal (magnitude=area)
    Scalar mobilityUpwindWeight_; //!< Upwind weight for mobility. Set to one for full upstream weighting

private:
    const FVElementGeometry* fvGeometryPtr_; //!< Information about the geometry of discretization
};

} // end namespace

#endif
