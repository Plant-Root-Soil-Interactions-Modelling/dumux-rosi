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
 * \brief This file contains the data which is required to calculate
 *        volume fluxes of fluid phases over a face of a finite volume by means
 *        of the Darcy approximation.
 *
 */
#ifndef DUMUX_MPFA_IMPLICIT_DARCY_FLUX_VARIABLES_HH
#define DUMUX_MPFA_IMPLICIT_DARCY_FLUX_VARIABLES_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

//#include <dumux/implicit/properties.hh>


namespace Dumux
{

namespace Properties
{
// forward declaration of properties
NEW_PROP_TAG(ImplicitMobilityUpwindWeight);
NEW_PROP_TAG(SpatialParams);
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(ProblemEnableGravity);
}

/*!
 * \ingroup ImplicitFluxVariables
 * \brief Evaluates the normal component of the Darcy velocity
 * on a (sub)control volume face.
 */
template <class TypeTag>
class MpfaImplicitDarcyFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, MpfaInteractionVolume) InteractionVolume;
    typedef typename GET_PROP_TYPE(TypeTag, MpfaBoundaryInteractionVolume) BoundaryInteractionVolume;
    typedef typename GET_PROP_TYPE(TypeTag, MpfaFluxCalculator) FluxCalculator;
    typedef typename GET_PROP_TYPE(TypeTag, MpfaBoundaryFluxCalculator) BoundaryFluxCalculator;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryLayers) BoundaryLayers;
    typedef typename GET_PROP(TypeTag, MpfaMethods) MpfaMethods;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename Dune::ReferenceElements<Scalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<Scalar, dim> ReferenceElement;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases)};

    //typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;
    //typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;
    //typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    //typedef Dune::FieldVector<Scalar, dim> DimVector;

    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

public:
    /*!
     * \brief The constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param fIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    MpfaImplicitDarcyFluxVariables(const Problem &problem,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const int fIdx,
                 const ElementVolumeVariables &elemVolVars,
                 const bool onBoundary = false)
    : fvGeometry_(fvGeometry), faceIdx_(fIdx), onBoundary_(onBoundary)
    {
        mobilityUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MobilityUpwindWeight);

        massUpwindWeight_ = 1; // GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);

        calculateFlux(problem, element, elemVolVars);
    }

public:

        /*!
     * \brief Return the volumetric flux over a face of a given phase.
     *
     *        This is the calculated velocity multiplied by the unit normal
     *        and the area of the face. face().normal has already the
     *        magnitude of the area.
     *
     * \param phaseIdx index of the phase
     */
    Scalar massFlux(const unsigned int phaseIdx) const
    { return massFlux_[phaseIdx]; }


    Scalar normalFlux(const unsigned int phaseIdx) const
    { return normalFlux_[phaseIdx]; }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase.
     *
     * \param phaseIdx index of the phase
     */
    const unsigned int downstreamIdx(const unsigned phaseIdx) const
    { return downstreamIdx_[phaseIdx]; }

    /*!
     * \brief Return the local index of the upstream control volume
     *        for a given phase.
     *
     * \param phaseIdx index of the phase
     */
    const unsigned int upstreamIdx(const unsigned phaseIdx) const
    { return upstreamIdx_[phaseIdx]; }

    /*!
     * \brief Return the SCV (sub-control-volume) face. This may be either
     *        a face within the element or a face on the element boundary,
     *        depending on the value of onBoundary_.
     */
    const SCVFace &face() const
    {
            return fvGeometry_.subContVolFace[faceIdx_];
    }



protected:

    /*!
     * \brief Calculation of the flux through an element face
     *
     * \param problem The problem
     * \param element The finite element
     * \param elemVolVars The volume variables of the current element
     */
    void calculateFlux(const Problem &problem,
                       const Element &element,
                       const ElementVolumeVariables &elemVolVars)
    {
        // Reference element of the current element
        const ReferenceElement& referenceElement = ReferenceElements::general(element.geometry().type());

        // loop over the phases
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
        {
            massFlux_[phaseIdx] = 0;
            normalFlux_[phaseIdx] = 0;

            // loop over the nodes of the face
            for (int faceNode = 0; faceNode < element.template subEntity<1>(faceIdx_).geometry().corners(); faceNode++)
            {
                // get the local node indices of the nodes with respect to element
                int localVertIdx = referenceElement.subEntity(faceIdx_, 1, faceNode, dim);

                // check whether we have to handle a boundary interaction volume or not
                int boundaryLayer = problem.model().interactionVolumeContainer().isBoundaryVolume(element, faceIdx_, localVertIdx);

                // calculate flux
                if (boundaryLayer == BoundaryLayers::interior && GET_PROP_VALUE(TypeTag, MpfaMethod) != MpfaMethods::oMethod)
                {
                    InteractionVolume& interactionVolume = *(problem.model().getInnerInteractionVolume(element, faceIdx_, localVertIdx));
                    int faceIndexInRegion = -1;
                    typename BoundaryInteractionVolume::DynamicVector neumannFluxes = BoundaryInteractionVolume::DynamicVector::Zero(interactionVolume.getNumberOfSubFaces());
                    Scalar tempVolumeFlux = FluxCalculator::calculateScvFaceFlux(problem, element, faceIdx_, interactionVolume, elemVolVars, fvGeometry_, phaseIdx, faceIndexInRegion, neumannFluxes);
                    normalFlux_[phaseIdx] += tempVolumeFlux + neumannFluxes(faceIndexInRegion);

                    // determine the upwind direction
                    if (tempVolumeFlux > 0)
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

                    tempVolumeFlux *= mobilityUpwindWeight_*upVolVars.mobility(phaseIdx)
                            + (1.0 - mobilityUpwindWeight_)*downVolVars.mobility(phaseIdx);

                    tempVolumeFlux *= massUpwindWeight_*upVolVars.fluidState().density(phaseIdx)
                            + (1.0 - massUpwindWeight_)*downVolVars.fluidState().density(phaseIdx);

                    massFlux_[phaseIdx] += (tempVolumeFlux + neumannFluxes(faceIndexInRegion))*elemVolVars[0].extrusionFactor();
                }
                else
                {
                    // initialize pointer to interaction volume
                    std::shared_ptr<BoundaryInteractionVolume> interactionVolume = std::make_shared<BoundaryInteractionVolume>();

                    if (boundaryLayer == BoundaryLayers::boundary)
                        interactionVolume = problem.model().getBoundaryInteractionVolume(element, faceIdx_, localVertIdx, phaseIdx);
                    else
                        // this gets any other interaction volume for the o-method
                        // and the volumes on the second half edge on a boundary for another method with full o-method boundary
                        interactionVolume = problem.model().getInsideInteractionVolume(element, faceIdx_, localVertIdx);

                    int faceIndexInRegion = -1;
                    typename BoundaryInteractionVolume::DynamicVector neumannFluxes = BoundaryInteractionVolume::DynamicVector::Zero(interactionVolume->getNumberOfSubFaces());
                    Scalar tempVolumeFlux = BoundaryFluxCalculator::calculateScvFaceFlux(problem, element, faceIdx_, *interactionVolume, elemVolVars, fvGeometry_, phaseIdx, faceIndexInRegion, neumannFluxes);
                    typename BoundaryInteractionVolume::SubVolumeFace& subFace = interactionVolume->getSubVolumeFace(faceIndexInRegion);
                    normalFlux_[phaseIdx] += tempVolumeFlux + neumannFluxes(faceIndexInRegion);

                    // initial guess of upstream and downstream volume variables
                    upstreamIdx_[phaseIdx] = face().i;
                    downstreamIdx_[phaseIdx] = face().j;
                    VolumeVariables upVolVarsTemp = elemVolVars[ upstreamIdx(phaseIdx) ];
                    VolumeVariables downVolVarsTemp = elemVolVars[ downstreamIdx(phaseIdx) ];
                    VolumeVariables &upVolVars = upVolVarsTemp;
                    VolumeVariables &downVolVars = downVolVarsTemp;

                    // set the flux
                    if (subFace.faceType == BoundaryInteractionVolume::FaceTypes::InternalDirichletFace)
                    {
                        PrimaryVariables interiorDirichletValues(0); VolumeVariables interiorDirichletVolVars;
                        problem.internalDirichlet(interiorDirichletValues, element, faceIdx_, *interactionVolume);
                        interiorDirichletVolVars.update(interiorDirichletValues, problem, element, fvGeometry_, 0, false);

                        // determine the upwind direction
                        if (tempVolumeFlux > 0)
                        {
                            // change downwind variables, because we are on an interior dirichlet face
                            if (GET_PROP_VALUE(TypeTag, FacetCoupling))
                                problem.getVolumeVarsOnFacet(downVolVars, element, faceIdx_, *interactionVolume);
                            else
                                downVolVars = interiorDirichletVolVars;
                        }
                        else
                        {
                            upstreamIdx_[phaseIdx] = face().j;
                            downstreamIdx_[phaseIdx] = face().i;

                            // change the upwind volume variables
                            if (GET_PROP_VALUE(TypeTag, FacetCoupling))
                                problem.getVolumeVarsOnFacet(upVolVars, element, faceIdx_, *interactionVolume);
                            else
                                upVolVars = interiorDirichletVolVars;
                            downVolVars = elemVolVars[ downstreamIdx(phaseIdx) ];
                        }

                        tempVolumeFlux *= mobilityUpwindWeight_*upVolVars.mobility(phaseIdx)
                                    + (1.0 - mobilityUpwindWeight_)*downVolVars.mobility(phaseIdx);
                        tempVolumeFlux *= massUpwindWeight_*upVolVars.fluidState().density(phaseIdx)
                                    + (1.0 - massUpwindWeight_)*downVolVars.fluidState().density(phaseIdx);
                    }
                    else if (subFace.faceType != BoundaryInteractionVolume::FaceTypes::NeumannFace
                                && subFace.faceType != BoundaryInteractionVolume::FaceTypes::InternalFluxFace)
                    {
                        // determine the upwind direction
                        if (tempVolumeFlux > 0)
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
                        upVolVars = elemVolVars[ upstreamIdx(phaseIdx) ];
                        downVolVars = elemVolVars[ downstreamIdx(phaseIdx) ];

                        tempVolumeFlux *= mobilityUpwindWeight_*upVolVars.mobility(phaseIdx)
                                + (1.0 - mobilityUpwindWeight_)*downVolVars.mobility(phaseIdx);

                        tempVolumeFlux *= massUpwindWeight_*upVolVars.fluidState().density(phaseIdx)
                                + (1.0 - massUpwindWeight_)*downVolVars.fluidState().density(phaseIdx);
                    }
                    else if (subFace.faceType == BoundaryInteractionVolume::FaceTypes::InternalFluxFace
                                && GET_PROP_VALUE(TypeTag, FacetCoupling))
                    {
                        // determine the upwind direction
                        if (tempVolumeFlux > 0)
                        {
                            // change downwind variables
                             problem.getVolumeVarsOnFacet(downVolVars, element, faceIdx_, *interactionVolume);
                        }
                        else
                        {
                            upstreamIdx_[phaseIdx] = face().j;
                            downstreamIdx_[phaseIdx] = face().i;

                            // change the volume variables
                            problem.getVolumeVarsOnFacet(upVolVars, element, faceIdx_, *interactionVolume);
                            downVolVars = elemVolVars[ downstreamIdx(phaseIdx) ];
                        }

                        tempVolumeFlux *= mobilityUpwindWeight_*upVolVars.mobility(phaseIdx)
                                    + (1.0 - mobilityUpwindWeight_)*downVolVars.mobility(phaseIdx);
                        tempVolumeFlux *= massUpwindWeight_*upVolVars.fluidState().density(phaseIdx)
                                    + (1.0 - massUpwindWeight_)*downVolVars.fluidState().density(phaseIdx);
                    }

                    // add the flux of the half face to the flux of the whole face
                    // in case of a neumann/internalFlux boundary, only the contribution of the neumann fluxes
                    // has to be considered. "tempVolumeFlux" for these half faces are (should be) zero (or rather eps) anyway
                   // if (subFace.faceType == BoundaryInteractionVolume::FaceTypes::NeumannFace
                   //             || subFace.faceType == BoundaryInteractionVolume::FaceTypes::InternalFluxFace)
                   //     massFlux_[phaseIdx] += neumannFluxes(faceIndexInRegion)*elemVolVars[0].extrusionFactor();
                   // else
                    massFlux_[phaseIdx] += (tempVolumeFlux + neumannFluxes(faceIndexInRegion))*elemVolVars[0].extrusionFactor();
                }
            }
        }
    }

    const FVElementGeometry &fvGeometry_;   	//!< Information about the geometry of discretization
    const unsigned int faceIdx_;            	//!< The index of the sub control volume face
    const bool      onBoundary_;                //!< Specifying whether we are currently on the boundary of the simulation domain
    unsigned int    upstreamIdx_[numPhases] , downstreamIdx_[numPhases]; //!< local index of the upstream / downstream vertex
    Scalar          massFlux_[numPhases] ;    //!< Mass flux over the face
    Scalar          normalFlux_[numPhases] ;    //!< Velocity multiplied with normal (magnitude=area)
    Scalar          mobilityUpwindWeight_;      //!< Upwind weight for mobility. Set to one for full upstream weighting
    Scalar          massUpwindWeight_;
};

} // end namespace

#endif // DUMUX_MPFA_IMPLICIT_DARCY_FLUX_VARIABLES_HH
