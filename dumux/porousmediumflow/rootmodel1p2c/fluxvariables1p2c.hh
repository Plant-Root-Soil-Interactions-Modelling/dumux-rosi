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
 *        all fluxes in a rootsystem over a face.
 *
 * This means pressure and temperature gradients, phase densities at
 * the integration point, etc.
 */
#ifndef DUMUX_ROOTSYSTEM_FLUX_VARIABLES_1P2C_HH
#define DUMUX_ROOTSYSTEM_FLUX_VARIABLES_1P2C_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>
#include <dumux/porousmediumflow/1d/rootsystem/fluxvariables.hh>
//#include <dumux/porousmediumflow/1p2c/implicit/fluxvariables.hh>

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
 * \ingroup RootSystemFluxVariables
 * \brief Evaluates the normal component of the Darcy velocity
 * on a (sub)control volume face.
 */
template <class TypeTag>
class RootsystemFluxOnePTwoCFluxVariables : public RootsystemFluxVariables<TypeTag>
//class RootsystemFluxOnePTwoCFluxVariables : public OnePTwoCFluxVariables<TypeTag>
{
    using ParentType = Dumux::RootsystemFluxVariables<TypeTag>;
    //using ParentType = Dumux::OnePTwoCFluxVariables<TypeTag>;
    friend class Dumux::RootsystemFluxVariables<TypeTag>;
    //friend class Dumux::OnePTwoCFluxVariables<TypeTag>;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
     typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases)} ;
    enum { phaseIdx = Indices::phaseIdx };
    enum { transportCompIdx = Indices::transportCompIdx };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

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

     calculateGradients_(problem, element, elemVolVars);
     calculateNormalVelocity_(problem, element, elemVolVars);
     calculateDiffusiveFlux_(problem, element, elemVolVars);

     }

     /*!
     * \brief Return the diffusive flux over a face of a given phase.
     *
     * \param compIdx index of the component
     */
    const Scalar diffusiveFlux(const unsigned int compIdx) const
    { return diffusiveFlux_; }

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
    { return volumeFlux_; }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase.
     *
     * \param phaseIdx index of the phase
     */
    const unsigned int downstreamIdx(const unsigned int phaseIdx = 0) const
    { return downstreamIdx_; }

    /*!
     * \brief Return the local index of the upstream control volume
     *        for a given phase.
     *
     * \param phaseIdx index of the phase
     */
    const unsigned int upstreamIdx(const unsigned int phaseIdx = 0) const
    { return upstreamIdx_; }

    /*!
     * \brief Return the SCV (sub-control-volume) face. This may be either
     *        a face within the element or a face on the element boundary,
     *        depending on the value of onBoundary_.
     */
    const SCVFace &face() const
    {
        if (onBoundary_)
            return fvGeometry_().boundaryFace[faceIdx_];
        else
            return fvGeometry_().subContVolFace[faceIdx_];
    }

    const Scalar transmissibility() const
    { return transmissibility_;}

protected:

    /*
     * \brief Calculation of the potential gradients
     *
     * \param problem The problem
     * \param element The finite element
     * \param elemVolVars The volume variables of the current element
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars)
    {
        transmissibility_ = 0.0;
        transmissibilityC_ = 0.0;
        density_ = 0.0;

        const VolumeVariables &volVarsI = elemVolVars[face().i];
        const VolumeVariables &volVarsJ = elemVolVars[face().j];
        const SpatialParams &spatialParams = problem.spatialParams();
        
        std::vector<Scalar> density(face().numFap);
        // to calcluate the integration point pressure we also need all the permeabilities
        //std::cout << faceIdx_ << std::endl;
        std::vector<Scalar> transmissibilitiesP(this->face().numFap);
        Scalar tPSum = 0.0;

        // to calcluate the integration point concentration we also need all the diffcoeff
        std::vector<Scalar> transmissibilitiesC(this->face().numFap);
        Scalar tCSum = 0.0;

        for (int fapIdx = 0; fapIdx < face().numFap; ++fapIdx)
        {
            {
                if(!isBox)
                {
                    if(onBoundary_)
                        {
                            transmissibilitiesP[fapIdx] = spatialParams.Kx(fvGeometry_().neighbors[0], fvGeometry_(), 0)/face().fapDistances[fapIdx];
                            transmissibilitiesC[fapIdx] = elemVolVars[0].diffCoeff()/this->face().fapDistances[fapIdx];
                        }
                    else
                        {
                            transmissibilitiesP[fapIdx] = spatialParams.Kx(fvGeometry_().neighbors[face().fapIndices[fapIdx]], fvGeometry_(), 0)/face().fapDistances[fapIdx];
                            transmissibilitiesC[fapIdx] = elemVolVars[face().fapIndices[fapIdx]].diffCoeff()/face().fapDistances[fapIdx];
                        }
                }
                else
                {
                    transmissibilitiesP[fapIdx] = spatialParams.Kx(element, fvGeometry_(), 0)
                                                                                      /face().fapDistances[fapIdx];
                    transmissibilitiesC[fapIdx] = elemVolVars[face().fapIndices[fapIdx]].diffCoeff()/face().fapDistances[fapIdx];
                }
            }
            tPSum += transmissibilitiesP[fapIdx];
            tCSum += transmissibilitiesC[fapIdx];
            int volVarsIdx = face().fapIndices[fapIdx];
            density[fapIdx] = elemVolVars[volVarsIdx].density();
        }

        if(onBoundary_)
        {
            transmissibility_ = transmissibilitiesP[0];
            transmissibilityC_ = transmissibilitiesC[0];
            density_ = density[0];
        }
        else
        {
            transmissibility_ = transmissibilitiesP[0]*transmissibilitiesP[1]/tPSum;
            transmissibilityC_ = transmissibilitiesC[0]*transmissibilitiesC[1]/tCSum;
            density_ = std::accumulate(density.begin(), density.end(), 0.0)/density.size();
        }

        // calculate the pressure gradient
        diffP_ = volVarsI.pressure() - volVarsJ.pressure();

        // correct the pressure diffrence by the gravitational acceleration
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
        {   // calculate the density
            const Element& elementI = fvGeometry_().neighbors[face().i];
            const Element& elementJ = fvGeometry_().neighbors[face().j];
            auto globalPosI = elementI.geometry().center();
            auto globalPosJ = elementJ.geometry().center();
            Scalar gravitationalTermI = density_*(problem.gravityAtPos(globalPosI)*globalPosI);
            Scalar gravitationalTermJ;
            if (onBoundary_)
            {
                gravitationalTermJ = density_*(problem.gravityAtPos(face().ipGlobal)*face().ipGlobal);
            }
            else
            {
                gravitationalTermJ = density_*(problem.gravityAtPos(globalPosJ)*globalPosJ);
            }
            diffP_ -= (gravitationalTermI - gravitationalTermJ);
        }

        // calculate the conc. gradient
        if(useMoles)
            diffC_ = volVarsI.moleFraction(transportCompIdx)*volVarsI.molarDensity() - volVarsJ.moleFraction(transportCompIdx)*volVarsJ.molarDensity();
        else
            diffC_ = volVarsI.massFraction(transportCompIdx)*volVarsI.density() - volVarsJ.massFraction(transportCompIdx)*volVarsJ.density();


    }

    /*!
    * \brief Calculation of the effective diffusion coefficient.
    *
    *        \param problem The considered problem file
    *        \param element The considered element of the grid
    *        \param elemVolVars The parameters stored in the considered element
    */
    void calculateDiffusiveFlux_(const Problem &problem,
                                   const Element &element,
                                   const ElementVolumeVariables &elemVolVars)
    {
        diffusiveFlux_ = transmissibilityC_*diffC_;
    }


    /*
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
        //get permeabily and distance from the neighbor element to calculate the velocity
        volumeFlux_ = transmissibility_*diffP_;
        // set the upstream and downstream vertices
        upstreamIdx_ = face().i;
        downstreamIdx_ = face().j;

        if (volumeFlux_ < 0)
            std::swap(upstreamIdx_, downstreamIdx_);
    }
      // return const reference to the fvGeometry
    const FVElementGeometry& fvGeometry_() const
    { return *fvGeometryPtr_; }


    unsigned int faceIdx_;            //!< The index of the sub control volume face
    bool onBoundary_;                //!< Specifying whether we are currently on the boundary of the simulation domain
    unsigned int upstreamIdx_ , downstreamIdx_; //!< local index of the upstream / downstream vertex
    Scalar volumeFlux_ ;    //!< Velocity multiplied with normal (magnitude=area)
    Scalar transmissibility_ ;                  //!< Transmissiblity i->j
    Scalar diffP_ ;                             //!< Pressure difference from cell i to cell j

    Scalar diffusiveFlux_;                      //!< Diffusive flux fo the transported component
    Scalar diffC_;                              //!< Concentration gradient
    Scalar transmissibilityC_ ;                 //!< Concentration on the face
    Scalar density_ ;

private:
    const FVElementGeometry* fvGeometryPtr_; //!< Information about the geometry of discretization

};

} // end namespace

#endif
