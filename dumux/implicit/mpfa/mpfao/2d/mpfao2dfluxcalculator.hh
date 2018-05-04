#ifndef DUMUX_MPFAO2DFLUXCALCULATOR_HH
#define DUMUX_MPFAO2DFLUXCALCULATOR_HH

#include <dumux/common/math.hh>
#include <eigen3/Eigen/Dense>

namespace Dumux{

template<class TypeTag>
class MpfaO2DFluxCalculator
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef Dumux::MpfaO2DInteractionVolume<TypeTag> InteractionVolume;

    typedef typename InteractionVolume::DynamicMatrix TransMatrix;
    typedef typename InteractionVolume::DynamicVector FluxVector;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Intersection Intersection;

    typedef typename GridView::Grid::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef typename Dune::ReferenceElements<Scalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<Scalar, dim> ReferenceElement;

public:

    // this routine is for compatibility reasons with other MPFA-methods, that use the mpfa-o method only for the boundaries
    static Scalar calculateScvFaceFlux(const Problem &problem, const Element &element, int localIsIndex,
                                        const InteractionVolume &interactionVolume, const ElementVolumeVariables &elemVolVars,
                                        const FVElementGeometry &fvGeometry, int phaseIdx, int &regionFaceIdx)
    {
        FluxVector neumannDummy = FluxVector::Zero(interactionVolume.getNumberOfSubFaces());
        return calculateScvFaceFlux(problem, element, localIsIndex, interactionVolume, elemVolVars, fvGeometry, phaseIdx, regionFaceIdx, neumannDummy);
    }

    static Scalar calculateScvFaceFlux(const Problem &problem, const Element &element, int localIsIndex,
                                        const InteractionVolume &interactionVolume, const ElementVolumeVariables &elemVolVars,
                                        const FVElementGeometry &fvGeometry, int phaseIdx, int &regionFaceIdx, FluxVector &neumannFluxes)
    {
        const TransMatrix &T = interactionVolume.getTransMatrix();
        // Set up vector with the primary variable and boundary condition values
        int noOfFaces = interactionVolume.getNumberOfSubFaces();
        int noOfSubVols = interactionVolume.getNumberOfSubVols();
        int noOfHeads = noOfSubVols + interactionVolume.getNumberOfBoundaryFaces();
        int noInteriorFaces = interactionVolume.getNumberOfInteriorFaces();
        FluxVector piezoHeads = FluxVector::Zero(noOfHeads);
        FluxVector fluxes = FluxVector::Zero(noOfFaces);

        // find the face index of the facet in the interaction region
        // and determine if the flux has to be multiplied by -1
        // boolean "switchSign" will be assigned with a value in function "getFaceIndexInRegion"
        bool switchSign;
        int faceIndexInRegion = interactionVolume.getFaceIndexInRegion(problem, element, localIsIndex, switchSign, regionFaceIdx);

        // set up piezometric head vector entries of the sub volumes
        for (int head = 0; head < noOfSubVols; head++)
        {
            const typename InteractionVolume::SubVolume& subVolume = interactionVolume.getSubVolume(head);
            int idxInStencil = fvGeometry.findElementInStencil(subVolume.globalIdx);
            Scalar pressure = elemVolVars[idxInStencil].fluidState().pressure(phaseIdx);

            // turn pressure into piezometric head
            if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            {
                // ask for the gravitational acceleration at the given point
                const GlobalPosition subVolCenter = subVolume.element.geometry().center();
                GlobalPosition g( problem.gravityAtPos(subVolCenter) );
                const Scalar density = elemVolVars[idxInStencil].fluidState().density(phaseIdx);
                // make gravity acceleration a force
                Scalar f = g * subVolCenter;
                f *= density;
                // calculate the final potential
                pressure -= f;
            }
            // set value in the piezometric heads vector
            piezoHeads(head) = pressure;
        }

        // set the entries in piezo head vector for the Dirichlet boundary conditions
        int dirichletFaceCounter = 0;
        const std::set<int> &dirichletFaceIndexSet = interactionVolume.getDirichletFaceIndexSet();
        std::set<int>::iterator it = dirichletFaceIndexSet.begin();
        for(; it != dirichletFaceIndexSet.end(); ++it)
        {
            const typename InteractionVolume::SubVolumeFace& subVolFace = interactionVolume.getSubVolumeFace(*it);

            // find corresponding intersection
            int posSubVolIndex = subVolFace.positiveSubVolIndex;
            int idxOnPosSubVol = subVolFace.localIdxOnPositive;
            const typename InteractionVolume::SubVolume& subVolume = interactionVolume.getSubVolume(posSubVolIndex);
            int facetIdx = subVolume.localIsIndex[idxOnPosSubVol];

            bool found = false;
            // find intersection
            IntersectionIterator it = problem.gridView().ibegin(subVolume.element);
            IntersectionIterator endIt = problem.gridView().iend(subVolume.element);
            for(; it != endIt; ++it)
                if (it->indexInInside() == facetIdx)
                    {found = true; break;}
            if (found == false)
                DUNE_THROW(Dune::InvalidStateException, "Error. Facet corresponding to index was not found!");

            // TODO!! Should the dirichlet boundary conditions be evaluated at the continuity points of the faces???
            Scalar pressure; PrimaryVariables values(0); typename VolumeVariables::FluidState fluidState;
            if (subVolFace.faceType == InteractionVolume::FaceTypes::DirichletFace)
            {
                problem.dirichlet(values, *it);
                VolumeVariables::completeFluidState(values, problem, element, fvGeometry, 0, fluidState);
                pressure = fluidState.pressure(phaseIdx);
            }
            else if ((subVolFace.faceType == InteractionVolume::FaceTypes::InternalDirichletFace))
            {
                problem.internalDirichlet(values, subVolume.element, facetIdx, interactionVolume);
                VolumeVariables::completeFluidState(values, problem, element, fvGeometry, 0, fluidState);
                pressure = fluidState.pressure(phaseIdx);
            }
            else
                DUNE_THROW(Dune::NotImplemented, "Dirichlet face has neither internal nor external dirichlet flag assigned to!");

            // evaluate gravity at the position of the continutiy point on the face
            if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            {
                // ask for the gravitational acceleration at the given point
                //const typename InteractionVolume::SubVolume& positiveSubVolume = interactionVolume.getSubVolume(subVolFace.positiveSubVolIndex);
                //const GlobalPosition contiPoint = positiveSubVolume.contiPoints[subVolFace.localIdxOnPositive];
                const GlobalPosition contiPoint = subVolume.contiPoints[idxOnPosSubVol];
                GlobalPosition g( problem.gravityAtPos(contiPoint) );
                const Scalar density = fluidState.density(phaseIdx);
                // make gravity acceleration a force
                Scalar f = g * contiPoint;
                f *= density;
                // calculate the final potential
                pressure -= f;
            }
            piezoHeads(noOfSubVols + dirichletFaceCounter) = pressure;
            dirichletFaceCounter++;
        }

        // set up the vector of eventual neumann fluxes
        int interiorFaceCounter = 0;
        FluxVector tmpNeumannFluxes = FluxVector::Zero(noInteriorFaces);
        const std::set<int> &interiorFaceIndexSet = interactionVolume.getInteriorFaceIndexSet();
        std::set<int>::iterator interiorIt = interiorFaceIndexSet.begin();
        for(; interiorIt != interiorFaceIndexSet.end(); ++interiorIt)
        {
            const typename InteractionVolume::SubVolumeFace& subVolFace = interactionVolume.getSubVolumeFace(*interiorIt);

            if (subVolFace.faceType != InteractionVolume::FaceTypes::NeumannFace
                    && subVolFace.faceType != InteractionVolume::FaceTypes::InternalFluxFace)
            {
                interiorFaceCounter++; continue;
            }
            else
            {
                // find corresponding intersection
                int posSubVolIndex = subVolFace.positiveSubVolIndex;
                int idxOnPosSubVol = subVolFace.localIdxOnPositive;
                const typename InteractionVolume::SubVolume& subVolume = interactionVolume.getSubVolume(posSubVolIndex);
                int facetIdx = subVolume.localIsIndex[idxOnPosSubVol];

                bool found = false;
                // find intersection
                IntersectionIterator it = problem.gridView().ibegin(subVolume.element);
                IntersectionIterator endIt = problem.gridView().iend(subVolume.element);
                for(; it != endIt; ++it)
                    if (it->indexInInside() == facetIdx)
                        {
                            found = true;
                            break;
                        }
                if (found == false)
                    DUNE_THROW(Dune::NotImplemented, "Error. Facet corresponding to index was not found!");

                if (subVolFace.faceType == InteractionVolume::FaceTypes::NeumannFace)
                {
                    PrimaryVariables values(0);
                    problem.neumann(values, it->inside(), fvGeometry, *it, /*dummy*/0, /*dummy*/0, interactionVolume);
                    values *= it->geometry().volume()/2;

                    int interiorFaceIdx = interactionVolume.getInteriorFaceIndex(*interiorIt);
                    tmpNeumannFluxes(interiorFaceIdx) = values[phaseIdx];
                    interiorFaceCounter++;
                }
                else if (subVolFace.faceType == InteractionVolume::FaceTypes::InternalFluxFace)
                {
                    PrimaryVariables values(0);
                    problem.internalFlux(values, it->inside(), facetIdx, interactionVolume);
                    values *= it->geometry().volume()/2;

                    int interiorFaceIdx = interactionVolume.getInteriorFaceIndex(*interiorIt);
                    tmpNeumannFluxes(interiorFaceIdx) = values[phaseIdx];
                    interiorFaceCounter++;
                }
            }
        }

        // Some error tracking... to be removed after testing!
        if (T.cols() != piezoHeads.size())
            DUNE_THROW(Dune::NotImplemented, "Error in mpfao2dfluxcalculator.hh, l. 123");
        if (T.rows() != fluxes.size())
            DUNE_THROW(Dune::NotImplemented, "T.rows() != fluxes.size()");
        if (interiorFaceCounter != noInteriorFaces)
            DUNE_THROW(Dune::NotImplemented, "interiorFaceCounter != noInteriorFaces");

        // calculate the fluxes over the sub control volume faces
        fluxes = T * piezoHeads;

        // more error tracking... to be removed after testing!
        if (interactionVolume.getCAInverse().rows() != noOfFaces)
            DUNE_THROW(Dune::NotImplemented, "interactionVolume.getCAInverse().rows() != noOfFaces");
        if (interactionVolume.getCAInverse().cols() != noInteriorFaces)
            DUNE_THROW(Dune::NotImplemented, "interactionVolume.getCAInverse().rows() != noInteriorFaces");

        // calculate the fluxes stemming from the neumann boundaries within the interaction volume
        if (noInteriorFaces > 0)
            neumannFluxes = interactionVolume.getCAInverse() * tmpNeumannFluxes;

        // return fluxes
        if (switchSign)
        {
            neumannFluxes *= -1;
            return -fluxes(faceIndexInRegion);
        }
        else
            return fluxes(faceIndexInRegion);
    }
};
}
#endif
